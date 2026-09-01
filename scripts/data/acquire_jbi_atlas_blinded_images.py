#!/usr/bin/env python3
"""Acquire one deterministic atlas shard and expose only renamed image objects."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import time
from typing import Any
from urllib.request import Request, urlopen

from PIL import Image


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.atlas_measurement import (
    measurement_shard,
    validate_execution_contract,
    validate_inference_contract,
)
from fcp_pipeline.flower_roi_v4_runtime import validate_scaleout_authorization
from scripts.data.validate_jbi_atlas_roi_v4_gate_evidence import (
    load_committed_locked_scaleout_result,
)


USER_AGENT = "FCP-image-first-atlas/3.0 (research image acquisition)"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def valid_image(path: Path) -> bool:
    try:
        with Image.open(path) as image:
            image.verify()
        return True
    except (OSError, ValueError):
        return False


def acquire(url: str, destination: Path, *, retries: int, timeout: float) -> str:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.is_file() and valid_image(destination):
        return sha256(destination)
    partial = destination.with_suffix(destination.suffix + ".partial")
    error: Exception | None = None
    for attempt in range(retries):
        try:
            request = Request(url, headers={"User-Agent": USER_AGENT})
            with urlopen(request, timeout=timeout) as response, partial.open("wb") as handle:
                if int(getattr(response, "status", 200)) != 200:
                    raise OSError(f"HTTP status {response.status}")
                while block := response.read(1024 * 1024):
                    handle.write(block)
            if not valid_image(partial):
                raise OSError("downloaded object is not a decodable image")
            partial.replace(destination)
            return sha256(destination)
        except Exception as exc:  # network and image decoders expose many subclasses
            error = exc
            if partial.exists():
                partial.unlink()
            if attempt + 1 < retries:
                time.sleep(min(2**attempt, 8))
    raise RuntimeError(str(error) if error else "image acquisition failed")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sealed-acquisition-key", type=Path, required=True)
    parser.add_argument("--firewall-manifest", type=Path, required=True)
    parser.add_argument("--roi-evidence-dir", type=Path, required=True)
    parser.add_argument(
        "--inference-contract",
        type=Path,
        default=Path("docs/supporting/jbi_image_first_atlas_inference_contract_v3.json"),
    )
    parser.add_argument(
        "--execution-contract",
        type=Path,
        default=Path("docs/supporting/jbi_atlas_scaleout_execution_contract_v1.json"),
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--shard-index", type=int, required=True)
    parser.add_argument("--shard-count", type=int, required=True)
    parser.add_argument("--retries", type=int, default=4)
    parser.add_argument("--timeout-seconds", type=float, default=60.0)
    parser.add_argument("--pause-seconds", type=float, default=0.05)
    return parser.parse_args()


def validate_firewall_for_acquisition(
    firewall: dict[str, Any], *, key_name: str, key_sha256: str
) -> None:
    if (
        firewall.get("status") != "pass_scaleout_measurement_firewall"
        or firewall.get("candidate_image_pixels_opened") is not False
        or firewall.get("coordinate_key_opened_by_measurement_worker") is not False
        or firewall.get("sealed_keys", {}).get(key_name) != key_sha256
        or not firewall.get("dated_source_gate_sha256")
        or not firewall.get("environmental_coverage_gate_sha256")
    ):
        raise RuntimeError("measurement firewall does not authorize this acquisition key")


def main() -> None:
    args = parse_args()
    if args.shard_index < 0 or args.shard_index >= args.shard_count:
        raise ValueError("shard_index must lie in [0, shard_count)")
    inference = json.loads(args.inference_contract.read_text(encoding="utf-8"))
    validate_inference_contract(inference)
    execution = json.loads(args.execution_contract.read_text(encoding="utf-8"))
    validate_execution_contract(execution)
    frozen = execution["acquisition"]
    if (
        args.shard_count != frozen["shard_count"]
        or args.retries != frozen["retries_per_image"]
        or args.timeout_seconds != frozen["timeout_seconds"]
        or args.pause_seconds != frozen["pause_seconds_after_each_terminal_image"]
    ):
        raise ValueError(
            "acquisition worker settings differ from the frozen execution contract"
        )
    roi = load_committed_locked_scaleout_result(args.roi_evidence_dir)
    validate_scaleout_authorization(roi)
    firewall = json.loads(args.firewall_manifest.read_text(encoding="utf-8"))
    validate_firewall_for_acquisition(
        firewall,
        key_name=args.sealed_acquisition_key.name,
        key_sha256=sha256(args.sealed_acquisition_key),
    )

    rows = read_csv(args.sealed_acquisition_key)
    expected = 60_000
    if len(rows) != expected or len({row["measurement_id"] for row in rows}) != expected:
        raise RuntimeError("sealed acquisition key does not match the 60,000-image denominator")
    selected = [
        row
        for row in rows
        if measurement_shard(row["measurement_id"], args.shard_count) == args.shard_index
    ]
    if not selected:
        raise RuntimeError("acquisition shard is empty; reduce shard_count")
    images_dir = args.output_dir / "images"
    records: list[dict[str, Any]] = []
    for index, row in enumerate(selected, start=1):
        measurement_id = row["measurement_id"]
        destination = images_dir / f"{measurement_id}.jpg"
        try:
            image_hash = acquire(
                row["photo_url_large"],
                destination,
                retries=args.retries,
                timeout=args.timeout_seconds,
            )
            status = "image_acquired"
            reason = ""
        except RuntimeError as exc:
            image_hash = ""
            status = "image_acquisition_failed"
            reason = str(exc)[:500]
        records.append(
            {
                "measurement_id": measurement_id,
                "image_filename": f"{measurement_id}.jpg",
                "acquisition_status": status,
                "image_sha256": image_hash,
                "failure_reason": reason,
            }
        )
        if args.pause_seconds > 0:
            time.sleep(args.pause_seconds)
        if index % 25 == 0 or index == len(selected):
            print(f"acquired_or_failed={index}/{len(selected)}", flush=True)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.output_dir / f"acquisition_shard_{args.shard_index:04d}.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)
    manifest = {
        "status": "complete_blinded_acquisition_shard",
        "protocol": inference["protocol"],
        "execution_protocol": execution["protocol"],
        "shard_index": args.shard_index,
        "shard_count": args.shard_count,
        "frozen_shard_denominator": len(selected),
        "images_acquired": sum(row["acquisition_status"] == "image_acquired" for row in records),
        "images_failed": sum(
            row["acquisition_status"] == "image_acquisition_failed" for row in records
        ),
        "measurement_worker_exposed_source_urls": False,
        "measurement_worker_exposed_coordinates": False,
        "result_sha256": sha256(csv_path),
    }
    (args.output_dir / f"acquisition_shard_{args.shard_index:04d}.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

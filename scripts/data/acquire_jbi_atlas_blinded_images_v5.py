#!/usr/bin/env python3
"""Acquire one terminal atlas shard only after the v5 pre-image firewall passes."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import time
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.atlas_acquisition_v5 import validate_v5_acquisition_firewall
from fcp_pipeline.atlas_measurement_v5 import measurement_shard, validate_measurement_execution_contract
from fcp_pipeline.flower_roi_v4_runtime import validate_scaleout_authorization
from scripts.data.acquire_jbi_atlas_blinded_images import acquire
from scripts.data.build_jbi_atlas_measurement_firewall_v5 import verify_repo_parent_blobs
from scripts.data.validate_jbi_atlas_roi_v4_gate_evidence import load_committed_locked_scaleout_result


DEFAULT_MEASUREMENT_CONTRACT = ROOT / "docs/supporting/jbi_atlas_measurement_execution_contract_v5.json"
DEFAULT_INFERENCE = ROOT / "docs/supporting/jbi_image_first_atlas_inference_contract_v5.json"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sealed-acquisition-key", type=Path, required=True)
    parser.add_argument("--firewall-manifest", type=Path, required=True)
    parser.add_argument("--roi-evidence-dir", type=Path, required=True)
    parser.add_argument("--measurement-contract", type=Path, default=DEFAULT_MEASUREMENT_CONTRACT)
    parser.add_argument("--inference-v5", type=Path, default=DEFAULT_INFERENCE)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--shard-index", type=int, required=True)
    parser.add_argument("--shard-count", type=int, default=16)
    parser.add_argument("--retries", type=int, default=4)
    parser.add_argument("--timeout-seconds", type=float, default=60.0)
    parser.add_argument("--pause-seconds", type=float, default=0.05)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.shard_index < 0 or args.shard_index >= args.shard_count:
        raise ValueError("shard_index must lie in [0, shard_count)")

    contract = load_json(args.measurement_contract)
    inference = load_json(args.inference_v5)
    validate_measurement_execution_contract(contract, inference)
    verify_repo_parent_blobs(contract)
    execution = contract["technical_execution"]
    if (
        args.shard_count != execution["acquisition_shard_count"]
        or args.retries != execution["acquisition_retries_per_image"]
        or args.timeout_seconds != execution["acquisition_timeout_seconds"]
        or args.pause_seconds != execution["acquisition_pause_seconds_after_each_terminal_image"]
    ):
        raise ValueError("v5 acquisition settings differ from the frozen contract")

    roi = load_committed_locked_scaleout_result(args.roi_evidence_dir)
    validate_scaleout_authorization(roi)

    firewall = load_json(args.firewall_manifest)
    key_hash = sha256(args.sealed_acquisition_key)
    validate_v5_acquisition_firewall(
        firewall,
        key_name=args.sealed_acquisition_key.name,
        key_sha256=key_hash,
        contract=contract,
        inference_v5=inference,
    )

    rows = read_csv(args.sealed_acquisition_key)
    if len(rows) != 60000 or len({row["measurement_id"] for row in rows}) != 60000:
        raise RuntimeError("v5 sealed acquisition key does not preserve the 60,000-image denominator")
    selected = [
        row for row in rows
        if measurement_shard(row["measurement_id"], args.shard_count) == args.shard_index
    ]
    if not selected:
        raise RuntimeError("v5 acquisition shard is empty")

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
    csv_path = args.output_dir / f"acquisition_v5_shard_{args.shard_index:04d}.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)

    manifest = {
        "status": "complete_blinded_acquisition_v5_shard",
        "protocol": contract["protocol"],
        "inference_version": inference["version"],
        "superseded_v3_ordered_inference_used": False,
        "shard_index": args.shard_index,
        "shard_count": args.shard_count,
        "frozen_shard_denominator": len(selected),
        "images_acquired": sum(row["acquisition_status"] == "image_acquired" for row in records),
        "images_failed": sum(row["acquisition_status"] == "image_acquisition_failed" for row in records),
        "measurement_worker_exposed_source_urls": False,
        "measurement_worker_exposed_coordinates": False,
        "firewall_manifest_sha256": sha256(args.firewall_manifest),
        "sealed_acquisition_key_sha256": key_hash,
        "result_sha256": sha256(csv_path),
    }
    (args.output_dir / f"acquisition_v5_shard_{args.shard_index:04d}.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

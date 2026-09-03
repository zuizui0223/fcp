#!/usr/bin/env python3
"""Acquire one fresh photo-first compute partition into blinded ephemeral files.

This is the only measurement-stage process that receives source URLs. Its durable
receipt omits URLs and source metadata. Successful images are written only to the
requested ephemeral directory using blinded filenames; download/decode failures
become terminal ``mixed_uncertain`` rows and are never replaced.
"""

from __future__ import annotations

import argparse
from io import BytesIO
import hashlib
import json
from pathlib import Path
import time
from urllib.request import Request, urlopen

import pandas as pd
from PIL import Image

from fcp_pipeline.photo_first_measurement_execution import (
    select_measurement_partition,
    validate_execution_contract,
    validate_terminal_partition_results,
)


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MEASUREMENT = ROOT / "docs/supporting/random_photo_first_measurement_contract_v1.json"
DEFAULT_EXECUTION = ROOT / "docs/supporting/random_photo_first_measurement_execution_v1.json"
USER_AGENT = "fcp-random-photo-first-measurement/1.0 (github.com/zuizui0223/fcp)"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--worker-manifest", type=Path, required=True)
    parser.add_argument("--acquisition-key", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--semantic-shard", type=int, required=True)
    parser.add_argument("--compute-partition", type=int, required=True)
    parser.add_argument("--measurement-contract", type=Path, default=DEFAULT_MEASUREMENT)
    parser.add_argument("--execution-contract", type=Path, default=DEFAULT_EXECUTION)
    return parser.parse_args()


def _download(url: str, *, retries: int, timeout: float) -> tuple[bytes | None, str | None]:
    last_error = None
    for attempt in range(retries + 1):
        try:
            request = Request(url, headers={"User-Agent": USER_AGENT, "Accept": "image/*"})
            with urlopen(request, timeout=timeout) as response:
                payload = response.read()
            if not payload:
                raise ValueError("empty_image_payload")
            with Image.open(BytesIO(payload)) as image:
                image.verify()
            return payload, None
        except Exception as exc:  # network/codec failures all become frozen terminal missingness
            last_error = f"{type(exc).__name__}:{str(exc)[:240]}"
            if attempt < retries:
                time.sleep(min(16.0, 2.0**attempt))
    return None, last_error or "unknown_acquisition_failure"


def main() -> int:
    args = parse_args()
    measurement_contract = json.loads(args.measurement_contract.read_text(encoding="utf-8"))
    execution_contract = json.loads(args.execution_contract.read_text(encoding="utf-8"))
    validate_execution_contract(measurement_contract, execution_contract)
    worker = pd.read_csv(args.worker_manifest, dtype=str).fillna("")
    acquisition = pd.read_csv(args.acquisition_key, dtype=str).fillna("")
    selected = select_measurement_partition(
        worker,
        acquisition,
        semantic_shard_index=args.semantic_shard,
        compute_partition_index=args.compute_partition,
        semantic_shards=int(execution_contract["partitioning"]["semantic_shards"]),
        compute_partitions_per_shard=int(
            execution_contract["partitioning"]["compute_partitions_per_semantic_shard"]
        ),
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    image_dir = args.output_dir / "ephemeral_images"
    image_dir.mkdir(parents=True, exist_ok=True)
    success_rows = []
    receipt_rows = []
    failures = []
    acquisition_rule = execution_contract["acquisition"]
    retries = int(acquisition_rule["retries_per_image"])
    timeout = float(acquisition_rule["timeout_seconds"])

    for position, row in enumerate(selected.itertuples(index=False), start=1):
        measurement = str(row.measurement_id)
        payload, error = _download(str(row.photo_url_large), retries=retries, timeout=timeout)
        if payload is None:
            failures.append(
                {
                    "measurement_id": measurement,
                    "morph": "mixed_uncertain",
                    "measurement_status": "image_acquisition_failed",
                    "roi_status": "not_opened_acquisition_failed",
                    "failure_reasons": error or "image_acquisition_failed",
                    "image_sha256": "",
                    "flower_effective_pixels": 0,
                }
            )
            receipt_rows.append(
                {
                    "measurement_id": measurement,
                    "acquisition_status": "image_acquisition_failed",
                    "image_sha256": "",
                    "image_bytes": 0,
                }
            )
        else:
            digest = hashlib.sha256(payload).hexdigest()
            destination = image_dir / str(row.image_filename)
            destination.write_bytes(payload)
            success_rows.append(
                {
                    "measurement_id": measurement,
                    "image_filename": str(row.image_filename),
                    "photo_license": str(row.photo_license),
                }
            )
            receipt_rows.append(
                {
                    "measurement_id": measurement,
                    "acquisition_status": "acquired_and_decode_verified",
                    "image_sha256": digest,
                    "image_bytes": len(payload),
                }
            )
        if position % 25 == 0 or position == len(selected):
            print(f"acquired_or_failed={position}/{len(selected)}", flush=True)

    success = pd.DataFrame(
        success_rows,
        columns=["measurement_id", "image_filename", "photo_license"],
    )
    failure = pd.DataFrame(
        failures,
        columns=[
            "measurement_id",
            "morph",
            "measurement_status",
            "roi_status",
            "failure_reasons",
            "image_sha256",
            "flower_effective_pixels",
        ],
    )
    receipt = pd.DataFrame(
        receipt_rows,
        columns=["measurement_id", "acquisition_status", "image_sha256", "image_bytes"],
    )
    selected_blind = selected[["measurement_id", "image_filename", "photo_license"]].copy()

    if len(failure):
        validate_terminal_partition_results(
            failure, failure["measurement_id"].astype(str).tolist()
        )
    success.to_csv(args.output_dir / "successful_worker_manifest.csv", index=False, lineterminator="\n")
    failure.to_csv(args.output_dir / "acquisition_failures.csv", index=False, lineterminator="\n")
    receipt.to_csv(args.output_dir / "acquisition_receipt.csv", index=False, lineterminator="\n")
    selected_blind.to_csv(args.output_dir / "selected_blind_manifest.csv", index=False, lineterminator="\n")

    manifest = {
        "protocol": execution_contract["protocol"],
        "status": "complete_blinded_partition_acquisition",
        "semantic_shard": int(args.semantic_shard),
        "compute_partition": int(args.compute_partition),
        "selected_rows": int(len(selected)),
        "acquired_rows": int(len(success)),
        "acquisition_failed_rows": int(len(failure)),
        "source_urls_persisted_in_receipt": False,
        "species_persisted_in_receipt": False,
        "coordinates_persisted_in_receipt": False,
        "replacement_after_failure": False,
        "candidate_images_are_ephemeral": True,
    }
    (args.output_dir / "acquisition_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

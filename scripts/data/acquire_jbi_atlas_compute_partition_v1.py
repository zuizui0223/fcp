#!/usr/bin/env python3
"""Acquire one deterministic terminal compute partition under the v5 preimage firewall."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import time
from typing import Any, Mapping, Sequence

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.atlas_acquisition_v5 import validate_v5_acquisition_firewall
from fcp_pipeline.atlas_compute_partition import (
    select_compute_partition,
    validate_compute_partition_contract,
)
from fcp_pipeline.atlas_measurement_v5 import validate_measurement_execution_contract
from fcp_pipeline.flower_roi_v4_runtime import validate_scaleout_authorization
from scripts.data.acquire_jbi_atlas_blinded_images import acquire
from scripts.data.build_jbi_atlas_measurement_firewall_v5 import verify_repo_parent_blobs
from scripts.data.validate_jbi_atlas_roi_v4_gate_evidence import load_committed_locked_scaleout_result

DEFAULT_COMPUTE = ROOT / "docs/supporting/jbi_atlas_compute_partition_amendment_v1.json"
DEFAULT_MEASUREMENT = ROOT / "docs/supporting/jbi_atlas_measurement_execution_contract_v5.json"
DEFAULT_INFERENCE = ROOT / "docs/supporting/jbi_image_first_atlas_inference_contract_v5.json"


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path} must contain one JSON object")
    return value


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: Sequence[Mapping[str, Any]], fields: Sequence[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        if not rows:
            raise ValueError(f"cannot infer fields for empty table: {path}")
        fields = list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_blob_sha(path: Path) -> str:
    payload = path.read_bytes()
    return hashlib.sha1(f"blob {len(payload)}\0".encode("ascii") + payload).hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sealed-acquisition-key", type=Path, required=True)
    parser.add_argument("--measurement-manifest", type=Path, required=True)
    parser.add_argument("--firewall-manifest", type=Path, required=True)
    parser.add_argument("--roi-evidence-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--semantic-shard-index", type=int, required=True)
    parser.add_argument("--compute-partition-index", type=int, required=True)
    parser.add_argument("--compute-contract", type=Path, default=DEFAULT_COMPUTE)
    parser.add_argument("--measurement-contract", type=Path, default=DEFAULT_MEASUREMENT)
    parser.add_argument("--inference-v5", type=Path, default=DEFAULT_INFERENCE)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    compute_contract = load_json(args.compute_contract)
    measurement_contract = load_json(args.measurement_contract)
    inference = load_json(args.inference_v5)
    validate_compute_partition_contract(compute_contract)
    validate_measurement_execution_contract(measurement_contract, inference)
    verify_repo_parent_blobs(measurement_contract)
    if compute_contract["immutable_parent"]["git_blob_sha"] != git_blob_sha(args.measurement_contract):
        raise RuntimeError("compute partition parent measurement contract changed")

    roi = load_committed_locked_scaleout_result(args.roi_evidence_dir)
    validate_scaleout_authorization(roi)
    firewall = load_json(args.firewall_manifest)
    key_hash = sha256(args.sealed_acquisition_key)
    validate_v5_acquisition_firewall(
        firewall,
        key_name=args.sealed_acquisition_key.name,
        key_sha256=key_hash,
        contract=measurement_contract,
        inference_v5=inference,
    )

    sealed_rows = read_csv(args.sealed_acquisition_key)
    worker_rows = read_csv(args.measurement_manifest)
    if len(sealed_rows) != 60000 or len(worker_rows) != 60000:
        raise RuntimeError("compute acquisition requires the exact 60,000-row preimage packet")
    if len({row["measurement_id"] for row in sealed_rows}) != 60000 or len({row["measurement_id"] for row in worker_rows}) != 60000:
        raise RuntimeError("preimage packet measurement IDs are not unique")

    private_partition = select_compute_partition(
        sealed_rows,
        semantic_shard_index=args.semantic_shard_index,
        compute_partition_index=args.compute_partition_index,
    )
    worker_partition = select_compute_partition(
        worker_rows,
        semantic_shard_index=args.semantic_shard_index,
        compute_partition_index=args.compute_partition_index,
    )
    private_ids = {row["measurement_id"] for row in private_partition}
    worker_ids = {row["measurement_id"] for row in worker_partition}
    if private_ids != worker_ids or not private_ids:
        raise RuntimeError("private and worker compute-partition membership differs or is empty")
    allowed_fields = list(measurement_contract["location_blind_measurement"]["measurement_worker_allowed_fields"])
    if any(set(row) != set(allowed_fields) for row in worker_partition):
        raise RuntimeError("location-blind compute worker packet fields changed")

    worker_path = args.output_dir / "worker_packet" / "measurement_manifest.csv"
    write_csv(worker_path, sorted(worker_partition, key=lambda row: row["measurement_id"]), allowed_fields)
    images_dir = args.output_dir / "ephemeral_images"
    images_dir.mkdir(parents=True, exist_ok=True)
    receipts: list[dict[str, Any]] = []
    execution = measurement_contract["technical_execution"]
    for index, row in enumerate(sorted(private_partition, key=lambda item: item["measurement_id"]), start=1):
        measurement_id = row["measurement_id"]
        destination = images_dir / f"{measurement_id}.jpg"
        try:
            image_hash = acquire(
                row["photo_url_large"],
                destination,
                retries=int(execution["acquisition_retries_per_image"]),
                timeout=float(execution["acquisition_timeout_seconds"]),
            )
            status = "image_acquired"
            reason = ""
        except RuntimeError as exc:
            image_hash = ""
            status = "image_acquisition_failed"
            reason = str(exc)[:500]
        receipts.append(
            {
                "measurement_id": measurement_id,
                "image_filename": f"{measurement_id}.jpg",
                "acquisition_status": status,
                "image_sha256": image_hash,
                "failure_reason": reason,
            }
        )
        pause = float(execution["acquisition_pause_seconds_after_each_terminal_image"])
        if pause > 0:
            time.sleep(pause)
        if index % 25 == 0 or index == len(private_partition):
            print(f"compute_partition_acquired_or_failed={index}/{len(private_partition)}", flush=True)

    receipt_path = args.output_dir / "acquisition_receipt.csv"
    write_csv(receipt_path, receipts)
    manifest = {
        "status": "complete_blinded_acquisition_v5_compute_partition",
        "protocol": measurement_contract["protocol"],
        "inference_version": inference["version"],
        "semantic_shard_index": args.semantic_shard_index,
        "compute_partition_index": args.compute_partition_index,
        "semantic_shard_count": 16,
        "compute_partitions_per_semantic_shard": 16,
        "frozen_partition_denominator": len(private_partition),
        "images_acquired": sum(row["acquisition_status"] == "image_acquired" for row in receipts),
        "images_failed": sum(row["acquisition_status"] == "image_acquisition_failed" for row in receipts),
        "worker_packet_sha256": sha256(worker_path),
        "acquisition_receipt_sha256": sha256(receipt_path),
        "compute_contract_git_blob_sha": git_blob_sha(args.compute_contract),
        "private_preimage_packet_must_be_deleted_before_measurement": True,
        "persist_after_partition_forbidden": ["ephemeral_images", "sealed acquisition key", "source URL", "coordinates", "taxon name"],
    }
    manifest_path = args.output_dir / "acquisition_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

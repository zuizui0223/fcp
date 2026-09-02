#!/usr/bin/env python3
"""Convert one v5 measurement worker output into immutable compute-partition evidence."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.atlas_compute_partition import (
    PARTITIONS_PER_SEMANTIC_SHARD,
    SEMANTIC_SHARDS,
    compute_partition_coordinates,
    validate_compute_partition_contract,
)

DEFAULT_COMPUTE = ROOT / "docs/supporting/jbi_atlas_compute_partition_amendment_v1.json"


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path} must contain one JSON object")
    return value


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


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
    parser.add_argument("--worker-result-csv", type=Path, required=True)
    parser.add_argument("--worker-result-manifest", type=Path, required=True)
    parser.add_argument("--semantic-shard-index", type=int, required=True)
    parser.add_argument("--compute-partition-index", type=int, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--compute-contract", type=Path, default=DEFAULT_COMPUTE)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    compute_contract = load_json(args.compute_contract)
    validate_compute_partition_contract(compute_contract)
    if not 0 <= args.semantic_shard_index < SEMANTIC_SHARDS or not 0 <= args.compute_partition_index < PARTITIONS_PER_SEMANTIC_SHARD:
        raise ValueError("compute partition coordinates must lie in 0..15")

    rows = read_csv(args.worker_result_csv)
    worker = load_json(args.worker_result_manifest)
    if not rows:
        raise RuntimeError("compute measurement result is empty")
    if (
        worker.get("status") != "complete_location_blind_roi_v4_measurement_v5_shard"
        or worker.get("shard_index") != args.semantic_shard_index
        or worker.get("shard_count") != SEMANTIC_SHARDS
        or worker.get("coordinates_opened") is not False
        or worker.get("taxon_names_opened") is not False
        or worker.get("terminal_records") != len(rows)
        or worker.get("result_sha256") != sha256(args.worker_result_csv)
    ):
        raise RuntimeError("location-blind worker result does not match the target semantic shard")
    ids = [str(row.get("measurement_id") or "") for row in rows]
    if len(set(ids)) != len(ids) or "" in ids:
        raise RuntimeError("compute measurement result IDs are empty or duplicated")
    wrong = [
        measurement_id
        for measurement_id in ids
        if compute_partition_coordinates(measurement_id)
        != (args.semantic_shard_index, args.compute_partition_index)
    ]
    if wrong:
        raise RuntimeError(f"worker result contains rows from another compute partition: {wrong[:10]}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    stem = f"measurement_v5_partition_s{args.semantic_shard_index:02d}_p{args.compute_partition_index:02d}"
    csv_path = args.output_dir / f"{stem}.csv"
    shutil.copyfile(args.worker_result_csv, csv_path)
    manifest = {
        "status": "complete_location_blind_roi_v4_measurement_v5_compute_partition",
        "protocol": worker["protocol"],
        "inference_version": worker["inference_version"],
        "semantic_shard_index": args.semantic_shard_index,
        "compute_partition_index": args.compute_partition_index,
        "semantic_shard_count": SEMANTIC_SHARDS,
        "compute_partitions_per_semantic_shard": PARTITIONS_PER_SEMANTIC_SHARD,
        "terminal_records": len(rows),
        "coordinates_opened": False,
        "taxon_names_opened": False,
        "result_sha256": sha256(csv_path),
        "model_id": worker["model_id"],
        "trained_weight_sha256": worker["trained_weight_sha256"],
        "roi_contract_sha256_lf_canonical_v1": worker["roi_contract_sha256_lf_canonical_v1"],
        "compute_contract_git_blob_sha": git_blob_sha(args.compute_contract),
        "worker_manifest_sha256": sha256(args.worker_result_manifest),
        "image_pixels_persisted": False,
        "source_urls_persisted": False,
        "coordinates_persisted": False,
    }
    json_path = args.output_dir / f"{stem}.json"
    json_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

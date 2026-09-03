#!/usr/bin/env python3
"""Build the location-blind worker packet after the fixed metadata capacity gate."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from fcp_pipeline.photo_first_measurement_execution import (
    build_measurement_firewall,
    compute_partition,
    semantic_shard,
)


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MEASUREMENT = ROOT / "docs/supporting/random_photo_first_measurement_contract_v1.json"
DEFAULT_EXECUTION = ROOT / "docs/supporting/random_photo_first_measurement_execution_v1.json"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-csv", type=Path, required=True)
    parser.add_argument("--candidate-manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--measurement-contract", type=Path, default=DEFAULT_MEASUREMENT)
    parser.add_argument("--execution-contract", type=Path, default=DEFAULT_EXECUTION)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    candidate = pd.read_csv(args.candidate_csv)
    candidate_manifest = json.loads(args.candidate_manifest.read_text(encoding="utf-8"))
    measurement_contract = json.loads(args.measurement_contract.read_text(encoding="utf-8"))
    execution_contract = json.loads(args.execution_contract.read_text(encoding="utf-8"))
    firewall = build_measurement_firewall(
        candidate,
        candidate_manifest,
        measurement_contract,
        execution_contract,
    )

    worker_dir = args.output_dir / "worker_packet"
    sealed_dir = args.output_dir / "sealed_keys"
    worker_dir.mkdir(parents=True, exist_ok=True)
    sealed_dir.mkdir(parents=True, exist_ok=True)
    worker_path = worker_dir / "measurement_manifest.csv"
    acquisition_path = sealed_dir / "acquisition_key.csv"
    metadata_path = sealed_dir / "metadata_join_key.csv"
    inventory_path = args.output_dir / "partition_inventory.csv"
    manifest_path = args.output_dir / "measurement_firewall_manifest.json"

    firewall.worker_manifest.to_csv(worker_path, index=False, lineterminator="\n")
    firewall.acquisition_key.to_csv(acquisition_path, index=False, lineterminator="\n")
    firewall.metadata_join_key.to_csv(metadata_path, index=False, lineterminator="\n")

    partition = execution_contract["partitioning"]
    n_shards = int(partition["semantic_shards"])
    n_parts = int(partition["compute_partitions_per_semantic_shard"])
    counts = {(shard, part): 0 for shard in range(n_shards) for part in range(n_parts)}
    for measurement in firewall.worker_manifest["measurement_id"].astype(str):
        key = (semantic_shard(measurement, n_shards), compute_partition(measurement, n_parts))
        counts[key] += 1
    inventory = pd.DataFrame(
        [
            {
                "semantic_shard": shard,
                "compute_partition": part,
                "candidate_rows": counts[(shard, part)],
            }
            for shard in range(n_shards)
            for part in range(n_parts)
        ]
    )
    inventory.to_csv(inventory_path, index=False, lineterminator="\n")

    manifest = dict(firewall.manifest)
    manifest.update(
        {
            "worker_manifest_sha256": sha256(worker_path),
            "acquisition_key_sha256": sha256(acquisition_path),
            "metadata_join_key_sha256": sha256(metadata_path),
            "partition_inventory_sha256": sha256(inventory_path),
            "total_compute_partitions": int(len(inventory)),
            "empty_compute_partitions": int((inventory["candidate_rows"] == 0).sum()),
            "minimum_partition_rows": int(inventory["candidate_rows"].min()),
            "maximum_partition_rows": int(inventory["candidate_rows"].max()),
            "candidate_pixels_opened": False,
        }
    )
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

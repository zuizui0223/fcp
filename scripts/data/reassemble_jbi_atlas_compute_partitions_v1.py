#!/usr/bin/env python3
"""Reassemble 256 location-blind compute partitions into the original 16 v5 shards."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
from typing import Any, Mapping, Sequence

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.atlas_compute_partition import (
    PARTITIONS_PER_SEMANTIC_SHARD,
    SEMANTIC_SHARDS,
    validate_compute_partition_contract,
    validate_partition_coverage,
)
from fcp_pipeline.atlas_measurement_v5 import validate_measurement_execution_contract
from scripts.data.build_jbi_atlas_measurement_firewall_v5 import verify_repo_parent_blobs

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


def write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        raise RuntimeError(f"cannot reassemble an empty semantic shard: {path.name}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    union = {key for row in rows for key in row}
    fields.extend(sorted(union - set(fields)))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
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
    parser.add_argument("--measurement-manifest", type=Path, required=True)
    parser.add_argument("--partition-results-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
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

    expected_rows = read_csv(args.measurement_manifest)
    expected_ids = [str(row["measurement_id"]) for row in expected_rows]
    if len(expected_ids) != 60000 or len(set(expected_ids)) != 60000:
        raise RuntimeError("reassembly requires the exact 60,000-row worker manifest")

    partition_rows: dict[tuple[int, int], list[dict[str, str]]] = {}
    partition_evidence: dict[tuple[int, int], dict[str, Any]] = {}
    model_ids: set[str] = set()
    weights: set[str] = set()
    roi_contracts: set[str] = set()
    for semantic in range(SEMANTIC_SHARDS):
        for partition in range(PARTITIONS_PER_SEMANTIC_SHARD):
            stem = f"measurement_v5_partition_s{semantic:02d}_p{partition:02d}"
            csv_path = args.partition_results_dir / f"{stem}.csv"
            json_path = args.partition_results_dir / f"{stem}.json"
            if not csv_path.is_file() or not json_path.is_file():
                raise RuntimeError(f"missing compute partition evidence: {stem}")
            rows = read_csv(csv_path)
            manifest = load_json(json_path)
            if (
                manifest.get("status") != "complete_location_blind_roi_v4_measurement_v5_compute_partition"
                or manifest.get("protocol") != measurement_contract["protocol"]
                or manifest.get("inference_version") != inference["version"]
                or manifest.get("semantic_shard_index") != semantic
                or manifest.get("compute_partition_index") != partition
                or manifest.get("semantic_shard_count") != SEMANTIC_SHARDS
                or manifest.get("compute_partitions_per_semantic_shard") != PARTITIONS_PER_SEMANTIC_SHARD
                or manifest.get("terminal_records") != len(rows)
                or manifest.get("coordinates_opened") is not False
                or manifest.get("taxon_names_opened") is not False
                or manifest.get("result_sha256") != sha256(csv_path)
                or manifest.get("compute_contract_git_blob_sha") != git_blob_sha(args.compute_contract)
            ):
                raise RuntimeError(f"compute partition evidence changed: {stem}")
            partition_rows[(semantic, partition)] = rows
            partition_evidence[(semantic, partition)] = {
                "csv": sha256(csv_path),
                "json": sha256(json_path),
            }
            model_ids.add(str(manifest.get("model_id") or ""))
            weights.add(str(manifest.get("trained_weight_sha256") or ""))
            roi_contracts.add(str(manifest.get("roi_contract_sha256_lf_canonical_v1") or ""))

    if any(len(values) != 1 or not next(iter(values)) for values in (model_ids, weights, roi_contracts)):
        raise RuntimeError("compute partitions disagree on frozen estimator identity")
    coverage = validate_partition_coverage(expected_ids, partition_rows)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    semantic_hashes: dict[str, dict[str, str]] = {}
    for semantic in range(SEMANTIC_SHARDS):
        combined = [
            row
            for partition in range(PARTITIONS_PER_SEMANTIC_SHARD)
            for row in partition_rows[(semantic, partition)]
        ]
        combined.sort(key=lambda row: str(row["measurement_id"]))
        csv_path = args.output_dir / f"measurement_v5_shard_{semantic:04d}.csv"
        write_csv(csv_path, combined)
        manifest = {
            "status": "complete_location_blind_roi_v4_measurement_v5_shard",
            "protocol": measurement_contract["protocol"],
            "inference_version": inference["version"],
            "superseded_v3_ordered_inference_used": False,
            "shard_index": semantic,
            "shard_count": SEMANTIC_SHARDS,
            "frozen_shard_denominator": len(combined),
            "terminal_records": len(combined),
            "coordinates_opened": False,
            "taxon_names_opened": False,
            "result_sha256": sha256(csv_path),
            "model_id": next(iter(model_ids)),
            "trained_weight_sha256": next(iter(weights)),
            "roi_contract_sha256_lf_canonical_v1": next(iter(roi_contracts)),
            "compute_partition_protocol": compute_contract["protocol"],
            "compute_contract_git_blob_sha": git_blob_sha(args.compute_contract),
            "compute_partitions_reassembled": PARTITIONS_PER_SEMANTIC_SHARD,
            "partition_evidence_sha256": {
                f"p{partition:02d}": partition_evidence[(semantic, partition)]
                for partition in range(PARTITIONS_PER_SEMANTIC_SHARD)
            },
        }
        json_path = args.output_dir / f"measurement_v5_shard_{semantic:04d}.json"
        json_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        semantic_hashes[f"shard_{semantic:04d}"] = {"csv": sha256(csv_path), "json": sha256(json_path)}

    receipt = {
        **coverage,
        "protocol": compute_contract["protocol"],
        "candidate_image_pixels_persisted": False,
        "coordinates_opened_by_measurement": False,
        "measurement_manifest_sha256": sha256(args.measurement_manifest),
        "compute_contract_git_blob_sha": git_blob_sha(args.compute_contract),
        "semantic_shards": semantic_hashes,
    }
    receipt_path = args.output_dir / "compute_partition_reassembly_v1.json"
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

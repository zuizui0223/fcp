#!/usr/bin/env python3
"""Reassemble all 128 fresh measurement partitions and only then open metadata join."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re

import pandas as pd

from fcp_pipeline.photo_first_measurement_execution import (
    compute_partition,
    reassemble_complete_measurement,
    semantic_shard,
    validate_terminal_partition_results,
)


PATTERN = re.compile(r"partition_s(?P<shard>\d{2})_p(?P<part>\d{2})\.csv$")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", type=Path, required=True)
    parser.add_argument("--worker-manifest", type=Path, required=True)
    parser.add_argument("--metadata-join-key", type=Path, required=True)
    parser.add_argument("--firewall-manifest", type=Path, required=True)
    parser.add_argument("--execution-contract", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--output-manifest", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    execution = json.loads(args.execution_contract.read_text(encoding="utf-8"))
    firewall = json.loads(args.firewall_manifest.read_text(encoding="utf-8"))
    worker = pd.read_csv(args.worker_manifest, dtype=str).fillna("")
    metadata = pd.read_csv(args.metadata_join_key, dtype=str).fillna("")
    n_shards = int(execution["partitioning"]["semantic_shards"])
    n_parts = int(execution["partitioning"]["compute_partitions_per_semantic_shard"])
    expected_pairs = {(shard, part) for shard in range(n_shards) for part in range(n_parts)}

    frames_by_pair: dict[tuple[int, int], pd.DataFrame] = {}
    receipts_by_pair: dict[tuple[int, int], dict] = {}
    for path in sorted(args.results_dir.glob("partition_s*_p*.csv")):
        match = PATTERN.match(path.name)
        if not match:
            continue
        pair = (int(match.group("shard")), int(match.group("part")))
        if pair in frames_by_pair:
            raise ValueError(f"duplicate terminal partition file: {pair}")
        frame = pd.read_csv(path, dtype={"measurement_id": str}).fillna("")
        frames_by_pair[pair] = frame
        receipt_path = path.with_suffix(".json")
        if not receipt_path.is_file():
            raise ValueError(f"partition receipt is missing: {receipt_path.name}")
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        if (
            receipt.get("status") != "complete_random_photo_first_terminal_partition"
            or int(receipt.get("semantic_shard", -1)) != pair[0]
            or int(receipt.get("compute_partition", -1)) != pair[1]
        ):
            raise ValueError(f"partition receipt failed identity check: {pair}")
        receipts_by_pair[pair] = receipt

    if set(frames_by_pair) != expected_pairs or set(receipts_by_pair) != expected_pairs:
        missing = sorted(expected_pairs.difference(frames_by_pair))
        extra = sorted(set(frames_by_pair).difference(expected_pairs))
        raise ValueError(
            "not_evaluable_incomplete_measurement_partitions: "
            f"missing={missing[:10]} extra={extra[:10]}"
        )

    worker_ids = worker["measurement_id"].astype(str)
    for pair in sorted(expected_pairs):
        expected_ids = [
            measurement
            for measurement in worker_ids
            if semantic_shard(measurement, n_shards) == pair[0]
            and compute_partition(measurement, n_parts) == pair[1]
        ]
        validate_terminal_partition_results(frames_by_pair[pair], expected_ids)
        if int(receipts_by_pair[pair]["terminal_rows"]) != len(expected_ids):
            raise ValueError(f"partition receipt denominator changed: {pair}")

    ordered_frames = [frames_by_pair[pair] for pair in sorted(expected_pairs)]
    reassembled = reassemble_complete_measurement(
        ordered_frames,
        worker,
        metadata,
        expected_partition_receipts=n_shards * n_parts,
    )
    joined = reassembled.joined_photos.copy()
    order_columns = [
        column
        for column in ("cell_id", "observation_id", "photo_id", "measurement_id")
        if column in joined.columns
    ]
    joined = joined.sort_values(order_columns).reset_index(drop=True)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    args.output_manifest.parent.mkdir(parents=True, exist_ok=True)
    joined.to_csv(args.output_csv, index=False, lineterminator="\n")

    manifest = dict(reassembled.result_manifest)
    manifest.update(
        {
            "measurement_table_sha256": sha256(args.output_csv),
            "worker_manifest_sha256": sha256(args.worker_manifest),
            "metadata_join_key_sha256": sha256(args.metadata_join_key),
            "firewall_manifest_sha256": sha256(args.firewall_manifest),
            "execution_contract_sha256": sha256(args.execution_contract),
            "firewall_candidate_table_sha256": firewall.get("candidate_table_sha256"),
            "all_partition_receipts_verified": True,
            "coordinate_colour_join_opened": True,
            "h1_run": False,
            "h2_run": False,
        }
    )
    args.output_manifest.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

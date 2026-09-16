#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from fcp_pipeline.photo_first_measurement_execution import (
    ACQUISITION_FIELDS,
    WORKER_FIELDS,
    compute_partition,
    semantic_shard,
)
from fcp_pipeline.third_cohort_prospective_execution_gate import (
    AUTHORIZED_METADATA_SHA256,
    AUTHORIZED_SPECIES_SHA256,
    EXPECTED_ROWS,
    EXPECTED_SPECIES,
    EXPECTED_TERMINAL_PARTITIONS,
    build_start_record,
    measurement_batch,
    measurement_id,
    validate_authorized_denominator,
    validate_candidate_metadata,
    validate_stage_transition,
)

ROOT = Path(__file__).resolve().parents[2]
PROTOCOL = ROOT / "docs" / "POLYMORPHISM_H2_THIRD_COHORT_PROSPECTIVE_MEASUREMENT_PROTOCOL_20260917.md"
METADATA = ROOT / "data" / "frozen" / "polymorphism_h2_third_cohort_authorized_metadata_v1.csv.gz"
AUTHORIZED_SPECIES = ROOT / "results" / "polymorphism_h2_third_cohort_candidate_metadata_20260916" / "authorized_species.csv"
PARENT_RESULT = ROOT / "results" / "polymorphism_h2_third_cohort_candidate_metadata_20260916" / "result.json"
TARGET_PER_SPECIES = 100
BLIND_BATCHES = 2
SEMANTIC_SHARDS = 32
COMPUTE_PARTITIONS = 4


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--head-sha", required=True)
    p.add_argument("--branch", required=True)
    p.add_argument("--known-output-paths-absent", action="store_true")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    for path in (PROTOCOL, METADATA, AUTHORIZED_SPECIES, PARENT_RESULT):
        if not path.exists():
            raise RuntimeError(f"missing frozen input: {path}")
    if sha256_file(METADATA) != AUTHORIZED_METADATA_SHA256:
        raise RuntimeError("third-cohort authorized metadata SHA drifted")
    if sha256_file(AUTHORIZED_SPECIES) != AUTHORIZED_SPECIES_SHA256:
        raise RuntimeError("third-cohort authorized species SHA drifted")

    parent = json.loads(PARENT_RESULT.read_text(encoding="utf-8"))
    validate_candidate_metadata(parent)
    species = pd.read_csv(AUTHORIZED_SPECIES)
    metadata = pd.read_csv(METADATA, low_memory=False)
    metadata = validate_authorized_denominator(
        metadata,
        species,
        expected_species=EXPECTED_SPECIES,
        target_per_species=TARGET_PER_SPECIES,
    )

    start = build_start_record(
        branch=args.branch,
        head_sha=args.head_sha,
        known_output_paths_absent=bool(args.known_output_paths_absent),
        working_tree_inputs_verified=True,
    )
    if start["opening_authorized"] is not True:
        raise RuntimeError("preopening start record did not authorize opening")

    metadata = metadata.reset_index(drop=True)
    metadata.insert(0, "measurement_id", [measurement_id(x) for x in metadata["photo_id"]])
    metadata.insert(1, "measurement_batch", [measurement_batch(x) for x in metadata["measurement_id"].astype(str)])
    if metadata["measurement_id"].nunique() != EXPECTED_ROWS:
        raise RuntimeError("third-cohort blinded measurement IDs are not unique")

    worker_all = pd.DataFrame({
        "measurement_id": metadata["measurement_id"].astype(str),
        "image_filename": metadata["measurement_id"].astype(str) + ".jpg",
        "photo_license": metadata["photo_license"].astype(str),
    })
    acquisition_all = pd.DataFrame({
        "measurement_id": worker_all["measurement_id"],
        "image_filename": worker_all["image_filename"],
        "photo_url_large": metadata["photo_url_large"].astype(str),
        "photo_license": worker_all["photo_license"],
    })
    if tuple(worker_all.columns) != WORKER_FIELDS:
        raise RuntimeError("third-cohort worker interface drifted")
    if tuple(acquisition_all.columns) != ACQUISITION_FIELDS:
        raise RuntimeError("third-cohort acquisition interface drifted")

    out = args.output_dir
    (out / "sealed_keys").mkdir(parents=True, exist_ok=True)
    metadata_join = metadata.drop(columns=["photo_url_large"]).copy()
    metadata_join.to_csv(out / "sealed_keys" / "metadata_join_key.csv", index=False, lineterminator="\n")

    assignments: list[pd.DataFrame] = []
    batch_counts: dict[str, int] = {}
    nonempty = 0
    max_partition_rows = 0
    for batch in range(BLIND_BATCHES):
        keep = metadata["measurement_batch"].astype(int).eq(batch).to_numpy()
        worker = worker_all.loc[keep].reset_index(drop=True)
        acquisition = acquisition_all.loc[keep].reset_index(drop=True)
        batch_dir = out / f"batch_{batch}"
        (batch_dir / "worker_packet").mkdir(parents=True, exist_ok=True)
        (batch_dir / "sealed_keys").mkdir(parents=True, exist_ok=True)
        worker.to_csv(batch_dir / "worker_packet" / "measurement_manifest.csv", index=False, lineterminator="\n")
        acquisition.to_csv(batch_dir / "sealed_keys" / "acquisition_key.csv", index=False, lineterminator="\n")
        batch_counts[str(batch)] = int(len(worker))
        local = pd.DataFrame({
            "measurement_id": worker["measurement_id"].astype(str),
            "measurement_batch": batch,
            "semantic_shard": [semantic_shard(x, SEMANTIC_SHARDS) for x in worker["measurement_id"].astype(str)],
            "compute_partition": [compute_partition(x, COMPUTE_PARTITIONS) for x in worker["measurement_id"].astype(str)],
        })
        assignments.append(local)
        counts = local.groupby(["semantic_shard", "compute_partition"], observed=True).size()
        nonempty += int(len(counts))
        if len(counts):
            max_partition_rows = max(max_partition_rows, int(counts.max()))

    assignment = pd.concat(assignments, ignore_index=True).sort_values("measurement_id", kind="mergesort")
    assignment.to_csv(out / "partition_assignments.csv", index=False, lineterminator="\n")
    if sum(batch_counts.values()) != EXPECTED_ROWS:
        raise RuntimeError("third-cohort blind batching changed denominator")

    firewall_stage = {
        "stage": "FIREWALL_FROZEN",
        "species": EXPECTED_SPECIES,
        "rows": EXPECTED_ROWS,
        "replacement_rows": 0,
        "replacement_species": 0,
        "measurement_ids": EXPECTED_ROWS,
        "terminal_partitions": EXPECTED_TERMINAL_PARTITIONS,
        "candidate_pixels_opened": False,
    }
    validate_stage_transition(start, firewall_stage)

    result = {
        "schema": "third_cohort_measurement_firewall_v1",
        "analysis": "polymorphism_h2_third_cohort_measurement_firewall",
        "status": "frozen_before_third_cohort_pixels",
        "stage": "FIREWALL_FROZEN",
        "protocol": str(PROTOCOL.relative_to(ROOT)),
        "parent_metadata_result": str(PARENT_RESULT.relative_to(ROOT)),
        "lineage": {
            "authorized_metadata_sha256": AUTHORIZED_METADATA_SHA256,
            "authorized_species_sha256": AUTHORIZED_SPECIES_SHA256,
            "protocol_sha256": sha256_file(PROTOCOL),
        },
        "frozen_species": EXPECTED_SPECIES,
        "frozen_rows": EXPECTED_ROWS,
        "target_rows_per_species": TARGET_PER_SPECIES,
        "blind_batches": BLIND_BATCHES,
        "rows_per_blind_batch": batch_counts,
        "semantic_shards_per_batch": SEMANTIC_SHARDS,
        "compute_partitions_per_semantic_shard": COMPUTE_PARTITIONS,
        "total_terminal_partitions": EXPECTED_TERMINAL_PARTITIONS,
        "nonempty_terminal_partitions": nonempty,
        "maximum_rows_in_terminal_partition": max_partition_rows,
        "worker_fields": list(worker_all.columns),
        "acquisition_fields": list(acquisition_all.columns),
        "worker_contains_species": False,
        "worker_contains_coordinates": False,
        "worker_contains_source_url": False,
        "coordinate_colour_join_opened": False,
        "candidate_pixels_opened": False,
        "H2_W_opened": False,
        "stage_receipts": {
            "start": start,
            "firewall": firewall_stage,
        },
    }
    (out / "measurement_firewall_manifest.json").write_text(
        json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

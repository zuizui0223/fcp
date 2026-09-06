#!/usr/bin/env python3
"""Run one deterministic subset of the frozen 999 RGFCA G1 nulls."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from fcp_pipeline.global_repeated_atlas import schedule_audit
from fcp_pipeline.global_rgfca_engine import COLOUR_COLUMNS, run_g1_shard

ROOT = Path(__file__).resolve().parents[2]
EXECUTION = ROOT / "docs/supporting/global_monte_carlo_inference_execution_contract_v1.json"
MEASUREMENT = ROOT / "docs/supporting/global_monte_carlo_measurement_result_v1.json"
MEASURED = ROOT / "data/derived/global_monte_carlo_measured_photos_v1.csv"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def array_sha256(*arrays: np.ndarray) -> str:
    h = hashlib.sha256()
    for array in arrays:
        a = np.ascontiguousarray(array)
        h.update(str(a.dtype).encode("ascii"))
        h.update(np.asarray(a.shape, dtype=np.int64).tobytes())
        h.update(a.tobytes())
    return h.hexdigest()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--shard-index", type=int, required=True)
    p.add_argument("--shard-count", type=int, default=9)
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    execution = json.loads(EXECUTION.read_text(encoding="utf-8"))
    measurement = json.loads(MEASUREMENT.read_text(encoding="utf-8"))
    if execution.get("status") != "frozen_before_capacity_outcome_before_candidate_pixels_and_before_global_colour_outcome":
        raise RuntimeError("global inference execution contract is not frozen pre-outcome")
    gate = execution["input_gate"]
    if measurement.get("status") != gate["measurement_status_required"]:
        raise RuntimeError("global measurement status does not authorize inference")
    post = measurement.get("postmeasurement_gate", {})
    if post.get("pass") is not True or post.get("decision") != gate["measurement_postgate_decision_required"]:
        raise RuntimeError("global measurement postgate did not authorize inference")
    if measurement.get("g1_g3_inference_run") is not False:
        raise RuntimeError("global measurement manifest says inference already ran")
    if sha256_file(MEASURED) != str(measurement.get("lineage", {}).get("measured_table_sha256") or ""):
        raise RuntimeError("measured table lineage differs from frozen measurement result")

    frame = pd.read_csv(MEASURED)
    required = {"photo_id", "species", "latitude", "longitude", "global_classifiable", *COLOUR_COLUMNS}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise RuntimeError(f"measured table lacks G1 inputs: {missing}")
    classifiable_flag = frame["global_classifiable"].astype(str).str.casefold().isin({"true", "1"})
    pool = frame.loc[classifiable_flag, ["photo_id", "species", "latitude", "longitude", *COLOUR_COLUMNS]].copy()
    minimum_pool = int(gate["minimum_classifiable_photos_per_species"])
    counts = pool.groupby("species", observed=True).size()
    eligible = counts[counts >= minimum_pool].index.astype(str)
    pool = pool.loc[pool["species"].astype(str).isin(set(eligible))].copy().reset_index(drop=True)
    eligible_species = int(pool["species"].nunique())
    if eligible_species != int(post["evaluable_species"]):
        raise RuntimeError("inference pool species count differs from measurement postgate")
    if eligible_species < int(gate["minimum_inferential_species"]):
        raise RuntimeError("inference pool is below frozen species gate")
    if pool["photo_id"].duplicated().any():
        raise RuntimeError("inference pool photo IDs are not unique")

    null_execution = execution["null_execution"]
    shard_count = int(args.shard_count)
    shard_index = int(args.shard_index)
    if shard_count != int(null_execution["deterministic_shards"]):
        raise RuntimeError("G1 null shard count differs from frozen execution")
    if shard_index < 0 or shard_index >= shard_count:
        raise RuntimeError("invalid G1 null shard index")
    null_indices = np.asarray([i for i in range(999) if i % shard_count == shard_index], dtype=np.int64)
    if len(null_indices) == 0:
        raise RuntimeError("G1 null shard unexpectedly empty")

    outer = execution["outer_schedule"]
    g1 = execution["g1_primary"]
    result = run_g1_shard(
        pool,
        null_indices=null_indices,
        n_outer=int(outer["observed_resamples"]),
        species_per_outer=int(outer["species_per_resample"]),
        photos_per_species=int(outer["photos_per_species"]),
        minimum_pool_photos_per_species=minimum_pool,
        k=int(g1["k"]),
        n_lon=36,
        n_sinlat=18,
        kernel_km=float(g1["kernel_bandwidth_km"]),
        cutoff_multiplier=float(g1["kernel_cutoff_multiplier"]),
        minimum_distinct_species=int(g1["minimum_distinct_species_support"]),
        species_seed=int(outer["species_seed"]),
        photo_master_seed=int(outer["photo_master_seed"]),
        null_master_seed=int(g1["null_master_seed"]),
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    npz_path = args.output_dir / f"global_rgfca_g1_shard_{shard_index:02d}.npz"
    np.savez_compressed(
        npz_path,
        observed_outer_fields=result.observed_outer_fields,
        observed_outer_opportunities=result.observed_outer_opportunities,
        observed_consensus=result.observed_consensus,
        observed_aggregate_opportunity=result.observed_aggregate_opportunity,
        observed_concentration=np.asarray([result.observed_concentration], dtype=float),
        null_indices=result.null_indices,
        null_consensus_fields=result.null_consensus_fields,
        null_concentrations=result.null_concentrations,
    )
    observed_digest = array_sha256(
        result.observed_outer_fields,
        result.observed_outer_opportunities,
        result.observed_consensus,
        result.observed_aggregate_opportunity,
        np.asarray([result.observed_concentration], dtype=np.float64),
    )
    manifest = {
        "protocol": execution["protocol"],
        "status": "complete_global_rgfca_g1_null_shard",
        "shard_index": shard_index,
        "shard_count": shard_count,
        "null_indices": [int(x) for x in result.null_indices],
        "null_count": int(len(result.null_indices)),
        "eligible_species": eligible_species,
        "classifiable_pool_rows": int(len(pool)),
        "observed_concentration": float(result.observed_concentration),
        "observed_digest_sha256": observed_digest,
        "null_unit": result.null_unit,
        "schedule_audit": schedule_audit(result.schedule),
        "lineage": {
            "execution_contract_sha256": sha256_file(EXECUTION),
            "measurement_manifest_sha256": sha256_file(MEASUREMENT),
            "measured_table_sha256": sha256_file(MEASURED),
            "npz_sha256": sha256_file(npz_path),
        },
    }
    manifest_path = args.output_dir / f"global_rgfca_g1_shard_{shard_index:02d}.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

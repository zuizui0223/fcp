#!/usr/bin/env python3
"""Execute only the G1 robustness diagnostics frozen before the RGFCA colour outcome.

This script cannot alter or rescue the primary G1 decision. It implements:
1) the two prespecified spatial-support sensitivities;
2) the prespecified minimum-distinct-species=10 sensitivity; and
3) the prespecified odd/even half-null excess diagnostics.

All modes reuse the exact measured pool, balanced outer schedule, colour estimator,
null indices 0..998, and frozen seeds from the primary G1 program.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from fcp_pipeline.global_barrier_field import equal_area_grid_centers
from fcp_pipeline.global_repeated_atlas import build_repeated_atlas_schedule, consensus_field
from fcp_pipeline.global_rgfca_engine import (
    COLOUR_COLUMNS,
    _field_from_rank,
    _field_matrix_from_rank,
    _rank_edges,
    _rank_edges_matrix,
    _raw_jsd_from_source_rows,
    _raw_jsd_matrix_from_source_rows,
    _weighted_concentration,
    build_pairwise_jsd_cache,
    canonical_colour_pool,
    null_source_row_matrix,
    prepare_sparse_outer_geometry,
    run_g1_shard,
)

ROOT = Path(__file__).resolve().parents[2]
EXECUTION = ROOT / "docs/supporting/global_monte_carlo_inference_execution_contract_v1.json"
BARRIER = ROOT / "docs/supporting/global_monte_carlo_barrier_field_contract_v1.json"
MEASUREMENT = ROOT / "docs/supporting/global_monte_carlo_measurement_result_v1.json"
PRIMARY = ROOT / "docs/supporting/global_rgfca_g1_result_v1.json"
MEASURED = ROOT / "data/derived/global_monte_carlo_measured_photos_v1.csv"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_inputs() -> tuple[pd.DataFrame, dict, dict, dict, dict]:
    execution = json.loads(EXECUTION.read_text())
    barrier = json.loads(BARRIER.read_text())
    measurement = json.loads(MEASUREMENT.read_text())
    primary = json.loads(PRIMARY.read_text())

    if barrier.get("status") != "frozen_before_new_global_species_discovery_outcome_and_before_new_colour_pixels":
        raise RuntimeError("barrier robustness contract status drift")
    if primary.get("status") != "complete_global_rgfca_g1_primary_inference":
        raise RuntimeError("primary G1 result is not frozen")
    if int(primary.get("null_permutations", 0)) != 999:
        raise RuntimeError("primary G1 null count drift")
    if primary.get("persistent_zone_extracted") is not False or primary.get("g4_overlay_run") is not False:
        raise RuntimeError("null primary G1 must not have produced zones/overlays")

    gate = execution["input_gate"]
    post = measurement.get("postmeasurement_gate", {})
    if post.get("pass") is not True:
        raise RuntimeError("measurement postgate no longer passes")
    if sha256_file(MEASURED) != str(measurement.get("lineage", {}).get("measured_table_sha256") or ""):
        raise RuntimeError("measured table lineage drift")

    frame = pd.read_csv(MEASURED)
    required = {"photo_id", "species", "latitude", "longitude", "global_classifiable", *COLOUR_COLUMNS}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise RuntimeError(f"measured table lacks inputs: {missing}")
    classifiable = frame["global_classifiable"].astype(str).str.casefold().isin({"true", "1"})
    pool = frame.loc[classifiable, ["photo_id", "species", "latitude", "longitude", *COLOUR_COLUMNS]].copy()
    minimum_pool = int(gate["minimum_classifiable_photos_per_species"])
    counts = pool.groupby("species", observed=True).size()
    eligible = set(counts[counts >= minimum_pool].index.astype(str))
    pool = pool.loc[pool["species"].astype(str).isin(eligible)].copy().reset_index(drop=True)
    if int(pool["species"].nunique()) != int(post["evaluable_species"]):
        raise RuntimeError("eligible species count drift")
    if len(pool) != int(primary["classifiable_pool_rows"]):
        raise RuntimeError("classifiable pool row count drift")
    return pool, execution, barrier, measurement, primary


def parse_grid(text: str) -> tuple[int, int]:
    left, right = str(text).lower().split("x", 1)
    return int(left), int(right)


def common_kwargs(execution: dict) -> dict[str, object]:
    outer = execution["outer_schedule"]
    g1 = execution["g1_primary"]
    return {
        "n_outer": int(outer["observed_resamples"]),
        "species_per_outer": int(outer["species_per_resample"]),
        "photos_per_species": int(outer["photos_per_species"]),
        "minimum_pool_photos_per_species": int(execution["input_gate"]["minimum_classifiable_photos_per_species"]),
        "k": int(g1["k"]),
        "cutoff_multiplier": float(g1["kernel_cutoff_multiplier"]),
        "species_seed": int(outer["species_seed"]),
        "photo_master_seed": int(outer["photo_master_seed"]),
        "null_master_seed": int(g1["null_master_seed"]),
    }


def run_spatial(name: str, output: Path) -> None:
    pool, execution, barrier, _, primary = load_inputs()
    primary_field = barrier["primary_field"]
    sensitivities = barrier["predeclared_spatial_support_sensitivities"]
    if len(sensitivities) != 2:
        raise RuntimeError("expected exactly two prespecified spatial sensitivities")

    if name == "coarse_24x12_1000km":
        spec = sensitivities[0]
        if spec != {
            "grid": "24x12",
            "kernel_bandwidth_km": 1000.0,
            "compact_support_cutoff_multiplier": 3.0,
            "minimum_distinct_species_support_per_cell": 5,
        }:
            raise RuntimeError("coarse sensitivity contract drift")
        n_lon, n_sinlat = parse_grid(spec["grid"])
        kernel = float(spec["kernel_bandwidth_km"])
        minimum_support = int(spec["minimum_distinct_species_support_per_cell"])
    elif name == "fine_72x36_250km":
        spec = sensitivities[1]
        if spec != {
            "grid": "72x36",
            "kernel_bandwidth_km": 250.0,
            "compact_support_cutoff_multiplier": 3.0,
            "minimum_distinct_species_support_per_cell": 5,
        }:
            raise RuntimeError("fine sensitivity contract drift")
        n_lon, n_sinlat = parse_grid(spec["grid"])
        kernel = float(spec["kernel_bandwidth_km"])
        minimum_support = int(spec["minimum_distinct_species_support_per_cell"])
    elif name == "support10_36x18_500km":
        threshold = barrier["support_threshold_sensitivity"]
        if int(threshold["primary_minimum_distinct_species"]) != 5 or int(threshold["secondary_minimum_distinct_species"]) != 10:
            raise RuntimeError("support threshold sensitivity contract drift")
        if threshold.get("secondary_cannot_replace_primary") is not True:
            raise RuntimeError("support10 is not allowed to replace primary")
        n_lon, n_sinlat = parse_grid(primary_field["grid"].split(" equal-area")[0])
        kernel = float(primary_field["kernel_bandwidth_km"])
        minimum_support = 10
    else:
        raise ValueError(f"unknown spatial robustness configuration: {name}")

    kwargs = common_kwargs(execution)
    result = run_g1_shard(
        pool,
        null_indices=np.arange(999, dtype=np.int64),
        n_lon=n_lon,
        n_sinlat=n_sinlat,
        kernel_km=kernel,
        minimum_distinct_species=minimum_support,
        **kwargs,
    )
    null_mean = float(np.mean(result.null_concentrations))
    observed = float(result.observed_concentration)
    p_upper = float((1 + np.count_nonzero(result.null_concentrations >= observed)) / 1000)
    payload = {
        "protocol": "global-rgfca-prespecified-robustness-execution-v1",
        "status": "complete_prespecified_g1_spatial_robustness",
        "configuration": name,
        "role": "robustness_only_cannot_rescue_primary_g1",
        "primary_g1_supported": bool(primary["g1_supported"]),
        "primary_g1_p_upper": float(primary["p_upper"]),
        "grid": f"{n_lon}x{n_sinlat}",
        "kernel_bandwidth_km": kernel,
        "minimum_distinct_species_support": minimum_support,
        "observed_concentration": observed,
        "null_mean_concentration": null_mean,
        "g1_excess": observed - null_mean,
        "g1_excess_ratio": observed / null_mean if null_mean > 0 else float("nan"),
        "positive_excess": bool(observed > null_mean),
        "descriptive_p_upper_not_primary_gate": p_upper,
        "null_permutations": 999,
        "null_indices_exact_0_998": bool(np.array_equal(result.null_indices, np.arange(999))),
        "scientific_rules_changed": False,
        "lineage": {
            "execution_contract_sha256": sha256_file(EXECUTION),
            "barrier_contract_sha256": sha256_file(BARRIER),
            "measurement_manifest_sha256": sha256_file(MEASUREMENT),
            "measured_table_sha256": sha256_file(MEASURED),
            "primary_g1_result_sha256": sha256_file(PRIMARY),
        },
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2), flush=True)


def _half_payload(observed, null_concentrations, opportunity, half_name: str) -> dict[str, object]:
    null_mean = float(np.mean(null_concentrations))
    obs = float(observed.concentration)
    return {
        "half": half_name,
        "observed_concentration": obs,
        "null_mean_concentration": null_mean,
        "g1_excess": obs - null_mean,
        "g1_excess_ratio": obs / null_mean if null_mean > 0 else float("nan"),
        "positive_excess": bool(obs > null_mean),
        "descriptive_p_upper_not_primary_gate": float((1 + np.count_nonzero(null_concentrations >= obs)) / 1000),
        "evaluable_cells": int(np.count_nonzero(opportunity > 0)),
    }


def run_halves(output: Path) -> None:
    pool_frame, execution, barrier, _, primary = load_inputs()
    kwargs = common_kwargs(execution)
    n_outer = int(kwargs["n_outer"])
    if n_outer != 200:
        raise RuntimeError("odd/even robustness expects frozen 200 outer realizations")

    pool = canonical_colour_pool(pool_frame)
    schedule = build_repeated_atlas_schedule(
        pool.photo_ids,
        pool.species,
        n_outer=n_outer,
        species_per_outer=int(kwargs["species_per_outer"]),
        photos_per_species=int(kwargs["photos_per_species"]),
        minimum_pool_photos_per_species=int(kwargs["minimum_pool_photos_per_species"]),
        species_seed=int(kwargs["species_seed"]),
        photo_master_seed=int(kwargs["photo_master_seed"]),
    )
    cache = build_pairwise_jsd_cache(pool)
    null_indices = np.arange(999, dtype=np.int64)
    source_null = null_source_row_matrix(pool, null_indices, master_seed=int(kwargs["null_master_seed"]))
    identity = np.arange(len(pool.photo_ids), dtype=np.int64)

    primary_field = barrier["primary_field"]
    n_lon, n_sinlat = parse_grid(primary_field["grid"].split(" equal-area")[0])
    grid = equal_area_grid_centers(n_lon, n_sinlat)
    kernel = float(primary_field["kernel_bandwidth_km"])
    minimum_support = int(primary_field["minimum_distinct_species_support_per_cell"])
    cutoff = float(primary_field["compact_support_cutoff_multiplier"])

    observed_fields = np.full((n_outer, grid.n_cells), np.nan, dtype=float)
    observed_opportunities = np.zeros((n_outer, grid.n_cells), dtype=float)
    odd_null_num = np.zeros((999, grid.n_cells), dtype=float)
    even_null_num = np.zeros((999, grid.n_cells), dtype=float)
    odd_opp = np.zeros(grid.n_cells, dtype=float)
    even_opp = np.zeros(grid.n_cells, dtype=float)

    for outer in range(n_outer):
        geometry = prepare_sparse_outer_geometry(
            pool,
            schedule.outer_photo_ids[outer],
            schedule.outer_species[outer],
            grid=grid,
            k=int(kwargs["k"]),
            kernel_km=kernel,
            cutoff_multiplier=cutoff,
            minimum_distinct_species=minimum_support,
        )
        raw_obs = _raw_jsd_from_source_rows(cache, geometry.edge_nodes, identity)
        rank_obs = _rank_edges(raw_obs, geometry.edge_species_slices)
        field_obs, _ = _field_from_rank(geometry, rank_obs)
        observed_fields[outer] = field_obs
        observed_opportunities[outer, geometry.evaluable] = geometry.opportunity[geometry.evaluable]

        raw_null = _raw_jsd_matrix_from_source_rows(cache, geometry.edge_nodes, source_null)
        rank_null = _rank_edges_matrix(raw_null, geometry.edge_species_slices)
        _, numerator_null = _field_matrix_from_rank(geometry, rank_null)
        if outer % 2 == 0:
            odd_null_num[:, geometry.evaluable] += numerator_null[:, geometry.evaluable]
            odd_opp[geometry.evaluable] += geometry.opportunity[geometry.evaluable]
        else:
            even_null_num[:, geometry.evaluable] += numerator_null[:, geometry.evaluable]
            even_opp[geometry.evaluable] += geometry.opportunity[geometry.evaluable]

    observed_odd = consensus_field(observed_fields[0::2], observed_opportunities[0::2])
    observed_even = consensus_field(observed_fields[1::2], observed_opportunities[1::2])
    if not np.allclose(observed_odd.aggregate_opportunity, odd_opp):
        raise RuntimeError("odd opportunity aggregation mismatch")
    if not np.allclose(observed_even.aggregate_opportunity, even_opp):
        raise RuntimeError("even opportunity aggregation mismatch")

    def null_concentrations(numerator: np.ndarray, opportunity: np.ndarray) -> np.ndarray:
        fields = np.full_like(numerator, np.nan, dtype=float)
        keep = opportunity > 0
        fields[:, keep] = numerator[:, keep] / opportunity[keep]
        return np.asarray([_weighted_concentration(row, opportunity) for row in fields], dtype=float)

    odd_null = null_concentrations(odd_null_num, odd_opp)
    even_null = null_concentrations(even_null_num, even_opp)
    common = np.isfinite(observed_odd.field) & np.isfinite(observed_even.field)
    odd_even_r = float(np.corrcoef(observed_odd.field[common], observed_even.field[common])[0, 1]) if int(common.sum()) >= 3 else float("nan")
    odd = _half_payload(observed_odd, odd_null, odd_opp, "odd_realizations_1_3_to_199")
    even = _half_payload(observed_even, even_null, even_opp, "even_realizations_2_4_to_200")
    threshold = barrier["stability_claim_thresholds"]
    payload = {
        "protocol": "global-rgfca-prespecified-robustness-execution-v1",
        "status": "complete_prespecified_g1_odd_even_half_null_robustness",
        "role": "robustness_only_cannot_rescue_primary_g1",
        "primary_g1_supported": bool(primary["g1_supported"]),
        "primary_g1_p_upper": float(primary["p_upper"]),
        "null_permutations_per_half": 999,
        "null_indices_exact_0_998": True,
        "odd": odd,
        "even": even,
        "odd_even_consensus_field_pearson_r": odd_even_r,
        "predeclared_r_minimum": float(threshold["odd_even_consensus_field_pearson_r_minimum"]),
        "r_threshold_pass": bool(odd_even_r >= float(threshold["odd_even_consensus_field_pearson_r_minimum"])),
        "both_half_excess_positive": bool(odd["positive_excess"] and even["positive_excess"]),
        "predeclared_half_excess_gate_pass": bool(odd["positive_excess"] and even["positive_excess"]),
        "scientific_rules_changed": False,
        "lineage": {
            "execution_contract_sha256": sha256_file(EXECUTION),
            "barrier_contract_sha256": sha256_file(BARRIER),
            "measurement_manifest_sha256": sha256_file(MEASUREMENT),
            "measured_table_sha256": sha256_file(MEASURED),
            "primary_g1_result_sha256": sha256_file(PRIMARY),
        },
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2), flush=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("spatial", "halves"), required=True)
    parser.add_argument("--configuration", default="")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.mode == "spatial":
        if not args.configuration:
            raise SystemExit("--configuration is required for spatial mode")
        run_spatial(args.configuration, args.output)
    else:
        if args.configuration:
            raise SystemExit("--configuration is not used for halves mode")
        run_halves(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

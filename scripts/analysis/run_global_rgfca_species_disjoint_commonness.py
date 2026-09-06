#!/usr/bin/env python3
"""Run the prospectively frozen species-disjoint RGFCA commonness test.

This test asks whether a discontinuity field learned from one set of species predicts
where discontinuities occur in entirely held-out species. It allows any number of
spatially disconnected recurrent zones and cannot reclassify the frozen primary G1.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from fcp_pipeline.global_barrier_field import equal_area_grid_centers
from fcp_pipeline.global_repeated_atlas import build_repeated_atlas_schedule
from fcp_pipeline.global_rgfca_engine import (
    COLOUR_COLUMNS,
    _field_from_rank,
    _field_matrix_from_rank,
    _rank_edges,
    _rank_edges_matrix,
    _raw_jsd_from_source_rows,
    _raw_jsd_matrix_from_source_rows,
    build_pairwise_jsd_cache,
    canonical_colour_pool,
    null_source_row_matrix,
    prepare_sparse_outer_geometry,
)

ROOT = Path(__file__).resolve().parents[2]
EXECUTION = ROOT / "docs/supporting/global_monte_carlo_inference_execution_contract_v1.json"
MEASUREMENT = ROOT / "docs/supporting/global_monte_carlo_measurement_result_v1.json"
PRIMARY = ROOT / "docs/supporting/global_rgfca_g1_result_v1.json"
COMMON = ROOT / "docs/supporting/global_rgfca_common_boundary_claim_contract_v1.json"
MAPPING = ROOT / "docs/supporting/global_rgfca_species_disjoint_execution_v1.json"
MEASURED = ROOT / "data/derived/global_monte_carlo_measured_photos_v1.csv"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def weighted_pearson(x: np.ndarray, y: np.ndarray, w: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    w = np.asarray(w, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y) & np.isfinite(w) & (w > 0)
    if int(keep.sum()) < 3:
        return float("nan")
    x = x[keep]
    y = y[keep]
    w = w[keep]
    sw = float(w.sum())
    if not np.isfinite(sw) or sw <= 0:
        return float("nan")
    mx = float(np.dot(w, x) / sw)
    my = float(np.dot(w, y) / sw)
    dx = x - mx
    dy = y - my
    vx = float(np.dot(w, dx * dx) / sw)
    vy = float(np.dot(w, dy * dy) / sw)
    if vx <= 0 or vy <= 0 or not np.isfinite(vx + vy):
        return float("nan")
    cov = float(np.dot(w, dx * dy) / sw)
    return float(cov / np.sqrt(vx * vy))


def deterministic_folds(species_labels: tuple[str, ...], *, seed: int, n_folds: int) -> np.ndarray:
    labels = np.asarray(species_labels, dtype=object)
    if len(labels) == 0 or len(np.unique(labels)) != len(labels):
        raise ValueError("species labels must be unique and non-empty")
    canonical_order = np.argsort(labels.astype(str), kind="stable")
    canonical_labels = labels[canonical_order]
    if not np.array_equal(canonical_labels.astype(str), np.sort(labels.astype(str))):
        raise RuntimeError("canonical species ordering failure")
    rng = np.random.default_rng(int(seed))
    permuted_positions = rng.permutation(len(canonical_labels))
    fold_for_canonical = np.empty(len(canonical_labels), dtype=np.int16)
    for rank, canonical_position in enumerate(permuted_positions):
        fold_for_canonical[int(canonical_position)] = int(rank % int(n_folds))
    out = np.empty(len(labels), dtype=np.int16)
    out[canonical_order] = fold_for_canonical
    return out


def load_inputs() -> tuple[pd.DataFrame, dict, dict, dict, dict]:
    execution = json.loads(EXECUTION.read_text())
    measurement = json.loads(MEASUREMENT.read_text())
    primary = json.loads(PRIMARY.read_text())
    common = json.loads(COMMON.read_text())
    mapping = json.loads(MAPPING.read_text())

    if primary.get("status") != "complete_global_rgfca_g1_primary_inference":
        raise RuntimeError("primary G1 result missing")
    if mapping.get("status") != "technical_execution_mapping_frozen_before_species_disjoint_commonness_result":
        raise RuntimeError("species-disjoint execution mapping is not frozen")
    if common.get("status") != "frozen_before_species_disjoint_commonness_result":
        raise RuntimeError("common-boundary claim contract status drift")
    if common["terminology"].get("single_universal_boundary_required") is not False:
        raise RuntimeError("single-boundary prohibition drift")
    if common["prospective_species_disjoint_test"]["claim_rule"].get("cannot_reclassify_primary_G1") is not True:
        raise RuntimeError("primary G1 rescue guard drift")
    if mapping["decision"].get("primary_g1_reclassification_allowed") is not False:
        raise RuntimeError("species-disjoint mapping permits primary rescue")
    if int(primary.get("null_permutations", 0)) != 999:
        raise RuntimeError("primary null count drift")
    if sha256_file(MEASURED) != str(measurement.get("lineage", {}).get("measured_table_sha256") or ""):
        raise RuntimeError("measured table hash drift")

    frame = pd.read_csv(MEASURED)
    required = {"photo_id", "species", "latitude", "longitude", "global_classifiable", *COLOUR_COLUMNS}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise RuntimeError(f"measured table lacks inputs: {missing}")
    classifiable = frame["global_classifiable"].astype(str).str.casefold().isin({"true", "1"})
    pool = frame.loc[classifiable, ["photo_id", "species", "latitude", "longitude", *COLOUR_COLUMNS]].copy()
    minimum_pool = int(execution["input_gate"]["minimum_classifiable_photos_per_species"])
    counts = pool.groupby("species", observed=True).size()
    eligible = set(counts[counts >= minimum_pool].index.astype(str))
    pool = pool.loc[pool["species"].astype(str).isin(eligible)].copy().reset_index(drop=True)
    expected_species = int(mapping["species_partition"]["eligible_species_count_expected"])
    if int(pool["species"].nunique()) != expected_species:
        raise RuntimeError("eligible species count drift")
    if len(pool) != int(primary["classifiable_pool_rows"]):
        raise RuntimeError("classifiable pool row count drift")
    return pool, execution, primary, common, mapping


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-result", type=Path, required=True)
    parser.add_argument("--output-folds", type=Path, required=True)
    parser.add_argument("--output-null", type=Path, required=True)
    args = parser.parse_args()

    frame, execution, primary, common, mapping = load_inputs()
    pool = canonical_colour_pool(frame)
    part = mapping["species_partition"]
    n_folds = int(part["folds"])
    fold_by_species = deterministic_folds(
        pool.species_labels,
        seed=int(part["seed"]),
        n_folds=n_folds,
    )
    fold_sizes = np.bincount(fold_by_species, minlength=n_folds).astype(int)
    if sorted(fold_sizes.tolist(), reverse=True) != sorted(part["balanced_fold_sizes_expected"], reverse=True):
        raise RuntimeError(f"fold-size drift: {fold_sizes.tolist()}")
    if set(np.unique(fold_by_species).tolist()) != set(range(n_folds)):
        raise RuntimeError("not all folds represented")

    fold_frame = pd.DataFrame({
        "species": list(pool.species_labels),
        "fold": fold_by_species.astype(int),
    }).sort_values(["fold", "species"], kind="stable")
    args.output_folds.parent.mkdir(parents=True, exist_ok=True)
    fold_frame.to_csv(args.output_folds, index=False)

    outer = execution["outer_schedule"]
    g1 = execution["g1_primary"]
    schedule = build_repeated_atlas_schedule(
        pool.photo_ids,
        pool.species,
        n_outer=int(outer["observed_resamples"]),
        species_per_outer=int(outer["species_per_resample"]),
        photos_per_species=int(outer["photos_per_species"]),
        minimum_pool_photos_per_species=int(execution["input_gate"]["minimum_classifiable_photos_per_species"]),
        species_seed=int(outer["species_seed"]),
        photo_master_seed=int(outer["photo_master_seed"]),
    )
    if schedule.n_outer != 200:
        raise RuntimeError("species-disjoint test requires frozen 200 outer schedule")

    grid = equal_area_grid_centers(36, 18)
    n_cells = grid.n_cells
    min_support = int(mapping["field"]["minimum_distinct_species_support_per_cell_per_subset_per_outer"])
    min_joint = int(mapping["cross_prediction"]["minimum_joint_cells_required"])
    if min_support != 5 or min_joint != 10:
        raise RuntimeError("species-disjoint technical thresholds drift")

    cache = build_pairwise_jsd_cache(pool)
    identity = np.arange(len(pool.photo_ids), dtype=np.int64)
    null_indices = np.arange(999, dtype=np.int64)
    source_null = null_source_row_matrix(pool, null_indices, master_seed=int(g1["null_master_seed"]))

    observed_train_num = np.zeros((n_folds, n_cells), dtype=float)
    observed_test_num = np.zeros((n_folds, n_cells), dtype=float)
    train_opp = np.zeros((n_folds, n_cells), dtype=float)
    test_opp = np.zeros((n_folds, n_cells), dtype=float)
    null_train_num = np.zeros((n_folds, 999, n_cells), dtype=float)
    null_test_num = np.zeros((n_folds, 999, n_cells), dtype=float)

    species_to_index = {str(label): i for i, label in enumerate(pool.species_labels)}

    for outer_index in range(schedule.n_outer):
        geometry = prepare_sparse_outer_geometry(
            pool,
            schedule.outer_photo_ids[outer_index],
            schedule.outer_species[outer_index],
            grid=grid,
            k=int(g1["k"]),
            kernel_km=float(g1["kernel_bandwidth_km"]),
            cutoff_multiplier=float(g1["kernel_cutoff_multiplier"]),
            minimum_distinct_species=1,
        )
        outer_folds = np.asarray(
            [fold_by_species[species_to_index[str(label)]] for label in schedule.outer_species[outer_index]],
            dtype=np.int16,
        )
        edge_folds = outer_folds[geometry.edge_species_index]

        raw_obs = _raw_jsd_from_source_rows(cache, geometry.edge_nodes, identity)
        rank_obs = _rank_edges(raw_obs, geometry.edge_species_slices)
        _, full_obs_num = _field_from_rank(geometry, rank_obs)

        raw_null = _raw_jsd_matrix_from_source_rows(cache, geometry.edge_nodes, source_null)
        rank_null = _rank_edges_matrix(raw_null, geometry.edge_species_slices)
        _, full_null_num = _field_matrix_from_rank(geometry, rank_null)

        weighted = geometry.weighted_kernel
        full_opp = geometry.opportunity
        full_support = geometry.distinct_species_support.astype(np.int32)

        for fold in range(n_folds):
            test_species_positions = np.flatnonzero(outer_folds == fold)
            if len(test_species_positions) == 0 or len(test_species_positions) == len(outer_folds):
                raise RuntimeError("outer realization has empty train or held-out fold")
            test_edge_mask = edge_folds == fold
            test_weighted = weighted[test_edge_mask]
            local_test_opp = np.asarray(test_weighted.sum(axis=0)).ravel()
            local_train_opp = np.maximum(full_opp - local_test_opp, 0.0)

            test_support = np.zeros(n_cells, dtype=np.int32)
            for species_position in test_species_positions:
                start, stop = geometry.edge_species_slices[int(species_position)]
                contribution = np.asarray(weighted[start:stop].sum(axis=0)).ravel()
                test_support += (contribution > 0).astype(np.int32)
            train_support = full_support - test_support
            if np.any(train_support < 0):
                raise RuntimeError("train support subtraction failure")
            test_eval = (test_support >= min_support) & (local_test_opp > 0)
            train_eval = (train_support >= min_support) & (local_train_opp > 0)

            test_obs_num = np.asarray(test_weighted.T.dot(rank_obs[test_edge_mask])).ravel()
            train_obs_num = full_obs_num - test_obs_num

            observed_test_num[fold, test_eval] += test_obs_num[test_eval]
            observed_train_num[fold, train_eval] += train_obs_num[train_eval]
            test_opp[fold, test_eval] += local_test_opp[test_eval]
            train_opp[fold, train_eval] += local_train_opp[train_eval]

            test_null_num = np.asarray(test_weighted.T.dot(rank_null[:, test_edge_mask].T).T)
            train_null = full_null_num - test_null_num
            null_test_num[fold][:, test_eval] += test_null_num[:, test_eval]
            null_train_num[fold][:, train_eval] += train_null[:, train_eval]

        if (outer_index + 1) % 20 == 0:
            print(json.dumps({"stage": "outer_progress", "completed": outer_index + 1}), flush=True)

    observed_fold_r = np.full(n_folds, np.nan, dtype=float)
    joint_cells = np.zeros(n_folds, dtype=int)
    null_fold_r = np.full((999, n_folds), np.nan, dtype=float)

    for fold in range(n_folds):
        train_field = np.full(n_cells, np.nan, dtype=float)
        test_field = np.full(n_cells, np.nan, dtype=float)
        train_keep = train_opp[fold] > 0
        test_keep = test_opp[fold] > 0
        train_field[train_keep] = observed_train_num[fold, train_keep] / train_opp[fold, train_keep]
        test_field[test_keep] = observed_test_num[fold, test_keep] / test_opp[fold, test_keep]
        joint = train_keep & test_keep & np.isfinite(train_field) & np.isfinite(test_field)
        joint_cells[fold] = int(joint.sum())
        if joint_cells[fold] < min_joint:
            raise RuntimeError(f"fold {fold} has only {joint_cells[fold]} jointly evaluable cells")
        weights = np.sqrt(train_opp[fold, joint] * test_opp[fold, joint])
        observed_fold_r[fold] = weighted_pearson(train_field[joint], test_field[joint], weights)
        if not np.isfinite(observed_fold_r[fold]):
            raise RuntimeError(f"fold {fold} observed correlation is not finite")

        null_train_field = null_train_num[fold][:, joint] / train_opp[fold, joint][None, :]
        null_test_field = null_test_num[fold][:, joint] / test_opp[fold, joint][None, :]
        for null_index in range(999):
            null_fold_r[null_index, fold] = weighted_pearson(
                null_train_field[null_index],
                null_test_field[null_index],
                weights,
            )
        if not np.isfinite(null_fold_r[:, fold]).all():
            raise RuntimeError(f"fold {fold} contains non-finite null correlation")

    observed_global = float(np.median(observed_fold_r))
    null_global = np.median(null_fold_r, axis=1)
    null_mean = float(np.mean(null_global))
    p_upper = float((1 + np.count_nonzero(null_global >= observed_global)) / 1000)
    positive_folds = int(np.count_nonzero(observed_fold_r > 0))
    supported = bool(p_upper < 0.05 and observed_global > null_mean and positive_folds >= 4)

    null_frame = pd.DataFrame({
        "null_index": null_indices,
        "global_median_fold_r": null_global,
    })
    for fold in range(n_folds):
        null_frame[f"fold_{fold}_r"] = null_fold_r[:, fold]
    args.output_null.parent.mkdir(parents=True, exist_ok=True)
    null_frame.to_csv(args.output_null, index=False)

    fold_results = []
    for fold in range(n_folds):
        fold_results.append({
            "fold": fold,
            "held_out_species": int(fold_sizes[fold]),
            "jointly_evaluable_cells": int(joint_cells[fold]),
            "observed_train_vs_heldout_weighted_r": float(observed_fold_r[fold]),
        })

    payload = {
        "protocol": "global-rgfca-species-disjoint-commonness-execution-v1",
        "status": "complete_species_disjoint_commonness_test",
        "question": "Do multiple high-discontinuity spatial zones learned from some species recur in entirely held-out species?",
        "single_boundary_required": False,
        "number_of_zones_fixed": False,
        "primary_g1": {
            "supported": bool(primary["g1_supported"]),
            "p_upper": float(primary["p_upper"]),
            "immutable": True,
            "reclassification_allowed": False
        },
        "species": {
            "eligible": len(pool.species_labels),
            "folds": n_folds,
            "fold_sizes": fold_sizes.tolist(),
            "partition_seed": int(part["seed"]),
            "each_species_held_out_exactly_once": True
        },
        "field": {
            "grid": "36x18",
            "kernel_bandwidth_km": float(g1["kernel_bandwidth_km"]),
            "minimum_distinct_species_support_per_subset_per_outer": min_support,
            "outer_realizations": schedule.n_outer
        },
        "fold_results": fold_results,
        "observed_global_median_fold_r": observed_global,
        "observed_positive_fold_count": positive_folds,
        "null_permutations": 999,
        "null_indices_exact_0_998": True,
        "null_global_mean": null_mean,
        "null_global_q025": float(np.quantile(null_global, 0.025)),
        "null_global_q975": float(np.quantile(null_global, 0.975)),
        "p_upper": p_upper,
        "alpha": 0.05,
        "support_gates": {
            "p_upper_lt_0_05": bool(p_upper < 0.05),
            "observed_above_null_mean": bool(observed_global > null_mean),
            "at_least_four_of_five_fold_correlations_positive": bool(positive_folds >= 4)
        },
        "common_across_species_supported": supported,
        "decision": (
            "support_common_multi_zone_flower_colour_boundary_architecture_across_species"
            if supported else
            "no_species_disjoint_support_for_common_boundary_architecture"
        ),
        "claim_ceiling": (
            "Common multi-zone flower-colour boundary architecture across species is supported; global/universal generalization still requires geographic and taxonomic blocking."
            if supported else
            "Retain recurrent spatial concentration and prespecified robustness only; do not claim common architecture across species."
        ),
        "lineage": {
            "execution_contract_sha256": sha256_file(EXECUTION),
            "measurement_manifest_sha256": sha256_file(MEASUREMENT),
            "measured_table_sha256": sha256_file(MEASURED),
            "primary_g1_result_sha256": sha256_file(PRIMARY),
            "common_boundary_claim_contract_sha256": sha256_file(COMMON),
            "species_disjoint_execution_mapping_sha256": sha256_file(MAPPING),
            "fold_mapping_sha256": sha256_file(args.output_folds),
            "null_result_sha256": sha256_file(args.output_null)
        }
    }
    args.output_result.parent.mkdir(parents=True, exist_ok=True)
    args.output_result.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Replicate the earlier H5 geographic-scale moderation in the balanced RGFCA frame.

This is explicitly post-outcome relative to the completed RGFCA WorldClim global
mean result.  It tests only the already-existing directional moderation hypothesis
and cannot rescue or reclassify any earlier result.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from fcp_pipeline.global_repeated_atlas import build_repeated_atlas_schedule
from fcp_pipeline.global_rgfca_engine import (
    _rank_edges,
    _rank_edges_matrix,
    _raw_jsd_from_source_rows,
    _raw_jsd_matrix_from_source_rows,
    build_pairwise_jsd_cache,
    canonical_colour_pool,
    null_source_row_matrix,
)
from fcp_pipeline.shared_transition_surface import EqualAreaGrid, equal_area_cell_centers, equal_area_cell_ids

ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "docs/supporting/global_rgfca_worldclim_scale_moderation_replication_contract_v1.json"
WORLDCLIM_RESULT = ROOT / "docs/supporting/global_rgfca_worldclim_edge_secondary_result_v1.json"
WORLDCLIM_SPECIES = ROOT / "data/derived/global_rgfca_worldclim_edge_secondary_species_v1.csv"
WORLDCLIM_SCRIPT = ROOT / "scripts/analysis/run_global_rgfca_worldclim_edge_secondary.py"


def load_worldclim_module():
    spec = importlib.util.spec_from_file_location("rgfca_worldclim_secondary", WORLDCLIM_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot import completed WorldClim secondary implementation")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def pearson(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.shape != y.shape or x.ndim != 1 or len(x) < 2:
        return float("nan")
    dx = x - float(np.mean(x))
    dy = y - float(np.mean(y))
    nx = float(np.sqrt(np.sum(dx * dx)))
    ny = float(np.sqrt(np.sum(dy * dy)))
    if nx <= 0.0 or ny <= 0.0:
        return float("nan")
    return float(np.sum(dx * dy) / (nx * ny))


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    return pearson(
        rankdata(np.asarray(x, dtype=float), method="average"),
        rankdata(np.asarray(y, dtype=float), method="average"),
    )


def great_circle_km(lat1, lon1, lat2, lon2) -> np.ndarray:
    lat1r = np.radians(np.asarray(lat1, dtype=float))
    lon1r = np.radians(np.asarray(lon1, dtype=float))
    lat2r = np.radians(np.asarray(lat2, dtype=float))
    lon2r = np.radians(np.asarray(lon2, dtype=float))
    dlat = lat2r - lat1r
    dlon = lon2r - lon1r
    a = np.sin(dlat / 2.0) ** 2 + np.cos(lat1r) * np.cos(lat2r) * np.sin(dlon / 2.0) ** 2
    return 6371.0088 * 2.0 * np.arcsin(np.minimum(1.0, np.sqrt(a)))


def species_span_by_index(pool, photo_complete: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    legacy_grid = EqualAreaGrid(18, 9)
    photo_cells = equal_area_cell_ids(pool.latitude, pool.longitude, legacy_grid)
    _, center_lat, center_lon = equal_area_cell_centers(legacy_grid)
    species_to_index = {str(label): i for i, label in enumerate(pool.species_labels)}
    occupied = [set() for _ in pool.species_labels]
    for row, complete in enumerate(np.asarray(photo_complete, dtype=bool)):
        if not complete:
            continue
        sid = species_to_index[str(pool.species[row])]
        occupied[sid].add(int(photo_cells[row]))
    spans = np.full(len(pool.species_labels), np.nan, dtype=float)
    cell_counts = np.zeros(len(pool.species_labels), dtype=np.int64)
    for sid, cells_set in enumerate(occupied):
        cells = np.asarray(sorted(cells_set), dtype=np.int64)
        cell_counts[sid] = len(cells)
        if len(cells) < 2:
            continue
        ii, jj = np.triu_indices(len(cells), k=1)
        d = great_circle_km(
            center_lat[cells[ii]], center_lon[cells[ii]],
            center_lat[cells[jj]], center_lon[cells[jj]],
        )
        spans[sid] = float(np.mean(d))
    return spans, cell_counts


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-result", type=Path, required=True)
    parser.add_argument("--output-species", type=Path, required=True)
    parser.add_argument("--output-null", type=Path, required=True)
    parser.add_argument("--null-batch-size", type=int, default=64)
    args = parser.parse_args()

    contract = json.loads(CONTRACT.read_text())
    completed = json.loads(WORLDCLIM_RESULT.read_text())
    if contract.get("status") != "postoutcome_replication_frozen_after_worldclim_global_result_before_any_rgfca_scale_moderation_result":
        raise RuntimeError("scale-moderation replication contract is not frozen")
    if completed.get("status") != "complete_global_rgfca_worldclim_edge_secondary":
        raise RuntimeError("completed WorldClim secondary result missing")
    if float(completed.get("p_upper")) != 0.699:
        raise RuntimeError("known WorldClim global result drift")
    if completed.get("secondary_supported") is not False:
        raise RuntimeError("WorldClim global support status drift")

    wc = load_worldclim_module()
    frame, execution, secondary, primary = wc.load_frozen_pool()
    pool = canonical_colour_pool(frame)
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
    blocks = wc.build_edge_blocks(pool, schedule, k=int(g1["k"]))

    lookup, grid_shape, _ = wc.load_worldclim_lookup(secondary)
    climate_grid = EqualAreaGrid(int(grid_shape[0]), int(grid_shape[1]))
    climate_cells = equal_area_cell_ids(pool.latitude, pool.longitude, climate_grid)
    photo_z = lookup[climate_cells]
    photo_complete = np.isfinite(photo_z).all(axis=1)
    spans, occupied_h1_cells = species_span_by_index(pool, photo_complete)

    species_parts = []
    distance_parts = []
    external_parts = []
    cursor = 0
    for block in blocks:
        left = photo_z[block.edge_nodes[:, 0]]
        right = photo_z[block.edge_nodes[:, 1]]
        complete_edge = np.isfinite(left).all(axis=1) & np.isfinite(right).all(axis=1)
        external = np.full(len(block.edge_nodes), np.nan, dtype=float)
        if np.any(complete_edge):
            delta = left[complete_edge] - right[complete_edge]
            external[complete_edge] = np.sqrt(np.mean(delta * delta, axis=1))
        keep = np.isfinite(external)
        block.keep_external = keep
        n_keep = int(np.count_nonzero(keep))
        species_parts.append(block.global_species_index[keep])
        distance_parts.append(block.edge_distance_km[keep])
        external_parts.append(external[keep])
        block.sorted_target_positions = np.arange(cursor, cursor + n_keep, dtype=np.int64)
        cursor += n_keep

    species_unsorted = np.concatenate(species_parts)
    distance_unsorted = np.concatenate(distance_parts)
    external_unsorted = np.concatenate(external_parts)
    sort_order = np.argsort(species_unsorted, kind="stable")
    inverse = np.empty(len(sort_order), dtype=np.int64)
    inverse[sort_order] = np.arange(len(sort_order), dtype=np.int64)
    species_sorted = species_unsorted[sort_order]
    distance_sorted = distance_unsorted[sort_order]
    external_sorted = external_unsorted[sort_order]
    species_counts = np.bincount(species_sorted, minlength=len(pool.species_labels)).astype(np.int64)

    cursor = 0
    for block in blocks:
        keep = block.keep_external
        assert keep is not None
        n_keep = int(np.count_nonzero(keep))
        block.sorted_target_positions = inverse[np.arange(cursor, cursor + n_keep, dtype=np.int64)]
        cursor += n_keep

    cache = build_pairwise_jsd_cache(pool)
    identity = np.arange(len(pool.photo_ids), dtype=np.int64)
    observed_scores = np.empty(len(species_sorted), dtype=float)
    for block in blocks:
        keep = block.keep_external
        target = block.sorted_target_positions
        assert keep is not None and target is not None
        raw = _raw_jsd_from_source_rows(cache, block.edge_nodes, identity)
        rank = _rank_edges(raw, block.edge_species_slices)
        observed_scores[target] = rank[keep]

    minimum_edges = int(secondary["test"]["minimum_edge_occurrences_per_species_after_external_filter"])
    observed_values, observed_species_ids = wc._species_partial_values(
        observed_scores, external_sorted, distance_sorted, species_counts, minimum_edges=minimum_edges
    )

    # Exact replay guard against the already-completed WorldClim species result.
    prior = pd.read_csv(WORLDCLIM_SPECIES).sort_values("species_index", kind="mergesort").reset_index(drop=True)
    replay = pd.DataFrame({
        "species_index": observed_species_ids,
        "partial_rho": observed_values,
    }).sort_values("species_index", kind="mergesort").reset_index(drop=True)
    if not np.array_equal(prior["species_index"].to_numpy(dtype=np.int64), replay["species_index"].to_numpy(dtype=np.int64)):
        raise RuntimeError("WorldClim species IDs do not replay completed result")
    if not np.allclose(
        prior["partial_rho_colour_vs_worldclim_given_distance"].to_numpy(dtype=float),
        replay["partial_rho"].to_numpy(dtype=float), atol=1e-12, rtol=0.0,
    ):
        raise RuntimeError("WorldClim species partial rhos do not replay completed result")

    min_cells = int(contract["moderator"]["minimum_occupied_H1_cells"])
    subset_mask = (
        (occupied_h1_cells[observed_species_ids] >= min_cells)
        & np.isfinite(spans[observed_species_ids])
    )
    eval_ids = observed_species_ids[subset_mask]
    eval_values = observed_values[subset_mask]
    eval_span = spans[eval_ids]
    if len(eval_ids) < int(contract["primary_test"]["minimum_evaluable_species"]):
        raise RuntimeError("fewer than frozen minimum species for scale moderation")
    x = np.log1p(eval_span)
    observed_mod = spearman(x, eval_values)
    if not np.isfinite(observed_mod):
        raise RuntimeError("observed scale-moderation statistic is not finite")
    print(json.dumps({
        "stage": "observed_scale_moderation",
        "evaluable_species": int(len(eval_ids)),
        "rho": float(observed_mod),
    }), flush=True)

    # Map observed species-id order to the fixed moderation subset once.
    pos = {int(sid): i for i, sid in enumerate(observed_species_ids)}
    eval_positions = np.asarray([pos[int(sid)] for sid in eval_ids], dtype=np.int64)

    batch_size = int(args.null_batch_size)
    if batch_size < 1 or batch_size > 256:
        raise ValueError("null batch size must lie in [1,256]")
    null_indices = np.arange(999, dtype=np.int64)
    null_mod = np.empty(999, dtype=float)
    for batch_start in range(0, 999, batch_size):
        batch_indices = null_indices[batch_start:batch_start + batch_size]
        source_rows = null_source_row_matrix(pool, batch_indices, master_seed=int(g1["null_master_seed"]))
        batch_scores = np.empty((len(batch_indices), len(species_sorted)), dtype=float)
        for block in blocks:
            keep = block.keep_external
            target = block.sorted_target_positions
            assert keep is not None and target is not None
            raw = _raw_jsd_matrix_from_source_rows(cache, block.edge_nodes, source_rows)
            rank = _rank_edges_matrix(raw, block.edge_species_slices)
            batch_scores[:, target] = rank[:, keep]
        for local_index, null_index in enumerate(batch_indices):
            values, ids = wc._species_partial_values(
                batch_scores[local_index], external_sorted, distance_sorted, species_counts,
                minimum_edges=minimum_edges,
            )
            if not np.array_equal(ids, observed_species_ids):
                raise RuntimeError(f"null {int(null_index)} changed evaluable species identity")
            value = spearman(x, values[eval_positions])
            if not np.isfinite(value):
                raise RuntimeError(f"null {int(null_index)} scale moderation is not finite")
            null_mod[int(null_index)] = value
        print(json.dumps({
            "stage": "null_batch_complete",
            "first_null": int(batch_indices[0]),
            "last_null": int(batch_indices[-1]),
            "completed": int(batch_indices[-1]) + 1,
            "total": 999,
        }), flush=True)

    p_lower = float((1 + np.count_nonzero(null_mod <= observed_mod)) / 1000)
    alpha = float(contract["null"]["alpha"])
    supported = bool(observed_mod < 0 and p_lower < alpha)

    args.output_result.parent.mkdir(parents=True, exist_ok=True)
    args.output_species.parent.mkdir(parents=True, exist_ok=True)
    args.output_null.parent.mkdir(parents=True, exist_ok=True)

    species_table = pd.DataFrame({
        "species": [pool.species_labels[int(i)] for i in eval_ids],
        "species_index": eval_ids,
        "occupied_legacy_H1_18x9_cells": occupied_h1_cells[eval_ids],
        "mean_pairwise_H1_cell_centroid_great_circle_km": eval_span,
        "log1p_mean_pairwise_H1_cell_centroid_great_circle_km": x,
        "partial_rho_colour_vs_worldclim_given_distance": eval_values,
    })
    species_table.to_csv(args.output_species, index=False, lineterminator="\n")
    pd.DataFrame({
        "null_index": null_indices,
        "scale_moderation_spearman_rho": null_mod,
    }).to_csv(args.output_null, index=False, lineterminator="\n")

    result = {
        "protocol": contract["protocol"],
        "status": "complete_global_rgfca_worldclim_scale_moderation_replication",
        "inferential_role": contract["inferential_role"],
        "response_source": "completed global-rgfca-worldclim-edge-secondary-execution-v1 species effects",
        "moderator": contract["moderator"]["name"],
        "moderator_transform": "log1p",
        "minimum_occupied_H1_cells": min_cells,
        "evaluable_species": int(len(eval_ids)),
        "observed_scale_moderation_spearman_rho": float(observed_mod),
        "null_permutations": 999,
        "null_indices_exact_0_998": True,
        "null_mean": float(np.mean(null_mod)),
        "null_median": float(np.median(null_mod)),
        "null_q025": float(np.quantile(null_mod, 0.025)),
        "null_q975": float(np.quantile(null_mod, 0.975)),
        "p_lower": p_lower,
        "alpha": alpha,
        "supported": supported,
        "claim_ceiling": contract["support_rule"]["claim_if_supported"] if supported else contract["support_rule"]["claim_if_not_supported"],
        "execution_audit": {
            "completed_worldclim_species_effects_replayed_within_1e_12": True,
            "moderator_uses_colour": False,
            "alternative_moderators_run": False,
            "sensitivities_run": False,
            "same_999_RGFCA_colour_nulls": True,
            "same_200_outer_edge_schedule": True,
        },
        "lineage": {
            "contract_sha256": wc.sha256_file(CONTRACT),
            "worldclim_result_sha256": wc.sha256_file(WORLDCLIM_RESULT),
            "worldclim_species_sha256": wc.sha256_file(WORLDCLIM_SPECIES),
            "measured_table_sha256": wc.sha256_file(wc.MEASURED),
        },
        "files": {
            "species": str(args.output_species),
            "null": str(args.output_null),
        },
    }
    args.output_result.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Run the frozen post-outcome WorldClim edge secondary validation.

This analysis cannot rescue or alter the frozen primary G1, species-disjoint
commonness, or RESOLVE results. It reuses exact pre-colour WorldClim bytes only
as a separately labelled secondary test of the shifted process-recurrence
hypothesis.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from fcp_pipeline.global_edge_mechanism_v2 import partial_spearman_distance
from fcp_pipeline.global_repeated_atlas import build_repeated_atlas_schedule
from fcp_pipeline.global_rgfca_engine import (
    COLOUR_COLUMNS,
    _rank_edges,
    _rank_edges_matrix,
    _raw_jsd_from_source_rows,
    _raw_jsd_matrix_from_source_rows,
    build_pairwise_jsd_cache,
    canonical_colour_pool,
    null_source_row_matrix,
    prepare_sparse_outer_geometry,
)
from fcp_pipeline.shared_transition_surface import (
    EqualAreaGrid,
    equal_area_cell_centers,
    equal_area_cell_ids,
)
from fcp_pipeline.spatial_graph import spherical_knn_edges

ROOT = Path(__file__).resolve().parents[2]
EXECUTION = ROOT / "docs/supporting/global_monte_carlo_inference_execution_contract_v1.json"
MEASUREMENT = ROOT / "docs/supporting/global_monte_carlo_measurement_result_v1.json"
MEASURED = ROOT / "data/derived/global_monte_carlo_measured_photos_v1.csv"
PRIMARY_G1 = ROOT / "docs/supporting/global_rgfca_g1_result_v1.json"
COMMONNESS = ROOT / "docs/supporting/global_rgfca_species_disjoint_commonness_result_v1.json"
RESOLVE_RESULT = ROOT / "docs/supporting/global_rgfca_resolve_edge_mechanism_result_v1.json"
SECONDARY = ROOT / "docs/supporting/global_rgfca_worldclim_edge_secondary_execution_v1.json"
SOURCE_MANIFEST = ROOT / "data/atlas/environment/climate_ecoregion_source_manifest.json"
OLD_H2_CONTRACT = ROOT / "docs/supporting/random_photo_first_h2_climate_contract_v1.json"
CLIMATE_GRID = ROOT / "data/atlas/environment/climate_ecoregion_grid_250km.csv"


@dataclass
class EdgeBlock:
    edge_nodes: np.ndarray
    edge_species_slices: tuple[tuple[int, int], ...]
    global_species_index: np.ndarray
    edge_distance_km: np.ndarray
    keep_external: np.ndarray | None = None
    sorted_target_positions: np.ndarray | None = None


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_frozen_pool() -> tuple[pd.DataFrame, dict, dict, dict]:
    execution = json.loads(EXECUTION.read_text())
    measurement = json.loads(MEASUREMENT.read_text())
    primary = json.loads(PRIMARY_G1.read_text())
    common = json.loads(COMMONNESS.read_text())
    resolve = json.loads(RESOLVE_RESULT.read_text())
    secondary = json.loads(SECONDARY.read_text())

    if secondary.get("status") != "postoutcome_secondary_execution_mapping_frozen_before_any_rgfca_worldclim_edge_alignment_result":
        raise RuntimeError("WorldClim secondary execution mapping is not frozen")
    if primary.get("status") != "complete_global_rgfca_g1_primary_inference":
        raise RuntimeError("primary G1 result missing")
    if primary.get("g1_supported") is not False or float(primary.get("p_upper")) != 0.07:
        raise RuntimeError("primary G1 result drift")
    if common.get("common_across_species_supported") is not False:
        raise RuntimeError("species-disjoint commonness result drift")
    if resolve.get("supported") is not False:
        raise RuntimeError("RESOLVE result drift")
    if secondary["inferential_role"].get("may_reclassify_primary_G1") is not False:
        raise RuntimeError("primary G1 rescue protection drift")

    expected_measured = measurement.get("lineage", {}).get("measured_table_sha256")
    if expected_measured != sha256_file(MEASURED):
        raise RuntimeError("measured table hash drift")
    if measurement.get("postmeasurement_gate", {}).get("pass") is not True:
        raise RuntimeError("postmeasurement gate no longer passes")

    frame = pd.read_csv(MEASURED)
    required = {"photo_id", "species", "latitude", "longitude", "global_classifiable", *COLOUR_COLUMNS}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise RuntimeError(f"measured table lacks required columns: {missing}")
    classifiable = frame["global_classifiable"].astype(str).str.casefold().isin({"true", "1"})
    pool = frame.loc[classifiable, ["photo_id", "species", "latitude", "longitude", *COLOUR_COLUMNS]].copy()
    minimum_pool = int(execution["input_gate"]["minimum_classifiable_photos_per_species"])
    counts = pool.groupby("species", observed=True).size()
    eligible = set(counts[counts >= minimum_pool].index.astype(str))
    pool = pool.loc[pool["species"].astype(str).isin(eligible)].copy().reset_index(drop=True)
    if int(pool["species"].nunique()) != int(measurement["postmeasurement_gate"]["evaluable_species"]):
        raise RuntimeError("eligible species count drift")
    if len(pool) != int(primary["classifiable_pool_rows"]):
        raise RuntimeError("classifiable pool row count drift")
    return pool, execution, secondary, primary


def load_worldclim_lookup(secondary: dict) -> tuple[np.ndarray, np.ndarray, dict]:
    manifest = json.loads(SOURCE_MANIFEST.read_text())
    old_h2 = json.loads(OLD_H2_CONTRACT.read_text())
    source = secondary["external_source"]

    expected_grid_sha = source["grid_content_sha256"]
    if sha256_file(CLIMATE_GRID) != expected_grid_sha:
        raise RuntimeError("frozen 250-km climate grid hash drift")
    if manifest["output_sha256"].get("climate_ecoregion_grid_250km.csv") != expected_grid_sha:
        raise RuntimeError("source manifest 250-km grid hash mismatch")
    if manifest["archive_sha256"].get("wc2.1_10m_bio.zip") != source["source_archive_sha256"]:
        raise RuntimeError("WorldClim source archive hash mismatch")
    if old_h2["environment_source"]["primary_grid_source"].get("content_sha256") != expected_grid_sha:
        raise RuntimeError("pre-outcome H2 grid lineage mismatch")
    if old_h2["environment_source"].get("source_archive_sha256") != source["source_archive_sha256"]:
        raise RuntimeError("pre-outcome H2 archive lineage mismatch")

    variables = list(source["variables"])
    if variables != ["bio1", "bio4", "bio12", "bio15"]:
        raise RuntimeError("WorldClim variable family drift")
    frame = pd.read_csv(CLIMATE_GRID)
    required = {"cell_id", "latitude", "longitude", *variables}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise RuntimeError(f"climate grid lacks required columns: {missing}")

    grid = EqualAreaGrid(int(source["grid_n_lon"]), int(source["grid_n_sinlat"]))
    if grid.n_cells != 8192:
        raise RuntimeError("250-km equal-area grid dimensions drift")
    cell_id = pd.to_numeric(frame["cell_id"], errors="raise").to_numpy(dtype=np.int64)
    if len(np.unique(cell_id)) != len(cell_id) or np.any((cell_id < 0) | (cell_id >= grid.n_cells)):
        raise RuntimeError("invalid or duplicated climate cell ids")
    all_ids, center_lat, center_lon = equal_area_cell_centers(grid)
    source_lat = pd.to_numeric(frame["latitude"], errors="raise").to_numpy(dtype=float)
    source_lon = pd.to_numeric(frame["longitude"], errors="raise").to_numpy(dtype=float)
    if np.max(np.abs(source_lat - center_lat[cell_id])) > 1e-9:
        raise RuntimeError("climate grid latitude centers do not match frozen equal-area geometry")
    if np.max(np.abs(source_lon - center_lon[cell_id])) > 1e-9:
        raise RuntimeError("climate grid longitude centers do not match frozen equal-area geometry")

    values = frame[variables].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    complete = np.isfinite(values).all(axis=1)
    if int(np.count_nonzero(complete)) < 2:
        raise RuntimeError("fewer than two complete climate source cells")
    mean = np.mean(values[complete], axis=0)
    sd = np.std(values[complete], axis=0, ddof=0)
    if not np.isfinite(sd).all() or np.any(sd <= 0):
        raise RuntimeError("cannot standardize one or more WorldClim variables")
    z = (values - mean[None, :]) / sd[None, :]
    lookup = np.full((grid.n_cells, len(variables)), np.nan, dtype=float)
    lookup[cell_id[complete]] = z[complete]
    audit = {
        "source_rows": int(len(frame)),
        "complete_source_rows": int(np.count_nonzero(complete)),
        "grid_cells": int(grid.n_cells),
        "variable_means_complete_source": {name: float(value) for name, value in zip(variables, mean)},
        "variable_sd_ddof0_complete_source": {name: float(value) for name, value in zip(variables, sd)},
    }
    return lookup, np.asarray([grid.n_lon, grid.n_sinlat], dtype=np.int64), audit


def build_edge_blocks(pool, schedule, *, k: int) -> list[EdgeBlock]:
    id_to_row = {int(pid): i for i, pid in enumerate(pool.photo_ids)}
    species_to_global = {str(label): i for i, label in enumerate(pool.species_labels)}
    blocks: list[EdgeBlock] = []
    for outer in range(schedule.n_outer):
        nodes: list[np.ndarray] = []
        global_species: list[np.ndarray] = []
        distances: list[np.ndarray] = []
        slices: list[tuple[int, int]] = []
        cursor = 0
        for label_raw, ids in zip(schedule.outer_species[outer], schedule.outer_photo_ids[outer]):
            label = str(label_raw)
            rows = np.asarray([id_to_row[int(pid)] for pid in ids], dtype=np.int64)
            if np.any(pool.species[rows] != label):
                raise RuntimeError("frozen outer schedule assigns a photo to the wrong species")
            local_edges, edge_distance = spherical_knn_edges(pool.latitude[rows], pool.longitude[rows], k=int(k))
            global_edges = rows[local_edges]
            nodes.append(global_edges)
            sid = int(species_to_global[label])
            global_species.append(np.full(len(global_edges), sid, dtype=np.int32))
            distances.append(np.asarray(edge_distance, dtype=float))
            slices.append((cursor, cursor + len(global_edges)))
            cursor += len(global_edges)
        blocks.append(EdgeBlock(
            edge_nodes=np.vstack(nodes),
            edge_species_slices=tuple(slices),
            global_species_index=np.concatenate(global_species),
            edge_distance_km=np.concatenate(distances),
        ))

    audit = prepare_sparse_outer_geometry(
        pool,
        schedule.outer_photo_ids[0],
        schedule.outer_species[0],
        k=int(k),
        kernel_km=500.0,
        cutoff_multiplier=3.0,
        minimum_distinct_species=5,
    )
    if not np.array_equal(blocks[0].edge_nodes, audit.edge_nodes):
        raise RuntimeError("secondary edge construction differs from exact G1 geometry")
    if blocks[0].edge_species_slices != audit.edge_species_slices:
        raise RuntimeError("secondary species edge slices differ from exact G1 geometry")
    return blocks


def _species_partial_values(
    colour_scores_sorted: np.ndarray,
    external_sorted: np.ndarray,
    distance_sorted: np.ndarray,
    species_counts: np.ndarray,
    *,
    minimum_edges: int,
) -> tuple[np.ndarray, np.ndarray]:
    values: list[float] = []
    species_ids: list[int] = []
    cursor = 0
    for sid, count_raw in enumerate(species_counts):
        count = int(count_raw)
        start, stop = cursor, cursor + count
        cursor = stop
        if count < int(minimum_edges):
            continue
        try:
            rho = partial_spearman_distance(
                colour_scores_sorted[start:stop],
                external_sorted[start:stop],
                distance_sorted[start:stop],
            )
        except ValueError:
            continue
        if np.isfinite(rho):
            values.append(float(rho))
            species_ids.append(int(sid))
    if cursor != len(colour_scores_sorted):
        raise RuntimeError("species slice accounting failure")
    return np.asarray(values, dtype=float), np.asarray(species_ids, dtype=np.int32)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-result", type=Path, required=True)
    parser.add_argument("--output-species", type=Path, required=True)
    parser.add_argument("--output-null", type=Path, required=True)
    parser.add_argument("--null-batch-size", type=int, default=64)
    args = parser.parse_args()

    frame, execution, secondary, primary = load_frozen_pool()
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
    blocks = build_edge_blocks(pool, schedule, k=int(g1["k"]))
    total_occurrences = int(sum(len(block.edge_nodes) for block in blocks))
    print(json.dumps({"stage": "edges_built", "outer": len(blocks), "edge_occurrences": total_occurrences}), flush=True)

    lookup, grid_shape, climate_audit = load_worldclim_lookup(secondary)
    grid = EqualAreaGrid(int(grid_shape[0]), int(grid_shape[1]))
    photo_cells = equal_area_cell_ids(pool.latitude, pool.longitude, grid)
    photo_z = lookup[photo_cells]
    photo_complete = np.isfinite(photo_z).all(axis=1)
    print(json.dumps({
        "stage": "climate_join_ready",
        "photos_with_complete_climate": int(np.count_nonzero(photo_complete)),
        "photos_total": int(len(photo_complete)),
        "photo_coverage_fraction": float(np.mean(photo_complete)),
    }), flush=True)

    species_parts: list[np.ndarray] = []
    distance_parts: list[np.ndarray] = []
    external_parts: list[np.ndarray] = []
    cursor = 0
    for block in blocks:
        left = photo_z[block.edge_nodes[:, 0]]
        right = photo_z[block.edge_nodes[:, 1]]
        complete = np.isfinite(left).all(axis=1) & np.isfinite(right).all(axis=1)
        external = np.full(len(block.edge_nodes), np.nan, dtype=float)
        if np.any(complete):
            delta = left[complete] - right[complete]
            external[complete] = np.sqrt(np.mean(delta * delta, axis=1))
        keep = np.isfinite(external)
        block.keep_external = keep
        n_keep = int(np.count_nonzero(keep))
        species_parts.append(block.global_species_index[keep])
        distance_parts.append(block.edge_distance_km[keep])
        external_parts.append(external[keep])
        block.sorted_target_positions = np.arange(cursor, cursor + n_keep, dtype=np.int64)
        cursor += n_keep
    if cursor == 0:
        raise RuntimeError("no edge occurrences have complete frozen WorldClim coverage")

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
    del species_parts, distance_parts, external_parts, species_unsorted, distance_unsorted, external_unsorted, sort_order, inverse

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

    test = secondary["test"]
    minimum_edges = int(test["minimum_edge_occurrences_per_species_after_external_filter"])
    minimum_species = int(test["minimum_evaluable_species"])
    observed_values, observed_species_ids = _species_partial_values(
        observed_scores, external_sorted, distance_sorted, species_counts, minimum_edges=minimum_edges
    )
    if len(observed_values) < minimum_species:
        raise RuntimeError(f"only {len(observed_values)} species evaluable; frozen minimum is {minimum_species}")
    observed_mean = float(np.mean(observed_values))
    observed_median = float(np.median(observed_values))
    observed_positive = float(np.mean(observed_values > 0))
    print(json.dumps({
        "stage": "observed_secondary",
        "n_species": int(len(observed_values)),
        "mean_partial_rho": observed_mean,
        "median_partial_rho": observed_median,
        "positive_fraction": observed_positive,
    }), flush=True)

    batch_size = int(args.null_batch_size)
    if batch_size < 1 or batch_size > 256:
        raise ValueError("null batch size must lie in [1,256]")
    null_indices = np.arange(999, dtype=np.int64)
    null_stat = np.empty(999, dtype=float)
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
            values, _ = _species_partial_values(
                batch_scores[local_index], external_sorted, distance_sorted, species_counts,
                minimum_edges=minimum_edges,
            )
            if len(values) != len(observed_values):
                raise RuntimeError(
                    f"null {int(null_index)} changed evaluable species count from {len(observed_values)} to {len(values)}"
                )
            null_stat[int(null_index)] = float(np.mean(values))
        print(json.dumps({
            "stage": "null_batch_complete",
            "first_null": int(batch_indices[0]),
            "last_null": int(batch_indices[-1]),
            "completed": int(batch_indices[-1]) + 1,
            "total": 999,
        }), flush=True)
        del source_rows, batch_scores

    p_upper = float((1 + np.count_nonzero(null_stat >= observed_mean)) / 1000)
    alpha = float(test["alpha"])
    supported = bool(observed_mean > 0 and p_upper < alpha)

    args.output_result.parent.mkdir(parents=True, exist_ok=True)
    args.output_species.parent.mkdir(parents=True, exist_ok=True)
    args.output_null.parent.mkdir(parents=True, exist_ok=True)

    species_table = pd.DataFrame({
        "species": [pool.species_labels[int(i)] for i in observed_species_ids],
        "species_index": observed_species_ids,
        "retained_edge_occurrences": species_counts[observed_species_ids],
        "partial_rho_colour_vs_worldclim_given_distance": observed_values,
    })
    species_table.to_csv(args.output_species, index=False, lineterminator="\n")
    null_table = pd.DataFrame({"null_index": null_indices, "global_mean_species_partial_rho": null_stat})
    null_table.to_csv(args.output_null, index=False, lineterminator="\n")

    payload = {
        "protocol": secondary["protocol"],
        "status": "complete_global_rgfca_worldclim_edge_secondary",
        "inferential_role": secondary["inferential_role"],
        "primary_g1_supported": bool(primary["g1_supported"]),
        "primary_g1_p_upper": float(primary["p_upper"]),
        "species_disjoint_common_architecture_supported": False,
        "resolve_boundary_alignment_supported": False,
        "predictor": secondary["external_predictor"]["name"],
        "predictor_source": "exact pre-colour-frozen WorldClim 2.1 250-km equal-area grid",
        "source_archive_sha256": secondary["external_source"]["source_archive_sha256"],
        "source_grid_sha256": sha256_file(CLIMATE_GRID),
        "statistic": test["statistic"],
        "observed_mean_species_partial_rho": observed_mean,
        "observed_median_species_partial_rho": observed_median,
        "observed_positive_species_fraction": observed_positive,
        "n_evaluable_species": int(len(observed_values)),
        "minimum_evaluable_species": minimum_species,
        "total_scheduled_edge_occurrences": total_occurrences,
        "climate_evaluable_edge_occurrences": int(len(species_sorted)),
        "climate_edge_occurrence_coverage_fraction": float(len(species_sorted) / total_occurrences),
        "photos_with_complete_climate": int(np.count_nonzero(photo_complete)),
        "photo_climate_coverage_fraction": float(np.mean(photo_complete)),
        "climate_source_audit": climate_audit,
        "null_permutations": 999,
        "null_indices_exact_0_998": True,
        "null_mean": float(np.mean(null_stat)),
        "null_q025": float(np.quantile(null_stat, 0.025)),
        "null_q975": float(np.quantile(null_stat, 0.975)),
        "p_upper": p_upper,
        "alpha": alpha,
        "secondary_supported": supported,
        "claim_ceiling": secondary["claim_if_supported"] if supported else secondary["claim_if_not_supported"],
        "causal_language_allowed": False,
        "execution_audit": {
            "all_200_outer_realizations_used": len(blocks) == 200,
            "same_G1_edge_geometry_audited": True,
            "repeated_edge_occurrences_deduplicated": False,
            "missing_climate_recoded_to_zero": False,
            "single_predictor_only": True,
            "BIO_axis_decomposition_run": False,
            "null_batching_changed_statistic": False
        },
        "lineage": {
            "secondary_execution_mapping_sha256": sha256_file(SECONDARY),
            "precolour_worldclim_source_manifest_sha256": sha256_file(SOURCE_MANIFEST),
            "preoutcome_old_H2_contract_sha256": sha256_file(OLD_H2_CONTRACT),
            "measurement_manifest_sha256": sha256_file(MEASUREMENT),
            "measured_table_sha256": sha256_file(MEASURED),
            "primary_g1_result_sha256": sha256_file(PRIMARY_G1),
            "species_disjoint_result_sha256": sha256_file(COMMONNESS),
            "resolve_result_sha256": sha256_file(RESOLVE_RESULT),
            "species_result_sha256": sha256_file(args.output_species),
            "null_result_sha256": sha256_file(args.output_null)
        },
        "files": {
            "species": str(args.output_species),
            "null": str(args.output_null)
        }
    }
    args.output_result.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

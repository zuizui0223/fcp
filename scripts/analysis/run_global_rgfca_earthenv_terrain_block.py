#!/usr/bin/env python3
"""Run the source-ready EarthEnv terrain member of the expanded environment panel.

The five-block panel and its nonselective execution order were frozen before any
expanded-panel block result.  This script computes only the terrain block raw
Monte Carlo result.  Final block support remains pending Holm adjustment across
all evaluable members of the fixed five-block panel.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio

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

ROOT = Path(__file__).resolve().parents[2]
PANEL_CONTRACT = ROOT / "docs/supporting/global_rgfca_expanded_environmental_process_panel_contract_v1.json"
PANEL_EXECUTION = ROOT / "docs/supporting/global_rgfca_expanded_environmental_panel_execution_v1.json"
TERRAIN_SOURCE = ROOT / "docs/supporting/global_rgfca_earthenv_terrain_source_verification_v1.json"
WORLDCLIM_SCRIPT = ROOT / "scripts/analysis/run_global_rgfca_worldclim_edge_secondary.py"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_worldclim_module():
    name = "rgfca_worldclim_secondary_for_terrain"
    spec = importlib.util.spec_from_file_location(name, WORLDCLIM_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot import frozen RGFCA edge implementation")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def load_terrain_z_at_photos(
    raster_paths: dict[str, Path], latitude: np.ndarray, longitude: np.ndarray
) -> tuple[np.ndarray, dict[str, object]]:
    fixed = ["elevation", "slope", "terrain_roughness_index", "vector_ruggedness_measure"]
    z_columns: list[np.ndarray] = []
    audit: dict[str, object] = {}
    reference_shape = None
    reference_transform = None
    reference_crs = None

    for name in fixed:
        path = raster_paths[name]
        with rasterio.open(path) as ds:
            arr = ds.read(1).astype(np.float64)
            nodata = ds.nodata
            valid = np.isfinite(arr)
            if nodata is not None:
                valid &= arr != float(nodata)
            values = arr[valid]
            if len(values) < 2:
                raise RuntimeError(f"terrain raster {name} has insufficient finite support")
            mean = float(np.mean(values))
            sd = float(np.std(values, ddof=0))
            if not np.isfinite(sd) or sd <= 0:
                raise RuntimeError(f"terrain raster {name} cannot be standardized")

            if reference_shape is None:
                reference_shape = (int(ds.height), int(ds.width))
                reference_transform = tuple(ds.transform)
                reference_crs = str(ds.crs)
            else:
                if (int(ds.height), int(ds.width)) != reference_shape:
                    raise RuntimeError("EarthEnv terrain raster grid shape mismatch")
                if not np.allclose(np.asarray(tuple(ds.transform)), np.asarray(reference_transform), atol=0, rtol=0):
                    raise RuntimeError("EarthEnv terrain raster transform mismatch")
                if str(ds.crs) != reference_crs:
                    raise RuntimeError("EarthEnv terrain raster CRS mismatch")

            sampled = np.full(len(latitude), np.nan, dtype=np.float64)
            rows, cols = rasterio.transform.rowcol(
                ds.transform,
                np.asarray(longitude, dtype=float),
                np.asarray(latitude, dtype=float),
            )
            rows = np.asarray(rows, dtype=np.int64)
            cols = np.asarray(cols, dtype=np.int64)
            inside = (
                (rows >= 0) & (rows < ds.height) & (cols >= 0) & (cols < ds.width)
            )
            if np.any(inside):
                raw = arr[rows[inside], cols[inside]]
                ok = np.isfinite(raw)
                if nodata is not None:
                    ok &= raw != float(nodata)
                idx = np.flatnonzero(inside)
                sampled[idx[ok]] = (raw[ok] - mean) / sd
            z_columns.append(sampled)
            audit[name] = {
                "global_finite_cells": int(np.count_nonzero(valid)),
                "global_mean": mean,
                "global_sd_ddof0": sd,
                "photo_complete_count": int(np.count_nonzero(np.isfinite(sampled))),
                "raster_sha256": sha256_file(path),
            }

    z = np.column_stack(z_columns)
    audit["joint_photo_complete_count"] = int(np.count_nonzero(np.isfinite(z).all(axis=1)))
    audit["photo_count"] = int(len(latitude))
    audit["joint_photo_complete_fraction"] = float(np.mean(np.isfinite(z).all(axis=1)))
    return z, audit


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--terrain-dir", type=Path, required=True)
    parser.add_argument("--output-result", type=Path, required=True)
    parser.add_argument("--output-species", type=Path, required=True)
    parser.add_argument("--output-null", type=Path, required=True)
    parser.add_argument("--null-batch-size", type=int, default=64)
    args = parser.parse_args()

    panel = json.loads(PANEL_CONTRACT.read_text())
    execution_map = json.loads(PANEL_EXECUTION.read_text())
    source = json.loads(TERRAIN_SOURCE.read_text())
    if panel.get("status") != "postoutcome_secondary_panel_frozen_before_any_expanded_environmental_panel_result_or_new_source_payload_is_opened":
        raise RuntimeError("expanded environmental panel contract drift")
    if execution_map.get("status") != "technical_execution_order_frozen_before_any_expanded_panel_block_colour_alignment_result":
        raise RuntimeError("expanded environmental execution-order drift")
    if execution_map.get("execution_order", [None])[0] != "terrain_structure":
        raise RuntimeError("terrain is no longer the frozen first source-ready block")
    if source.get("status") != "pass_exact_earthenv_terrain_source_acquisition_before_any_terrain_colour_alignment":
        raise RuntimeError("exact EarthEnv terrain source verification missing")
    if source.get("terrain_statistic_computed") is not False or source.get("colour_join_performed") is not False:
        raise RuntimeError("terrain source firewall drift")

    expected_files = {entry["id"]: entry for entry in source["files"]}
    expected_ids = ["elevation", "slope", "terrain_roughness_index", "vector_ruggedness_measure"]
    if sorted(expected_files) != sorted(expected_ids):
        raise RuntimeError("terrain source family drift")
    raster_paths: dict[str, Path] = {}
    for ident in expected_ids:
        entry = expected_files[ident]
        path = args.terrain_dir / entry["filename"]
        if not path.exists():
            raise RuntimeError(f"missing frozen terrain file: {path}")
        if sha256_file(path) != entry["sha256"]:
            raise RuntimeError(f"terrain source SHA drift: {ident}")
        raster_paths[ident] = path

    wc = load_worldclim_module()
    frame, inference, _secondary, primary = wc.load_frozen_pool()
    pool = canonical_colour_pool(frame)
    outer = inference["outer_schedule"]
    g1 = inference["g1_primary"]
    schedule = build_repeated_atlas_schedule(
        pool.photo_ids,
        pool.species,
        n_outer=int(outer["observed_resamples"]),
        species_per_outer=int(outer["species_per_resample"]),
        photos_per_species=int(outer["photos_per_species"]),
        minimum_pool_photos_per_species=int(inference["input_gate"]["minimum_classifiable_photos_per_species"]),
        species_seed=int(outer["species_seed"]),
        photo_master_seed=int(outer["photo_master_seed"]),
    )
    blocks = wc.build_edge_blocks(pool, schedule, k=int(g1["k"]))
    total_occurrences = int(sum(len(block.edge_nodes) for block in blocks))
    print(json.dumps({"stage":"edges_built","outer":len(blocks),"edge_occurrences":total_occurrences}), flush=True)

    photo_z, terrain_audit = load_terrain_z_at_photos(
        raster_paths, pool.latitude, pool.longitude
    )
    print(json.dumps({"stage":"terrain_join_ready", **terrain_audit}), flush=True)

    species_parts: list[np.ndarray] = []
    distance_parts: list[np.ndarray] = []
    external_parts: list[np.ndarray] = []
    cursor = 0
    for block in blocks:
        left = photo_z[block.edge_nodes[:, 0]]
        right = photo_z[block.edge_nodes[:, 1]]
        keep = np.isfinite(left).all(axis=1) & np.isfinite(right).all(axis=1)
        external = np.full(len(block.edge_nodes), np.nan, dtype=float)
        if np.any(keep):
            delta = left[keep] - right[keep]
            external[keep] = np.sqrt(np.mean(delta * delta, axis=1))
        block.keep_external = keep
        n_keep = int(np.count_nonzero(keep))
        species_parts.append(block.global_species_index[keep])
        distance_parts.append(block.edge_distance_km[keep])
        external_parts.append(external[keep])
        block.sorted_target_positions = np.arange(cursor, cursor + n_keep, dtype=np.int64)
        cursor += n_keep
    if cursor == 0:
        raise RuntimeError("no RGFCA edge occurrences have complete EarthEnv terrain coverage")

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

    minimum_edges = int(execution_map["terrain_execution_mapping"]["minimum_edges_per_species"])
    minimum_species = int(execution_map["terrain_execution_mapping"]["minimum_species"])
    observed_values, observed_species_ids = wc._species_partial_values(
        observed_scores, external_sorted, distance_sorted, species_counts,
        minimum_edges=minimum_edges,
    )
    if len(observed_values) < minimum_species:
        raise RuntimeError(f"only {len(observed_values)} terrain-evaluable species; frozen minimum {minimum_species}")
    observed_mean = float(np.mean(observed_values))
    observed_median = float(np.median(observed_values))
    observed_positive = float(np.mean(observed_values > 0))
    print(json.dumps({
        "stage":"observed_terrain_block",
        "n_species":int(len(observed_values)),
        "mean_partial_rho":observed_mean,
        "median_partial_rho":observed_median,
        "positive_fraction":observed_positive,
    }), flush=True)

    batch_size = int(args.null_batch_size)
    if batch_size < 1 or batch_size > 256:
        raise ValueError("null batch size must lie in [1,256]")
    null_indices = np.arange(999, dtype=np.int64)
    null_stat = np.empty(999, dtype=float)
    for batch_start in range(0, 999, batch_size):
        batch_indices = null_indices[batch_start:batch_start + batch_size]
        source_rows = null_source_row_matrix(
            pool, batch_indices,
            master_seed=int(execution_map["terrain_execution_mapping"]["null_master_seed"]),
        )
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
                raise RuntimeError(f"null {int(null_index)} changed terrain-evaluable species identity")
            null_stat[int(null_index)] = float(np.mean(values))
        print(json.dumps({
            "stage":"null_batch_complete",
            "first_null":int(batch_indices[0]),
            "last_null":int(batch_indices[-1]),
            "completed":int(batch_indices[-1])+1,
            "total":999,
        }), flush=True)

    p_upper = float((1 + np.count_nonzero(null_stat >= observed_mean)) / 1000)
    alpha = float(panel["primary_edge_test"]["alpha"])
    raw_positive = bool(observed_mean > 0 and p_upper < alpha)

    args.output_result.parent.mkdir(parents=True, exist_ok=True)
    args.output_species.parent.mkdir(parents=True, exist_ok=True)
    args.output_null.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({
        "species": [pool.species_labels[int(i)] for i in observed_species_ids],
        "species_index": observed_species_ids,
        "retained_edge_occurrences": species_counts[observed_species_ids],
        "partial_rho_colour_vs_terrain_given_distance": observed_values,
    }).to_csv(args.output_species, index=False, lineterminator="\n")
    pd.DataFrame({
        "null_index": null_indices,
        "global_mean_species_partial_rho": null_stat,
    }).to_csv(args.output_null, index=False, lineterminator="\n")

    payload = {
        "protocol": panel["protocol"],
        "status": "complete_terrain_structure_block_raw_pending_full_panel_holm",
        "block": "terrain_structure",
        "inferential_role": panel["inferential_role"],
        "final_panel_support_decision_available": False,
        "reason_final_support_pending": "Holm adjustment requires all evaluable fixed blocks in the nonselective five-block panel",
        "predictor_source": "exact frozen EarthEnv GMTED-derived 50-km median elevation+slope+TRI+VRM",
        "source_artifact": execution_map["terrain_execution_mapping"]["source_artifact_name"],
        "statistic": panel["primary_edge_test"]["global_statistic"],
        "observed_mean_species_partial_rho": observed_mean,
        "observed_median_species_partial_rho": observed_median,
        "observed_positive_species_fraction": observed_positive,
        "n_evaluable_species": int(len(observed_values)),
        "minimum_evaluable_species": minimum_species,
        "total_scheduled_edge_occurrences": total_occurrences,
        "terrain_evaluable_edge_occurrences": int(len(species_sorted)),
        "terrain_edge_occurrence_coverage_fraction": float(len(species_sorted) / total_occurrences),
        "terrain_source_audit": terrain_audit,
        "null_permutations": 999,
        "null_indices_exact_0_998": True,
        "null_mean": float(np.mean(null_stat)),
        "null_q025": float(np.quantile(null_stat, 0.025)),
        "null_q975": float(np.quantile(null_stat, 0.975)),
        "raw_p_upper": p_upper,
        "raw_nominal_positive_before_panel_Holm": raw_positive,
        "alpha": alpha,
        "p_holm": None,
        "panel_supported": None,
        "individual_variable_decomposition_open": False,
        "causal_language_allowed": False,
        "execution_audit": {
            "all_200_outer_realizations_used": len(blocks) == 200,
            "same_RGFCA_edge_geometry": True,
            "repeated_edge_occurrences_deduplicated": False,
            "missing_terrain_recoded_to_zero": False,
            "four_fixed_terrain_axes_only": True,
            "same_999_RGFCA_colour_nulls": True,
            "terrain_result_cannot_stop_later_fixed_blocks": True,
        },
        "lineage": {
            "panel_contract_sha256": sha256_file(PANEL_CONTRACT),
            "panel_execution_sha256": sha256_file(PANEL_EXECUTION),
            "terrain_source_verification_sha256": sha256_file(TERRAIN_SOURCE),
            "measurement_manifest_sha256": sha256_file(wc.MEASUREMENT),
            "measured_table_sha256": sha256_file(wc.MEASURED),
            "primary_g1_result_sha256": sha256_file(wc.PRIMARY_G1),
            "species_result_sha256": sha256_file(args.output_species),
            "null_result_sha256": sha256_file(args.output_null),
        },
        "files": {
            "species": str(args.output_species),
            "null": str(args.output_null),
        },
    }
    args.output_result.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

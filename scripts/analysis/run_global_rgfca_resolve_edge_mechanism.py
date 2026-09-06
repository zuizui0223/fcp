#!/usr/bin/env python3
"""Run the prospectively frozen RESOLVE biogeographic-boundary edge mechanism.

This is independent of primary G1 support and cannot rescue G1.  The only external
predictor evaluated here is the exact pre-colour RESOLVE 2017 payload, because the
other members of the fixed primary predictor family lack an exact pre-colour
payload freeze and must remain explicitly not evaluable.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from fcp_pipeline.global_edge_mechanism_v2 import holm_adjust, partial_spearman_distance
from fcp_pipeline.global_repeated_atlas import build_repeated_atlas_schedule
from fcp_pipeline.global_resolve_boundary import load_resolve_boundary_index, score_resolve_boundary_pairs
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
from fcp_pipeline.spatial_graph import spherical_knn_edges

ROOT = Path(__file__).resolve().parents[2]
EXECUTION = ROOT / "docs/supporting/global_monte_carlo_inference_execution_contract_v1.json"
MEASUREMENT = ROOT / "docs/supporting/global_monte_carlo_measurement_result_v1.json"
MEASURED = ROOT / "data/derived/global_monte_carlo_measured_photos_v1.csv"
PRIMARY_G1 = ROOT / "docs/supporting/global_rgfca_g1_result_v1.json"
MECHANISM_CONTRACT = ROOT / "docs/supporting/global_edge_mechanism_distance_control_amendment_v1.json"
OVERLAY_CONTRACT = ROOT / "docs/supporting/global_colour_biogeography_overlay_contract_v1.json"
RESOLVE_EXECUTION = ROOT / "docs/supporting/global_rgfca_resolve_edge_mechanism_execution_v1.json"
RESOLVE_VERIFICATION = ROOT / "docs/supporting/global_rgfca_resolve2017_source_verification_v1.json"


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
    mechanism = json.loads(RESOLVE_EXECUTION.read_text())
    source = json.loads(RESOLVE_VERIFICATION.read_text())
    old_mechanism = json.loads(MECHANISM_CONTRACT.read_text())
    overlay = json.loads(OVERLAY_CONTRACT.read_text())

    if mechanism.get("status") != "technical_execution_mapping_frozen_before_any_rgfca_colour_resolve_alignment_result":
        raise RuntimeError("RESOLVE edge execution mapping is not frozen")
    if source.get("status") != "pass_exact_precolour_resolve2017_payload_recovered":
        raise RuntimeError("exact pre-colour RESOLVE source is not verified")
    if source.get("archive_sha256") != mechanism["external_source_provenance"]["archive_sha256"]:
        raise RuntimeError("RESOLVE source hash drift")
    if primary.get("status") != "complete_global_rgfca_g1_primary_inference":
        raise RuntimeError("primary G1 result missing")
    if primary.get("g1_supported") is not False or mechanism["execution_policy"].get("primary_g1_may_be_rescued") is not False:
        raise RuntimeError("primary G1 immutability drift")
    if old_mechanism.get("status") != "frozen_before_any_global_monte_carlo_colour_field":
        raise RuntimeError("pre-outcome mechanism contract status drift")
    if overlay["inferential_levels"]["edge_mechanism_alignment"].get("requires_G1") is not False:
        raise RuntimeError("edge mechanism unexpectedly requires G1")
    if overlay["inferential_levels"]["edge_mechanism_alignment"].get("cannot_rescue_null_G1") is not True:
        raise RuntimeError("edge mechanism rescue guard drift")

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
    return pool, execution, mechanism, primary


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
            local_edges, edge_distance = spherical_knn_edges(
                pool.latitude[rows], pool.longitude[rows], k=int(k)
            )
            global_edges = rows[local_edges]
            nodes.append(global_edges)
            sid = int(species_to_global[label])
            global_species.append(np.full(len(global_edges), sid, dtype=np.int32))
            distances.append(np.asarray(edge_distance, dtype=float))
            slices.append((cursor, cursor + len(global_edges)))
            cursor += len(global_edges)
        block = EdgeBlock(
            edge_nodes=np.vstack(nodes),
            edge_species_slices=tuple(slices),
            global_species_index=np.concatenate(global_species),
            edge_distance_km=np.concatenate(distances),
        )
        blocks.append(block)

    # Audit one realization against the exact G1 geometry constructor. This is only
    # an identity check; the global grid/kernels are not used in the mechanism.
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
        raise RuntimeError("direct mechanism edge construction differs from G1 geometry")
    if blocks[0].edge_species_slices != audit.edge_species_slices:
        raise RuntimeError("direct mechanism species edge slices differ from G1 geometry")
    return blocks


def normalized_pairs(edge_nodes: np.ndarray) -> np.ndarray:
    nodes = np.asarray(edge_nodes, dtype=np.int64)
    return np.column_stack([np.minimum(nodes[:, 0], nodes[:, 1]), np.maximum(nodes[:, 0], nodes[:, 1])])


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
    parser.add_argument("--resolve-shapefile", type=Path, required=True)
    parser.add_argument("--output-result", type=Path, required=True)
    parser.add_argument("--output-species", type=Path, required=True)
    parser.add_argument("--output-null", type=Path, required=True)
    parser.add_argument("--null-batch-size", type=int, default=64)
    args = parser.parse_args()

    frame, execution, mechanism, primary = load_frozen_pool()
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

    all_pairs = np.vstack([normalized_pairs(block.edge_nodes) for block in blocks])
    unique_pairs = np.unique(all_pairs, axis=0)
    n_rows = len(pool.photo_ids)
    unique_keys = unique_pairs[:, 0].astype(np.int64) * np.int64(n_rows) + unique_pairs[:, 1].astype(np.int64)
    if len(unique_keys) and np.any(np.diff(unique_keys) <= 0):
        raise RuntimeError("unique pair key ordering failure")
    del all_pairs
    print(json.dumps({"stage": "unique_pairs", "unique_unordered_pairs": int(len(unique_pairs))}), flush=True)

    resolve_index = load_resolve_boundary_index(args.resolve_shapefile)
    pair_scores = score_resolve_boundary_pairs(
        unique_pairs,
        pool.latitude,
        pool.longitude,
        resolve_index,
        maximum_interval_km=float(mechanism["resolve_boundary_score"]["maximum_sample_interval_km"]),
    )
    finite_pair = np.isfinite(pair_scores)
    print(json.dumps({
        "stage": "resolve_scored",
        "resolve_polygons": int(resolve_index.polygon_count),
        "unique_pairs_evaluable": int(finite_pair.sum()),
        "unique_pairs_total": int(len(pair_scores)),
        "unique_pair_coverage_fraction": float(finite_pair.mean()) if len(pair_scores) else 0.0,
    }), flush=True)

    species_parts: list[np.ndarray] = []
    distance_parts: list[np.ndarray] = []
    external_parts: list[np.ndarray] = []
    cursor = 0
    for block in blocks:
        pairs = normalized_pairs(block.edge_nodes)
        keys = pairs[:, 0].astype(np.int64) * np.int64(n_rows) + pairs[:, 1].astype(np.int64)
        positions = np.searchsorted(unique_keys, keys)
        if np.any(positions < 0) or np.any(positions >= len(unique_keys)) or np.any(unique_keys[positions] != keys):
            raise RuntimeError("edge pair lookup failure")
        external = pair_scores[positions]
        keep = np.isfinite(external)
        block.keep_external = keep
        n_keep = int(keep.sum())
        species_parts.append(block.global_species_index[keep])
        distance_parts.append(block.edge_distance_km[keep])
        external_parts.append(external[keep])
        block.sorted_target_positions = np.arange(cursor, cursor + n_keep, dtype=np.int64)
        cursor += n_keep
    if cursor == 0:
        raise RuntimeError("no edge occurrences have complete terrestrial RESOLVE coverage")

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

    # Convert each block's temporary unsorted positions into the final species-sorted positions.
    cursor = 0
    for block in blocks:
        keep = block.keep_external
        assert keep is not None
        n_keep = int(keep.sum())
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

    minimum_edges = int(mechanism["edge_occurrence_estimand"]["minimum_edges_per_species_after_external_coverage_filter"])
    minimum_species = int(mechanism["edge_occurrence_estimand"]["minimum_evaluable_species"])
    observed_values, observed_species_ids = _species_partial_values(
        observed_scores, external_sorted, distance_sorted, species_counts, minimum_edges=minimum_edges
    )
    if len(observed_values) < minimum_species:
        raise RuntimeError(f"only {len(observed_values)} species evaluable; frozen minimum is {minimum_species}")
    observed_mean = float(np.mean(observed_values))
    observed_median = float(np.median(observed_values))
    observed_positive = float(np.mean(observed_values > 0))
    print(json.dumps({
        "stage": "observed_mechanism",
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
            # Match the fail-closed evaluability invariant in edge_mechanism_permutation_test.
            if len(values) != len(observed_values):
                raise RuntimeError(
                    f"null {int(null_index)} changed evaluable species count "
                    f"from {len(observed_values)} to {len(values)}"
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
    adjusted = holm_adjust({"biogeographic_boundary_crossing": p_upper})
    p_holm = float(adjusted["biogeographic_boundary_crossing"])
    alpha = float(mechanism["multiplicity"]["alpha"])
    supported = bool(observed_mean > 0 and p_holm < alpha)

    args.output_result.parent.mkdir(parents=True, exist_ok=True)
    args.output_species.parent.mkdir(parents=True, exist_ok=True)
    args.output_null.parent.mkdir(parents=True, exist_ok=True)

    species_table = pd.DataFrame({
        "species": [pool.species_labels[int(i)] for i in observed_species_ids],
        "species_index": observed_species_ids,
        "retained_edge_occurrences": species_counts[observed_species_ids],
        "partial_rho_colour_vs_resolve_given_distance": observed_values,
    })
    species_table.to_csv(args.output_species, index=False, lineterminator="\n")
    null_table = pd.DataFrame({"null_index": null_indices, "global_mean_species_partial_rho": null_stat})
    null_table.to_csv(args.output_null, index=False, lineterminator="\n")

    family_status = dict(mechanism["fixed_primary_predictor_family_evaluability"])
    payload = {
        "protocol": mechanism["protocol"],
        "status": "complete_global_rgfca_resolve_edge_mechanism",
        "g1_required": False,
        "cannot_rescue_primary_g1": True,
        "primary_g1_supported": bool(primary["g1_supported"]),
        "primary_g1_p_upper": float(primary["p_upper"]),
        "predictor": "biogeographic_boundary_crossing",
        "predictor_source": "exact pre-colour RESOLVE Ecoregions 2017",
        "predictor_archive_sha256": json.loads(RESOLVE_VERIFICATION.read_text())["archive_sha256"],
        "fixed_primary_predictor_family_evaluability": family_status,
        "multiplicity": "Holm across all evaluable predictors in the fixed primary family",
        "effective_evaluable_predictor_count": 1,
        "statistic": "equal-species mean within-species partial Spearman controlling ranked great-circle edge distance",
        "observed_mean_species_partial_rho": observed_mean,
        "observed_median_species_partial_rho": observed_median,
        "observed_positive_species_fraction": observed_positive,
        "n_evaluable_species": int(len(observed_values)),
        "minimum_evaluable_species": minimum_species,
        "total_scheduled_edge_occurrences": total_occurrences,
        "resolve_evaluable_edge_occurrences": int(len(species_sorted)),
        "resolve_edge_occurrence_coverage_fraction": float(len(species_sorted) / total_occurrences),
        "unique_unordered_photo_pairs": int(len(unique_pairs)),
        "resolve_evaluable_unique_pairs": int(finite_pair.sum()),
        "null_permutations": 999,
        "null_indices_exact_0_998": True,
        "null_mean": float(np.mean(null_stat)),
        "null_q025": float(np.quantile(null_stat, 0.025)),
        "null_q975": float(np.quantile(null_stat, 0.975)),
        "p_upper": p_upper,
        "p_holm": p_holm,
        "alpha": alpha,
        "supported": supported,
        "claim_ceiling": mechanism["claim_if_supported"] if supported else "No supported RESOLVE biogeographic-boundary concordance under the frozen distance-controlled edge mechanism test.",
        "causal_language_allowed": False,
        "execution_audit": {
            "all_200_outer_realizations_used": len(blocks) == 200,
            "repeated_edge_occurrences_deduplicated": False,
            "external_scores_memoized_only_by_unique_unordered_pair": True,
            "outside_RESOLVE_recode_to_zero": False,
            "null_batching_changed_statistic": False,
        },
        "lineage": {
            "execution_contract_sha256": sha256_file(EXECUTION),
            "mechanism_preoutcome_contract_sha256": sha256_file(MECHANISM_CONTRACT),
            "overlay_preoutcome_contract_sha256": sha256_file(OVERLAY_CONTRACT),
            "resolve_execution_mapping_sha256": sha256_file(RESOLVE_EXECUTION),
            "resolve_verification_sha256": sha256_file(RESOLVE_VERIFICATION),
            "measurement_manifest_sha256": sha256_file(MEASUREMENT),
            "measured_table_sha256": sha256_file(MEASURED),
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

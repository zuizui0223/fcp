#!/usr/bin/env python3
"""Run the single frozen exploratory thermal x atmospheric-dryness interaction test.

This does not reclassify the unsupported five-block environmental panel.  It asks
whether simultaneous thermal-regime and atmospheric-energy/dryness turnover is
associated with extra colour discontinuity beyond geographic distance and the two
additive environmental main effects.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
from dataclasses import dataclass
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

ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "docs/supporting/global_rgfca_thermal_dryness_interaction_contract_v1.json"
CHELSA_VERIFY = ROOT / "docs/supporting/global_rgfca_chelsa_expanded_source_verification_v1.json"
PANEL_FINAL = ROOT / "docs/supporting/global_rgfca_expanded_environmental_panel_final_result_v1.json"
CHELSA_RUNNER = ROOT / "scripts/analysis/run_global_rgfca_chelsa_block.py"
THERMAL = "thermal_regime"
DRYNESS = "atmospheric_energy_dryness"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_chelsa_runner():
    name = "rgfca_chelsa_parent_for_interaction"
    spec = importlib.util.spec_from_file_location(name, CHELSA_RUNNER)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot import frozen CHELSA parent runner")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def zrank_1d(values: np.ndarray) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    if x.ndim != 1 or len(x) < 2 or not np.all(np.isfinite(x)):
        raise ValueError("zrank requires a finite 1D vector")
    r = rankdata(x, method="average").astype(float)
    sd = float(np.std(r, ddof=0))
    if not np.isfinite(sd) or sd <= 1e-12:
        raise ValueError("ranked variable does not vary")
    return (r - float(np.mean(r))) / sd


def zrank_rows(values: np.ndarray) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    if x.ndim != 2 or not np.all(np.isfinite(x)):
        raise ValueError("row zrank requires a finite matrix")
    r = rankdata(x, axis=1, method="average").astype(float)
    means = np.mean(r, axis=1, keepdims=True)
    sds = np.std(r, axis=1, ddof=0, keepdims=True)
    if np.any(~np.isfinite(sds)) or np.any(sds <= 1e-12):
        raise ValueError("a null response became rank-constant")
    return (r - means) / sds


@dataclass(frozen=True)
class SpeciesInteractionDesign:
    species_index: int
    start: int
    stop: int
    interaction_weight: np.ndarray


def build_species_designs(
    thermal: np.ndarray,
    dryness: np.ndarray,
    distance: np.ndarray,
    colour_observed: np.ndarray,
    species_counts: np.ndarray,
    *,
    minimum_edges: int,
) -> tuple[list[SpeciesInteractionDesign], np.ndarray]:
    designs: list[SpeciesInteractionDesign] = []
    observed_beta: list[float] = []
    cursor = 0
    for species_id, count_raw in enumerate(species_counts):
        count = int(count_raw)
        start, stop = cursor, cursor + count
        cursor = stop
        if count < minimum_edges:
            continue
        t = thermal[start:stop]
        a = dryness[start:stop]
        g = distance[start:stop]
        y = colour_observed[start:stop]
        if not (
            np.all(np.isfinite(t))
            and np.all(np.isfinite(a))
            and np.all(np.isfinite(g))
            and np.all(np.isfinite(y))
        ):
            raise RuntimeError("sorted retained interaction arrays must be complete")
        try:
            zt = zrank_1d(t)
            za = zrank_1d(a)
            zg = zrank_1d(g)
            zy = zrank_1d(y)
        except ValueError:
            continue
        interaction = zt * za
        x = np.column_stack(
            [np.ones(count, dtype=float), zg, zt, za, interaction]
        )
        if int(np.linalg.matrix_rank(x)) != 5:
            continue
        pinv = np.linalg.pinv(x)
        weight = np.asarray(pinv[4], dtype=float)
        beta = float(weight @ zy)
        if not np.isfinite(beta):
            continue
        designs.append(
            SpeciesInteractionDesign(
                species_index=int(species_id),
                start=start,
                stop=stop,
                interaction_weight=weight,
            )
        )
        observed_beta.append(beta)
    if cursor != len(thermal):
        raise RuntimeError("species count partition does not cover retained edges")
    return designs, np.asarray(observed_beta, dtype=float)


def beta_matrix_for_designs(
    colour_batch: np.ndarray, designs: list[SpeciesInteractionDesign]
) -> np.ndarray:
    out = np.empty((colour_batch.shape[0], len(designs)), dtype=float)
    for j, design in enumerate(designs):
        zy = zrank_rows(colour_batch[:, design.start : design.stop])
        out[:, j] = zy @ design.interaction_weight
    if not np.all(np.isfinite(out)):
        raise RuntimeError("nonfinite interaction beta under null")
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-dir", type=Path, required=True)
    ap.add_argument("--output-result", type=Path, required=True)
    ap.add_argument("--output-species", type=Path, required=True)
    ap.add_argument("--output-null", type=Path, required=True)
    ap.add_argument("--null-batch-size", type=int, default=32)
    args = ap.parse_args()

    contract = json.loads(CONTRACT.read_text())
    verify = json.loads(CHELSA_VERIFY.read_text())
    panel_final = json.loads(PANEL_FINAL.read_text())
    if contract.get("status") != "postoutcome_exploratory_interaction_frozen_before_interaction_statistic_is_computed":
        raise RuntimeError("interaction contract drift")
    if verify.get("status") != "pass_exact_expanded_chelsa_source_acquisition_before_any_expanded_chelsa_colour_alignment":
        raise RuntimeError("CHELSA source verification missing")
    if panel_final.get("supported_blocks") != []:
        raise RuntimeError("parent five-block result drift")
    known = contract["known_parent_results_before_this_contract"]
    current_blocks = {x["block"]: x for x in panel_final["blocks"]}
    for block in (THERMAL, DRYNESS):
        expected = known[block]
        actual = current_blocks[block]
        for key in ("observed_mean_species_partial_rho", "raw_p_upper", "p_holm", "supported_after_five_block_holm"):
            if actual[key] != expected[key]:
                raise RuntimeError(f"known parent result drift: {block} {key}")

    entries_by_block: dict[str, list[dict]] = {}
    raster_paths_by_block: dict[str, list[Path]] = {}
    for block in (THERMAL, DRYNESS):
        entries = [x for x in verify["files"] if x["block"] == block]
        if len(entries) != 4:
            raise RuntimeError(f"expected four frozen axes for {block}")
        paths: list[Path] = []
        for entry in entries:
            p = args.source_dir / entry["filename"]
            if not p.exists():
                raise RuntimeError(f"missing CHELSA source: {p}")
            if sha256_file(p) != entry["sha256"]:
                raise RuntimeError(f"CHELSA source SHA drift: {entry['id']}")
            paths.append(p)
        entries_by_block[block] = entries
        raster_paths_by_block[block] = paths

    parent = load_chelsa_runner()
    wc = parent.load_worldclim_module()
    frame, inference, _secondary, _primary = wc.load_frozen_pool()
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
    if int(outer["species_seed"]) != 2026090401 or int(outer["photo_master_seed"]) != 2026090402:
        raise RuntimeError("frozen outer seeds drift")
    if int(g1["null_master_seed"]) != 2026090403:
        raise RuntimeError("frozen null seed drift")
    blocks = wc.build_edge_blocks(pool, schedule, k=int(g1["k"]))
    if len(blocks) != 200 or int(g1["k"]) != 3:
        raise RuntimeError("edge geometry drift")
    total_occurrences = int(sum(len(b.edge_nodes) for b in blocks))
    print(json.dumps({"stage": "edges_built", "edge_occurrences": total_occurrences}), flush=True)

    block_photo_z: dict[str, np.ndarray] = {}
    source_audit: dict[str, dict] = {}
    for block in (THERMAL, DRYNESS):
        columns: list[np.ndarray] = []
        for entry, path in zip(entries_by_block[block], raster_paths_by_block[block]):
            z, audit = parent.stream_global_z_and_sample(path, pool.latitude, pool.longitude)
            columns.append(z)
            source_audit[entry["id"]] = {**audit, "sha256": entry["sha256"], "block": block}
            print(json.dumps({"stage": "axis_ready", "block": block, "axis": entry["id"], "photo_complete": audit["photo_complete_count"]}), flush=True)
        block_photo_z[block] = np.column_stack(columns)

    joint_photo = np.isfinite(block_photo_z[THERMAL]).all(axis=1) & np.isfinite(block_photo_z[DRYNESS]).all(axis=1)
    print(json.dumps({"stage": "joint_photo_ready", "complete": int(np.count_nonzero(joint_photo)), "total": int(len(joint_photo)), "fraction": float(np.mean(joint_photo))}), flush=True)

    species_parts: list[np.ndarray] = []
    distance_parts: list[np.ndarray] = []
    thermal_parts: list[np.ndarray] = []
    dryness_parts: list[np.ndarray] = []
    cursor = 0
    for b in blocks:
        tleft = block_photo_z[THERMAL][b.edge_nodes[:, 0]]
        tright = block_photo_z[THERMAL][b.edge_nodes[:, 1]]
        aleft = block_photo_z[DRYNESS][b.edge_nodes[:, 0]]
        aright = block_photo_z[DRYNESS][b.edge_nodes[:, 1]]
        keep = (
            np.isfinite(tleft).all(axis=1)
            & np.isfinite(tright).all(axis=1)
            & np.isfinite(aleft).all(axis=1)
            & np.isfinite(aright).all(axis=1)
        )
        tscore = np.sqrt(np.mean((tleft[keep] - tright[keep]) ** 2, axis=1))
        ascore = np.sqrt(np.mean((aleft[keep] - aright[keep]) ** 2, axis=1))
        b.keep_external = keep
        n_keep = int(np.count_nonzero(keep))
        species_parts.append(b.global_species_index[keep])
        distance_parts.append(b.edge_distance_km[keep])
        thermal_parts.append(tscore)
        dryness_parts.append(ascore)
        b.sorted_target_positions = np.arange(cursor, cursor + n_keep, dtype=np.int64)
        cursor += n_keep
    if cursor == 0:
        raise RuntimeError("no jointly evaluable thermal-dryness edges")

    species_unsorted = np.concatenate(species_parts)
    distance_unsorted = np.concatenate(distance_parts)
    thermal_unsorted = np.concatenate(thermal_parts)
    dryness_unsorted = np.concatenate(dryness_parts)
    sort_order = np.argsort(species_unsorted, kind="stable")
    inverse = np.empty(len(sort_order), dtype=np.int64)
    inverse[sort_order] = np.arange(len(sort_order), dtype=np.int64)
    species_sorted = species_unsorted[sort_order]
    distance_sorted = distance_unsorted[sort_order]
    thermal_sorted = thermal_unsorted[sort_order]
    dryness_sorted = dryness_unsorted[sort_order]
    species_counts = np.bincount(species_sorted, minlength=len(pool.species_labels)).astype(np.int64)
    cursor = 0
    for b in blocks:
        keep = b.keep_external
        assert keep is not None
        n_keep = int(np.count_nonzero(keep))
        b.sorted_target_positions = inverse[np.arange(cursor, cursor + n_keep, dtype=np.int64)]
        cursor += n_keep

    cache = build_pairwise_jsd_cache(pool)
    identity = np.arange(len(pool.photo_ids), dtype=np.int64)
    observed_scores = np.empty(len(species_sorted), dtype=float)
    for b in blocks:
        keep = b.keep_external
        target = b.sorted_target_positions
        assert keep is not None and target is not None
        raw = _raw_jsd_from_source_rows(cache, b.edge_nodes, identity)
        rank = _rank_edges(raw, b.edge_species_slices)
        observed_scores[target] = rank[keep]

    minimum_edges = int(contract["primary_interaction_test"]["minimum_joint_edges_per_species"])
    minimum_species = int(contract["primary_interaction_test"]["minimum_evaluable_species"])
    designs, observed_beta = build_species_designs(
        thermal_sorted,
        dryness_sorted,
        distance_sorted,
        observed_scores,
        species_counts,
        minimum_edges=minimum_edges,
    )
    if len(designs) < minimum_species:
        raise RuntimeError(f"only {len(designs)} species evaluable for interaction")
    observed_mean = float(np.mean(observed_beta))
    observed_median = float(np.median(observed_beta))
    observed_positive = float(np.mean(observed_beta > 0))
    print(json.dumps({"stage": "observed_interaction", "n_species": len(designs), "mean_beta_interaction": observed_mean, "median_beta_interaction": observed_median, "positive_fraction": observed_positive}), flush=True)

    batch_size = int(args.null_batch_size)
    if batch_size < 1 or batch_size > 64:
        raise ValueError("null batch size must lie in [1,64]")
    null_indices = np.arange(999, dtype=np.int64)
    null_stat = np.empty(999, dtype=float)
    for batch_start in range(0, 999, batch_size):
        batch_indices = null_indices[batch_start : batch_start + batch_size]
        source_rows = null_source_row_matrix(pool, batch_indices, master_seed=int(g1["null_master_seed"]))
        batch_scores = np.empty((len(batch_indices), len(species_sorted)), dtype=float)
        for b in blocks:
            keep = b.keep_external
            target = b.sorted_target_positions
            assert keep is not None and target is not None
            raw = _raw_jsd_matrix_from_source_rows(cache, b.edge_nodes, source_rows)
            rank = _rank_edges_matrix(raw, b.edge_species_slices)
            batch_scores[:, target] = rank[:, keep]
        betas = beta_matrix_for_designs(batch_scores, designs)
        null_stat[batch_indices] = np.mean(betas, axis=1)
        print(json.dumps({"stage": "null_batch_complete", "first_null": int(batch_indices[0]), "last_null": int(batch_indices[-1]), "completed": int(batch_indices[-1]) + 1, "total": 999}), flush=True)

    p_upper = float((1 + np.count_nonzero(null_stat >= observed_mean)) / 1000)
    alpha = float(contract["null_and_decision"]["alpha"])
    supported = bool(observed_mean > 0 and p_upper < alpha)

    args.output_result.parent.mkdir(parents=True, exist_ok=True)
    args.output_species.parent.mkdir(parents=True, exist_ok=True)
    args.output_null.parent.mkdir(parents=True, exist_ok=True)
    species_ids = np.asarray([d.species_index for d in designs], dtype=np.int64)
    pd.DataFrame(
        {
            "species": [pool.species_labels[int(i)] for i in species_ids],
            "species_index": species_ids,
            "retained_joint_edge_occurrences": [d.stop - d.start for d in designs],
            "beta_thermal_x_atmospheric_dryness": observed_beta,
        }
    ).to_csv(args.output_species, index=False, lineterminator="\n")
    pd.DataFrame(
        {"null_index": null_indices, "global_mean_species_beta_thermal_x_atmospheric_dryness": null_stat}
    ).to_csv(args.output_null, index=False, lineterminator="\n")

    payload = {
        "protocol": contract["protocol"],
        "status": "complete_postoutcome_exploratory_thermal_dryness_interaction",
        "inferential_role": "postoutcome_exploratory_interaction_only",
        "parent_five_block_panel_remains_unsupported": True,
        "parent_supported_blocks_remain": [],
        "interaction_formula_precolour_preregistered": False,
        "observed_global_mean_species_beta_interaction": observed_mean,
        "observed_median_species_beta_interaction": observed_median,
        "observed_positive_species_fraction": observed_positive,
        "n_evaluable_species": int(len(designs)),
        "minimum_joint_edges_per_species": minimum_edges,
        "minimum_evaluable_species": minimum_species,
        "total_scheduled_edge_occurrences": total_occurrences,
        "joint_evaluable_edge_occurrences": int(len(species_sorted)),
        "joint_edge_occurrence_coverage_fraction": float(len(species_sorted) / total_occurrences),
        "joint_photo_complete_count": int(np.count_nonzero(joint_photo)),
        "joint_photo_complete_fraction": float(np.mean(joint_photo)),
        "source_audit": source_audit,
        "null_permutations": 999,
        "null_indices_exact_0_998": True,
        "null_mean": float(np.mean(null_stat)),
        "null_q025": float(np.quantile(null_stat, 0.025)),
        "null_q975": float(np.quantile(null_stat, 0.975)),
        "p_upper": p_upper,
        "alpha": alpha,
        "exploratory_interaction_supported": supported,
        "claim_ceiling": contract["claim_ceiling_if_supported"] if supported else contract["claim_if_not_supported"],
        "causal_language_allowed": False,
        "anti_fishing_guards_still_active": contract["anti_fishing_guards"],
        "execution_audit": {
            "all_200_outer_realizations_used": len(blocks) == 200,
            "same_RGFCA_edge_geometry": True,
            "repeated_edge_occurrences_deduplicated": False,
            "missing_environment_recoded_to_zero": False,
            "thermal_and_dryness_main_effects_included": True,
            "great_circle_distance_main_effect_included": True,
            "single_interaction_term_only": True,
            "same_999_RGFCA_colour_nulls": True,
        },
        "lineage": {
            "interaction_contract_sha256": sha256_file(CONTRACT),
            "chelsa_source_verification_sha256": sha256_file(CHELSA_VERIFY),
            "parent_five_block_final_result_sha256": sha256_file(PANEL_FINAL),
            "measurement_manifest_sha256": sha256_file(wc.MEASUREMENT),
            "measured_table_sha256": sha256_file(wc.MEASURED),
            "primary_g1_result_sha256": sha256_file(wc.PRIMARY_G1),
            "species_result_sha256": sha256_file(args.output_species),
            "null_result_sha256": sha256_file(args.output_null),
        },
        "files": {"species": str(args.output_species), "null": str(args.output_null)},
    }
    args.output_result.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

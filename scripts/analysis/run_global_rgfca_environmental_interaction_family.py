#!/usr/bin/env python3
"""Evaluate the frozen ten-pair environmental interaction family.

All ten unordered pairs among the five previously frozen environmental process
blocks are evaluated together.  The colour nulls and full RGFCA edge ranks are
computed once per batch and shared across pairs.  This analysis is explicitly
post-outcome exploratory and cannot reclassify the completed five-block panel.
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

from fcp_pipeline.global_edge_mechanism_v2 import holm_adjust
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
FAMILY = ROOT / "docs/supporting/global_rgfca_environmental_interaction_family_contract_v1.json"
PANEL_FINAL = ROOT / "docs/supporting/global_rgfca_expanded_environmental_panel_final_result_v1.json"
CHELSA_VERIFY = ROOT / "docs/supporting/global_rgfca_chelsa_expanded_source_verification_v1.json"
SOIL_VERIFY = ROOT / "docs/supporting/global_rgfca_soilgrids_expanded_source_verification_v1.json"
TERRAIN_VERIFY = ROOT / "docs/supporting/global_rgfca_earthenv_terrain_source_verification_v1.json"
CHELSA_RUNNER = ROOT / "scripts/analysis/run_global_rgfca_chelsa_block.py"
SOIL_RUNNER = ROOT / "scripts/analysis/run_global_rgfca_soilgrids_edaphic_block.py"
TERRAIN_RUNNER = ROOT / "scripts/analysis/run_global_rgfca_earthenv_terrain_block.py"
SINGLE_INTERACTION = ROOT / "scripts/analysis/run_global_rgfca_thermal_dryness_interaction.py"

CHELSA_BLOCKS = ("thermal_regime", "water_balance", "atmospheric_energy_dryness")
SOIL_BLOCK = "edaphic_regime"
TERRAIN_BLOCK = "terrain_structure"
ALL_BLOCKS = (*CHELSA_BLOCKS, SOIL_BLOCK, TERRAIN_BLOCK)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@dataclass
class PairPrepared:
    pair_id: str
    a: str
    b: str
    ecological_interpretation: str
    sorted_global_indices: np.ndarray
    species_counts: np.ndarray
    designs: list
    observed_beta: np.ndarray


def check_contracts() -> tuple[dict, dict, dict, dict, dict]:
    family = json.loads(FAMILY.read_text())
    panel = json.loads(PANEL_FINAL.read_text())
    chelsa = json.loads(CHELSA_VERIFY.read_text())
    soil = json.loads(SOIL_VERIFY.read_text())
    terrain = json.loads(TERRAIN_VERIFY.read_text())
    if family.get("status") != "postoutcome_exploratory_family_frozen_after_single_interaction_launch_but_before_any_interaction_outcome_is_opened":
        raise RuntimeError("interaction-family contract drift")
    if family["completeness_rule"].get("number_of_pairwise_interactions") != 10:
        raise RuntimeError("interaction-family size drift")
    if family["completeness_rule"].get("all_unordered_pairs_included") is not True:
        raise RuntimeError("interaction family is not complete")
    if panel.get("supported_blocks") != [] or panel.get("any_environmental_process_block_supported") is not False:
        raise RuntimeError("parent five-block result drift")
    if chelsa.get("status") != "pass_exact_expanded_chelsa_source_acquisition_before_any_expanded_chelsa_colour_alignment":
        raise RuntimeError("CHELSA verification drift")
    if soil.get("status") != "pass_exact_expanded_soilgrids_source_acquisition_before_any_edaphic_colour_alignment":
        raise RuntimeError("SoilGrids verification drift")
    if terrain.get("status") != "pass_exact_earthenv_terrain_source_acquisition_before_any_terrain_colour_alignment":
        raise RuntimeError("EarthEnv verification drift")
    return family, panel, chelsa, soil, terrain


def load_environment_at_photos(
    pool,
    *,
    chelsa_dir: Path,
    soil_dir: Path,
    terrain_dir: Path,
    chelsa_verify: dict,
    soil_verify: dict,
    terrain_verify: dict,
) -> tuple[dict[str, np.ndarray], dict[str, object]]:
    chelsa_runner = load_module("rgfca_chelsa_for_interaction_family", CHELSA_RUNNER)
    soil_runner = load_module("rgfca_soil_for_interaction_family", SOIL_RUNNER)
    terrain_runner = load_module("rgfca_terrain_for_interaction_family", TERRAIN_RUNNER)

    photo_blocks: dict[str, np.ndarray] = {}
    audit: dict[str, object] = {}

    for block in CHELSA_BLOCKS:
        entries = [e for e in chelsa_verify["files"] if e["block"] == block]
        if len(entries) != 4:
            raise RuntimeError(f"expected four CHELSA axes for {block}")
        cols: list[np.ndarray] = []
        block_audit: dict[str, object] = {}
        for entry in entries:
            path = chelsa_dir / entry["filename"]
            if not path.exists() or sha256_file(path) != entry["sha256"]:
                raise RuntimeError(f"CHELSA source SHA drift: {entry['id']}")
            z, a = chelsa_runner.stream_global_z_and_sample(path, pool.latitude, pool.longitude)
            cols.append(z)
            block_audit[entry["id"]] = {**a, "sha256": entry["sha256"]}
            print(json.dumps({"stage": "environment_axis_ready", "block": block, "axis": entry["id"], "photo_complete": a["photo_complete_count"]}), flush=True)
        matrix = np.column_stack(cols)
        photo_blocks[block] = matrix
        block_audit["joint_photo_complete_count"] = int(np.count_nonzero(np.isfinite(matrix).all(axis=1)))
        block_audit["joint_photo_complete_fraction"] = float(np.mean(np.isfinite(matrix).all(axis=1)))
        audit[block] = block_audit

    soil_by_key = {(e["property"], e["depth"]): e for e in soil_verify["files"]}
    soil_cols: list[np.ndarray] = []
    soil_audit: dict[str, object] = {}
    for prop in soil_runner.PROPERTIES:
        paths: list[Path] = []
        for depth, _weight in soil_runner.DEPTHS:
            entry = soil_by_key[(prop, depth)]
            p = soil_dir / prop / entry["filename"]
            if not p.exists() or sha256_file(p) != entry["sha256"]:
                raise RuntimeError(f"SoilGrids source SHA drift: {prop} {depth}")
            paths.append(p)
        z, a = soil_runner.derived_property_z_at_photos(paths, pool.latitude, pool.longitude)
        soil_cols.append(z)
        soil_audit[prop] = a
        print(json.dumps({"stage": "environment_axis_ready", "block": SOIL_BLOCK, "axis": prop, "photo_complete": a["photo_complete_count"]}), flush=True)
    soil_matrix = np.column_stack(soil_cols)
    photo_blocks[SOIL_BLOCK] = soil_matrix
    soil_audit["joint_photo_complete_count"] = int(np.count_nonzero(np.isfinite(soil_matrix).all(axis=1)))
    soil_audit["joint_photo_complete_fraction"] = float(np.mean(np.isfinite(soil_matrix).all(axis=1)))
    audit[SOIL_BLOCK] = soil_audit

    terrain_entries = {e["id"]: e for e in terrain_verify["files"]}
    terrain_paths: dict[str, Path] = {}
    for ident, entry in terrain_entries.items():
        p = terrain_dir / entry["filename"]
        if not p.exists() or sha256_file(p) != entry["sha256"]:
            raise RuntimeError(f"EarthEnv source SHA drift: {ident}")
        terrain_paths[ident] = p
    terrain_matrix, terrain_audit = terrain_runner.load_terrain_z_at_photos(
        terrain_paths, pool.latitude, pool.longitude
    )
    photo_blocks[TERRAIN_BLOCK] = terrain_matrix
    audit[TERRAIN_BLOCK] = terrain_audit
    print(json.dumps({"stage": "environment_block_ready", "block": TERRAIN_BLOCK, "joint_photo_complete": terrain_audit["joint_photo_complete_count"]}), flush=True)

    if set(photo_blocks) != set(ALL_BLOCKS):
        raise RuntimeError("not all five process blocks were built")
    return photo_blocks, audit


def build_full_edge_environment(blocks, photo_blocks: dict[str, np.ndarray]):
    species_parts: list[np.ndarray] = []
    distance_parts: list[np.ndarray] = []
    score_parts: dict[str, list[np.ndarray]] = {name: [] for name in ALL_BLOCKS}
    outer_slices: list[tuple[int, int]] = []
    cursor = 0
    for outer in blocks:
        n = len(outer.edge_nodes)
        species_parts.append(outer.global_species_index)
        distance_parts.append(outer.edge_distance_km)
        for name in ALL_BLOCKS:
            z = photo_blocks[name]
            left = z[outer.edge_nodes[:, 0]]
            right = z[outer.edge_nodes[:, 1]]
            keep = np.isfinite(left).all(axis=1) & np.isfinite(right).all(axis=1)
            score = np.full(n, np.nan, dtype=float)
            if np.any(keep):
                delta = left[keep] - right[keep]
                score[keep] = np.sqrt(np.mean(delta * delta, axis=1))
            score_parts[name].append(score)
        outer_slices.append((cursor, cursor + n))
        cursor += n
    species = np.concatenate(species_parts)
    distance = np.concatenate(distance_parts)
    scores = {name: np.concatenate(parts) for name, parts in score_parts.items()}
    return species, distance, scores, outer_slices


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chelsa-dir", type=Path, required=True)
    ap.add_argument("--soil-dir", type=Path, required=True)
    ap.add_argument("--terrain-dir", type=Path, required=True)
    ap.add_argument("--output-result", type=Path, required=True)
    ap.add_argument("--output-species", type=Path, required=True)
    ap.add_argument("--output-null", type=Path, required=True)
    ap.add_argument("--null-batch-size", type=int, default=4)
    args = ap.parse_args()

    family, panel, chelsa_verify, soil_verify, terrain_verify = check_contracts()
    chelsa_runner = load_module("rgfca_chelsa_parent_family", CHELSA_RUNNER)
    wc = chelsa_runner.load_worldclim_module()
    interaction_helper = load_module("rgfca_single_interaction_helpers", SINGLE_INTERACTION)

    frame, inference, _secondary, _primary = wc.load_frozen_pool()
    pool = canonical_colour_pool(frame)
    outer = inference["outer_schedule"]
    g1 = inference["g1_primary"]
    if int(outer["species_seed"]) != 2026090401 or int(outer["photo_master_seed"]) != 2026090402:
        raise RuntimeError("outer sampling seed drift")
    if int(g1["null_master_seed"]) != 2026090403:
        raise RuntimeError("null seed drift")
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
    if len(blocks) != 200 or int(g1["k"]) != 3:
        raise RuntimeError("RGFCA edge geometry drift")
    total_occurrences = int(sum(len(b.edge_nodes) for b in blocks))
    print(json.dumps({"stage": "edges_built", "outer": len(blocks), "edge_occurrences": total_occurrences}), flush=True)

    photo_blocks, source_audit = load_environment_at_photos(
        pool,
        chelsa_dir=args.chelsa_dir,
        soil_dir=args.soil_dir,
        terrain_dir=args.terrain_dir,
        chelsa_verify=chelsa_verify,
        soil_verify=soil_verify,
        terrain_verify=terrain_verify,
    )
    global_species, global_distance, block_scores, outer_slices = build_full_edge_environment(blocks, photo_blocks)
    if len(global_species) != total_occurrences or len(global_distance) != total_occurrences:
        raise RuntimeError("global edge concatenation drift")

    cache = build_pairwise_jsd_cache(pool)
    identity = np.arange(len(pool.photo_ids), dtype=np.int64)
    observed_full = np.empty(total_occurrences, dtype=float)
    for b, (start, stop) in zip(blocks, outer_slices):
        raw = _raw_jsd_from_source_rows(cache, b.edge_nodes, identity)
        observed_full[start:stop] = _rank_edges(raw, b.edge_species_slices)

    min_edges = int(family["pairwise_model"]["minimum_joint_edges_per_species"])
    min_species = int(family["pairwise_model"]["minimum_evaluable_species"])
    pair_prepared: list[PairPrepared] = []
    species_rows: list[dict[str, object]] = []

    for spec in family["interaction_family"]:
        pair_id, a_name, b_name = spec["id"], spec["a"], spec["b"]
        a_all, b_all = block_scores[a_name], block_scores[b_name]
        keep = np.isfinite(a_all) & np.isfinite(b_all) & np.isfinite(global_distance)
        kept_idx = np.flatnonzero(keep)
        species_unsorted = global_species[kept_idx]
        order_idx = np.argsort(species_unsorted, kind="stable")
        sorted_global_indices = kept_idx[order_idx]
        species_sorted = global_species[sorted_global_indices]
        distance_sorted = global_distance[sorted_global_indices]
        a_sorted = a_all[sorted_global_indices]
        b_sorted = b_all[sorted_global_indices]
        colour_sorted = observed_full[sorted_global_indices]
        species_counts = np.bincount(species_sorted, minlength=len(pool.species_labels)).astype(np.int64)
        designs, observed_beta = interaction_helper.build_species_designs(
            a_sorted,
            b_sorted,
            distance_sorted,
            colour_sorted,
            species_counts,
            minimum_edges=min_edges,
        )
        if len(designs) < min_species:
            raise RuntimeError(f"only {len(designs)} species evaluable for {pair_id}")
        if len(observed_beta) != len(designs):
            raise RuntimeError(f"observed beta/design mismatch for {pair_id}")
        for d, beta in zip(designs, observed_beta):
            species_rows.append({
                "interaction": pair_id,
                "block_a": a_name,
                "block_b": b_name,
                "species": pool.species_labels[int(d.species_index)],
                "species_index": int(d.species_index),
                "retained_edge_occurrences": int(species_counts[int(d.species_index)]),
                "beta_interaction": float(beta),
            })
        pair_prepared.append(PairPrepared(
            pair_id=pair_id,
            a=a_name,
            b=b_name,
            ecological_interpretation=spec["ecological_interpretation"],
            sorted_global_indices=sorted_global_indices,
            species_counts=species_counts,
            designs=designs,
            observed_beta=observed_beta,
        ))
        print(json.dumps({
            "stage": "observed_pair_ready",
            "interaction": pair_id,
            "joint_edge_occurrences": int(len(sorted_global_indices)),
            "coverage_fraction": float(len(sorted_global_indices) / total_occurrences),
            "n_species": int(len(designs)),
            "mean_beta": float(np.mean(observed_beta)),
            "median_beta": float(np.median(observed_beta)),
            "positive_fraction": float(np.mean(observed_beta > 0)),
        }), flush=True)

    batch_size = int(args.null_batch_size)
    if batch_size < 1 or batch_size > 8:
        raise ValueError("null batch size must lie in [1,8]")
    null_indices = np.arange(999, dtype=np.int64)
    null_stats: dict[str, np.ndarray] = {p.pair_id: np.empty(999, dtype=float) for p in pair_prepared}

    for batch_start in range(0, 999, batch_size):
        batch_indices = null_indices[batch_start:batch_start + batch_size]
        source_rows = null_source_row_matrix(pool, batch_indices, master_seed=int(g1["null_master_seed"]))
        full_batch = np.empty((len(batch_indices), total_occurrences), dtype=float)
        for b, (start, stop) in zip(blocks, outer_slices):
            raw = _raw_jsd_matrix_from_source_rows(cache, b.edge_nodes, source_rows)
            full_batch[:, start:stop] = _rank_edges_matrix(raw, b.edge_species_slices)
        for prep in pair_prepared:
            pair_batch = full_batch[:, prep.sorted_global_indices]
            betas = interaction_helper.beta_matrix_for_designs(pair_batch, prep.designs)
            if betas.shape != (len(batch_indices), len(prep.designs)):
                raise RuntimeError(f"null beta shape drift for {prep.pair_id}")
            null_stats[prep.pair_id][batch_indices] = np.mean(betas, axis=1)
        print(json.dumps({
            "stage": "null_batch_complete",
            "first_null": int(batch_indices[0]),
            "last_null": int(batch_indices[-1]),
            "completed": int(batch_indices[-1]) + 1,
            "total": 999,
        }), flush=True)

    raw_p: dict[str, float] = {}
    pair_payloads: list[dict[str, object]] = []
    for prep in pair_prepared:
        observed_mean = float(np.mean(prep.observed_beta))
        ns = null_stats[prep.pair_id]
        p_two = float((1 + np.count_nonzero(np.abs(ns) >= abs(observed_mean))) / 1000)
        raw_p[prep.pair_id] = p_two
    adjusted = holm_adjust(raw_p)
    alpha = float(family["null_and_multiplicity"]["alpha"])

    for prep in pair_prepared:
        observed_mean = float(np.mean(prep.observed_beta))
        p_two = raw_p[prep.pair_id]
        p_holm = float(adjusted[prep.pair_id])
        supported = bool(p_holm < alpha)
        sign = "synergistic_positive" if observed_mean > 0 else "antagonistic_or_buffering_negative" if observed_mean < 0 else "zero"
        pair_payloads.append({
            "interaction": prep.pair_id,
            "block_a": prep.a,
            "block_b": prep.b,
            "ecological_interpretation": prep.ecological_interpretation,
            "observed_mean_species_beta_interaction": observed_mean,
            "observed_median_species_beta_interaction": float(np.median(prep.observed_beta)),
            "observed_positive_species_fraction": float(np.mean(prep.observed_beta > 0)),
            "n_evaluable_species": int(len(prep.designs)),
            "joint_edge_occurrences": int(len(prep.sorted_global_indices)),
            "joint_edge_coverage_fraction": float(len(prep.sorted_global_indices) / total_occurrences),
            "null_mean": float(np.mean(null_stats[prep.pair_id])),
            "null_q025": float(np.quantile(null_stats[prep.pair_id], 0.025)),
            "null_q975": float(np.quantile(null_stats[prep.pair_id], 0.975)),
            "raw_p_two_sided": p_two,
            "p_holm": p_holm,
            "supported_after_ten_pair_holm": supported,
            "interaction_sign": sign,
            "variable_level_decomposition_open": supported,
        })

    supported_pairs = [x["interaction"] for x in pair_payloads if x["supported_after_ten_pair_holm"]]
    result = {
        "protocol": family["protocol"],
        "status": "complete_fixed_ten_pair_environmental_interaction_family",
        "inferential_role": family["inferential_role"],
        "family_size": 10,
        "multiplicity": "Holm across exactly all ten frozen block-level pairwise interactions",
        "pair_raw_p": "two-sided empirical absolute-tail permutation p",
        "alpha": alpha,
        "total_scheduled_edge_occurrences": total_occurrences,
        "null_permutations": 999,
        "null_indices_exact_0_998": True,
        "null_master_seed": int(g1["null_master_seed"]),
        "pairs": pair_payloads,
        "supported_pairs": supported_pairs,
        "any_pair_supported": bool(supported_pairs),
        "variable_level_decomposition_open_for": supported_pairs,
        "standalone_thermal_dryness_result_used_for_family_support": False,
        "higher_order_interactions_open": False,
        "cannot_rescue_or_reclassify_five_block_main_effect_panel": True,
        "cannot_rescue_or_reclassify_G1": True,
        "causal_or_local_adaptation_language_allowed": False,
        "source_audit": source_audit,
        "lineage": {
            "family_contract_sha256": sha256_file(FAMILY),
            "parent_panel_final_sha256": sha256_file(PANEL_FINAL),
            "chelsa_verification_sha256": sha256_file(CHELSA_VERIFY),
            "soilgrids_verification_sha256": sha256_file(SOIL_VERIFY),
            "earthenv_verification_sha256": sha256_file(TERRAIN_VERIFY),
        },
    }

    args.output_result.parent.mkdir(parents=True, exist_ok=True)
    args.output_species.parent.mkdir(parents=True, exist_ok=True)
    args.output_null.parent.mkdir(parents=True, exist_ok=True)
    args.output_result.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    pd.DataFrame(species_rows).to_csv(args.output_species, index=False, lineterminator="\n")
    null_frame = pd.DataFrame({"null_index": null_indices})
    for prep in pair_prepared:
        null_frame[prep.pair_id] = null_stats[prep.pair_id]
    null_frame.to_csv(args.output_null, index=False, lineterminator="\n")
    print(json.dumps({
        "stage": "complete",
        "supported_pairs": supported_pairs,
        "pairs": [{"interaction": x["interaction"], "mean_beta": x["observed_mean_species_beta_interaction"], "raw_p_two_sided": x["raw_p_two_sided"], "p_holm": x["p_holm"]} for x in pair_payloads],
    }), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

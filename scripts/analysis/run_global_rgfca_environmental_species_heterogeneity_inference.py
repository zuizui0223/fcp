#!/usr/bin/env python3
"""Exact shared-null inference for the frozen 5+10 environmental species-effect matrix."""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
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
CONTRACT = ROOT / "docs/supporting/global_rgfca_environmental_species_heterogeneity_inference_contract_v1.json"
PANEL_CONTRACT = ROOT / "docs/supporting/global_rgfca_expanded_environmental_process_panel_contract_v1.json"
PANEL_RESULT = ROOT / "docs/supporting/global_rgfca_expanded_environmental_panel_final_result_v1.json"
INTERACTION_RESULT = ROOT / "docs/supporting/global_rgfca_environmental_interaction_family_verification_v1.json"
INTERACTION_RUNNER = ROOT / "scripts/analysis/run_global_rgfca_environmental_interaction_family.py"
CHELSA_RUNNER = ROOT / "scripts/analysis/run_global_rgfca_chelsa_block.py"
SINGLE_INTERACTION = ROOT / "scripts/analysis/run_global_rgfca_thermal_dryness_interaction.py"


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def pc1_fraction(matrix: np.ndarray) -> float:
    x = np.asarray(matrix, dtype=float)
    if x.ndim != 2 or x.shape[1] < 2 or not np.isfinite(x).all():
        raise RuntimeError("PC1 matrix must be complete")
    sd = np.std(x, axis=0, ddof=1)
    if np.any(~np.isfinite(sd)) or np.any(sd <= 0):
        raise RuntimeError("degenerate species-effect feature")
    corr = np.corrcoef(x, rowvar=False)
    eig = np.linalg.eigvalsh(corr)
    return float(eig[-1] / x.shape[1])


def prepare_main(name, scores, global_species, global_distance, observed_full, wc, n_species, minimum_edges):
    keep = np.isfinite(scores) & np.isfinite(global_distance)
    kept = np.flatnonzero(keep)
    order = np.argsort(global_species[kept], kind="stable")
    idx = kept[order]
    sid = global_species[idx]
    counts = np.bincount(sid, minlength=n_species).astype(np.int64)
    values, ids = wc._species_partial_values(
        observed_full[idx], scores[idx], global_distance[idx], counts, minimum_edges=minimum_edges
    )
    if len(ids) != n_species or not np.array_equal(ids, np.arange(n_species, dtype=ids.dtype)):
        raise RuntimeError(f"{name}: expected all {n_species} species")
    return {"name": name, "idx": idx, "counts": counts, "observed": values}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chelsa-dir", type=Path, required=True)
    ap.add_argument("--soil-dir", type=Path, required=True)
    ap.add_argument("--terrain-dir", type=Path, required=True)
    ap.add_argument("--output-result", type=Path, required=True)
    ap.add_argument("--output-null", type=Path, required=True)
    ap.add_argument("--output-species", type=Path, required=True)
    ap.add_argument("--null-batch-size", type=int, default=4)
    args = ap.parse_args()

    contract = json.loads(CONTRACT.read_text())
    panel_contract = json.loads(PANEL_CONTRACT.read_text())
    panel_result = json.loads(PANEL_RESULT.read_text())
    interaction_verification = json.loads(INTERACTION_RESULT.read_text())
    if contract.get("status") != "postoutcome_exploratory_frozen_before_any_species_variance_or_syndrome_null_result_is_computed":
        raise RuntimeError("heterogeneity inference contract drift")
    if panel_result.get("supported_blocks") != [] or panel_result.get("any_environmental_process_block_supported") is not False:
        raise RuntimeError("parent main-effect panel drift")
    if interaction_verification.get("family", {}).get("supported_after_holm_count") != 0:
        raise RuntimeError("parent interaction family drift")
    if contract["null"]["master_seed"] != 2026090403 or contract["null"]["permutations"] != 999:
        raise RuntimeError("shared-null contract drift")

    env = load_module("rgfca_env_family_for_heterogeneity", INTERACTION_RUNNER)
    chelsa = load_module("rgfca_chelsa_for_heterogeneity", CHELSA_RUNNER)
    helper = load_module("rgfca_interaction_helper_for_heterogeneity", SINGLE_INTERACTION)
    family, _panel, chelsa_verify, soil_verify, terrain_verify = env.check_contracts()
    wc = chelsa.load_worldclim_module()
    frame, inference, _secondary, _primary = wc.load_frozen_pool()
    pool = canonical_colour_pool(frame)
    n_species = len(pool.species_labels)
    if n_species != 369:
        raise RuntimeError(f"species frame drift: {n_species}")

    outer = inference["outer_schedule"]
    g1 = inference["g1_primary"]
    schedule = build_repeated_atlas_schedule(
        pool.photo_ids, pool.species,
        n_outer=int(outer["observed_resamples"]),
        species_per_outer=int(outer["species_per_resample"]),
        photos_per_species=int(outer["photos_per_species"]),
        minimum_pool_photos_per_species=int(inference["input_gate"]["minimum_classifiable_photos_per_species"]),
        species_seed=int(outer["species_seed"]),
        photo_master_seed=int(outer["photo_master_seed"]),
    )
    blocks = wc.build_edge_blocks(pool, schedule, k=int(g1["k"]))
    total_edges = int(sum(len(b.edge_nodes) for b in blocks))
    photo_blocks, source_audit = env.load_environment_at_photos(
        pool,
        chelsa_dir=args.chelsa_dir,
        soil_dir=args.soil_dir,
        terrain_dir=args.terrain_dir,
        chelsa_verify=chelsa_verify,
        soil_verify=soil_verify,
        terrain_verify=terrain_verify,
    )
    global_species, global_distance, block_scores, outer_slices = env.build_full_edge_environment(blocks, photo_blocks)
    if len(global_species) != total_edges:
        raise RuntimeError("edge concatenation drift")

    cache = build_pairwise_jsd_cache(pool)
    identity = np.arange(len(pool.photo_ids), dtype=np.int64)
    observed_full = np.empty(total_edges, dtype=float)
    for b, (start, stop) in zip(blocks, outer_slices):
        observed_full[start:stop] = _rank_edges(
            _raw_jsd_from_source_rows(cache, b.edge_nodes, identity), b.edge_species_slices
        )

    main_names = list(contract["fixed_features"]["main_effect_family"])
    interaction_names = list(contract["fixed_features"]["interaction_family"])
    if main_names != list(env.ALL_BLOCKS):
        raise RuntimeError("main feature order drift")
    if interaction_names != [x["id"] for x in family["interaction_family"]]:
        raise RuntimeError("interaction feature order drift")

    min_main_edges = int(panel_contract["primary_edge_test"]["minimum_edges_per_species"])
    mains = [prepare_main(n, block_scores[n], global_species, global_distance, observed_full, wc, n_species, min_main_edges) for n in main_names]

    min_interaction_edges = int(family["pairwise_model"]["minimum_joint_edges_per_species"])
    interactions = []
    for spec in family["interaction_family"]:
        a, b = block_scores[spec["a"]], block_scores[spec["b"]]
        keep = np.isfinite(a) & np.isfinite(b) & np.isfinite(global_distance)
        kept = np.flatnonzero(keep)
        order = np.argsort(global_species[kept], kind="stable")
        idx = kept[order]
        sid = global_species[idx]
        counts = np.bincount(sid, minlength=n_species).astype(np.int64)
        designs, beta = helper.build_species_designs(
            a[idx], b[idx], global_distance[idx], observed_full[idx], counts, minimum_edges=min_interaction_edges
        )
        ids = np.asarray([int(d.species_index) for d in designs], dtype=np.int64)
        if len(ids) != n_species or not np.array_equal(ids, np.arange(n_species)):
            raise RuntimeError(f"{spec['id']}: expected all {n_species} species")
        interactions.append({"name": spec["id"], "idx": idx, "designs": designs, "observed": np.asarray(beta, dtype=float)})

    observed_matrix = np.column_stack([x["observed"] for x in mains] + [x["observed"] for x in interactions])
    if observed_matrix.shape != (369, 15) or not np.isfinite(observed_matrix).all():
        raise RuntimeError("observed 369x15 matrix drift")
    observed_main_var = np.var(observed_matrix[:, :5], axis=0, ddof=1)
    observed_interaction_var = np.var(observed_matrix[:, 5:], axis=0, ddof=1)
    observed_pc = np.asarray([pc1_fraction(observed_matrix[:, :5]), pc1_fraction(observed_matrix)], dtype=float)

    null_main_var = np.empty((999, 5), dtype=float)
    null_interaction_var = np.empty((999, 10), dtype=float)
    null_pc = np.empty((999, 2), dtype=float)
    batch_size = int(args.null_batch_size)
    if batch_size < 1 or batch_size > 8:
        raise ValueError("null batch size must lie in [1,8]")

    for batch_start in range(0, 999, batch_size):
        batch_indices = np.arange(batch_start, min(999, batch_start + batch_size), dtype=np.int64)
        source_rows = null_source_row_matrix(pool, batch_indices, master_seed=2026090403)
        full_batch = np.empty((len(batch_indices), total_edges), dtype=float)
        for edge_block, (start, stop) in zip(blocks, outer_slices):
            raw = _raw_jsd_matrix_from_source_rows(cache, edge_block.edge_nodes, source_rows)
            full_batch[:, start:stop] = _rank_edges_matrix(raw, edge_block.edge_species_slices)

        matrices = np.empty((len(batch_indices), 369, 15), dtype=float)
        for j, prep in enumerate(mains):
            scores = full_batch[:, prep["idx"]]
            for local in range(len(batch_indices)):
                values, ids = wc._species_partial_values(
                    scores[local], block_scores[prep["name"]][prep["idx"]], global_distance[prep["idx"]], prep["counts"], minimum_edges=min_main_edges
                )
                if len(ids) != 369 or not np.array_equal(ids, np.arange(369, dtype=ids.dtype)):
                    raise RuntimeError("null main species identity drift")
                matrices[local, :, j] = values
        for j, prep in enumerate(interactions):
            beta = helper.beta_matrix_for_designs(full_batch[:, prep["idx"]], prep["designs"])
            if beta.shape != (len(batch_indices), 369):
                raise RuntimeError("null interaction beta shape drift")
            matrices[:, :, 5 + j] = beta

        for local, null_index in enumerate(batch_indices):
            m = matrices[local]
            null_main_var[null_index] = np.var(m[:, :5], axis=0, ddof=1)
            null_interaction_var[null_index] = np.var(m[:, 5:], axis=0, ddof=1)
            null_pc[null_index, 0] = pc1_fraction(m[:, :5])
            null_pc[null_index, 1] = pc1_fraction(m)
        print(json.dumps({"stage":"null_batch_complete","first":int(batch_indices[0]),"last":int(batch_indices[-1])}), flush=True)

    main_raw = {n: float((1 + np.count_nonzero(null_main_var[:, j] >= observed_main_var[j])) / 1000) for j, n in enumerate(main_names)}
    interaction_raw = {n: float((1 + np.count_nonzero(null_interaction_var[:, j] >= observed_interaction_var[j])) / 1000) for j, n in enumerate(interaction_names)}
    pc_names = ["main5_pc1_fraction", "all15_pc1_fraction"]
    pc_raw = {n: float((1 + np.count_nonzero(null_pc[:, j] >= observed_pc[j])) / 1000) for j, n in enumerate(pc_names)}
    main_holm = holm_adjust(main_raw)
    interaction_holm = holm_adjust(interaction_raw)
    pc_holm = holm_adjust(pc_raw)
    alpha = 0.05

    main_rows = [{"feature":n,"observed_variance":float(observed_main_var[j]),"null_mean_variance":float(np.mean(null_main_var[:,j])),"raw_p_upper":main_raw[n],"p_holm":float(main_holm[n]),"supported_excess_species_heterogeneity":bool(main_holm[n] < alpha)} for j,n in enumerate(main_names)]
    interaction_rows = [{"feature":n,"observed_variance":float(observed_interaction_var[j]),"null_mean_variance":float(np.mean(null_interaction_var[:,j])),"raw_p_upper":interaction_raw[n],"p_holm":float(interaction_holm[n]),"supported_excess_species_heterogeneity":bool(interaction_holm[n] < alpha)} for j,n in enumerate(interaction_names)]
    pc_rows = [{"statistic":n,"observed_pc1_fraction":float(observed_pc[j]),"null_mean":float(np.mean(null_pc[:,j])),"null_q025":float(np.quantile(null_pc[:,j],0.025)),"null_q975":float(np.quantile(null_pc[:,j],0.975)),"raw_p_upper":pc_raw[n],"p_holm":float(pc_holm[n]),"supported_coordinated_syndrome":bool(pc_holm[n] < alpha)} for j,n in enumerate(pc_names)]

    payload = {
        "protocol": contract["protocol"],
        "status": "complete_exact_shared_null_environmental_species_heterogeneity_inference",
        "inferential_role": "postoutcome_exploratory_heterogeneity_not_mean_effect_rescue",
        "n_species": 369,
        "features": 15,
        "null_permutations": 999,
        "null_indices_exact_0_998": True,
        "null_master_seed": 2026090403,
        "main_effect_variance_family": main_rows,
        "interaction_variance_family": interaction_rows,
        "coordinated_syndrome_omnibus": pc_rows,
        "supported_main_variance_features": [x["feature"] for x in main_rows if x["supported_excess_species_heterogeneity"]],
        "supported_interaction_variance_features": [x["feature"] for x in interaction_rows if x["supported_excess_species_heterogeneity"]],
        "supported_syndrome_statistics": [x["statistic"] for x in pc_rows if x["supported_coordinated_syndrome"]],
        "parent_mean_effects_reclassified": False,
        "individual_variable_decomposition_opened": False,
        "causal_or_local_adaptation_claim_allowed": False,
        "source_audit": source_audit,
    }

    args.output_result.parent.mkdir(parents=True, exist_ok=True)
    args.output_null.parent.mkdir(parents=True, exist_ok=True)
    args.output_species.parent.mkdir(parents=True, exist_ok=True)
    args.output_result.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n")
    species_frame = pd.DataFrame({"species": pool.species_labels})
    for j,n in enumerate(main_names): species_frame[f"main__{n}"] = observed_matrix[:,j]
    for j,n in enumerate(interaction_names): species_frame[f"interaction__{n}"] = observed_matrix[:,5+j]
    species_frame.to_csv(args.output_species, index=False, lineterminator="\n")
    nf = pd.DataFrame({"null_index": np.arange(999,dtype=int)})
    for j,n in enumerate(main_names): nf[f"variance__main__{n}"] = null_main_var[:,j]
    for j,n in enumerate(interaction_names): nf[f"variance__interaction__{n}"] = null_interaction_var[:,j]
    nf["pc1_fraction__main5"] = null_pc[:,0]
    nf["pc1_fraction__all15"] = null_pc[:,1]
    nf.to_csv(args.output_null, index=False, lineterminator="\n")
    print(json.dumps({"stage":"complete","supported_main_variance":payload["supported_main_variance_features"],"supported_interaction_variance":payload["supported_interaction_variance_features"],"supported_syndrome":payload["supported_syndrome_statistics"]}), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

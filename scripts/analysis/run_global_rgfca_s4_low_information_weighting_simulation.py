#!/usr/bin/env python3
"""Run the frozen S4 low-information weighting methodological simulation."""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
EXECUTION = ROOT / "docs/supporting/global_rgfca_s4_low_information_weighting_execution_v1.json"
PARENT = ROOT / "docs/supporting/global_rgfca_design_power_identifiability_simulation_contract_v1.json"
H6B_CONTRACT = ROOT / "docs/supporting/random_photo_first_h6b_photo_level_null_diagnostic_contract_v1.json"
H6B_RUNNER = ROOT / "scripts/analysis/run_random_photo_first_h6b_photo_level_null_diagnostic.py"
H6B_RESULT = ROOT / "docs/supporting/random_photo_first_h6b_photo_level_null_diagnostic_result_v1.json"
H6_RESULT = ROOT / "docs/supporting/random_photo_first_h6_species_specific_spatial_structure_result_v1.json"


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def summarize(x: np.ndarray, weights: np.ndarray) -> dict[str, float]:
    return {
        "simulation_mean": float(np.mean(x)),
        "simulation_sd": float(np.std(x, ddof=1)),
        "q005": float(np.quantile(x, 0.005)),
        "q025": float(np.quantile(x, 0.025)),
        "q975": float(np.quantile(x, 0.975)),
        "q995": float(np.quantile(x, 0.995)),
        "q95_absolute_mean": float(np.quantile(np.abs(x), 0.95)),
        "probability_abs_mean_ge_0_025": float(np.mean(np.abs(x) >= 0.025)),
        "probability_abs_mean_ge_0_05": float(np.mean(np.abs(x) >= 0.05)),
        "probability_abs_mean_ge_0_10": float(np.mean(np.abs(x) >= 0.10)),
        "effective_species_number_1_over_sum_w2": float(1.0 / np.sum(weights * weights)),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output-result", type=Path, required=True)
    ap.add_argument("--output-species", type=Path, required=True)
    ap.add_argument("--output-draws", type=Path, required=True)
    ap.add_argument("--batch-size", type=int, default=5000)
    args = ap.parse_args()

    execution = json.loads(EXECUTION.read_text())
    parent = json.loads(PARENT.read_text())
    h6b_contract = json.loads(H6B_CONTRACT.read_text())
    h6b_result = json.loads(H6B_RESULT.read_text())
    h6_result = json.loads(H6_RESULT.read_text())
    if execution.get("status") != "postoutcome_methodological_execution_frozen_before_any_s4_simulation_result":
        raise RuntimeError("S4 execution contract drift")
    if parent["S4_low_information_weighting_simulation"]["seed"] != 2026090604 or parent["S4_low_information_weighting_simulation"]["replicates"] != 100000:
        raise RuntimeError("S4 parent contract drift")
    if h6b_contract["stricter_null"]["seed"] != 20260908 or h6b_contract["stricter_null"]["permutations"] != 999:
        raise RuntimeError("H6b null drift")
    if h6b_result.get("status") != "complete_h6b_diagnostic_evaluable":
        raise RuntimeError("H6b completion drift")
    if h6_result.get("decision") != "no_support_general_species_specific_spatial_structuring":
        raise RuntimeError("H6 decision drift")

    h6b = load_module("rgfca_h6b_for_s4", H6B_RUNNER)
    measured = pd.read_csv(h6b.MEASURED)
    required = {"measurement_id", "species", "cell_id", "morph", *h6b.PALETTE}
    if not required.issubset(measured.columns):
        raise RuntimeError("measurement input lacks H6b fields")
    photos = measured.loc[measured["morph"].astype(str).isin(h6b.BIOLOGICAL_MORPHS)].copy()
    photos["species"] = photos["species"].astype(str)
    photos = h6b.add_soft_colour(photos)

    capacity = pd.read_csv(h6b.CAPACITY)
    capacity["species"] = capacity["species"].astype(str)
    target_species = sorted(capacity.loc[
        (capacity["classifiable_photos"] >= 5) &
        (capacity["occupied_h1_cells"] >= 3) &
        (capacity["morph_levels"] >= 2), "species"
    ].astype(str))
    reliable_species = set(capacity.loc[
        (capacity["classifiable_photos"] >= 10) &
        (capacity["occupied_h1_cells"] >= 5) &
        (capacity["morph_levels"] >= 2), "species"
    ].astype(str))
    if len(target_species) != 188 or len(reliable_species) != 74:
        raise RuntimeError("S4 sparse/reliable frame drift")

    records = []
    matrix = np.empty((188, 999), dtype=np.float64)
    for i, species in enumerate(target_species):
        rec, null = h6b.prepare_species_photo_null(
            species,
            photos.loc[photos["species"] == species],
            permutations=999,
            seed=20260908,
        )
        centered = np.asarray(null, dtype=np.float64) - float(np.mean(null))
        var = float(np.var(centered, ddof=1))
        if not np.isfinite(var) or var <= 0:
            raise RuntimeError(f"zero/nonfinite species null variance: {species}")
        matrix[i] = centered
        records.append({
            "species": species,
            "classifiable_photos": int(rec["classifiable_photos"]),
            "occupied_h1_cells": int(rec["occupied_h1_cells"]),
            "cell_pairs": int(rec["cell_pairs"]),
            "null_mean_before_centering": float(np.mean(null)),
            "null_sd_ddof1": float(np.sqrt(var)),
            "reliable_74": bool(species in reliable_species),
        })
        if (i + 1) % 20 == 0 or i + 1 == 188:
            print(json.dumps({"stage":"species_nulls","completed":i+1,"total":188}), flush=True)

    table = pd.DataFrame(records)
    pair = table["cell_pairs"].to_numpy(dtype=np.float64)
    w_equal = np.full(188, 1.0 / 188.0)
    w_pair = pair / np.sum(pair)
    invvar = 1.0 / np.square(table["null_sd_ddof1"].to_numpy(dtype=np.float64))
    w_inv = invvar / np.sum(invvar)
    reliable_mask = table["reliable_74"].to_numpy(dtype=bool)
    if int(np.count_nonzero(reliable_mask)) != 74:
        raise RuntimeError("reliable subset drift")
    w_rel = np.full(74, 1.0 / 74.0)
    table["weight_equal188"] = w_equal
    table["weight_paircount188"] = w_pair
    table["weight_inverse_null_variance188"] = w_inv

    reps = 100000
    batch = int(args.batch_size)
    if batch < 100 or batch > 20000:
        raise ValueError("batch size must be in [100,20000]")
    rng = np.random.default_rng(2026090604)
    out_equal = np.empty(reps, dtype=np.float64)
    out_pair = np.empty(reps, dtype=np.float64)
    out_inv = np.empty(reps, dtype=np.float64)
    out_rel = np.empty(reps, dtype=np.float64)
    rows = np.arange(188, dtype=np.int64)[None, :]
    cursor = 0
    while cursor < reps:
        stop = min(reps, cursor + batch)
        n = stop - cursor
        pick = rng.integers(0, 999, size=(n, 188), endpoint=False)
        draw = matrix[rows, pick]
        out_equal[cursor:stop] = draw @ w_equal
        out_pair[cursor:stop] = draw @ w_pair
        out_inv[cursor:stop] = draw @ w_inv
        out_rel[cursor:stop] = draw[:, reliable_mask] @ w_rel
        cursor = stop
        print(json.dumps({"stage":"simulation","completed":cursor,"total":reps}), flush=True)

    estimators = {
        "equal_species_188": (out_equal, w_equal),
        "pair_count_weighted_188": (out_pair, w_pair),
        "inverse_null_variance_188": (out_inv, w_inv),
        "reliable_74_equal_species": (out_rel, w_rel),
    }
    summaries = {name: summarize(values, weights) for name, (values, weights) in estimators.items()}
    var_equal = float(np.var(out_equal, ddof=1)); var_pair = float(np.var(out_pair, ddof=1)); var_inv = float(np.var(out_inv, ddof=1)); var_rel = float(np.var(out_rel, ddof=1))
    rho, rho_p = spearmanr(np.log1p(pair), table["null_sd_ddof1"].to_numpy(dtype=float))
    comparisons = {
        "variance_ratio_equal188_to_inverse_variance188": var_equal / var_inv,
        "variance_ratio_paircount188_to_inverse_variance188": var_pair / var_inv,
        "variance_ratio_equal188_to_reliable74": var_equal / var_rel,
        "spearman_log1p_cell_pairs_vs_species_null_sd": float(rho),
        "spearman_p_descriptive_only": float(rho_p),
    }

    result = {
        "protocol": execution["protocol"],
        "status": "complete_frozen_s4_low_information_weighting_simulation",
        "role": execution["role"],
        "target_species": 188,
        "reliable_species": 74,
        "h6b_null_per_species": 999,
        "h6b_seed": 20260908,
        "simulation_replicates": reps,
        "simulation_seed": 2026090604,
        "observed_species_rho_used_to_set_simulation": False,
        "estimators": summaries,
        "comparisons": comparisons,
        "H6_or_H6b_reclassified": False,
        "claim_ceiling": execution["claim_ceiling"],
    }

    args.output_result.parent.mkdir(parents=True, exist_ok=True)
    args.output_species.parent.mkdir(parents=True, exist_ok=True)
    args.output_draws.parent.mkdir(parents=True, exist_ok=True)
    args.output_result.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    table.to_csv(args.output_species, index=False, lineterminator="\n")
    pd.DataFrame({
        "replicate": np.arange(reps, dtype=int),
        "equal_species_188": out_equal,
        "pair_count_weighted_188": out_pair,
        "inverse_null_variance_188": out_inv,
        "reliable_74_equal_species": out_rel,
    }).to_csv(args.output_draws, index=False, lineterminator="\n")
    print(json.dumps({"stage":"complete","comparisons":comparisons}), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

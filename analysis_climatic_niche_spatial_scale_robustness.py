#!/usr/bin/env python3
"""Robustness checks for within- versus among-population climatic niche contrasts."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

METRICS = ["pca_hull_area", "moisture_breadth"]


def zscore(x: pd.Series) -> pd.Series:
    x = pd.to_numeric(x, errors="coerce")
    sd = x.std(ddof=0)
    return (x - x.mean()) / sd if np.isfinite(sd) and sd > 0 else pd.Series(np.nan, index=x.index)


def fit_beta(d: pd.DataFrame, metric: str) -> float:
    x = d.copy()
    x["among"] = (x["spatial_scale"] == "among_population").astype(int)
    x["metric_z"] = zscore(x[metric])
    x["effort_z"] = zscore(np.log1p(pd.to_numeric(x["n_climate_cells"], errors="coerce")))
    x = x.dropna(subset=["among", "metric_z", "effort_z"])
    if len(x) < 20 or x["among"].nunique() < 2:
        return np.nan
    model = sm.GLM(x["among"], sm.add_constant(x[["metric_z", "effort_z"]], has_constant="add"), family=sm.families.Binomial())
    return float(model.fit().params["metric_z"])


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--dataset", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--permutations", type=int, default=9999)
    p.add_argument("--seed", type=int, default=20260719)
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(args.dataset)
    required = {"canonical_name", "family", "spatial_scale", "n_climate_cells", *METRICS}
    missing = sorted(required - set(data.columns))
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    d = data.loc[data["n_climate_cells"] >= 20].copy()
    rng = np.random.default_rng(args.seed)
    rows = []
    loo_rows = []

    for metric in METRICS:
        observed = fit_beta(d, metric)
        permuted = []
        for _ in range(args.permutations):
            xp = d.copy()
            xp["spatial_scale"] = rng.permutation(xp["spatial_scale"].to_numpy())
            permuted.append(fit_beta(xp, metric))
        permuted = np.asarray(permuted, dtype=float)
        valid = permuted[np.isfinite(permuted)]
        p_two = float((1 + np.sum(np.abs(valid) >= abs(observed))) / (1 + len(valid))) if np.isfinite(observed) else np.nan

        families = sorted(d["family"].dropna().astype(str).unique())
        loo = []
        for fam in families:
            beta = fit_beta(d.loc[d["family"].astype(str) != fam], metric)
            loo.append(beta)
            loo_rows.append({"metric": metric, "omitted_family": fam, "estimate": beta, "odds_ratio": float(np.exp(beta)) if np.isfinite(beta) else np.nan})
        loo_arr = np.asarray(loo, dtype=float)
        valid_loo = loo_arr[np.isfinite(loo_arr)]
        rows.append({
            "metric": metric,
            "n_species": int(len(d)),
            "n_within": int((d["spatial_scale"] == "within_population").sum()),
            "n_among": int((d["spatial_scale"] == "among_population").sum()),
            "estimate": observed,
            "odds_ratio": float(np.exp(observed)),
            "permutation_p_two_sided": p_two,
            "permutations_valid": int(len(valid)),
            "loo_min_odds_ratio": float(np.exp(np.nanmin(valid_loo))) if len(valid_loo) else np.nan,
            "loo_max_odds_ratio": float(np.exp(np.nanmax(valid_loo))) if len(valid_loo) else np.nan,
            "loo_same_direction_fraction": float(np.mean(np.sign(valid_loo) == np.sign(observed))) if len(valid_loo) else np.nan,
        })

    summary = pd.DataFrame(rows)
    loo_table = pd.DataFrame(loo_rows)
    summary.to_csv(outdir / "climatic_niche_spatial_scale_robustness.csv", index=False)
    loo_table.to_csv(outdir / "climatic_niche_spatial_scale_leave_one_family_out.csv", index=False)
    manifest = {
        "min_cells": 20,
        "metrics": METRICS,
        "permutations": args.permutations,
        "results": summary.to_dict("records"),
        "interpretation_guard": "Exploratory robustness analysis of two preselected near-threshold metrics; permutation p-values and leave-one-family-out stability do not establish causation or physiological niche breadth.",
    }
    (outdir / "climatic_niche_spatial_scale_robustness_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

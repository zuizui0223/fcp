#!/usr/bin/env python3
"""Run the frozen 34-species robustness design across all five niche metrics.

This extends the existing focal two-metric robustness analysis without changing
classification, occurrence data, model formula, uncertainty estimator, or seed.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

import analysis_climatic_niche_spatial_scale_robustness as base

METRICS = [
    "temperature_breadth",
    "moisture_breadth",
    "climatic_heterogeneity",
    "pca_dispersion",
    "pca_hull_area",
]


def adjust_within_set(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    out["wald_p_holm_five_metrics"] = np.nan
    out["permutation_p_holm_five_metrics"] = np.nan
    for analysis_set, idx in out.groupby("analysis_set").groups.items():
        rows = list(idx)
        complete = [i for i in rows if out.loc[i, "analysis_status"] == "complete"]
        if not complete:
            continue
        wald = pd.to_numeric(out.loc[complete, "wald_p_value_clustered"], errors="coerce")
        perm = pd.to_numeric(out.loc[complete, "permutation_p_two_sided"], errors="coerce")
        valid_w = wald.notna()
        valid_p = perm.notna()
        if valid_w.any():
            out.loc[wald.index[valid_w], "wald_p_holm_five_metrics"] = multipletests(
                wald[valid_w], method="holm"
            )[1]
        if valid_p.any():
            out.loc[perm.index[valid_p], "permutation_p_holm_five_metrics"] = multipletests(
                perm[valid_p], method="holm"
            )[1]
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--permutations", type=int, default=9999)
    parser.add_argument("--seed", type=int, default=20260719)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(args.dataset)
    required = {
        "canonical_name", "family", "spatial_scale", "classification_source",
        "n_climate_cells", *METRICS,
    }
    missing = sorted(required - set(data.columns))
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    eligible = data.loc[data["n_climate_cells"] >= 20].copy()
    baseline = eligible.loc[
        eligible["classification_source"] == "baseline_unambiguous"
    ].copy()
    if len(baseline) != 34:
        raise ValueError(f"Frozen baseline must contain 34 species; found {len(baseline)}")

    rng = np.random.default_rng(args.seed)
    rows, loo_rows = base.analyse_set(
        baseline,
        "baseline_unambiguous_34",
        METRICS,
        args.permutations,
        rng,
    )
    results = adjust_within_set(pd.DataFrame(rows))
    loo = pd.DataFrame(loo_rows)

    results.to_csv(outdir / "environmental_niche_five_metric_models.csv", index=False)
    loo.to_csv(outdir / "environmental_niche_five_metric_leave_one_family_out.csv", index=False)
    baseline.sort_values("canonical_name").to_csv(
        outdir / "environmental_niche_five_metric_model_dataset.csv", index=False
    )

    manifest = {
        "status": "complete",
        "dataset_role": "frozen_34_species_baseline_only",
        "n_species": int(len(baseline)),
        "n_families": int(baseline["family"].nunique()),
        "n_within": int((baseline["spatial_scale"] == "within_population").sum()),
        "n_among": int((baseline["spatial_scale"] == "among_population").sum()),
        "metrics": METRICS,
        "model_formula": base.MODEL_FORMULA,
        "covariance": "family-clustered sandwich",
        "permutations": args.permutations,
        "seed": args.seed,
        "multiplicity": "Holm correction across the five niche metrics within the frozen baseline",
        "results": results.to_dict("records"),
        "interpretation_guard": (
            "All five metrics are evaluated symmetrically. Moisture breadth remains a reported focal "
            "association identified within this five-metric family rather than a preregistered endpoint. "
            "Occupied climatic breadth is not physiological tolerance and associations are not causal."
        ),
    }
    (outdir / "environmental_niche_five_metric_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

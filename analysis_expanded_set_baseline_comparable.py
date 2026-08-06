#!/usr/bin/env python3
"""Apply the published 34-species focal model unchanged to the expanded unreviewed set.

This intentionally mirrors analysis_climatic_niche_spatial_scale_robustness.py:
- metrics: moisture_breadth and pca_hull_area
- minimum 20 occupied climate cells
- among ~ metric_z + effort_z
- effort = z(log1p(n_climate_cells))
- family-clustered sandwich SE
- 9,999 label permutations
- leave-one-family-out sensitivity

The expanded labels remain unreviewed and this script never creates a frozen manifest.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels
import statsmodels.api as sm

METRICS = ["pca_hull_area", "moisture_breadth"]
MIN_MODEL_SPECIES = 20
MIN_CLIMATE_CELLS = 20
MODEL_FORMULA = "among ~ metric_z + effort_z"
CI_METHOD = "Wald 95% CI from family-clustered sandwich standard errors"


def zscore(x: pd.Series) -> pd.Series:
    x = pd.to_numeric(x, errors="coerce")
    sd = x.std(ddof=0)
    return (x - x.mean()) / sd if np.isfinite(sd) and sd > 0 else pd.Series(np.nan, index=x.index)


def prepare_model_data(d: pd.DataFrame, metric: str) -> pd.DataFrame:
    x = d.copy()
    x["among"] = (x["spatial_scale"] == "among_population").astype(int)
    x["metric_z"] = zscore(x[metric])
    x["effort_z"] = zscore(np.log1p(pd.to_numeric(x["n_climate_cells"], errors="coerce")))
    return x.dropna(subset=["among", "metric_z", "effort_z", "family"])


def fit_model(d: pd.DataFrame, metric: str, clustered: bool = True):
    x = prepare_model_data(d, metric)
    if len(x) < MIN_MODEL_SPECIES or x["among"].nunique() < 2 or x["family"].nunique() < 2:
        return None, x
    model = sm.GLM(
        x["among"],
        sm.add_constant(x[["metric_z", "effort_z"]], has_constant="add"),
        family=sm.families.Binomial(),
    )
    if clustered:
        return model.fit(cov_type="cluster", cov_kwds={"groups": x["family"]}), x
    return model.fit(), x


def fit_beta(d: pd.DataFrame, metric: str) -> float:
    fit, _ = fit_model(d, metric, clustered=False)
    return float(fit.params["metric_z"]) if fit is not None else np.nan


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--classification", required=True)
    p.add_argument("--metrics", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--permutations", type=int, default=9999)
    p.add_argument("--seed", type=int, default=20260719)
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    classification = pd.read_csv(args.classification)
    metrics = pd.read_csv(args.metrics)

    class_required = {"canonical_name", "family", "spatial_scale"}
    metric_required = {"canonical_name", "n_climate_cells", *METRICS}
    missing = sorted((class_required - set(classification.columns)) | (metric_required - set(metrics.columns)))
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    classification = classification.loc[
        classification["spatial_scale"].isin(["within_population", "among_population"]),
        ["canonical_name", "family", "spatial_scale"],
    ].drop_duplicates("canonical_name")
    metric_cols = ["canonical_name", "n_climate_cells", *METRICS]
    data = classification.merge(metrics[metric_cols], on="canonical_name", how="inner", validate="one_to_one")
    data = data.loc[pd.to_numeric(data["n_climate_cells"], errors="coerce") >= MIN_CLIMATE_CELLS].copy()
    data.to_csv(outdir / "expanded_baseline_comparable_model_dataset.csv", index=False)

    rng = np.random.default_rng(args.seed)
    rows: list[dict] = []
    loo_rows: list[dict] = []

    for metric in METRICS:
        fit, model_data = fit_model(data, metric, clustered=True)
        estimable = fit is not None
        observed = float(fit.params["metric_z"]) if estimable else np.nan
        se = float(fit.bse["metric_z"]) if estimable else np.nan

        permuted: list[float] = []
        if estimable:
            for _ in range(args.permutations):
                xp = model_data.copy()
                xp["spatial_scale"] = rng.permutation(xp["spatial_scale"].to_numpy())
                permuted.append(fit_beta(xp, metric))
        valid = np.asarray(permuted, dtype=float)
        valid = valid[np.isfinite(valid)]
        p_two = (
            float((1 + np.sum(np.abs(valid) >= abs(observed))) / (1 + len(valid)))
            if np.isfinite(observed) else np.nan
        )

        loo: list[float] = []
        if estimable:
            for family in sorted(model_data["family"].dropna().astype(str).unique()):
                beta = fit_beta(model_data.loc[model_data["family"].astype(str) != family], metric)
                loo.append(beta)
                loo_rows.append({
                    "analysis_set": "expanded_unreviewed_exact_34spec_model",
                    "metric": metric,
                    "omitted_family": family,
                    "estimate": beta,
                    "odds_ratio": float(np.exp(beta)) if np.isfinite(beta) else np.nan,
                })
        valid_loo = np.asarray(loo, dtype=float)
        valid_loo = valid_loo[np.isfinite(valid_loo)]
        ci_low = observed - 1.96 * se if estimable else np.nan
        ci_high = observed + 1.96 * se if estimable else np.nan
        fit_history = getattr(fit, "fit_history", {}) if estimable else {}

        rows.append({
            "analysis_set": "expanded_unreviewed_exact_34spec_model",
            "metric": metric,
            "analysis_status": "complete" if estimable else "not_estimable",
            "model_formula": MODEL_FORMULA,
            "estimator": "statsmodels GLM Binomial(logit)",
            "covariance": "family-clustered sandwich",
            "ci_method": CI_METHOD,
            "statsmodels_version": statsmodels.__version__,
            "n_species": int(len(model_data)),
            "n_families": int(model_data["family"].nunique()),
            "n_within": int((model_data["spatial_scale"] == "within_population").sum()),
            "n_among": int((model_data["spatial_scale"] == "among_population").sum()),
            "estimate": observed,
            "std_error_clustered": se,
            "estimate_ci_low": ci_low,
            "estimate_ci_high": ci_high,
            "odds_ratio": float(np.exp(observed)) if np.isfinite(observed) else np.nan,
            "odds_ratio_ci_low": float(np.exp(ci_low)) if np.isfinite(ci_low) else np.nan,
            "odds_ratio_ci_high": float(np.exp(ci_high)) if np.isfinite(ci_high) else np.nan,
            "wald_p_value_clustered": float(fit.pvalues["metric_z"]) if estimable else np.nan,
            "permutation_p_two_sided": p_two,
            "permutations_requested": int(args.permutations),
            "permutations_valid": int(len(valid)),
            "converged": bool(getattr(fit, "converged", False)) if estimable else False,
            "iterations": int(fit_history.get("iteration", -1)) if estimable else -1,
            "loo_min_odds_ratio": float(np.exp(np.nanmin(valid_loo))) if len(valid_loo) else np.nan,
            "loo_max_odds_ratio": float(np.exp(np.nanmax(valid_loo))) if len(valid_loo) else np.nan,
            "loo_same_direction_fraction": (
                float(np.mean(np.sign(valid_loo) == np.sign(observed)))
                if len(valid_loo) and np.isfinite(observed) else np.nan
            ),
        })

    summary = pd.DataFrame(rows)
    loo = pd.DataFrame(loo_rows)
    summary.to_csv(outdir / "expanded_exact_34species_model_results.csv", index=False)
    loo.to_csv(outdir / "expanded_exact_34species_leave_one_family_out.csv", index=False)

    manifest = {
        "status": "complete",
        "review_status": "unreviewed",
        "eligible_for_freeze": False,
        "comparison_target": "analysis_climatic_niche_spatial_scale_robustness.py on the frozen 34-species baseline",
        "min_cells": MIN_CLIMATE_CELLS,
        "minimum_model_species": MIN_MODEL_SPECIES,
        "model_formula": MODEL_FORMULA,
        "estimator": "statsmodels GLM Binomial(logit)",
        "covariance": "family-clustered sandwich",
        "confidence_interval_method": CI_METHOD,
        "metrics": METRICS,
        "permutations": args.permutations,
        "seed": args.seed,
        "results": summary.to_dict("records"),
        "interpretation_guard": (
            "This is a model-specification comparison only. Expanded spatial labels remain unreviewed; "
            "the result cannot replace the frozen 34-species ecological analysis before adjudication."
        ),
    }
    (outdir / "expanded_exact_34species_model_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

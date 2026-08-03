#!/usr/bin/env python3
"""Fit climate-breadth models adjusted for geographic range and discovery effort."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.special import expit
from statsmodels.stats.multitest import multipletests

COMPONENTS = [
    "bio1_q95q05", "bio5_q95q05", "bio6_q95q05", "bio7_q95q05",
    "bio12_q95q05", "bio14_q95q05", "bio15_q95q05", "bio17_q95q05",
]
COVARIATES = [
    "geographic_hull_area_km2",
    "occupied_geographic_grid_cells",
    "n_occurrence_records",
    "n_supporting_p1_records",
]


def zscore(series: pd.Series) -> np.ndarray:
    values = series.to_numpy(float)
    sd = np.nanstd(values, ddof=0)
    if not np.isfinite(sd) or sd <= 0:
        return np.full(len(values), np.nan)
    return (values - np.nanmean(values)) / sd


def batch_betas(x: np.ndarray, y: np.ndarray, max_iter: int = 60) -> np.ndarray:
    beta = np.zeros((y.shape[0], x.shape[1]), dtype=float)
    for _ in range(max_iter):
        mu = expit(beta @ x.T)
        weights = np.clip(mu * (1 - mu), 1e-9, None)
        score = (y - mu) @ x
        hessian = np.einsum("ni,mn,nj->mij", x, weights, x, optimize=True)
        delta = np.vstack([
            np.linalg.lstsq(hessian[i], score[i], rcond=None)[0]
            for i in range(len(hessian))
        ])
        beta += delta
        if np.nanmax(np.abs(delta)) < 1e-9:
            break
    return beta


def make_design(data: pd.DataFrame, metric: str):
    y = (data.spatial_scale == "among_population").astype(int).to_numpy()
    term_names = ["intercept", metric, *COVARIATES]
    arrays = [np.ones(len(data)), zscore(data[metric])]
    arrays.extend(zscore(np.log1p(data[column].clip(lower=0))) for column in COVARIATES)
    return np.column_stack(arrays), y, term_names


def analyse(data: pd.DataFrame, metric: str, permutations: int, seed: int):
    x, y, term_names = make_design(data, metric)
    finite = np.isfinite(x).all(axis=1)
    x, y = x[finite], y[finite]
    subset = data.loc[finite].copy()
    if len(subset) < 20 or subset.spatial_scale.nunique() < 2:
        raise ValueError(f"{metric}: insufficient finite data after covariate adjustment")
    fit = sm.GLM(y, x, family=sm.families.Binomial()).fit(
        cov_type="cluster", cov_kwds={"groups": subset.family.astype(str).to_numpy()}
    )
    observed = float(fit.params[1])
    rng = np.random.default_rng(seed)
    permutation_indices = np.vstack([rng.permutation(len(y)) for _ in range(permutations)])
    permutation_betas = batch_betas(x, y[permutation_indices])[:, 1]
    valid = permutation_betas[np.isfinite(permutation_betas)]
    permutation_p = float((1 + np.sum(np.abs(valid) >= abs(observed))) / (1 + len(valid)))

    loo_rows, loo_betas = [], []
    families = subset.family.astype(str).to_numpy()
    for family in sorted(np.unique(families)):
        keep = families != family
        if keep.sum() < x.shape[1] + 2 or len(np.unique(y[keep])) < 2:
            continue
        loo_fit = sm.GLM(y[keep], x[keep], family=sm.families.Binomial()).fit()
        beta = float(loo_fit.params[1])
        loo_betas.append(beta)
        loo_rows.append({"metric": metric, "omitted_family": family, "estimate": beta, "odds_ratio": float(np.exp(beta))})
    loo = np.asarray(loo_betas, dtype=float)
    se = float(fit.bse[1])
    result = {
        "metric": metric,
        "n_species": int(len(subset)),
        "n_families": int(subset.family.nunique()),
        "n_within": int((subset.spatial_scale == "within_population").sum()),
        "n_among": int((subset.spatial_scale == "among_population").sum()),
        "estimate": observed,
        "std_error_family_clustered": se,
        "odds_ratio": float(np.exp(observed)),
        "odds_ratio_ci_low": float(np.exp(observed - 1.96 * se)),
        "odds_ratio_ci_high": float(np.exp(observed + 1.96 * se)),
        "wald_p_value_family_clustered": float(fit.pvalues[1]),
        "permutation_p_two_sided": permutation_p,
        "permutations_valid": int(len(valid)),
        "loo_same_direction_fraction": float(np.mean(np.sign(loo) == np.sign(observed))) if len(loo) else np.nan,
        "loo_min_odds_ratio": float(np.exp(loo.min())) if len(loo) else np.nan,
        "loo_max_odds_ratio": float(np.exp(loo.max())) if len(loo) else np.nan,
        "model_formula": "among ~ climate_component_z + log_hull_area_z + log_geographic_cells_z + log_occurrence_records_z + log_supporting_P1_records_z",
        "classification_source": "exploratory_automated_preserved",
    }
    coefficients = [
        {
            "metric": metric,
            "term": name,
            "estimate": float(fit.params[i]),
            "std_error_family_clustered": float(fit.bse[i]),
            "p_value_family_clustered": float(fit.pvalues[i]),
            "odds_ratio": float(np.exp(fit.params[i])),
        }
        for i, name in enumerate(term_names)
    ]
    return result, loo_rows, coefficients


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--classification", required=True)
    parser.add_argument("--metrics", required=True)
    parser.add_argument("--extent", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--min-cells", type=int, default=20)
    parser.add_argument("--permutations", type=int, default=9999)
    parser.add_argument("--seed", type=int, default=20260804)
    args = parser.parse_args()

    classification = pd.read_csv(args.classification)
    metrics = pd.read_csv(args.metrics)
    extent = pd.read_csv(args.extent)
    data = classification.merge(metrics, on="canonical_name", how="inner", suffixes=("_class", "_metric"))
    data = data.merge(extent, on="canonical_name", how="inner")
    if "family_class" in data:
        data["family"] = data.family_class
    if "classification_source" not in data or set(data.classification_source) != {"exploratory_automated"}:
        raise ValueError("Automated classification provenance was altered before ecological modelling")
    data = data.loc[
        data.spatial_scale.isin(["within_population", "among_population"])
        & (data.metric_status == "complete")
        & (data.n_climate_cells >= args.min_cells)
    ].copy()
    required = {"canonical_name", "family", "spatial_scale", "classification_source", *COMPONENTS, *COVARIATES}
    missing = sorted(required - set(data.columns))
    if missing:
        raise ValueError(f"Missing adjusted-model columns: {missing}")

    model_rows, loo_rows, coefficient_rows = [], [], []
    for offset, metric in enumerate(COMPONENTS):
        result, loo, coefficients = analyse(data, metric, args.permutations, args.seed + offset)
        model_rows.append(result)
        loo_rows.extend(loo)
        coefficient_rows.extend(coefficients)
    models = pd.DataFrame(model_rows)
    models["wald_p_holm_eight_components"] = multipletests(models.wald_p_value_family_clustered, method="holm")[1]
    models["permutation_p_holm_eight_components"] = multipletests(models.permutation_p_two_sided, method="holm")[1]

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    data.to_csv(outdir / "range_adjusted_model_dataset.csv", index=False)
    models.to_csv(outdir / "range_adjusted_climate_component_models.csv", index=False)
    pd.DataFrame(loo_rows).to_csv(outdir / "range_adjusted_leave_one_family_out.csv", index=False)
    pd.DataFrame(coefficient_rows).to_csv(outdir / "range_adjusted_all_coefficients.csv", index=False)
    manifest = {
        "status": "complete",
        "n_species": int(models.n_species.max()),
        "n_within": int(models.n_within.max()),
        "n_among": int(models.n_among.max()),
        "n_families": int(models.n_families.max()),
        "classification_source": "exploratory_automated",
        "focal_components": COMPONENTS,
        "range_and_effort_covariates": COVARIATES,
        "results": models.to_dict("records"),
        "interpretation_guard": "Climate breadth is estimated conditional on sampled range size, GBIF occurrence effort and P1 literature support. Residual confounding and automated-label error remain possible; associations are not causal.",
    }
    (outdir / "range_adjusted_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Fit ecological models that separate climatic breadth from range size and study effort.

Response: among-population versus within-population automated spatial organisation.
For each BIO breadth component, the focal coefficient is estimated while adjusting for
geographic hull area, occupied geographic grid cells, retained occurrence count and
number of supporting P1 sources.  Family-clustered uncertainty, common label
permutations and leave-one-family-out diagnostics are reported.
"""
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
    "n_distinct_supporting_sources",
]


def zscore(series: pd.Series) -> np.ndarray:
    values = series.to_numpy(float)
    sd = np.nanstd(values, ddof=0)
    if not np.isfinite(sd) or sd <= 0:
        return np.full(len(values), np.nan)
    return (values - np.nanmean(values)) / sd


def batch_logit_betas(x: np.ndarray, y: np.ndarray, max_iter: int = 60) -> np.ndarray:
    beta = np.zeros((y.shape[0], x.shape[1]), dtype=float)
    for _ in range(max_iter):
        eta = beta @ x.T
        mu = expit(eta)
        weights = np.clip(mu * (1.0 - mu), 1e-9, None)
        score = (y - mu) @ x
        hessian = np.einsum("ni,mn,nj->mij", x, weights, x, optimize=True)
        delta = np.vstack([
            np.linalg.lstsq(hessian[index], score[index], rcond=None)[0]
            for index in range(len(hessian))
        ])
        beta += delta
        if np.nanmax(np.abs(delta)) < 1e-9:
            break
    return beta


def design(data: pd.DataFrame, metric: str) -> tuple[np.ndarray, np.ndarray, list[str]]:
    y = (data["spatial_scale"] == "among_population").astype(int).to_numpy()
    columns = ["intercept", metric, *COVARIATES]
    arrays = [np.ones(len(data)), zscore(data[metric])]
    for column in COVARIATES:
        arrays.append(zscore(np.log1p(data[column].clip(lower=0))))
    return np.column_stack(arrays), y, columns


def fit_metric(data: pd.DataFrame, metric: str, permutations: int, seed: int):
    x, y, names = design(data, metric)
    finite = np.isfinite(x).all(axis=1)
    x = x[finite]
    y = y[finite]
    subset = data.loc[finite].copy()
    fit = sm.GLM(y, x, family=sm.families.Binomial()).fit(
        cov_type="cluster",
        cov_kwds={"groups": subset["family"].astype(str).to_numpy()},
    )
    observed = float(fit.params[1])
    rng = np.random.default_rng(seed)
    indices = np.vstack([rng.permutation(len(y)) for _ in range(permutations)])
    permuted_y = y[indices]
    permutation_betas = batch_logit_betas(x, permuted_y)[:, 1]
    valid = permutation_betas[np.isfinite(permutation_betas)]
    permutation_p = float((1 + np.sum(np.abs(valid) >= abs(observed))) / (1 + len(valid)))

    loo: list[dict[str, object]] = []
    loo_betas: list[float] = []
    for family in sorted(subset["family"].astype(str).unique()):
        keep = subset["family"].astype(str).to_numpy() != family
        if keep.sum() < x.shape[1] + 2 or len(np.unique(y[keep])) < 2:
            continue
        loo_fit = sm.GLM(y[keep], x[keep], family=sm.families.Binomial()).fit()
        beta = float(loo_fit.params[1])
        loo_betas.append(beta)
        loo.append({
            "metric": metric,
            "omitted_family": family,
            "estimate": beta,
            "odds_ratio": float(np.exp(beta)),
        })
    loo_array = np.asarray(loo_betas, dtype=float)
    se = float(fit.bse[1])
    row = {
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
        "loo_same_direction_fraction": float(np.mean(np.sign(loo_array) == np.sign(observed))) if len(loo_array) else np.nan,
        "loo_min_odds_ratio": float(np.exp(loo_array.min())) if len(loo_array) else np.nan,
        "loo_max_odds_ratio": float(np.exp(loo_array.max())) if len(loo_array) else np.nan,
        "model_formula": "among ~ climate_component_z + log_geographic_hull_z + log_geographic_grid_cells_z + log_occurrence_records_z + log_supporting_sources_z",
        "classification_source": "exploratory_automated_preserved",
    }
    coefficients = []
    for index, name in enumerate(names):
        coefficients.append({
            "metric": metric,
            "term": name,
            "estimate": float(fit.params[index]),
            "std_error_family_clustered": float(fit.bse[index]),
            "p_value_family_clustered": float(fit.pvalues[index]),
            "odds_ratio": float(np.exp(fit.params[index])),
        })
    return row, loo, coefficients


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
        data["family"] = data["family_class"]
    elif "family" not in data:
        raise ValueError("Family column missing")
    if "classification_source" not in data:
        raise ValueError("classification_source missing")
    if not (data["classification_source"] == "exploratory_automated").all():
        raise ValueError("Automated classification provenance was altered before ecological modelling")
    if "n_distinct_supporting_sources" not in data:
        data["n_distinct_supporting_sources"] = data.get("n_supporting_p1_records", 1)
    data = data.loc[
        data["spatial_scale"].isin(["within_population", "among_population"])
        & (data["metric_status"] == "complete")
        & (data["n_climate_cells"] >= args.min_cells)
    ].copy()

    required = {"canonical_name", "family", "spatial_scale", "classification_source", *COMPONENTS, *COVARIATES}
    missing = sorted(required - set(data.columns))
    if missing:
        raise ValueError(f"Missing adjusted-model columns: {missing}")
    if len(data) < 20 or data.spatial_scale.nunique() < 2:
        raise ValueError("Insufficient complete species for adjusted two-class model")

    rows = []
    loo_rows = []
    coefficient_rows = []
    for offset, metric in enumerate(COMPONENTS):
        row, loo, coefficients = fit_metric(data, metric, args.permutations, args.seed + offset)
        rows.append(row)
        loo_rows.extend(loo)
        coefficient_rows.extend(coefficients)

    result = pd.DataFrame(rows)
    result["wald_p_holm_eight_components"] = multipletests(
        result["wald_p_value_family_clustered"], method="holm"
    )[1]
    result["permutation_p_holm_eight_components"] = multipletests(
        result["permutation_p_two_sided"], method="holm"
    )[1]

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    data.to_csv(outdir / "range_adjusted_model_dataset.csv", index=False)
    result.to_csv(outdir / "range_adjusted_climate_component_models.csv", index=False)
    pd.DataFrame(loo_rows).to_csv(outdir / "range_adjusted_leave_one_family_out.csv", index=False)
    pd.DataFrame(coefficient_rows).to_csv(outdir / "range_adjusted_all_coefficients.csv", index=False)
    manifest = {
        "status": "complete",
        "n_species": int(len(data)),
        "n_within": int((data.spatial_scale == "within_population").sum()),
        "n_among": int((data.spatial_scale == "among_population").sum()),
        "n_families": int(data.family.nunique()),
        "classification_source": "exploratory_automated",
        "focal_components": COMPONENTS,
        "range_and_effort_covariates": COVARIATES,
        "results": result.to_dict("records"),
        "interpretation_guard": "Focal effects estimate climate breadth conditional on sampled geographic extent, occurrence effort and literature support. Residual confounding and automated-label error remain possible; associations are not causal.",
    }
    (outdir / "range_adjusted_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

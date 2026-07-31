#!/usr/bin/env python3
"""Decompose the manuscript moisture-breadth metric into BIO12/14/15/17.

Each component is analysed separately with the same response, occurrence-effort
covariate, family-clustered covariance, common label-permutation schedule and
leave-one-family-out procedure. The existing arithmetic-mean moisture breadth is
retained only as a comparator. No component is combined into a new index.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels
import statsmodels.api as sm
from scipy.special import expit
from statsmodels.stats.multitest import multipletests

COMPONENTS = ["bio12_q95q05", "bio14_q95q05", "bio15_q95q05", "bio17_q95q05"]
COMPARATOR = "moisture_breadth"
LABELS = {
    "bio12_q95q05": "BIO12 annual precipitation breadth",
    "bio14_q95q05": "BIO14 driest-month precipitation breadth",
    "bio15_q95q05": "BIO15 precipitation-seasonality breadth",
    "bio17_q95q05": "BIO17 driest-quarter precipitation breadth",
    "moisture_breadth": "Legacy arithmetic-mean moisture breadth",
}


def zscore(values: pd.Series | np.ndarray) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    sd = np.nanstd(x, ddof=0)
    if not np.isfinite(sd) or sd <= 0:
        return np.full_like(x, np.nan)
    return (x - np.nanmean(x)) / sd


def batch_logit_betas(x: np.ndarray, y: np.ndarray, max_iter: int = 60, tol: float = 1e-10) -> np.ndarray:
    """Fit many unclustered logistic models with a shared design matrix."""
    beta = np.zeros((y.shape[0], x.shape[1]), dtype=float)
    for _ in range(max_iter):
        eta = beta @ x.T
        mu = expit(eta)
        weights = np.clip(mu * (1 - mu), 1e-9, None)
        score = (y - mu) @ x
        hessian = np.einsum("ni,mn,nj->mij", x, weights, x, optimize=True)
        try:
            delta = np.linalg.solve(hessian, score[..., None])[..., 0]
        except np.linalg.LinAlgError:
            delta = np.vstack([
                np.linalg.lstsq(hessian[i], score[i], rcond=None)[0]
                for i in range(len(hessian))
            ])
        beta_new = beta + delta
        if np.nanmax(np.abs(delta)) < tol:
            beta = beta_new
            break
        beta = beta_new
    return beta


def merge_and_validate(model_path: str, metrics_path: str, dataset_label: str) -> pd.DataFrame:
    model = pd.read_csv(model_path)
    metrics = pd.read_csv(metrics_path)
    required_model = {
        "canonical_name", "family", "spatial_scale", "classification_source",
        "n_climate_cells", COMPARATOR,
    }
    required_metrics = {"canonical_name", "n_climate_cells", COMPARATOR, *COMPONENTS}
    missing_model = sorted(required_model - set(model.columns))
    missing_metrics = sorted(required_metrics - set(metrics.columns))
    if missing_model or missing_metrics:
        raise ValueError({"missing_model": missing_model, "missing_metrics": missing_metrics})

    keep = ["canonical_name", "n_climate_cells", COMPARATOR, *COMPONENTS]
    merged = model.merge(metrics[keep], on="canonical_name", suffixes=("_frozen", "_recomputed"), validate="one_to_one")
    if len(merged) != len(model):
        raise ValueError(f"{dataset_label}: model/metrics merge lost rows")

    n_diff = np.abs(merged["n_climate_cells_frozen"] - merged["n_climate_cells_recomputed"])
    moisture_diff = np.abs(merged[f"{COMPARATOR}_frozen"] - merged[f"{COMPARATOR}_recomputed"])
    if float(n_diff.max()) != 0:
        raise ValueError(f"{dataset_label}: climate-cell counts do not match frozen model data")
    if float(moisture_diff.max()) > 1e-6:
        raise ValueError(
            f"{dataset_label}: recomputed moisture breadth differs from frozen data; max={moisture_diff.max()}"
        )

    merged = merged.rename(columns={
        "n_climate_cells_frozen": "n_climate_cells",
        f"{COMPARATOR}_frozen": COMPARATOR,
    })
    merged = merged.drop(columns=["n_climate_cells_recomputed", f"{COMPARATOR}_recomputed"])
    merged["occurrence_dataset"] = dataset_label
    return merged


def fit_observed(data: pd.DataFrame, metric: str):
    y = (data["spatial_scale"] == "among_population").astype(int).to_numpy()
    metric_z = zscore(data[metric])
    effort_z = zscore(np.log1p(data["n_climate_cells"]))
    x = np.column_stack([np.ones(len(data)), metric_z, effort_z])
    fit = sm.GLM(y, x, family=sm.families.Binomial()).fit(
        cov_type="cluster", cov_kwds={"groups": data["family"].to_numpy()}
    )
    beta_check = batch_logit_betas(x, y[None, :])[0]
    if not np.allclose(beta_check, fit.params, atol=1e-7, rtol=1e-7):
        raise RuntimeError({"metric": metric, "statsmodels": fit.params.tolist(), "batch": beta_check.tolist()})
    return fit, x, y


def analyse_dataset(data: pd.DataFrame, dataset_label: str, permutations: int, seed: int):
    data = data.loc[
        (data["classification_source"] == "baseline_unambiguous")
        & data["spatial_scale"].isin(["within_population", "among_population"])
        & (data["n_climate_cells"] >= 20)
    ].copy()
    if len(data) != 34:
        raise ValueError(f"{dataset_label}: expected 34 baseline species, found {len(data)}")

    rng = np.random.default_rng(seed)
    permutation_indices = np.vstack([rng.permutation(len(data)) for _ in range(permutations)])
    base_labels = data["spatial_scale"].to_numpy()
    permuted_y = (base_labels[permutation_indices] == "among_population").astype(float)

    rows: list[dict] = []
    loo_rows: list[dict] = []
    for metric in [*COMPONENTS, COMPARATOR]:
        fit, x, y = fit_observed(data, metric)
        observed = float(fit.params[1])
        se = float(fit.bse[1])
        permutation_betas = batch_logit_betas(x, permuted_y)[:, 1]
        valid = permutation_betas[np.isfinite(permutation_betas)]
        permutation_p = float((1 + np.sum(np.abs(valid) >= abs(observed))) / (1 + len(valid)))

        loo_betas: list[float] = []
        for family in sorted(data["family"].astype(str).unique()):
            keep = data["family"].astype(str).to_numpy() != family
            loo_fit = sm.GLM(y[keep], x[keep], family=sm.families.Binomial()).fit()
            beta = float(loo_fit.params[1])
            loo_betas.append(beta)
            loo_rows.append({
                "occurrence_dataset": dataset_label,
                "metric": metric,
                "metric_label": LABELS[metric],
                "omitted_family": family,
                "estimate": beta,
                "odds_ratio": float(np.exp(beta)),
            })
        loo_array = np.asarray(loo_betas)
        low_beta = observed - 1.96 * se
        high_beta = observed + 1.96 * se
        rows.append({
            "occurrence_dataset": dataset_label,
            "metric": metric,
            "metric_label": LABELS[metric],
            "is_component": metric in COMPONENTS,
            "n_species": int(len(data)),
            "n_families": int(data["family"].nunique()),
            "n_within": int((data["spatial_scale"] == "within_population").sum()),
            "n_among": int((data["spatial_scale"] == "among_population").sum()),
            "estimate": observed,
            "std_error_clustered": se,
            "odds_ratio": float(np.exp(observed)),
            "odds_ratio_ci_low": float(np.exp(low_beta)),
            "odds_ratio_ci_high": float(np.exp(high_beta)),
            "wald_p_value_clustered": float(fit.pvalues[1]),
            "permutation_p_two_sided": permutation_p,
            "permutations_requested": int(permutations),
            "permutations_valid": int(len(valid)),
            "loo_min_odds_ratio": float(np.exp(loo_array.min())),
            "loo_max_odds_ratio": float(np.exp(loo_array.max())),
            "loo_same_direction_fraction": float(np.mean(np.sign(loo_array) == np.sign(observed))),
            "metric_raw_mean": float(data[metric].mean()),
            "metric_raw_sd": float(data[metric].std(ddof=0)),
            "model_formula": "among ~ component_z + effort_z",
            "covariance": "family-clustered sandwich",
            "permutation_schedule": "same label permutations shared across all five metrics within occurrence dataset",
            "seed": int(seed),
        })

    result = pd.DataFrame(rows)
    component_mask = result["is_component"]
    for p_column, adjusted_column in [
        ("wald_p_value_clustered", "wald_p_holm_four_components"),
        ("permutation_p_two_sided", "permutation_p_holm_four_components"),
    ]:
        adjusted = multipletests(result.loc[component_mask, p_column], method="holm")[1]
        result[adjusted_column] = np.nan
        result.loc[component_mask, adjusted_column] = adjusted

    metric_columns = [*COMPONENTS, COMPARATOR]
    correlations = data[metric_columns].corr(method="pearson")
    correlation_rows = []
    for left in metric_columns:
        for right in metric_columns:
            correlation_rows.append({
                "occurrence_dataset": dataset_label,
                "metric_1": left,
                "metric_2": right,
                "pearson_r": float(correlations.loc[left, right]),
                "r_squared": float(correlations.loc[left, right] ** 2),
            })
    return data, result, pd.DataFrame(loo_rows), pd.DataFrame(correlation_rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--primary-model", required=True)
    parser.add_argument("--primary-metrics", required=True)
    parser.add_argument("--paginated-model", required=True)
    parser.add_argument("--paginated-metrics", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--permutations", type=int, default=9999)
    parser.add_argument("--seed", type=int, default=20260801)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    primary = merge_and_validate(args.primary_model, args.primary_metrics, "primary")
    paginated = merge_and_validate(args.paginated_model, args.paginated_metrics, "paginated_quality_filtered")

    outputs = []
    model_tables = []
    loo_tables = []
    correlation_tables = []
    for offset, (label, frame) in enumerate([
        ("primary", primary),
        ("paginated_quality_filtered", paginated),
    ]):
        model_data, models, loo, correlations = analyse_dataset(
            frame, label, args.permutations, args.seed + offset
        )
        model_data.to_csv(outdir / f"precipitation_component_model_dataset_{label}.csv", index=False)
        model_tables.append(models)
        loo_tables.append(loo)
        correlation_tables.append(correlations)
        outputs.append({
            "occurrence_dataset": label,
            "n_species": int(len(model_data)),
            "composite_component_correlations": {
                component: float(model_data[[component, COMPARATOR]].corr().iloc[0, 1])
                for component in COMPONENTS
            },
        })

    models = pd.concat(model_tables, ignore_index=True)
    loo = pd.concat(loo_tables, ignore_index=True)
    correlations = pd.concat(correlation_tables, ignore_index=True)
    models.to_csv(outdir / "precipitation_component_breadth_models.csv", index=False)
    loo.to_csv(outdir / "precipitation_component_breadth_leave_one_family_out.csv", index=False)
    correlations.to_csv(outdir / "precipitation_component_breadth_correlations.csv", index=False)

    manifest = {
        "status": "complete",
        "components": COMPONENTS,
        "legacy_comparator": COMPARATOR,
        "component_definitions": {
            "bio12_q95q05": "95th minus 5th percentile of annual precipitation across occupied climate cells (mm)",
            "bio14_q95q05": "95th minus 5th percentile of precipitation of the driest month (mm)",
            "bio15_q95q05": "95th minus 5th percentile of precipitation seasonality (coefficient of variation; not mm)",
            "bio17_q95q05": "95th minus 5th percentile of precipitation of the driest quarter (mm)",
        },
        "legacy_metric_warning": (
            "The arithmetic-mean moisture breadth combines three precipitation amounts in mm with BIO15 precipitation seasonality, a coefficient of variation. Raw-scale averaging can therefore be dominated by BIO12 and should not be interpreted as a unit-homogeneous moisture index."
        ),
        "model_formula": "among ~ separately_standardised_component + standardised_log1p_n_climate_cells",
        "multiple_testing": "Holm adjustment across the four component tests within each occurrence dataset",
        "permutations": int(args.permutations),
        "seed_primary": int(args.seed),
        "seed_paginated": int(args.seed + 1),
        "statsmodels_version": statsmodels.__version__,
        "datasets": outputs,
        "results": models.to_dict("records"),
        "interpretation_guard": (
            "These models identify which occupied-climate precipitation breadth variable carries the association. They do not establish a causal moisture mechanism, morph-specific tolerance or local adaptation. Separate component models are correlated tests and should be interpreted with multiplicity and collinearity in view."
        ),
    }
    (outdir / "precipitation_component_breadth_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

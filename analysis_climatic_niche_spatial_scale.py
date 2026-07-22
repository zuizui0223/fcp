#!/usr/bin/env python3
"""Compare occupied climatic niches across documented spatial organizations of flower-colour variation."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels
import statsmodels.api as sm

METRICS = [
    "pca_dispersion",
    "climatic_heterogeneity",
    "pca_hull_area",
    "temperature_breadth",
    "moisture_breadth",
]
MODEL_FORMULA = "among ~ metric_z + effort_z"
CI_METHOD = "Wald 95% CI from family-clustered sandwich standard errors"


def zscore(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    sd = values.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(np.nan, index=series.index)
    return (values - values.mean()) / sd


def fit_one(data: pd.DataFrame, metric: str, min_cells: int) -> dict:
    d = data.loc[data["n_climate_cells"] >= min_cells].copy()
    d["among"] = (d["spatial_scale"] == "among_population").astype(int)
    d["metric_z"] = zscore(d[metric])
    d["effort_z"] = zscore(np.log1p(d["n_climate_cells"]))
    d = d.dropna(subset=["among", "metric_z", "effort_z", "family"])
    result = {
        "min_cells": min_cells,
        "metric": metric,
        "model_formula": MODEL_FORMULA,
        "estimator": "statsmodels GLM Binomial(logit)",
        "covariance": "family-clustered sandwich",
        "ci_method": CI_METHOD,
        "statsmodels_version": statsmodels.__version__,
        "n_species": int(len(d)),
        "n_families": int(d["family"].nunique()),
        "n_within": int((d["among"] == 0).sum()),
        "n_among": int((d["among"] == 1).sum()),
        "status": "insufficient",
        "analysis_reason": "",
    }
    if len(d) < 20 or d["among"].nunique() < 2 or d["family"].nunique() < 2:
        result["analysis_reason"] = "fewer than 20 species, fewer than two families, or only one response class"
        return result
    try:
        X = sm.add_constant(d[["metric_z", "effort_z"]], has_constant="add")
        model = sm.GLM(d["among"], X, family=sm.families.Binomial())
        fit = model.fit(cov_type="cluster", cov_kwds={"groups": d["family"]})
        beta = float(fit.params["metric_z"])
        se = float(fit.bse["metric_z"])
        ci_low_beta = beta - 1.96 * se
        ci_high_beta = beta + 1.96 * se
        fit_history = getattr(fit, "fit_history", {})
        result.update({
            "status": "complete",
            "estimate": beta,
            "std_error_clustered": se,
            "estimate_ci_low": ci_low_beta,
            "estimate_ci_high": ci_high_beta,
            "odds_ratio": float(np.exp(beta)),
            "odds_ratio_ci_low": float(np.exp(ci_low_beta)),
            "odds_ratio_ci_high": float(np.exp(ci_high_beta)),
            "p_value_clustered": float(fit.pvalues["metric_z"]),
            "effort_estimate": float(fit.params["effort_z"]),
            "effort_odds_ratio": float(np.exp(fit.params["effort_z"])),
            "converged": bool(getattr(fit, "converged", False)),
            "iterations": int(fit_history.get("iteration", -1)),
            "predicted_probability_min": float(np.min(fit.fittedvalues)),
            "predicted_probability_max": float(np.max(fit.fittedvalues)),
        })
    except Exception as exc:
        result.update({"status": "failed", "analysis_reason": f"{type(exc).__name__}: {exc}"})
    return result


def resolve_family_column(dataset: pd.DataFrame) -> pd.Series:
    """Resolve family robustly regardless of whether merge suffixes were applied."""
    candidates = [c for c in ("family_class", "family", "family_metric") if c in dataset.columns]
    if not candidates:
        raise ValueError("No family column remained after merging classification and niche metrics")
    family = dataset[candidates[0]].copy()
    for column in candidates[1:]:
        family = family.fillna(dataset[column])
    return family


def resolve_scale_column(classification: pd.DataFrame) -> str:
    """Return the integrated spatial-scale column emitted by supported classifiers."""
    for column in ("enriched_scale", "current_scale", "spatial_scale"):
        if column in classification.columns:
            return column
    raise ValueError(
        "Classification input has no supported spatial-scale column; expected one of "
        "enriched_scale, current_scale, or spatial_scale"
    )


def coalesce_text(frame: pd.DataFrame, columns: tuple[str, ...]) -> pd.Series:
    result = pd.Series("", index=frame.index, dtype="object")
    for column in columns:
        if column in frame.columns:
            values = frame[column].fillna("").astype(str).str.strip()
            result = result.mask(result.eq(""), values)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--classification", required=True)
    parser.add_argument("--metrics", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    classification = pd.read_csv(args.classification)
    metrics = pd.read_csv(args.metrics)

    scale_column = resolve_scale_column(classification)
    required_class = {"canonical_name", "family", scale_column}
    required_metrics = {"canonical_name", "role", "n_climate_cells", *METRICS}
    missing = sorted((required_class - set(classification.columns)) | (required_metrics - set(metrics.columns)))
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    class_columns = ["canonical_name", "family", scale_column]
    passthrough_columns = [
        c
        for c in (
            "classification_source",
            "baseline_scale",
            "best_title",
            "best_doi",
            "best_openalex_id",
            "best_match_evidence",
            "review_reason",
            "review_status",
            "evidence_source",
            "source_id",
            "decision_note",
        )
        if c in classification.columns
    ]
    class_columns.extend(passthrough_columns)
    classified = classification.loc[
        classification[scale_column].isin(["within_population", "among_population"]),
        class_columns,
    ].drop_duplicates("canonical_name")
    classified = classified.rename(columns={scale_column: "spatial_scale"})

    focal_metrics = metrics.loc[metrics["role"] == "focal"].drop_duplicates("canonical_name")
    dataset = classified.merge(focal_metrics, on="canonical_name", how="inner", suffixes=("_class", "_metric"))
    dataset["family"] = resolve_family_column(dataset)
    dataset["evidence_source"] = coalesce_text(dataset, ("evidence_source", "best_title"))
    dataset["source_id"] = coalesce_text(dataset, ("source_id", "best_doi", "best_openalex_id"))
    dataset["decision_note"] = coalesce_text(dataset, ("decision_note", "best_match_evidence", "review_reason"))
    dataset.to_csv(outdir / "climatic_niche_spatial_scale_dataset.csv", index=False)

    rows = [fit_one(dataset, metric, threshold) for threshold in (10, 20, 30, 50) for metric in METRICS]
    results = pd.DataFrame(rows)
    results.to_csv(outdir / "climatic_niche_spatial_scale_models.csv", index=False)
    results.to_csv(outdir / "climatic_niche_metric_by_threshold_table.csv", index=False)

    primary = results.loc[(results["min_cells"] == 20) & (results["status"] == "complete")].to_dict("records")
    manifest = {
        "classification_scale_column": scale_column,
        "classified_species": int(len(classified)),
        "classified_species_with_metrics": int(len(dataset)),
        "within_with_metrics": int((dataset["spatial_scale"] == "within_population").sum()),
        "among_with_metrics": int((dataset["spatial_scale"] == "among_population").sum()),
        "evidence_columns_preserved": passthrough_columns,
        "normalized_evidence_columns": ["evidence_source", "source_id", "decision_note"],
        "model_formula": MODEL_FORMULA,
        "estimator": "statsmodels GLM Binomial(logit)",
        "covariance": "family-clustered sandwich",
        "confidence_interval_method": CI_METHOD,
        "statsmodels_version": statsmodels.__version__,
        "metrics_evaluated": METRICS,
        "thresholds_evaluated": [10, 20, 30, 50],
        "specifications": int(len(results)),
        "complete_specifications": int((results["status"] == "complete").sum()),
        "primary_results": primary,
        "interpretation_rule": "Odds ratios above one indicate that broader occupied climatic niche is associated with among-population rather than within-population colour variation.",
        "semantic_guard": "This compares evidence-classified spatial organizations among documented cases; mixed and unclear cases are excluded, occupied niche is not physiological tolerance, multiple metrics were evaluated, and associations are not causal.",
    }
    (outdir / "climatic_niche_spatial_scale_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

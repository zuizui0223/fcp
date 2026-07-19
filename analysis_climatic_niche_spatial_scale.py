#!/usr/bin/env python3
"""Test whether occupied climatic niche breadth differs between within- and among-population FCP."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

METRICS = [
    "pca_dispersion",
    "climatic_heterogeneity",
    "pca_hull_area",
    "temperature_breadth",
    "moisture_breadth",
]


def zscore(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    sd = values.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(np.nan, index=series.index)
    return (values - values.mean()) / sd


def fit_one(data: pd.DataFrame, metric: str, min_cells: int) -> dict:
    d = data.loc[data["n_climate_cells"] >= min_cells].copy()
    d["among"] = (d["current_scale"] == "among_population").astype(int)
    d["metric_z"] = zscore(d[metric])
    d["effort_z"] = zscore(np.log1p(d["n_climate_cells"]))
    d = d.dropna(subset=["among", "metric_z", "effort_z", "family"])
    result = {
        "min_cells": min_cells,
        "metric": metric,
        "n_species": int(len(d)),
        "n_within": int((d["among"] == 0).sum()),
        "n_among": int((d["among"] == 1).sum()),
        "status": "insufficient",
    }
    if len(d) < 20 or d["among"].nunique() < 2 or d["family"].nunique() < 2:
        return result
    try:
        X = sm.add_constant(d[["metric_z", "effort_z"]], has_constant="add")
        model = sm.GLM(d["among"], X, family=sm.families.Binomial())
        fit = model.fit(cov_type="cluster", cov_kwds={"groups": d["family"]})
        beta = float(fit.params["metric_z"])
        se = float(fit.bse["metric_z"])
        result.update({
            "status": "complete",
            "estimate": beta,
            "std_error": se,
            "odds_ratio": float(np.exp(beta)),
            "ci_low": float(np.exp(beta - 1.96 * se)),
            "ci_high": float(np.exp(beta + 1.96 * se)),
            "p_value": float(fit.pvalues["metric_z"]),
            "effort_odds_ratio": float(np.exp(fit.params["effort_z"])),
        })
    except Exception as exc:
        result.update({"status": "failed", "error": str(exc)})
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

    required_class = {"canonical_name", "family", "current_scale"}
    required_metrics = {"canonical_name", "n_climate_cells", *METRICS}
    missing = sorted((required_class - set(classification.columns)) | (required_metrics - set(metrics.columns)))
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    classified = classification.loc[
        classification["current_scale"].isin(["within_population", "among_population"]),
        ["canonical_name", "family", "current_scale", "classification_source"],
    ].drop_duplicates("canonical_name")
    focal_metrics = metrics.loc[metrics["role"] == "focal"].drop_duplicates("canonical_name")
    dataset = classified.merge(focal_metrics, on="canonical_name", how="inner", suffixes=("_class", "_metric"))
    dataset["family"] = dataset["family_class"].fillna(dataset.get("family_metric"))
    dataset.to_csv(outdir / "climatic_niche_spatial_scale_dataset.csv", index=False)

    rows = [fit_one(dataset, metric, threshold) for threshold in (10, 20, 30, 50) for metric in METRICS]
    results = pd.DataFrame(rows)
    results.to_csv(outdir / "climatic_niche_spatial_scale_models.csv", index=False)

    primary = results.loc[(results["min_cells"] == 20) & (results["status"] == "complete")].to_dict("records")
    manifest = {
        "classified_species": int(len(classified)),
        "classified_species_with_metrics": int(len(dataset)),
        "within_with_metrics": int((dataset["current_scale"] == "within_population").sum()),
        "among_with_metrics": int((dataset["current_scale"] == "among_population").sum()),
        "specifications": int(len(results)),
        "complete_specifications": int((results["status"] == "complete").sum()),
        "primary_results": primary,
        "interpretation_rule": "Odds ratios above one indicate that broader occupied climatic niche is associated with among-population rather than within-population colour variation.",
        "semantic_guard": "This compares evidence-classified spatial scales among verified FCP species; mixed and unclear cases are excluded, occupied niche is not physiological tolerance, and associations are not causal.",
    }
    (outdir / "climatic_niche_spatial_scale_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

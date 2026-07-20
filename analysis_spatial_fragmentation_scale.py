#!/usr/bin/env python3
"""Test whether observed GBIF distribution fragmentation distinguishes FCP spatial scales.

The analysis is deliberately parsimonious because the high-confidence spatial-scale sample is
small. Models retain occupied moisture-niche breadth and sampling effort, then add one observed
spatial-connectivity metric at a time. Connected components describe the sampled GBIF point cloud;
they are not inferred biological populations.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm


FRAGMENTATION_METRICS = [
    "median_nearest_neighbor_km",
    "spatial_radius_q95_km",
    "connected_components_50km",
    "connected_components_100km",
    "largest_component_fraction_100km",
    "occupied_grid_cells_1deg",
]


def zscore(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    sd = values.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(np.nan, index=series.index, dtype=float)
    return (values - values.mean()) / sd


def fit_model(data: pd.DataFrame, name: str, predictors: list[str]) -> dict:
    needed = ["among", "family", *predictors]
    d = data.dropna(subset=needed).copy()
    result = {
        "model": name,
        "predictors": predictors,
        "n_species": int(len(d)),
        "n_within": int((d["among"] == 0).sum()),
        "n_among": int((d["among"] == 1).sum()),
        "n_families": int(d["family"].nunique()),
        "status": "insufficient",
    }
    if len(d) < 20 or d["among"].nunique() < 2 or d["family"].nunique() < 2:
        return result
    try:
        X = sm.add_constant(d[predictors], has_constant="add")
        model = sm.GLM(d["among"], X, family=sm.families.Binomial())
        fit = model.fit(cov_type="cluster", cov_kwds={"groups": d["family"]})
        terms = {}
        for term in predictors:
            beta = float(fit.params[term])
            se = float(fit.bse[term])
            terms[term] = {
                "estimate": beta,
                "std_error": se,
                "odds_ratio": float(np.exp(beta)),
                "ci_low": float(np.exp(beta - 1.96 * se)),
                "ci_high": float(np.exp(beta + 1.96 * se)),
                "p_value": float(fit.pvalues[term]),
            }
        result.update({
            "status": "complete",
            "aic": float(fit.aic),
            "terms": terms,
        })
    except Exception as exc:
        result.update({"status": "failed", "error": f"{type(exc).__name__}: {exc}"})
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scale-dataset", required=True)
    parser.add_argument("--geography", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--min-records", type=int, default=20)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    scale = pd.read_csv(args.scale_dataset)
    geography = pd.read_csv(args.geography)

    required_scale = {"canonical_name", "family", "spatial_scale", "moisture_breadth", "n_climate_cells"}
    required_geo = {"canonical_name", "gbif_coordinate_records_sampled", *FRAGMENTATION_METRICS}
    missing = sorted((required_scale - set(scale.columns)) | (required_geo - set(geography.columns)))
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    geo_columns = ["canonical_name", "gbif_coordinate_records_sampled", *FRAGMENTATION_METRICS]
    dataset = scale.merge(geography[geo_columns].drop_duplicates("canonical_name"), on="canonical_name", how="inner")
    dataset["among"] = (dataset["spatial_scale"] == "among_population").astype(int)
    dataset = dataset.loc[pd.to_numeric(dataset["gbif_coordinate_records_sampled"], errors="coerce") >= args.min_records].copy()

    dataset["moisture_z"] = zscore(dataset["moisture_breadth"])
    dataset["climate_effort_z"] = zscore(np.log1p(pd.to_numeric(dataset["n_climate_cells"], errors="coerce")))
    dataset["gbif_effort_z"] = zscore(np.log1p(pd.to_numeric(dataset["gbif_coordinate_records_sampled"], errors="coerce")))

    transformed: dict[str, str] = {}
    for metric in FRAGMENTATION_METRICS:
        values = pd.to_numeric(dataset[metric], errors="coerce")
        if metric == "largest_component_fraction_100km":
            transformed_name = f"{metric}_z"
            dataset[transformed_name] = zscore(values)
        else:
            transformed_name = f"log1p_{metric}_z"
            dataset[transformed_name] = zscore(np.log1p(values.clip(lower=0)))
        transformed[metric] = transformed_name

    dataset.to_csv(outdir / "spatial_fragmentation_scale_dataset.csv", index=False)

    base = ["moisture_z", "climate_effort_z", "gbif_effort_z"]
    models = [fit_model(dataset, "baseline_moisture", base)]
    for metric, transformed_name in transformed.items():
        models.append(fit_model(dataset, f"add_{metric}", [*base, transformed_name]))

    integrated = [
        *base,
        transformed["connected_components_100km"],
        transformed["largest_component_fraction_100km"],
    ]
    models.append(fit_model(dataset, "integrated_100km_connectivity", integrated))

    rows = []
    for model in models:
        if model["status"] != "complete":
            rows.append({k: v for k, v in model.items() if k != "terms"})
            continue
        for term, values in model["terms"].items():
            rows.append({
                "model": model["model"],
                "n_species": model["n_species"],
                "n_within": model["n_within"],
                "n_among": model["n_among"],
                "n_families": model["n_families"],
                "aic": model["aic"],
                "status": model["status"],
                "term": term,
                **values,
            })
    pd.DataFrame(rows).to_csv(outdir / "spatial_fragmentation_scale_models.csv", index=False)

    manifest = {
        "min_gbif_coordinate_records": args.min_records,
        "species_available": int(len(dataset)),
        "within_available": int((dataset["among"] == 0).sum()),
        "among_available": int((dataset["among"] == 1).sum()),
        "models_complete": int(sum(m["status"] == "complete" for m in models)),
        "models": models,
        "interpretation_rule": "Odds ratios above one for component counts or nearest-neighbour distance indicate that a more fragmented observed GBIF point cloud is associated with among-population rather than within-population colour variation.",
        "semantic_guard": "GBIF connected components and distance summaries describe spatial connectivity of sampled occurrence points, not biological populations, dispersal barriers, gene flow, or causal mechanisms.",
    }
    (outdir / "spatial_fragmentation_scale_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

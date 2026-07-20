#!/usr/bin/env python3
"""Sensitivity analysis for environmental turnover across spatial clustering choices."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analysis_environmental_turnover_scale import CLIMATE, fit_model, summarize_species, zscore


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--climate-cells", required=True)
    parser.add_argument("--scale-dataset", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--thresholds-km", default="50,100,200")
    parser.add_argument("--min-component-records", default="3,5")
    args = parser.parse_args()

    thresholds = [float(x) for x in args.thresholds_km.split(",") if x.strip()]
    minimums = [int(x) for x in args.min_component_records.split(",") if x.strip()]
    if not thresholds or not minimums:
        raise ValueError("At least one threshold and one minimum component size are required")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    cells = pd.read_csv(args.climate_cells)
    scale = pd.read_csv(args.scale_dataset)

    required_cells = {"canonical_name", "role", "decimalLatitude", "decimalLongitude", *CLIMATE}
    required_scale = {"canonical_name", "family", "spatial_scale", "moisture_breadth", "n_climate_cells"}
    missing = sorted((required_cells - set(cells.columns)) | (required_scale - set(scale.columns)))
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    focal_names = set(scale["canonical_name"].astype(str))
    focal_cells = cells[(cells["role"] == "focal") & cells["canonical_name"].astype(str).isin(focal_names)].copy()

    model_rows: list[dict[str, object]] = []
    species_rows: list[dict[str, object]] = []
    for threshold in thresholds:
        for minimum in minimums:
            summaries = [
                summarize_species(group, threshold, minimum)
                for _, group in focal_cells.groupby("canonical_name", sort=True)
            ]
            turnover = pd.DataFrame(summaries)
            turnover["specification"] = f"{threshold:g}km_min{minimum}"
            species_rows.extend(turnover.to_dict("records"))

            dataset = scale.merge(turnover.drop(columns=["specification"]), on="canonical_name", how="left")
            dataset["among"] = (dataset["spatial_scale"] == "among_population").astype(int)
            dataset["moisture_z"] = zscore(dataset["moisture_breadth"])
            dataset["climate_effort_z"] = zscore(np.log1p(pd.to_numeric(dataset["n_climate_cells"], errors="coerce")))
            complete = dataset["status"] == "complete"
            dataset["turnover_z"] = np.nan
            dataset.loc[complete, "turnover_z"] = zscore(dataset.loc[complete, "turnover_weighted_mean"])
            model = fit_model(dataset, "turnover_weighted_mean")
            turnover_term = (model.get("terms") or {}).get("turnover_z", {})
            moisture_term = (model.get("terms") or {}).get("moisture_z", {})
            model_rows.append({
                "specification": f"{threshold:g}km_min{minimum}",
                "threshold_km": threshold,
                "min_component_records": minimum,
                "species_with_two_components": int((turnover["status"] == "complete").sum()),
                "status": model.get("status"),
                "n_species": model.get("n_species"),
                "n_within": model.get("n_within"),
                "n_among": model.get("n_among"),
                "n_families": model.get("n_families"),
                "aic": model.get("aic"),
                "turnover_odds_ratio": turnover_term.get("odds_ratio"),
                "turnover_ci_low": turnover_term.get("ci_low"),
                "turnover_ci_high": turnover_term.get("ci_high"),
                "turnover_p_value": turnover_term.get("p_value"),
                "moisture_odds_ratio": moisture_term.get("odds_ratio"),
                "moisture_p_value": moisture_term.get("p_value"),
                "error": model.get("error"),
            })

    models = pd.DataFrame(model_rows)
    species = pd.DataFrame(species_rows)
    models.to_csv(outdir / "environmental_turnover_robustness_models.csv", index=False)
    species.to_csv(outdir / "environmental_turnover_robustness_species.csv", index=False)

    complete = models[models["status"] == "complete"].copy()
    manifest = {
        "specifications": int(len(models)),
        "complete_specifications": int(len(complete)),
        "thresholds_km": thresholds,
        "minimum_component_records": minimums,
        "turnover_direction_positive": int((complete["turnover_odds_ratio"] > 1).sum()) if len(complete) else 0,
        "turnover_direction_negative": int((complete["turnover_odds_ratio"] < 1).sum()) if len(complete) else 0,
        "turnover_p_below_0_05": int((complete["turnover_p_value"] < 0.05).sum()) if len(complete) else 0,
        "turnover_or_range": [float(complete["turnover_odds_ratio"].min()), float(complete["turnover_odds_ratio"].max())] if len(complete) else None,
        "semantic_guard": "Sensitivity specifications vary distance thresholds and minimum sampled-component sizes; occurrence clusters remain observational GBIF point-cloud components, not biological populations.",
    }
    (outdir / "environmental_turnover_robustness_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

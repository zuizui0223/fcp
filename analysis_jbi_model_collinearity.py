#!/usr/bin/env python3
"""Explicit multicollinearity diagnostics for the JBI spatial-scale models.

The script diagnoses the models that are actually interpreted in the manuscript:
(1) the 34-species baseline moisture-breadth model, and
(2) each coarse occurrence-cloud sensitivity model. It does not use automated
variable selection; biologically prespecified covariates are retained and range
metrics are assessed one at a time, except for the explicitly labelled integrated
100-km connectivity sensitivity model.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.outliers_influence import variance_inflation_factor


def zscore(x: pd.Series) -> pd.Series:
    values = pd.to_numeric(x, errors="coerce")
    sd = values.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(np.nan, index=x.index, dtype=float)
    return (values - values.mean()) / sd


def diagnose(model_id: str, frame: pd.DataFrame, predictors: list[str]) -> tuple[list[dict], dict]:
    d = frame[predictors].apply(pd.to_numeric, errors="coerce").dropna().copy()
    if len(d) < 3:
        return [], {
            "model_id": model_id,
            "n_species": int(len(d)),
            "n_predictors": len(predictors),
            "condition_number": np.nan,
            "max_vif": np.nan,
            "status": "insufficient",
        }
    design = sm.add_constant(d, has_constant="add")
    matrix = design.to_numpy(dtype=float)
    rows = []
    for i, predictor in enumerate(predictors, start=1):
        vif = float(variance_inflation_factor(matrix, i))
        rows.append({
            "model_id": model_id,
            "predictor": predictor,
            "n_species": int(len(d)),
            "vif": vif,
        })
    condition_number = float(np.linalg.cond(matrix))
    return rows, {
        "model_id": model_id,
        "n_species": int(len(d)),
        "n_predictors": len(predictors),
        "condition_number": condition_number,
        "max_vif": max(r["vif"] for r in rows),
        "status": "complete",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scale-dataset", required=True)
    parser.add_argument("--fragmentation-dataset", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    scale = pd.read_csv(args.scale_dataset)
    fragmentation = pd.read_csv(args.fragmentation_dataset)

    required_scale = {
        "classification_source", "n_climate_cells", "moisture_breadth",
    }
    missing_scale = sorted(required_scale - set(scale.columns))
    if missing_scale:
        raise ValueError(f"Missing scale-dataset columns: {missing_scale}")

    baseline = scale.loc[
        (scale["classification_source"] == "baseline_unambiguous")
        & (pd.to_numeric(scale["n_climate_cells"], errors="coerce") >= 20)
    ].copy()
    baseline["moisture_z"] = zscore(baseline["moisture_breadth"])
    baseline["climate_effort_z"] = zscore(
        np.log1p(pd.to_numeric(baseline["n_climate_cells"], errors="coerce"))
    )

    required_fragmentation = {
        "moisture_z", "climate_effort_z", "gbif_effort_z",
        "log1p_median_nearest_neighbor_km_z",
        "log1p_spatial_extent_q95_km_z",
        "log1p_components_50km_z",
        "log1p_components_100km_z",
        "largest_component_fraction_100km_z",
        "log1p_occupied_grid_cells_1deg_z",
    }
    missing_fragmentation = sorted(required_fragmentation - set(fragmentation.columns))
    if missing_fragmentation:
        raise ValueError(f"Missing fragmentation-dataset columns: {missing_fragmentation}")

    specifications: list[tuple[str, pd.DataFrame, list[str]]] = [
        (
            "baseline34_moisture",
            baseline,
            ["moisture_z", "climate_effort_z"],
        )
    ]
    base = ["moisture_z", "climate_effort_z", "gbif_effort_z"]
    one_at_a_time = {
        "add_median_nearest_neighbor": "log1p_median_nearest_neighbor_km_z",
        "add_spatial_extent_q95": "log1p_spatial_extent_q95_km_z",
        "add_components_50km": "log1p_components_50km_z",
        "add_components_100km": "log1p_components_100km_z",
        "add_largest_component_fraction_100km": "largest_component_fraction_100km_z",
        "add_occupied_grid_cells_1deg": "log1p_occupied_grid_cells_1deg_z",
    }
    specifications.append(("fragmentation_baseline", fragmentation, base))
    for model_id, metric in one_at_a_time.items():
        specifications.append((model_id, fragmentation, [*base, metric]))
    specifications.append((
        "integrated_100km_connectivity",
        fragmentation,
        [*base, "log1p_components_100km_z", "largest_component_fraction_100km_z"],
    ))

    vif_rows: list[dict] = []
    model_rows: list[dict] = []
    for model_id, frame, predictors in specifications:
        rows, summary = diagnose(model_id, frame, predictors)
        vif_rows.extend(rows)
        model_rows.append({**summary, "predictors": " + ".join(predictors)})

    vif_table = pd.DataFrame(vif_rows)
    model_table = pd.DataFrame(model_rows)
    vif_table.to_csv(outdir / "jbi_model_vif_diagnostics.csv", index=False)
    model_table.to_csv(outdir / "jbi_model_condition_diagnostics.csv", index=False)

    baseline_row = model_table.loc[model_table["model_id"] == "baseline34_moisture"].iloc[0]
    sensitivity = model_table.loc[model_table["model_id"] != "baseline34_moisture"]
    manifest = {
        "standardization_note": (
            "All diagnosed predictors are the standardized variables used by their fitted models; "
            "standardization changes scale but does not remove collinearity."
        ),
        "selection_rule": (
            "No VIF-based stepwise selection is used. The focal climate and effort terms are retained "
            "a priori, and coarse range/connectivity metrics are added one at a time except for the "
            "explicit integrated 100-km sensitivity model."
        ),
        "baseline34": {
            "n_species": int(baseline_row["n_species"]),
            "max_vif": float(baseline_row["max_vif"]),
            "condition_number": float(baseline_row["condition_number"]),
        },
        "sensitivity_models": {
            "n_models": int(len(sensitivity)),
            "max_vif_across_models": float(sensitivity["max_vif"].max()),
            "max_condition_number": float(sensitivity["condition_number"].max()),
            "models_with_max_vif_gt_5": sensitivity.loc[
                sensitivity["max_vif"] > 5, "model_id"
            ].astype(str).tolist(),
        },
    }
    (outdir / "jbi_model_collinearity_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

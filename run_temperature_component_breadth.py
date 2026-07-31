#!/usr/bin/env python3
"""Run the temperature-component analysis with model/metric merge parity.

The frozen S8/S9 model tables contain classifications and effort, while the legacy
temperature breadth and its four components are sourced from the matching climate
metric artifacts. This wrapper keeps that provenance explicit.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import analysis_temperature_component_breadth as analysis


def merge_and_validate(model_path: str, metrics_path: str, dataset_label: str) -> pd.DataFrame:
    model = pd.read_csv(model_path).drop(columns=[analysis.COMPARATOR], errors="ignore")
    metrics = pd.read_csv(metrics_path)
    required_model = {
        "canonical_name", "family", "spatial_scale", "classification_source",
        "n_climate_cells",
    }
    required_metrics = {
        "canonical_name", "n_climate_cells", analysis.COMPARATOR,
        *analysis.COMPONENTS,
    }
    missing_model = sorted(required_model - set(model.columns))
    missing_metrics = sorted(required_metrics - set(metrics.columns))
    if missing_model or missing_metrics:
        raise ValueError({"missing_model": missing_model, "missing_metrics": missing_metrics})

    keep = ["canonical_name", "n_climate_cells", analysis.COMPARATOR, *analysis.COMPONENTS]
    merged = model.merge(
        metrics[keep], on="canonical_name", suffixes=("_frozen", "_metrics"),
        validate="one_to_one"
    )
    if len(merged) != len(model):
        raise ValueError(f"{dataset_label}: model/metrics merge lost rows")
    n_diff = np.abs(merged["n_climate_cells_frozen"] - merged["n_climate_cells_metrics"])
    if float(n_diff.max()) != 0:
        raise ValueError(f"{dataset_label}: climate-cell counts do not match frozen model data")
    merged = merged.rename(columns={"n_climate_cells_frozen": "n_climate_cells"})
    merged = merged.drop(columns=["n_climate_cells_metrics"])
    merged["occurrence_dataset"] = dataset_label
    return merged


analysis.merge_and_validate = merge_and_validate
analysis.main()

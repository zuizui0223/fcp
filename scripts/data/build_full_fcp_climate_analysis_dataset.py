#!/usr/bin/env python3
"""Join the independently constructed all-species FCP universe to climate/geography.

No fixed species count is assumed. Climate eligibility is applied only after species
membership and C/S evidence states have been constructed from literature.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

METRICS = [
    "temperature_breadth",
    "moisture_breadth",
    "climatic_heterogeneity",
    "pca_dispersion",
    "pca_hull_area",
]
GEO = [
    "n_geographic_records",
    "geographic_radius_median_km",
    "geographic_radius_95_km",
    "geographic_dispersion_rms_km",
    "geographic_bbox_area_km2_approx",
    "latitudinal_span_deg",
    "minimal_longitude_span_deg",
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--universe", required=True)
    p.add_argument("--metrics", required=True)
    p.add_argument("--geography", required=True)
    p.add_argument("--out-full", required=True)
    p.add_argument("--out-analysis", required=True)
    p.add_argument("--summary-out", required=True)
    p.add_argument("--min-cells", type=int, default=20)
    p.add_argument("--worldclim-resolution", default="10 arc-minute")
    p.add_argument("--worldclim-url", default="")
    args = p.parse_args()

    universe_path = Path(args.universe)
    metrics_path = Path(args.metrics)
    geography_path = Path(args.geography)
    universe = pd.read_csv(universe_path)
    metrics = pd.read_csv(metrics_path)
    geography = pd.read_csv(geography_path)

    for label, df in [("universe", universe), ("climate", metrics), ("geography", geography)]:
        if df.canonical_name.duplicated().any():
            raise SystemExit(f"Duplicate species in {label} table")

    required_universe = {
        "canonical_name", "family", "organization_state",
        "C_local_coexistence_documented", "S_spatial_segregation_documented",
        "n_FCP_eligible_sources",
    }
    missing = required_universe - set(universe.columns)
    if missing:
        raise SystemExit(f"Universe missing columns: {sorted(missing)}")

    keep_metrics = ["canonical_name", "n_climate_cells", "metric_status", *METRICS]
    missing_metrics = set(keep_metrics) - set(metrics.columns)
    if missing_metrics:
        raise SystemExit(f"Climate table missing columns: {sorted(missing_metrics)}")
    missing_geo = {"canonical_name", *GEO} - set(geography.columns)
    if missing_geo:
        raise SystemExit(f"Geography table missing columns: {sorted(missing_geo)}")

    merged = universe.merge(metrics[keep_metrics], how="left", on="canonical_name", validate="one_to_one")
    merged = merged.merge(geography[["canonical_name", *GEO]], how="left", on="canonical_name", validate="one_to_one")
    merged["climate_eligible"] = (
        merged["metric_status"].eq("complete")
        & merged["n_climate_cells"].fillna(0).ge(args.min_cells)
        & merged[METRICS].notna().all(axis=1)
        & merged["geographic_radius_95_km"].notna()
    )
    merged = merged.sort_values("canonical_name").reset_index(drop=True)
    analysis = merged.loc[merged.climate_eligible].copy()

    out_full = Path(args.out_full)
    out_analysis = Path(args.out_analysis)
    out_full.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out_full, index=False)
    analysis.to_csv(out_analysis, index=False)

    states = [
        "local_coexistence_only", "spatial_segregation_only",
        "coexistence_and_segregation", "organization_unresolved",
    ]
    summary = {
        "status": "complete",
        "universe_sha256": sha256(universe_path),
        "climate_metrics_sha256": sha256(metrics_path),
        "geographic_metrics_sha256": sha256(geography_path),
        "full_join_sha256": sha256(out_full),
        "analysis_dataset_sha256": sha256(out_analysis),
        "universe_species": int(len(universe)),
        "climate_metrics_species": int(len(metrics)),
        "geographic_metrics_species": int(len(geography)),
        "analysis_species": int(len(analysis)),
        "excluded_for_climate_or_geography": int(len(universe) - len(analysis)),
        "min_occupied_climate_cells": int(args.min_cells),
        "worldclim_version": "2.1",
        "worldclim_resolution": args.worldclim_resolution,
        "worldclim_bioclim_variables": [1,4,5,6,7,12,14,15,17],
        "worldclim_url": args.worldclim_url,
        "state_counts_universe": {s:int((universe.organization_state==s).sum()) for s in states},
        "state_counts_analysis": {s:int((analysis.organization_state==s).sum()) for s in states},
        "C_positive_analysis_species": int(analysis.C_local_coexistence_documented.sum()),
        "S_positive_analysis_species": int(analysis.S_spatial_segregation_documented.sum()),
        "organization_unresolved_analysis_species": int((analysis.organization_state=="organization_unresolved").sum()),
        "excluded_species": merged.loc[~merged.climate_eligible, ["canonical_name","n_climate_cells","metric_status","geographic_radius_95_km"]].fillna("").to_dict(orient="records"),
        "semantic_guard": (
            "Species membership is independent of C/S positivity and of climate. C=0 or S=0 means not documented, not biological absence. "
            "Geographic extent is retained as a confounder rather than folded into climatic niche breadth."
        ),
    }
    Path(args.summary_out).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))

    if len(universe) <= 34:
        raise SystemExit(f"Universe unexpectedly collapsed to historical scale: {len(universe)}")
    if len(analysis) < 0.5 * len(universe):
        raise SystemExit(f"More than half of universe lost after climate/geography QC: {len(analysis)}/{len(universe)}")


if __name__ == "__main__":
    main()

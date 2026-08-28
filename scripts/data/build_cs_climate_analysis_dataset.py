#!/usr/bin/env python3
"""Join the frozen coexistence/segregation evidence set to rebuilt climate metrics.

The evidence freeze is immutable input. Species failing the >=20 occupied-climate-cell
criterion remain visible in the full joined table but are not forced into the primary
analysis dataset.
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


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--evidence-freeze", required=True)
    p.add_argument("--metrics", required=True)
    p.add_argument("--out-full", required=True)
    p.add_argument("--out-analysis", required=True)
    p.add_argument("--summary-out", required=True)
    p.add_argument("--min-cells", type=int, default=20)
    p.add_argument("--worldclim-resolution", default="10 arc-minute")
    p.add_argument("--worldclim-url", default="")
    args = p.parse_args()

    evidence_path = Path(args.evidence_freeze)
    metrics_path = Path(args.metrics)
    evidence = pd.read_csv(evidence_path)
    metrics = pd.read_csv(metrics_path)

    if evidence.canonical_name.duplicated().any():
        raise SystemExit("Duplicate species in C/S evidence freeze")
    if metrics.canonical_name.duplicated().any():
        raise SystemExit("Duplicate species in climate metrics table")

    required_evidence = {
        "canonical_name",
        "family",
        "organization_state",
        "C_local_coexistence_documented",
        "S_spatial_segregation_documented",
    }
    missing = required_evidence - set(evidence.columns)
    if missing:
        raise SystemExit(f"Evidence freeze missing columns: {sorted(missing)}")

    keep_metrics = [
        "canonical_name",
        "n_climate_cells",
        "metric_status",
        *METRICS,
    ]
    missing_metrics = set(keep_metrics) - set(metrics.columns)
    if missing_metrics:
        raise SystemExit(f"Climate table missing columns: {sorted(missing_metrics)}")

    merged = evidence.merge(metrics[keep_metrics], how="left", on="canonical_name", validate="one_to_one")
    merged["climate_eligible"] = (
        merged["metric_status"].eq("complete")
        & merged["n_climate_cells"].fillna(0).ge(args.min_cells)
        & merged[METRICS].notna().all(axis=1)
    )
    merged = merged.sort_values("canonical_name").reset_index(drop=True)
    analysis = merged.loc[merged.climate_eligible].copy()

    out_full = Path(args.out_full)
    out_analysis = Path(args.out_analysis)
    out_full.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out_full, index=False)
    analysis.to_csv(out_analysis, index=False)

    state_order = [
        "local_coexistence_only",
        "spatial_segregation_only",
        "coexistence_and_segregation",
    ]
    summary = {
        "evidence_freeze_sha256": sha256(evidence_path),
        "climate_metrics_sha256": sha256(metrics_path),
        "full_join_sha256": sha256(out_full),
        "analysis_dataset_sha256": sha256(out_analysis),
        "evidence_species": int(len(evidence)),
        "climate_metrics_species": int(len(metrics)),
        "analysis_species": int(len(analysis)),
        "excluded_for_climate_cells_or_missing_metrics": int(len(evidence) - len(analysis)),
        "min_occupied_climate_cells": int(args.min_cells),
        "worldclim_version": "2.1",
        "worldclim_resolution": args.worldclim_resolution,
        "worldclim_bioclim_variables": [1, 4, 5, 6, 7, 12, 14, 15, 17],
        "worldclim_url": args.worldclim_url,
        "state_counts_evidence_freeze": {
            state: int((evidence.organization_state == state).sum()) for state in state_order
        },
        "state_counts_analysis": {
            state: int((analysis.organization_state == state).sum()) for state in state_order
        },
        "C_positive_analysis_species": int(analysis.C_local_coexistence_documented.sum()),
        "S_positive_analysis_species": int(analysis.S_spatial_segregation_documented.sum()),
        "excluded_species": merged.loc[~merged.climate_eligible, ["canonical_name", "n_climate_cells", "metric_status"]]
        .fillna("")
        .to_dict(orient="records"),
        "semantic_guard": (
            "C and S are documented-evidence outcomes. Climate eligibility is applied after evidence freezing; "
            "species are never reclassified to satisfy the climate-cell threshold."
        ),
    }
    Path(args.summary_out).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))

    if len(evidence) != 34:
        raise SystemExit(f"Expected 34 evidence-freeze species, found {len(evidence)}")
    if len(analysis) < 10:
        raise SystemExit(f"Unexpectedly few climate-eligible species: {len(analysis)}")


if __name__ == "__main__":
    main()

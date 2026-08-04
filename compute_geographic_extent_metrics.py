#!/usr/bin/env python3
"""Compute species geographic-extent covariates from filtered GBIF occurrences.

These variables are intended to separate climatic heterogeneity from the trivial effect
of occupying a geographically large range.  Distances and areas are descriptive
properties of the retained occurrence sample, not formal range maps.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull

EARTH_KM_PER_DEGREE = 111.195


def projected_xy(frame: pd.DataFrame) -> np.ndarray:
    mean_lat = np.deg2rad(frame["decimalLatitude"].mean())
    x = frame["decimalLongitude"].to_numpy(float) * np.cos(mean_lat) * EARTH_KM_PER_DEGREE
    y = frame["decimalLatitude"].to_numpy(float) * EARTH_KM_PER_DEGREE
    return np.column_stack([x, y])


def max_distance_approx(points: np.ndarray) -> float:
    if len(points) < 2:
        return 0.0
    if len(points) > 500:
        indices = np.linspace(0, len(points) - 1, 500, dtype=int)
        points = points[indices]
    delta = points[:, None, :] - points[None, :, :]
    return float(np.sqrt(np.sum(delta * delta, axis=2)).max())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--occurrences", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--qc-out", required=True)
    parser.add_argument("--grid-degrees", type=float, default=1.0)
    args = parser.parse_args()

    data = pd.read_csv(args.occurrences)
    required = {"canonical_name", "decimalLatitude", "decimalLongitude"}
    missing = sorted(required - set(data.columns))
    if missing:
        raise ValueError(f"Missing occurrence columns: {missing}")

    rows: list[dict[str, object]] = []
    for name, group in data.groupby("canonical_name", sort=True):
        group = group.dropna(subset=["decimalLatitude", "decimalLongitude"]).copy()
        points = projected_xy(group)
        hull_area = 0.0
        if len(points) >= 3:
            try:
                hull_area = float(ConvexHull(points).volume)
            except Exception:
                hull_area = 0.0
        grid_lat = np.floor((group["decimalLatitude"].to_numpy(float) + 90.0) / args.grid_degrees)
        grid_lon = np.floor((group["decimalLongitude"].to_numpy(float) + 180.0) / args.grid_degrees)
        occupied_grid = len(set(zip(grid_lat.astype(int), grid_lon.astype(int))))
        rows.append({
            "canonical_name": name,
            "n_occurrence_records": int(len(group)),
            "occupied_geographic_grid_cells": int(occupied_grid),
            "geographic_hull_area_km2": hull_area,
            "maximum_pairwise_distance_km_approx": max_distance_approx(points),
            "latitude_span_degrees": float(group["decimalLatitude"].max() - group["decimalLatitude"].min()),
            "longitude_span_degrees": float(group["decimalLongitude"].max() - group["decimalLongitude"].min()),
        })

    output = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.out, index=False)
    qc = {
        "status": "complete",
        "species": int(output.canonical_name.nunique()),
        "grid_degrees": args.grid_degrees,
        "covariates": [
            "n_occurrence_records",
            "occupied_geographic_grid_cells",
            "geographic_hull_area_km2",
            "maximum_pairwise_distance_km_approx",
        ],
        "interpretation_guard": "Extent metrics describe the retained GBIF sample and are included to reduce range-size circularity; they are not authoritative range maps.",
    }
    Path(args.qc_out).write_text(json.dumps(qc, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(qc, indent=2))


if __name__ == "__main__":
    main()

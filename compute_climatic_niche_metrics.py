#!/usr/bin/env python3
"""Extract WorldClim values and compute comparable occupied climatic niche metrics."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from scipy.spatial import ConvexHull
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

SELECT = [1, 4, 5, 6, 7, 12, 14, 15, 17]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--occurrences", required=True)
    parser.add_argument("--worldclim-dir", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--qc-out", required=True)
    parser.add_argument("--cell-out")
    parser.add_argument("--min-cells", type=int, default=10)
    args = parser.parse_args()

    data = pd.read_csv(args.occurrences)
    coordinates = list(zip(data.decimalLongitude, data.decimalLatitude))
    extracted = []
    for bio in SELECT:
        candidates = list(Path(args.worldclim_dir).glob(f"*bio_{bio}.tif")) + list(
            Path(args.worldclim_dir).glob(f"*bio{bio}.tif")
        )
        if not candidates:
            raise FileNotFoundError(f"BIO{bio} raster not found")
        with rasterio.open(candidates[0]) as source:
            extracted.append(
                [value[0] if value and value[0] != source.nodata else np.nan for value in source.sample(coordinates)]
            )

    for bio, values in zip(SELECT, extracted):
        data[f"bio{bio}"] = values
    climate = [f"bio{bio}" for bio in SELECT]
    data = data.dropna(subset=climate).copy()
    data = data.drop_duplicates(["canonical_name", *climate])

    scaler = StandardScaler()
    standardized = scaler.fit_transform(data[climate])
    pca = PCA(n_components=3).fit(standardized)
    pcs = pca.transform(standardized)
    for index, column in enumerate(climate):
        data[f"z_{column}"] = standardized[:, index]
    for index in range(3):
        data[f"pc{index + 1}"] = pcs[:, index]

    if args.cell_out:
        cell_path = Path(args.cell_out)
        cell_path.parent.mkdir(parents=True, exist_ok=True)
        data.to_csv(cell_path, index=False)

    result = []
    group_columns = ["canonical_name", "family", "role", "focal_species", "match_level"]
    for keys, group in data.groupby(group_columns, dropna=False):
        n = len(group)
        row = dict(zip(group_columns, keys))
        row["n_climate_cells"] = n
        if n < args.min_cells:
            row["metric_status"] = "insufficient_cells"
            result.append(row)
            continue
        quantiles = group[climate].quantile([0.05, 0.95])
        row.update(
            {
                f"{column}_q95q05": float(quantiles.loc[0.95, column] - quantiles.loc[0.05, column])
                for column in climate
            }
        )
        row["temperature_breadth"] = float(np.mean([row[f"bio{bio}_q95q05"] for bio in [1, 5, 6, 7]]))
        row["moisture_breadth"] = float(np.mean([row[f"bio{bio}_q95q05"] for bio in [12, 14, 15, 17]]))
        row["climatic_heterogeneity"] = float(
            np.mean([group[f"z_{column}"].std(ddof=0) for column in climate])
        )
        points = group[["pc1", "pc2", "pc3"]].to_numpy()
        centre = points.mean(axis=0)
        row["pca_dispersion"] = float(np.mean(np.linalg.norm(points - centre, axis=1)))
        row["pc1_q95q05"] = float(group.pc1.quantile(0.95) - group.pc1.quantile(0.05))
        row["pc2_q95q05"] = float(group.pc2.quantile(0.95) - group.pc2.quantile(0.05))
        try:
            row["pca_hull_area"] = float(ConvexHull(group[["pc1", "pc2"]].to_numpy()).volume) if n >= 3 else np.nan
        except Exception:
            row["pca_hull_area"] = np.nan
        row["metric_status"] = "complete"
        result.append(row)

    output = pd.DataFrame(result)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.out, index=False)
    qc = {
        "occurrence_rows_with_climate": int(len(data)),
        "species_total": int(output.canonical_name.nunique()),
        "species_complete": int((output.metric_status == "complete").sum()),
        "min_cells": args.min_cells,
        "worldclim_bio_variables": SELECT,
        "global_pca_variance_explained": pca.explained_variance_ratio_.tolist(),
        "cell_output_written": bool(args.cell_out),
        "semantic_guard": "metrics describe occupied climatic breadth and heterogeneity, not physiological tolerance or causal niche expansion",
    }
    Path(args.qc_out).write_text(json.dumps(qc, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(qc))


if __name__ == "__main__":
    main()

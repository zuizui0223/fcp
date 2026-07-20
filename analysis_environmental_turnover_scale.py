#!/usr/bin/env python3
"""Test whether climate turnover among observed occurrence clusters differs by FCP scale.

Occurrences are linked into distance-threshold components. Climate turnover is measured from
WorldClim-standardized component centroids. Components describe the sampled GBIF point cloud and
must not be interpreted as biological populations or barriers to gene flow.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

EARTH_RADIUS_KM = 6371.0088
CLIMATE = ["z_bio1", "z_bio4", "z_bio5", "z_bio6", "z_bio7", "z_bio12", "z_bio14", "z_bio15", "z_bio17"]


def zscore(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    sd = values.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(np.nan, index=series.index, dtype=float)
    return (values - values.mean()) / sd


def pairwise_haversine(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    lat_r = np.radians(lat)
    lon_r = np.radians(lon)
    dlat = lat_r[:, None] - lat_r[None, :]
    dlon = lon_r[:, None] - lon_r[None, :]
    a = np.sin(dlat / 2) ** 2 + np.cos(lat_r[:, None]) * np.cos(lat_r[None, :]) * np.sin(dlon / 2) ** 2
    return 2 * EARTH_RADIUS_KM * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def component_labels(distances: np.ndarray, threshold: float) -> np.ndarray:
    n = len(distances)
    parent = list(range(n))
    size = [1] * n

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i: int, j: int) -> None:
        ri, rj = find(i), find(j)
        if ri == rj:
            return
        if size[ri] < size[rj]:
            ri, rj = rj, ri
        parent[rj] = ri
        size[ri] += size[rj]

    rows, cols = np.where(np.triu(distances <= threshold, k=1))
    for i, j in zip(rows.tolist(), cols.tolist()):
        union(i, j)
    roots = [find(i) for i in range(n)]
    mapping = {root: index for index, root in enumerate(sorted(set(roots)))}
    return np.asarray([mapping[root] for root in roots], dtype=int)


def summarize_species(group: pd.DataFrame, threshold: float, min_component_records: int) -> dict[str, object]:
    g = group.dropna(subset=["decimalLatitude", "decimalLongitude", *CLIMATE]).copy()
    result: dict[str, object] = {
        "canonical_name": group["canonical_name"].iloc[0],
        "n_occurrences": int(len(g)),
        "threshold_km": threshold,
        "min_component_records": min_component_records,
        "status": "insufficient_occurrences",
    }
    if len(g) < max(10, min_component_records * 2):
        return result
    distances = pairwise_haversine(g["decimalLatitude"].to_numpy(float), g["decimalLongitude"].to_numpy(float))
    labels = component_labels(distances, threshold)
    g["component"] = labels
    sizes = g.groupby("component").size()
    retained = sizes[sizes >= min_component_records].index
    g = g[g["component"].isin(retained)].copy()
    result["components_total"] = int(len(sizes))
    result["components_retained"] = int(len(retained))
    result["records_retained"] = int(len(g))
    if len(retained) < 2:
        result["status"] = "single_retained_component"
        return result

    centroids = g.groupby("component")[CLIMATE].mean()
    weights = g.groupby("component").size().reindex(centroids.index).to_numpy(float)
    climate_distances = np.linalg.norm(centroids.to_numpy()[:, None, :] - centroids.to_numpy()[None, :, :], axis=2)
    upper_i, upper_j = np.triu_indices(len(centroids), k=1)
    pair_values = climate_distances[upper_i, upper_j]
    pair_weights = weights[upper_i] * weights[upper_j]
    weighted_mean = float(np.average(pair_values, weights=pair_weights))

    result.update({
        "status": "complete",
        "turnover_weighted_mean": weighted_mean,
        "turnover_median": float(np.median(pair_values)),
        "turnover_max": float(np.max(pair_values)),
        "largest_component_fraction_retained": float(weights.max() / weights.sum()),
    })
    return result


def fit_model(data: pd.DataFrame, turnover_column: str) -> dict[str, object]:
    d = data.dropna(subset=["among", "family", "moisture_z", "climate_effort_z", "turnover_z"]).copy()
    result: dict[str, object] = {
        "turnover_metric": turnover_column,
        "n_species": int(len(d)),
        "n_within": int((d["among"] == 0).sum()),
        "n_among": int((d["among"] == 1).sum()),
        "n_families": int(d["family"].nunique()),
        "status": "insufficient",
    }
    if len(d) < 20 or d["among"].nunique() < 2 or d["family"].nunique() < 2:
        return result
    try:
        predictors = ["moisture_z", "climate_effort_z", "turnover_z"]
        X = sm.add_constant(d[predictors], has_constant="add")
        fit = sm.GLM(d["among"], X, family=sm.families.Binomial()).fit(
            cov_type="cluster", cov_kwds={"groups": d["family"]}
        )
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
        result.update({"status": "complete", "aic": float(fit.aic), "terms": terms})
    except Exception as exc:
        result.update({"status": "failed", "error": f"{type(exc).__name__}: {exc}"})
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--climate-cells", required=True)
    parser.add_argument("--scale-dataset", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--threshold-km", type=float, default=100.0)
    parser.add_argument("--min-component-records", type=int, default=3)
    args = parser.parse_args()

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
    summaries = [
        summarize_species(group, args.threshold_km, args.min_component_records)
        for _, group in focal_cells.groupby("canonical_name", sort=True)
    ]
    turnover = pd.DataFrame(summaries)
    turnover.to_csv(outdir / "environmental_turnover_species.csv", index=False)

    dataset = scale.merge(turnover, on="canonical_name", how="left")
    dataset["among"] = (dataset["spatial_scale"] == "among_population").astype(int)
    dataset["moisture_z"] = zscore(dataset["moisture_breadth"])
    dataset["climate_effort_z"] = zscore(np.log1p(pd.to_numeric(dataset["n_climate_cells"], errors="coerce")))
    complete = dataset["status"] == "complete"
    dataset["turnover_z"] = np.nan
    dataset.loc[complete, "turnover_z"] = zscore(dataset.loc[complete, "turnover_weighted_mean"])
    dataset.to_csv(outdir / "environmental_turnover_scale_dataset.csv", index=False)

    model = fit_model(dataset, "turnover_weighted_mean")
    manifest = {
        "threshold_km": args.threshold_km,
        "min_component_records": args.min_component_records,
        "species_summarized": int(len(turnover)),
        "species_with_two_retained_components": int((turnover.get("status") == "complete").sum()),
        "model": model,
        "interpretation_rule": "An odds ratio above one for turnover_z indicates that greater climate separation among observed GBIF occurrence components is associated with among-population rather than within-population colour variation.",
        "semantic_guard": "Occurrence components are distance-threshold clusters of sampled GBIF records, not biological populations; climate-centroid turnover is observational and does not establish selection, barriers, gene flow, or causation.",
    }
    (outdir / "environmental_turnover_scale_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

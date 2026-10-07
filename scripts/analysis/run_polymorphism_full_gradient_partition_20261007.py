#!/usr/bin/env python3
"""Full FCP geographic–climate–soil gradient decomposition.

Fits climate, soil, and block-balanced climate+soil PCA representations on
DISCOVERY environmental data only, transports them unchanged to validation and
third cohort, then decomposes continuous flower-colour turnover into unique
partial-rank components.

Primary component model:
  great-circle distance + latitude difference + elevation difference
  + climate-PC95 distance + soil-PC95 distance

Primary combined model:
  great-circle distance + latitude difference + elevation difference
  + block-balanced abiotic-PC95 distance

Inference uses matched complete-colour-vector vertex permutations.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio.warp import transform as rio_transform
from scipy.stats import rankdata

EARTH_RADIUS_KM = 6371.0088
MORPHS = ["white", "yellow_orange", "red_pink", "blue_purple"]
BIO9 = ["white", "yellow", "orange", "red", "pink", "magenta", "purple", "blue", "bronze"]
LEGACY_COLS = [f"palette_count_{x}" for x in BIO9]
THIRD_COLS = [f"flower_fraction_{x}" for x in BIO9]

CLIMATE = [f"bio{i}" for i in range(1, 20)] + ["srad_mean"]
SOIL_RAW = ["phh2o", "soc", "nitrogen", "cec", "bdod", "cfvo", "clay", "sand", "wv0033", "wv1500"]
SOIL = ["phh2o", "soc", "nitrogen", "cec", "bdod", "cfvo", "clay", "sand", "awc"]
SOIL_LOG1P = {"soc", "nitrogen", "cec", "cfvo"}
SOIL_CONVERSION = {
    "phh2o": 10.0,
    "soc": 10.0,
    "nitrogen": 100.0,
    "cec": 10.0,
    "bdod": 100.0,
    "cfvo": 10.0,
    "clay": 10.0,
    "sand": 10.0,
    "wv0033": 10.0,
    "wv1500": 10.0,
}
SOIL_DEPTHS = [("0-5cm", 5.0), ("5-15cm", 10.0), ("15-30cm", 15.0)]

PERMUTATIONS = 199
SIGNFLIP = 9999
MASTER_SEED = 2026100757
SHA = {
    "discovery": "ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4",
    "validation": "0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6",
    "third": "57630fc9f281bce94a0c40a70aaf7bce879dde93d6154175adcd021e8f5c1186",
}


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def stable_seed(*parts: object) -> int:
    payload = "|".join(map(str, (MASTER_SEED, *parts))).encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little")


def as_bool(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False).astype(bool)
    return s.fillna("").astype(str).str.strip().str.casefold().isin({"true", "1", "yes", "y"})


def sample_raster(path: Path, lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    out = np.full(len(lon), np.nan, float)
    good = np.flatnonzero(np.isfinite(lon) & np.isfinite(lat))
    if len(good) == 0:
        return out
    with rasterio.open(path) as src:
        x = lon[good].astype(float)
        y = lat[good].astype(float)
        if src.crs is not None and src.crs.to_string().upper() not in {"EPSG:4326", "OGC:CRS84"}:
            xx, yy = rio_transform("EPSG:4326", src.crs, x.tolist(), y.tolist())
            pts = list(zip(xx, yy))
        else:
            pts = list(zip(x.tolist(), y.tolist()))
        vals = list(src.sample(pts))
        nodata = src.nodata
        for idx, raw in zip(good, vals):
            value = float(raw[0])
            if nodata is not None and np.isclose(value, nodata):
                value = np.nan
            out[idx] = value
    return out


def pairwise_geo(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    latr = np.deg2rad(np.asarray(lat, float))
    lonr = np.deg2rad(np.asarray(lon, float))
    c = np.cos(latr)
    xyz = np.column_stack([c * np.cos(lonr), c * np.sin(lonr), np.sin(latr)])
    dot = np.clip(xyz @ xyz.T, -1.0, 1.0)
    return np.arccos(dot) * EARTH_RADIUS_KM


def pairwise_abs(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, float)
    return np.abs(x[:, None] - x[None, :])


def pairwise_euclidean(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, float)
    gram = x @ x.T
    sq = np.clip(np.diag(gram)[:, None] + np.diag(gram)[None, :] - 2.0 * gram, 0.0, None)
    return np.sqrt(sq)


def pairwise_jsd(prob: np.ndarray) -> np.ndarray:
    p = np.asarray(prob, float)
    mass = p.sum(axis=1)
    if p.ndim != 2 or np.any(~np.isfinite(p)) or np.any(p < 0) or np.any(mass <= 0):
        raise ValueError("invalid colour composition")
    p = p / mass[:, None]
    a = p[:, None, :]
    b = p[None, :, :]
    m = 0.5 * (a + b)
    with np.errstate(divide="ignore", invalid="ignore"):
        ka = np.where(a > 0, a * np.log2(a / m), 0).sum(axis=2)
        kb = np.where(b > 0, b * np.log2(b / m), 0).sum(axis=2)
    return np.clip(0.5 * (ka + kb), 0.0, 1.0)


def load_measured(path: Path, cohort: str) -> pd.DataFrame:
    if file_sha256(path) != SHA[cohort]:
        raise RuntimeError(f"{cohort} measured table SHA256 mismatch")
    d = pd.read_csv(path, low_memory=False)
    colour_cols = THIRD_COLS if cohort == "third" else LEGACY_COLS
    req = {
        "inat_taxon_id", "species", "photo_id", "latitude", "longitude",
        "morph", "global_classifiable", *colour_cols,
    }
    missing = sorted(req - set(d.columns))
    if missing:
        raise RuntimeError(f"{cohort}: missing {missing}")
    keep = as_bool(d["global_classifiable"]) & d["morph"].astype(str).isin(MORPHS)
    d = d.loc[keep].copy()
    counts = d.groupby("inat_taxon_id").size()
    d = d.loc[d["inat_taxon_id"].isin(counts[counts >= 40].index)].copy()
    d["cohort"] = cohort
    return d.sort_values(["inat_taxon_id", "photo_id"], kind="stable").reset_index(drop=True)


def attach_climate_elevation(
    d: pd.DataFrame,
    bio_dir: Path,
    srad_dir: Path,
    elev: Path,
) -> pd.DataFrame:
    d = d.copy()
    lon = pd.to_numeric(d["longitude"], errors="coerce").to_numpy(float)
    lat = pd.to_numeric(d["latitude"], errors="coerce").to_numpy(float)
    for i in range(1, 20):
        d[f"bio{i}"] = sample_raster(bio_dir / f"wc2.1_10m_bio_{i}.tif", lon, lat)
    srad = np.column_stack([
        sample_raster(srad_dir / f"wc2.1_10m_srad_{m:02d}.tif", lon, lat)
        for m in range(1, 13)
    ])
    with np.errstate(invalid="ignore"):
        d["srad_mean"] = np.nanmean(srad, axis=1)
    d["elevation_m"] = sample_raster(elev, lon, lat)
    return d


def attach_soil(d: pd.DataFrame, soil_dir: Path) -> pd.DataFrame:
    d = d.copy()
    lon = pd.to_numeric(d["longitude"], errors="coerce").to_numpy(float)
    lat = pd.to_numeric(d["latitude"], errors="coerce").to_numpy(float)
    raw_values = {}
    for prop in SOIL_RAW:
        layers = []
        weights = []
        for depth, thickness in SOIL_DEPTHS:
            path = soil_dir / prop / f"{prop}_{depth}_mean_5000.tif"
            layers.append(sample_raster(path, lon, lat))
            weights.append(thickness)
        arr = np.column_stack(layers)
        w = np.asarray(weights, float)
        good = np.isfinite(arr).all(axis=1)
        value = np.full(len(d), np.nan, float)
        value[good] = (arr[good] @ w) / w.sum()
        value = value / SOIL_CONVERSION[prop]
        if prop in SOIL_LOG1P:
            value = np.where(np.isfinite(value) & (value >= 0), np.log1p(value), np.nan)
        raw_values[prop] = value
    for prop in ["phh2o", "soc", "nitrogen", "cec", "bdod", "cfvo", "clay", "sand"]:
        d[f"soil_{prop}"] = raw_values[prop]
    awc = raw_values["wv0033"] - raw_values["wv1500"]
    d["soil_awc"] = np.where(np.isfinite(awc) & (awc >= 0), awc, np.nan)
    return d


def fit_standardized_pca(x: np.ndarray, names: list[str], *, block_scale: np.ndarray | None = None) -> dict:
    x = np.asarray(x, float)
    good = np.isfinite(x).all(axis=1)
    x = x[good]
    if len(x) < 1000:
        raise RuntimeError("too few complete discovery rows for PCA")
    mean = x.mean(axis=0)
    sd = x.std(axis=0, ddof=0)
    if np.any(~np.isfinite(sd)) or np.any(sd <= 0):
        raise RuntimeError("degenerate variable in PCA block")
    z = (x - mean) / sd
    if block_scale is not None:
        z = z * block_scale[None, :]
    _, s, vt = np.linalg.svd(z, full_matrices=False)
    eig = (s * s) / (len(z) - 1)
    ratio = eig / eig.sum()
    cum = np.cumsum(ratio)
    return {
        "names": list(names),
        "mean": mean,
        "sd": sd,
        "block_scale": np.ones(len(names), float) if block_scale is None else block_scale,
        "loadings": vt.T,
        "explained": ratio,
        "cumulative": cum,
        "k80": int(np.searchsorted(cum, 0.80) + 1),
        "k90": int(np.searchsorted(cum, 0.90) + 1),
        "k95": int(np.searchsorted(cum, 0.95) + 1),
        "n_fit_rows": int(len(z)),
    }


def add_pca_scores(d: pd.DataFrame, columns: list[str], pca: dict, prefix: str) -> pd.DataFrame:
    d = d.copy()
    x = d[columns].to_numpy(float)
    good = np.isfinite(x).all(axis=1)
    z = np.full_like(x, np.nan, dtype=float)
    z[good] = ((x[good] - pca["mean"]) / pca["sd"]) * pca["block_scale"]
    scores = np.full((len(d), len(columns)), np.nan, float)
    scores[good] = z[good] @ pca["loadings"]
    for j in range(len(columns)):
        d[f"{prefix}_pc_{j+1}"] = scores[:, j]
    return d


def prepare_pca_blocks(discovery: pd.DataFrame) -> tuple[dict, dict, dict]:
    soil_cols = [f"soil_{x}" for x in SOIL]
    climate_pca = fit_standardized_pca(discovery[CLIMATE].to_numpy(float), CLIMATE)
    soil_pca = fit_standardized_pca(discovery[soil_cols].to_numpy(float), soil_cols)

    combined_cols = CLIMATE + soil_cols
    scale = np.concatenate([
        np.repeat(1.0 / np.sqrt(len(CLIMATE)), len(CLIMATE)),
        np.repeat(1.0 / np.sqrt(len(soil_cols)), len(soil_cols)),
    ])
    combined_pca = fit_standardized_pca(discovery[combined_cols].to_numpy(float), combined_cols, block_scale=scale)
    return climate_pca, soil_pca, combined_pca


def residual_basis(covariates: list[np.ndarray]) -> np.ndarray:
    if not covariates:
        return np.ones((len(covariates[0]), 1), float)
    X = np.column_stack([np.ones(len(covariates[0]), float), *covariates])
    Q, _ = np.linalg.qr(X, mode="reduced")
    return Q


def partial_rank_components(y: np.ndarray, predictors: dict[str, np.ndarray]) -> tuple[dict[str, float], dict[str, tuple[np.ndarray, np.ndarray, float]]]:
    yr = rankdata(y, method="average")
    pred_rank = {k: rankdata(v, method="average") for k, v in predictors.items()}
    observed = {}
    state = {}
    for focal, xr in pred_rank.items():
        cov = [v for k, v in pred_rank.items() if k != focal]
        X = np.column_stack([np.ones(len(yr), float), *cov])
        Q, _ = np.linalg.qr(X, mode="reduced")
        xres = xr - Q @ (Q.T @ xr)
        yres = yr - Q @ (Q.T @ yr)
        xn = float(np.linalg.norm(xres))
        yn = float(np.linalg.norm(yres))
        rho = 0.0 if xn <= 1e-14 or yn <= 1e-14 else float(np.dot(xres, yres) / (xn * yn))
        observed[focal] = rho
        state[focal] = (Q, xres, xn)
    return observed, state


def permuted_component_nulls(
    rank_matrix: np.ndarray,
    u: np.ndarray,
    v: np.ndarray,
    n: int,
    state: dict[str, tuple[np.ndarray, np.ndarray, float]],
    *,
    cohort: str,
    taxon: int,
    model: str,
) -> dict[str, np.ndarray]:
    out = {k: np.empty(PERMUTATIONS, float) for k in state}
    rng = np.random.default_rng(stable_seed(cohort, taxon, model))
    batch = 24
    for start in range(0, PERMUTATIONS, batch):
        stop = min(PERMUTATIONS, start + batch)
        pp = np.stack([rng.permutation(n) for _ in range(stop - start)])
        vals = rank_matrix[pp[:, u], pp[:, v]]
        for name, (Q, xres, xn) in state.items():
            proj = (vals @ Q) @ Q.T
            yres = vals - proj
            yn = np.linalg.norm(yres, axis=1)
            num = yres @ xres
            out[name][start:stop] = np.where((yn <= 1e-14) | (xn <= 1e-14), 0.0, num / (yn * xn))
    return out


def species_analysis(g: pd.DataFrame, cohort: str, climate_pca: dict, soil_pca: dict, combined_pca: dict) -> tuple[dict, dict, dict]:
    colour_cols = THIRD_COLS if cohort == "third" else LEGACY_COLS
    required = [
        "latitude", "longitude", "elevation_m", *colour_cols,
        *[f"clim_pc_{j+1}" for j in range(climate_pca["k95"])],
        *[f"soil_pc_{j+1}" for j in range(soil_pca["k95"])],
        *[f"abiotic_pc_{j+1}" for j in range(combined_pca["k95"])],
    ]
    g = g.dropna(subset=required).copy()
    if len(g) < 40:
        raise ValueError("below full-gradient complete-case eligibility")
    g = g.sort_values("photo_id", kind="stable").reset_index(drop=True)
    n = len(g)
    u, v = np.triu_indices(n, k=1)

    lat = g["latitude"].to_numpy(float)
    lon = g["longitude"].to_numpy(float)
    elev = g["elevation_m"].to_numpy(float)

    geo = pairwise_geo(lat, lon)[u, v]
    latdiff = pairwise_abs(lat)[u, v]
    elevdiff = pairwise_abs(elev)[u, v]

    climate_scores = g[[f"clim_pc_{j+1}" for j in range(climate_pca["k95"])]].to_numpy(float)
    soil_scores = g[[f"soil_pc_{j+1}" for j in range(soil_pca["k95"])]].to_numpy(float)
    abiotic_scores = g[[f"abiotic_pc_{j+1}" for j in range(combined_pca["k95"])]].to_numpy(float)

    clim = pairwise_euclidean(climate_scores)[u, v]
    soil = pairwise_euclidean(soil_scores)[u, v]
    abiotic = pairwise_euclidean(abiotic_scores)[u, v]

    colour = pairwise_jsd(g[colour_cols].to_numpy(float))
    y = colour[u, v]
    if np.ptp(y) <= 1e-15:
        raise ValueError("degenerate continuous colour distance")

    model_a = {
        "distance": geo,
        "latitude": latdiff,
        "elevation": elevdiff,
        "climate": clim,
        "soil": soil,
    }
    model_b = {
        "distance": geo,
        "latitude": latdiff,
        "elevation": elevdiff,
        "abiotic": abiotic,
    }
    obs_a, state_a = partial_rank_components(y, model_a)
    obs_b, state_b = partial_rank_components(y, model_b)

    yr = rankdata(y, method="average")
    rank_matrix = np.zeros((n, n), float)
    rank_matrix[u, v] = yr
    rank_matrix[v, u] = yr

    taxon = int(g["inat_taxon_id"].iloc[0])
    null_a = permuted_component_nulls(rank_matrix, u, v, n, state_a, cohort=cohort, taxon=taxon, model="A")
    null_b = permuted_component_nulls(rank_matrix, u, v, n, state_b, cohort=cohort, taxon=taxon, model="B")

    row = {
        "cohort": cohort,
        "inat_taxon_id": taxon,
        "species": str(g["species"].iloc[0]),
        "n": int(n),
    }
    for k, val in obs_a.items():
        row[f"A_{k}"] = val
    for k, val in obs_b.items():
        row[f"B_{k}"] = val
    return row, null_a, null_b


def holm(raw: dict[str, float]) -> dict[str, float]:
    order = sorted(raw, key=lambda k: raw[k])
    out = {}
    running = 0.0
    m = len(order)
    for i, name in enumerate(order):
        adjusted = min(1.0, (m - i) * raw[name])
        running = max(running, adjusted)
        out[name] = running
    return out


def summarize_model(df: pd.DataFrame, null_list: list[dict[str, np.ndarray]], prefix: str, cohort: str) -> dict:
    names = [c[len(prefix)+1:] for c in df.columns if c.startswith(prefix + "_")]
    names = sorted(names)
    raw_p = {}
    results = {}
    mats = {}
    for name in names:
        observed = df[f"{prefix}_{name}"].to_numpy(float)
        mat = np.vstack([x[name] for x in null_list])
        mats[name] = mat
        null_mean = mat.mean(axis=0)
        obs_mean = float(np.mean(observed))
        p = float((1 + np.count_nonzero(null_mean >= obs_mean)) / (PERMUTATIONS + 1))
        raw_p[name] = p
        results[name] = {
            "mean_partial_rho": obs_mean,
            "median_partial_rho": float(np.median(observed)),
            "positive_species_fraction": float(np.mean(observed > 0)),
            "raw_matched_p_upper": p,
            "null_mean": float(np.mean(null_mean)),
            "null_q025": float(np.quantile(null_mean, 0.025)),
            "null_q975": float(np.quantile(null_mean, 0.975)),
        }
    adj = holm(raw_p)
    for name in names:
        results[name]["holm_p"] = float(adj[name])
        results[name]["supported_holm"] = bool(results[name]["mean_partial_rho"] > 0 and adj[name] <= 0.05)

    if prefix == "B":
        env = df["B_abiotic"].to_numpy(float)
        geo = df["B_distance"].to_numpy(float)
        delta = env - geo
        null_delta = (mats["abiotic"] - mats["distance"]).mean(axis=0)
        obs = float(np.mean(delta))
        p_perm = float((1 + np.count_nonzero(null_delta >= obs)) / (PERMUTATIONS + 1))
        rng = np.random.default_rng(stable_seed("signflip", cohort, "B"))
        sf = np.empty(SIGNFLIP, float)
        for i in range(SIGNFLIP):
            signs = rng.choice(np.array([-1.0, 1.0]), size=len(delta))
            sf[i] = float(np.mean(delta * signs))
        p_sign = float((1 + np.count_nonzero(sf >= obs)) / (SIGNFLIP + 1))
        results["abiotic_vs_distance"] = {
            "mean_abiotic_minus_distance": obs,
            "median_abiotic_minus_distance": float(np.median(delta)),
            "fraction_abiotic_gt_distance": float(np.mean(delta > 0)),
            "matched_null_p_upper": p_perm,
            "species_signflip_p_upper": p_sign,
            "abiotic_stronger_supported": bool(obs > 0 and p_perm < 0.05 and p_sign < 0.05),
        }
    return results


def cohort_analysis(d: pd.DataFrame, cohort: str, climate_pca: dict, soil_pca: dict, combined_pca: dict) -> tuple[dict, pd.DataFrame]:
    rows = []
    null_a = []
    null_b = []
    skipped = []
    for tid, g in d.groupby("inat_taxon_id", sort=True):
        try:
            row, na, nb = species_analysis(g, cohort, climate_pca, soil_pca, combined_pca)
        except ValueError as exc:
            skipped.append({"inat_taxon_id": int(tid), "reason": str(exc)})
            continue
        rows.append(row)
        null_a.append(na)
        null_b.append(nb)
    if not rows:
        raise RuntimeError(f"{cohort}: no complete species")
    df = pd.DataFrame(rows).sort_values("inat_taxon_id").reset_index(drop=True)
    return {
        "n_species": int(len(df)),
        "skipped_species": int(len(skipped)),
        "model_A_unique_components": summarize_model(df, null_a, "A", cohort),
        "model_B_combined_abiotic": summarize_model(df, null_b, "B", cohort),
    }, df


def pca_summary(pca: dict) -> dict:
    k = pca["k95"]
    pcs = []
    for j in range(k):
        load = pca["loadings"][:, j]
        idx = np.argsort(np.abs(load))[::-1][:6]
        pcs.append({
            "pc": j + 1,
            "explained_ratio": float(pca["explained"][j]),
            "top_loadings": [{"variable": pca["names"][i], "loading": float(load[i])} for i in idx],
        })
    return {
        "n_variables": len(pca["names"]),
        "n_fit_rows": pca["n_fit_rows"],
        "k80": pca["k80"],
        "k90": pca["k90"],
        "k95": pca["k95"],
        "pc95_cumulative": float(pca["cumulative"][pca["k95"] - 1]),
        "primary_pc_loadings": pcs,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--discovery", type=Path, required=True)
    ap.add_argument("--validation", type=Path, required=True)
    ap.add_argument("--third", type=Path, required=True)
    ap.add_argument("--bio-dir", type=Path, required=True)
    ap.add_argument("--srad-dir", type=Path, required=True)
    ap.add_argument("--elevation", type=Path, required=True)
    ap.add_argument("--soil-dir", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()

    raw = {
        "discovery": load_measured(args.discovery, "discovery"),
        "validation": load_measured(args.validation, "validation"),
        "third": load_measured(args.third, "third"),
    }
    enriched = {}
    for cohort, d in raw.items():
        d = attach_climate_elevation(d, args.bio_dir, args.srad_dir, args.elevation)
        d = attach_soil(d, args.soil_dir)
        enriched[cohort] = d

    climate_pca, soil_pca, combined_pca = prepare_pca_blocks(enriched["discovery"])
    soil_cols = [f"soil_{x}" for x in SOIL]
    combined_cols = CLIMATE + soil_cols

    for cohort in enriched:
        d = add_pca_scores(enriched[cohort], CLIMATE, climate_pca, "clim")
        d = add_pca_scores(d, soil_cols, soil_pca, "soil")
        d = add_pca_scores(d, combined_cols, combined_pca, "abiotic")
        enriched[cohort] = d

    disc, ddf = cohort_analysis(enriched["discovery"], "discovery", climate_pca, soil_pca, combined_pca)
    val, vdf = cohort_analysis(enriched["validation"], "validation", climate_pca, soil_pca, combined_pca)
    third, tdf = cohort_analysis(enriched["third"], "third", climate_pca, soil_pca, combined_pca)

    A_names = ["distance", "latitude", "elevation", "climate", "soil"]
    B_names = ["distance", "latitude", "elevation", "abiotic"]
    result = {
        "schema": "fcp_full_gradient_partition_posthoc_v1",
        "date_jst": "2026-10-07",
        "status": "complete_full_geographic_climate_soil_decomposition",
        "confirmatory_decisions_changed": False,
        "environmental_sources": {
            "climate": "WorldClim 2.1 BIO1-BIO19 + mean monthly solar radiation",
            "elevation": "WorldClim 2.1 10-minute elevation",
            "soil": "SoilGrids 2.0 5-km aggregated mean predictions, thickness-weighted 0-30 cm",
            "soil_source_properties": SOIL_RAW,
            "soil_PCA_features": SOIL,
            "soil_available_water_proxy": "awc = thickness-weighted wv0033 - thickness-weighted wv1500",
            "soil_log1p_properties": sorted(SOIL_LOG1P),
        },
        "PCA": {
            "climate": pca_summary(climate_pca),
            "soil": pca_summary(soil_pca),
            "combined_block_balanced_abiotic": pca_summary(combined_pca),
            "combined_block_weighting": {
                "climate_variable_multiplier": float(1.0 / np.sqrt(len(CLIMATE))),
                "soil_variable_multiplier": float(1.0 / np.sqrt(len(SOIL))),
            },
        },
        "discovery": disc,
        "validation": val,
        "third": third,
        "replication": {
            "model_A_500_plus_500": {
                name: bool(
                    disc["model_A_unique_components"][name]["supported_holm"]
                    and val["model_A_unique_components"][name]["supported_holm"]
                )
                for name in A_names
            },
            "model_A_all_three": {
                name: bool(
                    all(x["model_A_unique_components"][name]["supported_holm"] for x in (disc, val, third))
                )
                for name in A_names
            },
            "model_B_500_plus_500": {
                name: bool(
                    disc["model_B_combined_abiotic"][name]["supported_holm"]
                    and val["model_B_combined_abiotic"][name]["supported_holm"]
                )
                for name in B_names
            },
            "model_B_all_three": {
                name: bool(
                    all(x["model_B_combined_abiotic"][name]["supported_holm"] for x in (disc, val, third))
                )
                for name in B_names
            },
            "abiotic_stronger_than_distance_500_plus_500": bool(
                disc["model_B_combined_abiotic"]["abiotic_vs_distance"]["abiotic_stronger_supported"]
                and val["model_B_combined_abiotic"]["abiotic_vs_distance"]["abiotic_stronger_supported"]
            ),
            "abiotic_stronger_than_distance_all_three": bool(
                all(x["model_B_combined_abiotic"]["abiotic_vs_distance"]["abiotic_stronger_supported"] for x in (disc, val, third))
            ),
        },
        "hard_nonclaims": [
            "post hoc decomposition",
            "SoilGrids is modelled macro-scale soil, not plant-level soil measurement",
            "macroclimate plus soil is not total environment",
            "does not include pollinators, herbivory, direct UV-B, land use, demography or neutral genetic structure",
            "not genetic isolation by environment",
            "does not establish local adaptation or causal selection",
            "cannot distinguish genetic differentiation from phenotypic plasticity",
        ],
    }

    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)
    pd.concat([ddf, vdf, tdf], ignore_index=True).to_csv(out / "species_full_gradient_metrics.csv", index=False)
    (out / "result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Post hoc FCP IBD versus multivariate macroclimate IBE-like decomposition.

The environmental representation is constructed without flower-colour outcomes:
WorldClim BIO1--BIO19 plus mean monthly solar radiation are standardized and
PCA-fitted on discovery observations only, then transported unchanged to
validation and third-cohort observations.

Primary environmental distance uses the smallest number of discovery PCs
explaining >=95% of environmental variance. Sensitivities use 80%, 90%, and
full standardized 20-dimensional environmental distance.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from scipy.stats import rankdata, spearmanr

EARTH_RADIUS_KM = 6371.0088
MORPHS = ["white", "yellow_orange", "red_pink", "blue_purple"]
BIO9 = ["white", "yellow", "orange", "red", "pink", "magenta", "purple", "blue", "bronze"]
LEGACY_COLS = [f"palette_count_{x}" for x in BIO9]
THIRD_COLS = [f"flower_fraction_{x}" for x in BIO9]
ENV_NAMES = [f"bio{i}" for i in range(1, 20)] + ["srad_mean"]
THRESHOLDS = {"pc80": 0.80, "pc90": 0.90, "pc95": 0.95}
PERMUTATIONS = 199
SIGNFLIP = 9999
MASTER_SEED = 2026100751
SHA = {
    "discovery": "ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4",
    "validation": "0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6",
    "third": "57630fc9f281bce94a0c40a70aaf7bce879dde93d6154175adcd021e8f5c1186",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def stable_seed(*parts: object) -> int:
    b = "|".join(map(str, (MASTER_SEED, *parts))).encode()
    return int.from_bytes(hashlib.sha256(b).digest()[:8], "little")


def as_bool(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False).astype(bool)
    return s.fillna("").astype(str).str.strip().str.casefold().isin({"true", "1", "yes", "y"})


def sample_raster(path: Path, lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    out = np.full(len(lon), np.nan, float)
    good = [i for i, (x, y) in enumerate(zip(lon, lat)) if np.isfinite(x) and np.isfinite(y)]
    if not good:
        return out
    points = [(float(lon[i]), float(lat[i])) for i in good]
    with rasterio.open(path) as src:
        vals = list(src.sample(points))
        for i, v in zip(good, vals):
            x = float(v[0])
            if src.nodata is not None and np.isclose(x, src.nodata):
                x = np.nan
            out[i] = x
    return out


def pairwise_geo(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    latr = np.deg2rad(np.asarray(lat, float))
    lonr = np.deg2rad(np.asarray(lon, float))
    c = np.cos(latr)
    xyz = np.column_stack([c * np.cos(lonr), c * np.sin(lonr), np.sin(latr)])
    dot = np.clip(xyz @ xyz.T, -1.0, 1.0)
    return np.arccos(dot) * EARTH_RADIUS_KM


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
    return np.clip(0.5 * (ka + kb), 0, 1)


def partial_from_centered(yc: np.ndarray, xc: np.ndarray, zc: np.ndarray) -> float:
    z2 = float(np.dot(zc, zc))
    yr = yc - (0.0 if z2 <= 1e-15 else float(np.dot(zc, yc) / z2)) * zc
    xr = xc - (0.0 if z2 <= 1e-15 else float(np.dot(zc, xc) / z2)) * zc
    den = float(np.linalg.norm(yr) * np.linalg.norm(xr))
    return 0.0 if den <= 1e-14 else float(np.dot(yr, xr) / den)


def r2(y: np.ndarray, X: np.ndarray) -> float:
    y = np.asarray(y, float)
    X = np.asarray(X, float)
    A = np.column_stack([np.ones(len(y)), X])
    b, *_ = np.linalg.lstsq(A, y, rcond=None)
    fit = A @ b
    ss = float(np.dot(y - y.mean(), y - y.mean()))
    if ss <= 1e-15:
        return 0.0
    return max(0.0, min(1.0, 1 - float(np.dot(y - fit, y - fit)) / ss))


def load_measured(path: Path, cohort: str) -> pd.DataFrame:
    if sha256(path) != SHA[cohort]:
        raise RuntimeError(f"{cohort} SHA256 mismatch")
    d = pd.read_csv(path, low_memory=False)
    colour_cols = THIRD_COLS if cohort == "third" else LEGACY_COLS
    required = {
        "inat_taxon_id", "species", "photo_id", "latitude", "longitude",
        "morph", "global_classifiable", *colour_cols,
    }
    missing = sorted(required - set(d.columns))
    if missing:
        raise RuntimeError(f"{cohort} missing columns: {missing}")
    keep = as_bool(d["global_classifiable"]) & d["morph"].astype(str).isin(MORPHS)
    d = d.loc[keep].copy()
    counts = d.groupby("inat_taxon_id").size()
    d = d.loc[d["inat_taxon_id"].isin(counts[counts >= 40].index)].copy()
    d["cohort"] = cohort
    return d.sort_values(["inat_taxon_id", "photo_id"], kind="stable").reset_index(drop=True)


def attach_environment(d: pd.DataFrame, bio_dir: Path, srad_dir: Path) -> pd.DataFrame:
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
    return d


def fit_discovery_pca(d: pd.DataFrame) -> dict:
    x = d[ENV_NAMES].to_numpy(float)
    good = np.isfinite(x).all(axis=1)
    x = x[good]
    if len(x) < 1000:
        raise RuntimeError("too few discovery rows for environmental PCA")
    mean = x.mean(axis=0)
    sd = x.std(axis=0, ddof=0)
    if np.any(~np.isfinite(sd)) or np.any(sd <= 0):
        raise RuntimeError("degenerate environmental variable in discovery")
    z = (x - mean) / sd
    _, s, vt = np.linalg.svd(z, full_matrices=False)
    eigen = (s * s) / (len(z) - 1)
    ratio = eigen / eigen.sum()
    cumulative = np.cumsum(ratio)
    k = {name: int(np.searchsorted(cumulative, threshold) + 1) for name, threshold in THRESHOLDS.items()}
    return {
        "mean": mean,
        "sd": sd,
        "loadings": vt.T,
        "eigenvalues": eigen,
        "explained_ratio": ratio,
        "cumulative": cumulative,
        "k": k,
        "n_fit_rows": int(len(z)),
    }


def transform_environment(d: pd.DataFrame, pca: dict) -> pd.DataFrame:
    d = d.copy()
    x = d[ENV_NAMES].to_numpy(float)
    good = np.isfinite(x).all(axis=1)
    z = np.full_like(x, np.nan, dtype=float)
    z[good] = (x[good] - pca["mean"]) / pca["sd"]
    scores = np.full((len(d), len(ENV_NAMES)), np.nan, float)
    scores[good] = z[good] @ pca["loadings"]
    for j in range(len(ENV_NAMES)):
        d[f"env_z_{j+1}"] = z[:, j]
        d[f"env_pc_{j+1}"] = scores[:, j]
    return d


def species_stats(
    g: pd.DataFrame,
    cohort: str,
    pca: dict,
) -> tuple[dict, np.ndarray, np.ndarray]:
    colour_cols = THIRD_COLS if cohort == "third" else LEGACY_COLS
    required = ["latitude", "longitude", *colour_cols] + [f"env_z_{j+1}" for j in range(len(ENV_NAMES))]
    g = g.dropna(subset=required).copy()
    if len(g) < 40:
        raise ValueError("below complete-environment eligibility")
    g = g.sort_values("photo_id", kind="stable").reset_index(drop=True)
    n = len(g)
    u, v = np.triu_indices(n, k=1)

    gp = pairwise_geo(g["latitude"].to_numpy(float), g["longitude"].to_numpy(float))[u, v]
    colour = pairwise_jsd(g[colour_cols].to_numpy(float))
    yp = colour[u, v]
    if np.ptp(yp) <= 1e-15 or np.ptp(gp) <= 1e-15:
        raise ValueError("degenerate colour or geography")

    env_dist = {}
    for name in ("pc80", "pc90", "pc95"):
        k = pca["k"][name]
        scores = g[[f"env_pc_{j+1}" for j in range(k)]].to_numpy(float)
        env_dist[name] = pairwise_euclidean(scores)[u, v]
    z = g[[f"env_z_{j+1}" for j in range(len(ENV_NAMES))]].to_numpy(float)
    env_dist["full20"] = pairwise_euclidean(z)[u, v]

    gr = rankdata(gp, method="average")
    yr = rankdata(yp, method="average")
    gc = gr - gr.mean()
    yc = yr - yr.mean()

    observed = {}
    env_centered = {}
    for name, ep in env_dist.items():
        er = rankdata(ep, method="average")
        ec = er - er.mean()
        env_centered[name] = ec
        ibe = partial_from_centered(yc, ec, gc)
        ibd = partial_from_centered(yc, gc, ec)
        Rg = r2(yr, gr[:, None])
        Re = r2(yr, er[:, None])
        Rge = r2(yr, np.column_stack([gr, er]))
        observed[name] = {
            "rho_ENV": ibe,
            "rho_GEO": ibd,
            "delta_ENV_minus_GEO": ibe - ibd,
            "R2_geo": Rg,
            "R2_env": Re,
            "R2_both": Rge,
            "unique_geo": Rge - Re,
            "unique_env": Rge - Rg,
            "shared_geo_env": Rg + Re - Rge,
            "rho_geo_env": float(spearmanr(gp, ep).statistic),
        }

    # Primary matched null: PC95 only.
    primary = "pc95"
    ec = env_centered[primary]
    rank_matrix = np.zeros((n, n), float)
    rank_matrix[u, v] = yr
    rank_matrix[v, u] = yr
    rng = np.random.default_rng(stable_seed(cohort, int(g["inat_taxon_id"].iloc[0]), primary))
    env_null = np.empty(PERMUTATIONS, float)
    geo_null = np.empty(PERMUTATIONS, float)
    batch = 32

    g2 = float(np.dot(gc, gc))
    e2 = float(np.dot(ec, ec))
    e_res = ec - (0.0 if g2 <= 1e-15 else float(np.dot(gc, ec) / g2)) * gc
    g_res = gc - (0.0 if e2 <= 1e-15 else float(np.dot(ec, gc) / e2)) * ec
    en = float(np.linalg.norm(e_res))
    gn = float(np.linalg.norm(g_res))

    for start in range(0, PERMUTATIONS, batch):
        stop = min(PERMUTATIONS, start + batch)
        pp = np.stack([rng.permutation(n) for _ in range(stop - start)])
        vals = rank_matrix[pp[:, u], pp[:, v]]
        vals = vals - vals.mean(axis=1, keepdims=True)

        beta_g = np.zeros(stop - start) if g2 <= 1e-15 else (vals @ gc) / g2
        y_rg = vals - beta_g[:, None] * gc
        yn = np.linalg.norm(y_rg, axis=1)
        env_null[start:stop] = np.where((yn <= 1e-14) | (en <= 1e-14), 0.0, (y_rg @ e_res) / (yn * en))

        beta_e = np.zeros(stop - start) if e2 <= 1e-15 else (vals @ ec) / e2
        y_re = vals - beta_e[:, None] * ec
        yn2 = np.linalg.norm(y_re, axis=1)
        geo_null[start:stop] = np.where((yn2 <= 1e-14) | (gn <= 1e-14), 0.0, (y_re @ g_res) / (yn2 * gn))

    row = {
        "cohort": cohort,
        "inat_taxon_id": int(g["inat_taxon_id"].iloc[0]),
        "species": str(g["species"].iloc[0]),
        "n": int(n),
    }
    for name, vals in observed.items():
        for key, value in vals.items():
            row[f"{name}_{key}"] = value
    return row, env_null, geo_null


def summarize(df: pd.DataFrame, env_nulls: list[np.ndarray], geo_nulls: list[np.ndarray], cohort: str) -> dict:
    primary = "pc95"
    env = df[f"{primary}_rho_ENV"].to_numpy(float)
    geo = df[f"{primary}_rho_GEO"].to_numpy(float)
    env_mat = np.vstack(env_nulls)
    geo_mat = np.vstack(geo_nulls)
    env_null = env_mat.mean(axis=0)
    geo_null = geo_mat.mean(axis=0)

    obs_env = float(np.mean(env))
    obs_geo = float(np.mean(geo))
    p_env = float((1 + np.count_nonzero(env_null >= obs_env)) / (PERMUTATIONS + 1))
    p_geo = float((1 + np.count_nonzero(geo_null >= obs_geo)) / (PERMUTATIONS + 1))

    delta = env - geo
    obs_delta = float(np.mean(delta))
    null_delta = (env_mat - geo_mat).mean(axis=0)
    p_delta_perm = float((1 + np.count_nonzero(null_delta >= obs_delta)) / (PERMUTATIONS + 1))

    rng = np.random.default_rng(stable_seed("signflip", cohort, primary))
    sf = np.empty(SIGNFLIP, float)
    for i in range(SIGNFLIP):
        signs = rng.choice(np.array([-1.0, 1.0]), size=len(delta))
        sf[i] = float(np.mean(delta * signs))
    p_delta_sign = float((1 + np.count_nonzero(sf >= obs_delta)) / (SIGNFLIP + 1))

    sensitivities = {}
    for name in ("pc80", "pc90", "pc95", "full20"):
        e = df[f"{name}_rho_ENV"].to_numpy(float)
        g = df[f"{name}_rho_GEO"].to_numpy(float)
        sensitivities[name] = {
            "mean_rho_ENV": float(np.mean(e)),
            "median_rho_ENV": float(np.median(e)),
            "positive_ENV_fraction": float(np.mean(e > 0)),
            "mean_rho_GEO": float(np.mean(g)),
            "median_rho_GEO": float(np.median(g)),
            "positive_GEO_fraction": float(np.mean(g > 0)),
            "mean_ENV_minus_GEO": float(np.mean(e - g)),
            "mean_unique_env_R2": float(df[f"{name}_unique_env"].mean()),
            "mean_unique_geo_R2": float(df[f"{name}_unique_geo"].mean()),
            "mean_shared_R2": float(df[f"{name}_shared_geo_env"].mean()),
            "mean_geo_env_rho": float(df[f"{name}_rho_geo_env"].mean()),
        }

    return {
        "n_species": int(len(df)),
        "primary_pc95": {
            "ENV_IBE_like": {
                "mean_partial_rho": obs_env,
                "median_partial_rho": float(np.median(env)),
                "positive_species_fraction": float(np.mean(env > 0)),
                "matched_vertex_p_upper": p_env,
                "null_mean": float(np.mean(env_null)),
                "null_q025": float(np.quantile(env_null, 0.025)),
                "null_q975": float(np.quantile(env_null, 0.975)),
                "supported": bool(obs_env > 0 and p_env < 0.05),
            },
            "GEO_IBD_like": {
                "mean_partial_rho": obs_geo,
                "median_partial_rho": float(np.median(geo)),
                "positive_species_fraction": float(np.mean(geo > 0)),
                "matched_vertex_p_upper": p_geo,
                "null_mean": float(np.mean(geo_null)),
                "null_q025": float(np.quantile(geo_null, 0.025)),
                "null_q975": float(np.quantile(geo_null, 0.975)),
                "supported": bool(obs_geo > 0 and p_geo < 0.05),
            },
            "relative_balance": {
                "mean_ENV_minus_GEO": obs_delta,
                "median_ENV_minus_GEO": float(np.median(delta)),
                "fraction_ENV_gt_GEO": float(np.mean(delta > 0)),
                "matched_null_p_upper_for_positive_delta": p_delta_perm,
                "species_signflip_p_upper_for_positive_delta": p_delta_sign,
                "ENV_stronger_supported": bool(obs_delta > 0 and p_delta_perm < 0.05 and p_delta_sign < 0.05),
            },
        },
        "sensitivity": sensitivities,
    }


def run_cohort(d: pd.DataFrame, cohort: str, pca: dict) -> tuple[dict, pd.DataFrame]:
    rows = []
    env_nulls = []
    geo_nulls = []
    skipped = []
    for tid, g in d.groupby("inat_taxon_id", sort=True):
        try:
            row, en, gn = species_stats(g, cohort, pca)
        except ValueError as exc:
            skipped.append({"inat_taxon_id": int(tid), "reason": str(exc)})
            continue
        rows.append(row)
        env_nulls.append(en)
        geo_nulls.append(gn)
    if not rows:
        raise RuntimeError(f"{cohort}: no evaluable species")
    df = pd.DataFrame(rows).sort_values("inat_taxon_id").reset_index(drop=True)
    summary = summarize(df, env_nulls, geo_nulls, cohort)
    summary["skipped_species"] = int(len(skipped))
    return summary, df


def top_loadings(pca: dict, n_pc: int, n_var: int = 6) -> list[dict]:
    out = []
    for j in range(n_pc):
        load = pca["loadings"][:, j]
        order = np.argsort(np.abs(load))[::-1][:n_var]
        out.append({
            "pc": j + 1,
            "explained_ratio": float(pca["explained_ratio"][j]),
            "top_loadings": [
                {"variable": ENV_NAMES[i], "loading": float(load[i])}
                for i in order
            ],
        })
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--discovery", type=Path, required=True)
    ap.add_argument("--validation", type=Path, required=True)
    ap.add_argument("--third", type=Path, required=True)
    ap.add_argument("--bio-dir", type=Path, required=True)
    ap.add_argument("--srad-dir", type=Path, required=True)
    ap.add_argument("--bio5-result", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    a = ap.parse_args()

    raw = {
        "discovery": load_measured(a.discovery, "discovery"),
        "validation": load_measured(a.validation, "validation"),
        "third": load_measured(a.third, "third"),
    }
    env = {name: attach_environment(d, a.bio_dir, a.srad_dir) for name, d in raw.items()}

    pca = fit_discovery_pca(env["discovery"])
    transformed = {name: transform_environment(d, pca) for name, d in env.items()}

    disc, ddf = run_cohort(transformed["discovery"], "discovery", pca)
    val, vdf = run_cohort(transformed["validation"], "validation", pca)
    third, tdf = run_cohort(transformed["third"], "third", pca)

    bio5 = json.loads(a.bio5_result.read_text(encoding="utf-8"))
    old = bio5["cohorts"]

    k95 = pca["k"]["pc95"]
    result = {
        "schema": "fcp_macroclimate_PC_IBD_IBE_posthoc_v1",
        "date_jst": "2026-10-07",
        "status": "complete_fixed_multivariate_macroclimate_decomposition",
        "role": "posthoc_outcome_independent_macroclimate_representation_after_BIO5_only_analysis",
        "confirmatory_decisions_changed": False,
        "environment": {
            "variables": ENV_NAMES,
            "n_variables": len(ENV_NAMES),
            "fit_cohort": "discovery only",
            "fit_rows": pca["n_fit_rows"],
            "standardization": "discovery mean and population SD transported unchanged",
            "primary_rule": "smallest K with >=95% cumulative discovery environmental variance",
            "k_pc80": pca["k"]["pc80"],
            "k_pc90": pca["k"]["pc90"],
            "k_pc95": k95,
            "primary_explained_variance": float(pca["cumulative"][k95 - 1]),
            "explained_variance_ratio": [float(x) for x in pca["explained_ratio"]],
            "means": {name: float(x) for name, x in zip(ENV_NAMES, pca["mean"])},
            "sds": {name: float(x) for name, x in zip(ENV_NAMES, pca["sd"])},
            "primary_PC_loadings": top_loadings(pca, k95),
        },
        "discovery": disc,
        "validation": val,
        "third": third,
        "replication": {
            "macroclimate_ENV_500_plus_500": bool(
                disc["primary_pc95"]["ENV_IBE_like"]["supported"]
                and val["primary_pc95"]["ENV_IBE_like"]["supported"]
            ),
            "macroclimate_ENV_all_three": bool(
                all(x["primary_pc95"]["ENV_IBE_like"]["supported"] for x in (disc, val, third))
            ),
            "GEO_500_plus_500": bool(
                disc["primary_pc95"]["GEO_IBD_like"]["supported"]
                and val["primary_pc95"]["GEO_IBD_like"]["supported"]
            ),
            "GEO_all_three": bool(
                all(x["primary_pc95"]["GEO_IBD_like"]["supported"] for x in (disc, val, third))
            ),
            "macroclimate_ENV_stronger_than_GEO_500_plus_500": bool(
                disc["primary_pc95"]["relative_balance"]["ENV_stronger_supported"]
                and val["primary_pc95"]["relative_balance"]["ENV_stronger_supported"]
            ),
            "macroclimate_ENV_stronger_than_GEO_all_three": bool(
                all(x["primary_pc95"]["relative_balance"]["ENV_stronger_supported"] for x in (disc, val, third))
            ),
        },
        "BIO5_comparison": {
            cohort: {
                "BIO5_ENV_mean_partial_rho": float(old[cohort]["IBE_like"]["mean_partial_rho"]),
                "BIO5_GEO_mean_partial_rho": float(old[cohort]["IBD_like"]["mean_partial_rho"]),
                "macroclimate_PC95_ENV_mean_partial_rho": float(cur["primary_pc95"]["ENV_IBE_like"]["mean_partial_rho"]),
                "macroclimate_PC95_GEO_mean_partial_rho": float(cur["primary_pc95"]["GEO_IBD_like"]["mean_partial_rho"]),
                "macroclimate_minus_BIO5_ENV": float(
                    cur["primary_pc95"]["ENV_IBE_like"]["mean_partial_rho"] - old[cohort]["IBE_like"]["mean_partial_rho"]
                ),
            }
            for cohort, cur in (("discovery", disc), ("validation", val), ("third", third))
        },
        "interpretation": {
            "hard_nonclaims": [
                "macroclimate is not total environment",
                "does not include soil chemistry, direct UV-B, pollinators, herbivores, demography or genetic structure",
                "not genetic isolation by environment",
                "does not establish causal climate selection",
                "does not measure fitness or local adaptation",
                "cannot distinguish genetic differentiation from plasticity",
                "post hoc analysis cannot alter frozen confirmatory decisions",
            ]
        },
    }

    out = a.outdir
    out.mkdir(parents=True, exist_ok=True)
    pd.concat([ddf, vdf, tdf], ignore_index=True).to_csv(out / "species_macroclimate_PC_IBD_IBE_metrics.csv", index=False)
    (out / "result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

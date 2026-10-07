#!/usr/bin/env python3
"""Post hoc adaptive-sorting analysis of FCP discovery/validation spatial structure.

Two questions are tested with the frozen 500+500 high-depth cohorts:

1) Do species sampled across more heterogeneous climates have greater
   species-wide sampled four-state colour diversity D after accounting for
   sampled geographic span, classifiable-row count and observer count?

2) Within species, are flower-colour differences larger between
   environmentally different locations even after geographic distance is
   partialled out?

The second question uses matched within-species vertex permutations of complete
colour vectors, preserving each species' coordinates, environmental values,
colour-vector multiset and D. Discovery and validation are reported separately.

This is post-outcome ecological diagnosis. It cannot change frozen H1/H2 or
D-spatial decisions and cannot establish fitness-mediated local adaptation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from scipy.stats import rankdata, spearmanr, wilcoxon

EARTH_RADIUS_KM = 6371.0088
MORPHS = ["white", "yellow_orange", "red_pink", "blue_purple"]
COLOURS = ["colour_white", "colour_yellow_orange", "colour_red_pink", "colour_blue_purple"]
ENV_NAMES = ["bio5", "bio14", "srad"]
SORT_METRICS = ["multivariate", "bio5", "bio14", "srad"]
DISC_SHA256 = "ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4"
VAL_SHA256 = "0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6"
PERMUTATIONS = 199
HET_PERMUTATIONS = 999
MASTER_SEED = 2026100717


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def as_bool(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False).astype(bool)
    return s.fillna("").astype(str).str.strip().str.casefold().isin({"true", "1", "yes", "y"})


def sample_raster(path: Path, lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    pts = [(float(x), float(y)) for x, y in zip(lon, lat)]
    out = np.full(len(pts), np.nan, float)
    good = [i for i, (x, y) in enumerate(pts) if np.isfinite(x) and np.isfinite(y)]
    if not good:
        return out
    with rasterio.open(path) as src:
        values = list(src.sample([pts[i] for i in good]))
        for i, raw in zip(good, values):
            value = float(raw[0])
            if src.nodata is not None and np.isclose(value, src.nodata):
                value = np.nan
            out[i] = value
    return out


def sample_srad_mean(paths: list[Path], lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    stack = np.column_stack([sample_raster(p, lon, lat) for p in paths])
    with np.errstate(invalid="ignore"):
        return np.nanmean(stack, axis=1)


def pairwise_geo_km(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    latr = np.deg2rad(np.asarray(lat, float))
    lonr = np.deg2rad(np.asarray(lon, float))
    c = np.cos(latr)
    xyz = np.column_stack([c * np.cos(lonr), c * np.sin(lonr), np.sin(latr)])
    dot = np.clip(xyz @ xyz.T, -1.0, 1.0)
    return np.arccos(dot) * EARTH_RADIUS_KM


def pairwise_euclidean(x: np.ndarray) -> np.ndarray:
    d = x[:, None, :] - x[None, :, :]
    return np.sqrt(np.sum(d * d, axis=2))


def pairwise_abs(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, float)
    return np.abs(x[:, None] - x[None, :])


def pairwise_jsd_matrix(prob: np.ndarray) -> np.ndarray:
    p = np.asarray(prob, float)
    if p.ndim != 2 or np.any(~np.isfinite(p)) or np.any(p < 0):
        raise ValueError("invalid colour probability rows")
    mass = p.sum(axis=1)
    if np.any(mass <= 0):
        raise ValueError("zero-mass colour row")
    p = p / mass[:, None]
    a = p[:, None, :]
    b = p[None, :, :]
    m = 0.5 * (a + b)
    with np.errstate(divide="ignore", invalid="ignore"):
        ka = np.where(a > 0, a * np.log2(a / m), 0.0).sum(axis=2)
        kb = np.where(b > 0, b * np.log2(b / m), 0.0).sum(axis=2)
    return np.clip(0.5 * (ka + kb), 0.0, 1.0)


def residualize_rank(y: np.ndarray, covariates: list[np.ndarray]) -> np.ndarray:
    ry = rankdata(np.asarray(y, float), method="average")
    cols = [np.ones(len(ry), float)]
    for c in covariates:
        cols.append(rankdata(np.asarray(c, float), method="average"))
    X = np.column_stack(cols)
    beta, *_ = np.linalg.lstsq(X, ry, rcond=None)
    return ry - X @ beta


def partial_rank(x: np.ndarray, y: np.ndarray, covariates: list[np.ndarray]) -> float:
    rx = residualize_rank(x, covariates)
    ry = residualize_rank(y, covariates)
    nx = float(np.linalg.norm(rx))
    ny = float(np.linalg.norm(ry))
    if nx <= 1e-14 or ny <= 1e-14:
        return 0.0
    return float(np.dot(rx, ry) / (nx * ny))


def stable_seed(*parts: object) -> int:
    payload = "|".join(map(str, (MASTER_SEED, *parts))).encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little", signed=False)


def D_from_morphs(morph: pd.Series) -> float:
    counts = morph.astype(str).value_counts()
    p = np.array([counts.get(m, 0) / len(morph) for m in MORPHS], float)
    return float(1.0 - np.sum(p * p))


def prepare(path: Path, cohort: str, bio5: Path, bio14: Path, srad_paths: list[Path]) -> pd.DataFrame:
    d = pd.read_csv(path, low_memory=False)
    required = {
        "inat_taxon_id", "species", "photo_id", "observer_id", "latitude", "longitude",
        "morph", "global_classifiable", *COLOURS,
    }
    missing = sorted(required - set(d.columns))
    if missing:
        raise RuntimeError(f"{cohort}: missing columns {missing}")

    keep = as_bool(d["global_classifiable"]) & d["morph"].isin(MORPHS)
    d = d.loc[keep].copy()
    counts = d.groupby("inat_taxon_id").size()
    d = d.loc[d["inat_taxon_id"].isin(counts[counts >= 40].index)].copy()
    d = d.sort_values(["inat_taxon_id", "photo_id"], kind="stable").reset_index(drop=True)

    lon = pd.to_numeric(d["longitude"], errors="coerce").to_numpy(float)
    lat = pd.to_numeric(d["latitude"], errors="coerce").to_numpy(float)
    d["bio5"] = sample_raster(bio5, lon, lat)
    d["bio14"] = sample_raster(bio14, lon, lat)
    d["srad"] = sample_srad_mean(srad_paths, lon, lat)
    d["cohort"] = cohort
    return d


def discovery_scaling(d: pd.DataFrame) -> dict:
    out = {}
    for name in ENV_NAMES:
        x = pd.to_numeric(d[name], errors="coerce").to_numpy(float)
        x = x[np.isfinite(x)]
        mu = float(np.mean(x))
        sd = float(np.std(x, ddof=0))
        if not np.isfinite(sd) or sd <= 0:
            raise RuntimeError(f"discovery scaling degenerate for {name}")
        out[name] = {"mean": mu, "sd": sd}
    return out


def apply_scaling(d: pd.DataFrame, scaling: dict) -> pd.DataFrame:
    x = d.copy()
    for name in ENV_NAMES:
        x[f"z_{name}"] = (pd.to_numeric(x[name], errors="coerce") - scaling[name]["mean"]) / scaling[name]["sd"]
    return x


def species_environment_table(d: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    rows = []
    sorting_cache = {}
    for taxon, g in d.groupby("inat_taxon_id", sort=True):
        g = g.dropna(subset=["z_bio5", "z_bio14", "z_srad", "latitude", "longitude", *COLOURS]).copy()
        if len(g) < 40:
            continue
        g = g.sort_values("photo_id", kind="stable").reset_index(drop=True)
        env = g[["z_bio5", "z_bio14", "z_srad"]].to_numpy(float)
        lat = g["latitude"].to_numpy(float)
        lon = g["longitude"].to_numpy(float)
        geo = pairwise_geo_km(lat, lon)
        env_multi = pairwise_euclidean(env)
        u, v = np.triu_indices(len(g), k=1)
        geo_pair = geo[u, v]
        env_multi_pair = env_multi[u, v]

        row = {
            "cohort": str(g["cohort"].iloc[0]),
            "inat_taxon_id": int(taxon),
            "species": str(g["species"].iloc[0]),
            "n_classifiable": int(len(g)),
            "n_observers": int(g["observer_id"].astype(str).nunique()),
            "D": D_from_morphs(g["morph"]),
            "sampled_span_km": float(np.max(geo_pair)),
            "env_heterogeneity_multivariate": float(np.median(env_multi_pair)),
            "env_heterogeneity_bio5": float(np.std(env[:, 0], ddof=0)),
            "env_heterogeneity_bio14": float(np.std(env[:, 1], ddof=0)),
            "env_heterogeneity_srad": float(np.std(env[:, 2], ddof=0)),
        }

        colours = g[COLOURS].to_numpy(float)
        colour_jsd = pairwise_jsd_matrix(colours)
        colour_pair = colour_jsd[u, v]
        colour_nonconstant = bool(np.ptp(colour_pair) > 1e-15)
        row["colour_distance_nondegenerate"] = colour_nonconstant

        env_mats = {
            "multivariate": env_multi,
            "bio5": pairwise_abs(env[:, 0]),
            "bio14": pairwise_abs(env[:, 1]),
            "srad": pairwise_abs(env[:, 2]),
        }

        if colour_nonconstant:
            colour_rank_matrix = np.zeros_like(colour_jsd)
            ranks = rankdata(colour_pair, method="average")
            colour_rank_matrix[u, v] = ranks
            colour_rank_matrix[v, u] = ranks
            geo_rank = rankdata(geo_pair, method="average")
            gc = geo_rank - geo_rank.mean()
            gnorm2 = float(np.dot(gc, gc))

            metric_state = {}
            for metric, mat in env_mats.items():
                ep = mat[u, v]
                if np.ptp(ep) <= 1e-15:
                    row[f"sorting_partial_rho_{metric}"] = np.nan
                    continue
                er = rankdata(ep, method="average")
                ec = er - er.mean()
                beta = 0.0 if gnorm2 <= 1e-15 else float(np.dot(gc, ec) / gnorm2)
                eres = ec - beta * gc
                enorm = float(np.linalg.norm(eres))
                if enorm <= 1e-14:
                    row[f"sorting_partial_rho_{metric}"] = np.nan
                    continue

                yc = ranks - ranks.mean()
                ybeta = 0.0 if gnorm2 <= 1e-15 else float(np.dot(gc, yc) / gnorm2)
                yres = yc - ybeta * gc
                ynorm = float(np.linalg.norm(yres))
                obs = 0.0 if ynorm <= 1e-14 else float(np.dot(eres, yres) / (enorm * ynorm))
                row[f"sorting_partial_rho_{metric}"] = obs
                metric_state[metric] = {"eres": eres, "enorm": enorm}

            if metric_state:
                sorting_cache[int(taxon)] = {
                    "species": row["species"],
                    "n": len(g),
                    "u": u,
                    "v": v,
                    "colour_rank_matrix": colour_rank_matrix,
                    "geo_centered_rank": gc,
                    "geo_norm2": gnorm2,
                    "metrics": metric_state,
                }
        rows.append(row)
    return pd.DataFrame(rows).sort_values("inat_taxon_id").reset_index(drop=True), sorting_cache


def heterogeneity_test(species: pd.DataFrame, column: str, cohort: str) -> dict:
    x = species.dropna(subset=["D", column, "sampled_span_km", "n_classifiable", "n_observers"]).copy()
    cov = [
        np.log1p(x["sampled_span_km"].to_numpy(float)),
        x["n_classifiable"].to_numpy(float),
        x["n_observers"].to_numpy(float),
    ]
    dres = residualize_rank(x["D"].to_numpy(float), cov)
    eres = residualize_rank(x[column].to_numpy(float), cov)
    denom = float(np.linalg.norm(dres) * np.linalg.norm(eres))
    observed = 0.0 if denom <= 1e-14 else float(np.dot(dres, eres) / denom)

    null = np.empty(HET_PERMUTATIONS, float)
    rng = np.random.default_rng(stable_seed("heterogeneity", cohort, column))
    for i in range(HET_PERMUTATIONS):
        p = rng.permutation(len(eres))
        ep = eres[p]
        den = float(np.linalg.norm(dres) * np.linalg.norm(ep))
        null[i] = 0.0 if den <= 1e-14 else float(np.dot(dres, ep) / den)
    p_upper = float((1 + np.count_nonzero(null >= observed)) / (HET_PERMUTATIONS + 1))
    return {
        "n_species": int(len(x)),
        "partial_spearman": observed,
        "p_upper_permutation": p_upper,
        "null_mean": float(np.mean(null)),
        "null_q025": float(np.quantile(null, 0.025)),
        "null_q975": float(np.quantile(null, 0.975)),
        "positive_supported_at_0_05": bool(observed > 0 and p_upper < 0.05),
    }


def sorting_test(species: pd.DataFrame, cache: dict, metric: str, cohort: str) -> tuple[dict, pd.DataFrame]:
    eligible = species.dropna(subset=[f"sorting_partial_rho_{metric}"]).copy()
    tids = eligible["inat_taxon_id"].astype(int).tolist()
    observed_values = eligible[f"sorting_partial_rho_{metric}"].to_numpy(float)
    if not tids:
        return {"n_species": 0, "evaluable": False}, pd.DataFrame()

    null_sum = np.zeros(PERMUTATIONS, float)
    for tid in tids:
        state = cache[tid]
        m = state["metrics"][metric]
        rng = np.random.default_rng(stable_seed("sorting", cohort, metric, tid))
        n = state["n"]
        u, v = state["u"], state["v"]
        gc = state["geo_centered_rank"]
        gnorm2 = state["geo_norm2"]
        crm = state["colour_rank_matrix"]
        eres = m["eres"]
        enorm = m["enorm"]
        per = np.empty(PERMUTATIONS, float)
        for i in range(PERMUTATIONS):
            p = rng.permutation(n)
            yr = crm[p[u], p[v]]
            yc = yr - yr.mean()
            beta = 0.0 if gnorm2 <= 1e-15 else float(np.dot(gc, yc) / gnorm2)
            yres = yc - beta * gc
            ynorm = float(np.linalg.norm(yres))
            per[i] = 0.0 if ynorm <= 1e-14 else float(np.dot(eres, yres) / (enorm * ynorm))
        null_sum += per

    null_mean = null_sum / len(tids)
    observed_mean = float(np.mean(observed_values))
    p_upper = float((1 + np.count_nonzero(null_mean >= observed_mean)) / (PERMUTATIONS + 1))
    try:
        w = wilcoxon(observed_values, zero_method="wilcox", alternative="greater")
        wilcoxon_p = float(w.pvalue)
    except ValueError:
        wilcoxon_p = 1.0

    d_assoc = float(spearmanr(
        eligible["D"].to_numpy(float),
        observed_values,
    ).statistic)

    result = {
        "n_species": int(len(tids)),
        "observed_equal_species_mean_partial_rho": observed_mean,
        "observed_median_species_partial_rho": float(np.median(observed_values)),
        "positive_species_fraction": float(np.mean(observed_values > 0)),
        "wilcoxon_greater_p": wilcoxon_p,
        "matched_vertex_null_p_upper": p_upper,
        "null_mean": float(np.mean(null_mean)),
        "null_q025": float(np.quantile(null_mean, 0.025)),
        "null_q975": float(np.quantile(null_mean, 0.975)),
        "rho_D_vs_species_sorting_strength": d_assoc,
        "positive_supported_at_0_05": bool(observed_mean > 0 and p_upper < 0.05),
    }
    null_df = pd.DataFrame({
        "cohort": cohort,
        "metric": metric,
        "permutation_index": np.arange(PERMUTATIONS),
        "equal_species_mean_partial_rho": null_mean,
    })
    return result, null_df


def analyze_cohort(d: pd.DataFrame, cohort: str) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    species, cache = species_environment_table(d)
    heterogeneity = {
        "multivariate": heterogeneity_test(species, "env_heterogeneity_multivariate", cohort),
        "bio5": heterogeneity_test(species, "env_heterogeneity_bio5", cohort),
        "bio14": heterogeneity_test(species, "env_heterogeneity_bio14", cohort),
        "srad": heterogeneity_test(species, "env_heterogeneity_srad", cohort),
    }
    sorting = {}
    nulls = []
    for metric in SORT_METRICS:
        res, ndf = sorting_test(species, cache, metric, cohort)
        sorting[metric] = res
        if not ndf.empty:
            nulls.append(ndf)
    return {
        "species": int(len(species)),
        "environmental_heterogeneity_vs_D": heterogeneity,
        "environment_colour_sorting_given_geography": sorting,
    }, species, (pd.concat(nulls, ignore_index=True) if nulls else pd.DataFrame())


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--discovery", type=Path, required=True)
    ap.add_argument("--validation", type=Path, required=True)
    ap.add_argument("--bio5", type=Path, required=True)
    ap.add_argument("--bio14", type=Path, required=True)
    ap.add_argument("--srad-dir", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()

    if file_sha256(args.discovery) != DISC_SHA256:
        raise RuntimeError("discovery measured table SHA256 mismatch")
    if file_sha256(args.validation) != VAL_SHA256:
        raise RuntimeError("validation measured table SHA256 mismatch")

    srad_paths = [args.srad_dir / f"wc2.1_10m_srad_{m:02d}.tif" for m in range(1, 13)]
    if not all(p.exists() for p in srad_paths):
        raise RuntimeError("missing monthly SRAD rasters")

    disc0 = prepare(args.discovery, "discovery", args.bio5, args.bio14, srad_paths)
    val0 = prepare(args.validation, "validation", args.bio5, args.bio14, srad_paths)
    scaling = discovery_scaling(disc0)
    disc = apply_scaling(disc0, scaling)
    val = apply_scaling(val0, scaling)

    disc_result, disc_species, disc_null = analyze_cohort(disc, "discovery")
    val_result, val_species, val_null = analyze_cohort(val, "validation")

    het_rep = bool(
        disc_result["environmental_heterogeneity_vs_D"]["multivariate"]["positive_supported_at_0_05"]
        and val_result["environmental_heterogeneity_vs_D"]["multivariate"]["positive_supported_at_0_05"]
    )
    sort_rep = bool(
        disc_result["environment_colour_sorting_given_geography"]["multivariate"]["positive_supported_at_0_05"]
        and val_result["environment_colour_sorting_given_geography"]["multivariate"]["positive_supported_at_0_05"]
    )

    if het_rep and sort_rep:
        verdict = "POSTHOC_ENVIRONMENTAL_SORTING_PATTERN_REPLICATED_ON_BOTH_TESTS"
    elif sort_rep:
        verdict = "POSTHOC_WITHIN_SPECIES_ENVIRONMENTAL_SORTING_REPLICATED_BUT_HETEROGENEITY_D_NOT_REPLICATED"
    elif het_rep:
        verdict = "POSTHOC_ENVIRONMENTAL_HETEROGENEITY_D_REPLICATED_BUT_WITHIN_SPECIES_SORTING_NOT_REPLICATED"
    else:
        verdict = "POSTHOC_ADAPTIVE_SORTING_PATTERN_NOT_JOINTLY_REPLICATED"

    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)
    pd.concat([disc_species, val_species], ignore_index=True).to_csv(out / "species_environment_spatial_metrics.csv", index=False)
    pd.concat([disc_null, val_null], ignore_index=True).to_csv(out / "sorting_vertex_nulls.csv", index=False)

    result = {
        "schema": "fcp_adaptive_spatial_sorting_posthoc_v1",
        "date_jst": "2026-10-07",
        "status": "complete_posthoc_ecological_diagnostic",
        "role": "post_outcome_test_of_environmental_sorting_consistent_with_adaptive_or_plastic_spatial_response",
        "confirmatory_decisions_changed": False,
        "source_commit": "5142f7951af0dde5364bb047a566d67e8c479e51",
        "discovery_sha256": DISC_SHA256,
        "validation_sha256": VAL_SHA256,
        "environmental_variables": [
            "WorldClim 2.1 BIO5 maximum temperature of warmest month",
            "WorldClim 2.1 BIO14 precipitation of driest month",
            "WorldClim 2.1 mean monthly solar radiation",
        ],
        "environment_scaling": {
            "rule": "mean/SD estimated once from discovery classifiable rows and transported unchanged to validation",
            "values": scaling,
        },
        "permutations": {
            "within_species_colour_vertex": PERMUTATIONS,
            "cross_species_heterogeneity_partial_rank": HET_PERMUTATIONS,
        },
        "discovery": disc_result,
        "validation": val_result,
        "joint_posthoc_verdict": verdict,
        "interpretation": {
            "if_supported": "Replicated positive environmental-distance/colour-distance coupling after geographic distance control would support ecological sorting of intraspecific flower-colour variation. A replicated environmental-heterogeneity/D association would further support the hypothesis that heterogeneous environments contribute to species-wide colour diversity.",
            "hard_nonclaims": [
                "does not measure fitness",
                "does not demonstrate local adaptation",
                "does not distinguish genetic differentiation from phenotypic plasticity",
                "does not identify the causal environmental variable when predictors covary",
                "does not establish universality beyond the same iNaturalist measurement system",
                "post hoc analysis cannot upgrade frozen prospective H2",
            ],
        },
    }
    (out / "result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

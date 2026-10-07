#!/usr/bin/env python3
"""Fixed third-cohort transport of two selected FCP environmental signatures."""
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
COARSE = ["colour_white", "colour_yellow_orange", "colour_red_pink", "colour_blue_purple"]
THIRD_SHA256 = "57630fc9f281bce94a0c40a70aaf7bce879dde93d6154175adcd021e8f5c1186"
VERTEX_PERMUTATIONS = 199
HET_PERMUTATIONS = 999
MASTER_SEED = 2026100721


def sha256(path: Path) -> str:
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
    with rasterio.open(path) as src:
        vals = list(src.sample([pts[i] for i in good]))
        for i, raw in zip(good, vals):
            v = float(raw[0])
            if src.nodata is not None and np.isclose(v, src.nodata):
                v = np.nan
            out[i] = v
    return out


def sample_srad_mean(paths: list[Path], lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    a = np.column_stack([sample_raster(p, lon, lat) for p in paths])
    with np.errstate(invalid="ignore"):
        return np.nanmean(a, axis=1)


def pairwise_geo_km(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    latr = np.deg2rad(np.asarray(lat, float))
    lonr = np.deg2rad(np.asarray(lon, float))
    c = np.cos(latr)
    xyz = np.column_stack([c * np.cos(lonr), c * np.sin(lonr), np.sin(latr)])
    dot = np.clip(xyz @ xyz.T, -1.0, 1.0)
    return np.arccos(dot) * EARTH_RADIUS_KM


def pairwise_jsd(prob: np.ndarray) -> np.ndarray:
    p = np.asarray(prob, float)
    if p.ndim != 2 or np.any(~np.isfinite(p)) or np.any(p < 0):
        raise ValueError("invalid composition")
    mass = p.sum(axis=1)
    if np.any(mass <= 0):
        raise ValueError("zero-mass row")
    p = p / mass[:, None]
    a, b = p[:, None, :], p[None, :, :]
    m = 0.5 * (a + b)
    with np.errstate(divide="ignore", invalid="ignore"):
        ka = np.where(a > 0, a * np.log2(a / m), 0.0).sum(axis=2)
        kb = np.where(b > 0, b * np.log2(b / m), 0.0).sum(axis=2)
    return np.clip(0.5 * (ka + kb), 0.0, 1.0)


def D_from_morphs(morph: pd.Series) -> float:
    c = morph.astype(str).value_counts()
    p = np.array([c.get(m, 0) / len(morph) for m in MORPHS], float)
    return float(1.0 - np.sum(p * p))


def residualize_rank(y: np.ndarray, covariates: list[np.ndarray]) -> np.ndarray:
    ry = rankdata(np.asarray(y, float), method="average")
    X = np.column_stack([np.ones(len(ry)), *[
        rankdata(np.asarray(c, float), method="average") for c in covariates
    ]])
    beta, *_ = np.linalg.lstsq(X, ry, rcond=None)
    return ry - X @ beta


def stable_seed(*parts: object) -> int:
    s = "|".join(map(str, (MASTER_SEED, *parts))).encode()
    return int.from_bytes(hashlib.sha256(s).digest()[:8], "little")


def holm_two(p1: float, p2: float) -> tuple[float, float]:
    if p1 <= p2:
        return min(1.0, 2 * p1), max(min(1.0, 2 * p1), p2)
    a2, a1 = holm_two(p2, p1)
    return a1, a2


def target_a(species: pd.DataFrame) -> dict:
    x = species.dropna(subset=["D", "srad_heterogeneity", "sampled_span_km", "n_classifiable", "n_observers"]).copy()
    cov = [
        np.log1p(x["sampled_span_km"].to_numpy(float)),
        x["n_classifiable"].to_numpy(float),
        x["n_observers"].to_numpy(float),
    ]
    dres = residualize_rank(x["D"].to_numpy(float), cov)
    sres = residualize_rank(x["srad_heterogeneity"].to_numpy(float), cov)
    denom = float(np.linalg.norm(dres) * np.linalg.norm(sres))
    obs = 0.0 if denom <= 1e-14 else float(np.dot(dres, sres) / denom)
    rng = np.random.default_rng(stable_seed("target_a"))
    null = np.empty(HET_PERMUTATIONS, float)
    for i in range(HET_PERMUTATIONS):
        sp = sres[rng.permutation(len(sres))]
        den = float(np.linalg.norm(dres) * np.linalg.norm(sp))
        null[i] = 0.0 if den <= 1e-14 else float(np.dot(dres, sp) / den)
    p = float((1 + np.count_nonzero(null >= obs)) / (HET_PERMUTATIONS + 1))
    return {
        "n_species": int(len(x)),
        "partial_spearman": obs,
        "p_upper": p,
        "null_mean": float(null.mean()),
        "null_q025": float(np.quantile(null, .025)),
        "null_q975": float(np.quantile(null, .975)),
        "supported_raw_0_05": bool(obs > 0 and p < .05),
    }


def partial_sorting_with_null(distance: np.ndarray, geo_pair: np.ndarray, bio_pair: np.ndarray,
                              *, seed_parts: tuple[object, ...]) -> tuple[float, np.ndarray]:
    n = distance.shape[0]
    u, v = np.triu_indices(n, k=1)
    y = distance[u, v]
    if np.ptp(y) <= 1e-15 or np.ptp(bio_pair) <= 1e-15:
        return np.nan, np.full(VERTEX_PERMUTATIONS, np.nan)

    gr = rankdata(geo_pair, method="average")
    gc = gr - gr.mean()
    g2 = float(np.dot(gc, gc))

    xr = rankdata(bio_pair, method="average")
    xc = xr - xr.mean()
    xb = 0.0 if g2 <= 1e-15 else float(np.dot(gc, xc) / g2)
    xres = xc - xb * gc
    xnorm = float(np.linalg.norm(xres))
    if xnorm <= 1e-14:
        return np.nan, np.full(VERTEX_PERMUTATIONS, np.nan)

    yr = rankdata(y, method="average")
    yc = yr - yr.mean()
    yb = 0.0 if g2 <= 1e-15 else float(np.dot(gc, yc) / g2)
    yres = yc - yb * gc
    ynorm = float(np.linalg.norm(yres))
    obs = 0.0 if ynorm <= 1e-14 else float(np.dot(xres, yres) / (xnorm * ynorm))

    rank_matrix = np.zeros_like(distance, dtype=float)
    rank_matrix[u, v] = yr
    rank_matrix[v, u] = yr

    rng = np.random.default_rng(stable_seed(*seed_parts))
    perms = np.stack([rng.permutation(n) for _ in range(VERTEX_PERMUTATIONS)])
    out = np.empty(VERTEX_PERMUTATIONS, float)
    batch = 32
    for start in range(0, VERTEX_PERMUTATIONS, batch):
        stop = min(VERTEX_PERMUTATIONS, start + batch)
        pp = perms[start:stop]
        vals = rank_matrix[pp[:, u], pp[:, v]]
        center = vals - vals.mean(axis=1, keepdims=True)
        beta = np.zeros(stop - start) if g2 <= 1e-15 else (center @ gc) / g2
        resid = center - beta[:, None] * gc
        norm = np.linalg.norm(resid, axis=1)
        num = resid @ xres
        out[start:stop] = np.where(norm <= 1e-14, 0.0, num / (norm * xnorm))
    return obs, out


def target_b(species_rows: list[dict], nulls: dict[str, list[np.ndarray]]) -> dict:
    out = {}
    for rep in ["coarse_four_state", "flower_palette", "flower_minus_background"]:
        vals = np.array([r.get(rep, np.nan) for r in species_rows], float)
        good = np.isfinite(vals)
        arrays = [a for a, ok in zip(nulls[rep], good) if ok and np.isfinite(a).all()]
        if good.sum() == 0 or not arrays:
            out[rep] = {"evaluable": False, "n_species": 0}
            continue
        nmat = np.vstack(arrays)
        null_mean = nmat.mean(axis=0)
        obs = float(vals[good].mean())
        p = float((1 + np.count_nonzero(null_mean >= obs)) / (VERTEX_PERMUTATIONS + 1))
        out[rep] = {
            "evaluable": True,
            "n_species": int(good.sum()),
            "observed_equal_species_mean_partial_rho": obs,
            "observed_median_species_partial_rho": float(np.median(vals[good])),
            "positive_species_fraction": float(np.mean(vals[good] > 0)),
            "matched_vertex_p_upper": p,
            "null_mean": float(null_mean.mean()),
            "null_q025": float(np.quantile(null_mean, .025)),
            "null_q975": float(np.quantile(null_mean, .975)),
            "supported_raw_0_05": bool(obs > 0 and p < .05),
        }
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--third", type=Path, required=True)
    ap.add_argument("--bio5", type=Path, required=True)
    ap.add_argument("--srad-dir", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()

    if sha256(args.third) != THIRD_SHA256:
        raise RuntimeError("third-cohort measured table SHA256 mismatch")

    d = pd.read_csv(args.third, low_memory=False)
    required = {"inat_taxon_id", "species", "photo_id", "observer_id", "latitude", "longitude",
                "morph", "global_classifiable", *COARSE}
    miss = sorted(required - set(d.columns))
    if miss:
        raise RuntimeError(f"missing required columns: {miss}")

    keep = as_bool(d["global_classifiable"]) & d["morph"].isin(MORPHS)
    d = d.loc[keep].copy()
    counts = d.groupby("inat_taxon_id").size()
    d = d.loc[d["inat_taxon_id"].isin(counts[counts >= 40].index)].copy()
    d = d.sort_values(["inat_taxon_id", "photo_id"], kind="stable").reset_index(drop=True)

    lon = pd.to_numeric(d["longitude"], errors="coerce").to_numpy(float)
    lat = pd.to_numeric(d["latitude"], errors="coerce").to_numpy(float)
    d["bio5"] = sample_raster(args.bio5, lon, lat)
    srad_paths = [args.srad_dir / f"wc2.1_10m_srad_{m:02d}.tif" for m in range(1, 13)]
    d["srad"] = sample_srad_mean(srad_paths, lon, lat)

    flower_cols = sorted([c for c in d.columns if c.startswith("palette_count_")])
    bg_cols = [f"background_{c}" for c in flower_cols]
    have_background = bool(flower_cols) and all(c in d.columns for c in bg_cols)

    species_table = []
    sorting_rows = []
    nulls = {k: [] for k in ["coarse_four_state", "flower_palette", "flower_minus_background"]}

    for taxon, g in d.groupby("inat_taxon_id", sort=True):
        g = g.dropna(subset=["bio5", "srad", "latitude", "longitude", *COARSE]).copy()
        if len(g) < 40:
            continue
        g = g.sort_values("photo_id", kind="stable").reset_index(drop=True)
        lat = g["latitude"].to_numpy(float)
        lon = g["longitude"].to_numpy(float)
        geo = pairwise_geo_km(lat, lon)
        u, v = np.triu_indices(len(g), k=1)
        geo_pair = geo[u, v]
        bio = g["bio5"].to_numpy(float)
        bio_pair = np.abs(bio[u] - bio[v])

        row = {
            "inat_taxon_id": int(taxon),
            "species": str(g["species"].iloc[0]),
            "n_classifiable": int(len(g)),
            "n_observers": int(g["observer_id"].astype(str).nunique()),
            "D": D_from_morphs(g["morph"]),
            "sampled_span_km": float(np.max(geo_pair)),
            "srad_heterogeneity": float(np.std(g["srad"].to_numpy(float), ddof=0)),
        }

        srow = {"inat_taxon_id": int(taxon), "species": row["species"]}

        coarse = pairwise_jsd(g[COARSE].to_numpy(float))
        obs, nul = partial_sorting_with_null(coarse, geo_pair, bio_pair, seed_parts=("coarse", int(taxon)))
        srow["coarse_four_state"] = obs
        nulls["coarse_four_state"].append(nul)

        if flower_cols:
            flower = pairwise_jsd(g[flower_cols].to_numpy(float))
            obs, nul = partial_sorting_with_null(flower, geo_pair, bio_pair, seed_parts=("flower", int(taxon)))
            srow["flower_palette"] = obs
            nulls["flower_palette"].append(nul)
        else:
            srow["flower_palette"] = np.nan
            nulls["flower_palette"].append(np.full(VERTEX_PERMUTATIONS, np.nan))

        if have_background:
            flower = pairwise_jsd(g[flower_cols].to_numpy(float))
            back = pairwise_jsd(g[bg_cols].to_numpy(float))
            diff = flower - back
            obs, nul = partial_sorting_with_null(diff, geo_pair, bio_pair, seed_parts=("diff", int(taxon)))
            srow["flower_minus_background"] = obs
            nulls["flower_minus_background"].append(nul)
        else:
            srow["flower_minus_background"] = np.nan
            nulls["flower_minus_background"].append(np.full(VERTEX_PERMUTATIONS, np.nan))

        species_table.append(row)
        sorting_rows.append(srow)

    species = pd.DataFrame(species_table)
    A = target_a(species)
    B = target_b(sorting_rows, nulls)

    pA = float(A["p_upper"])
    pB = float(B["coarse_four_state"]["matched_vertex_p_upper"])
    hA, hB = holm_two(pA, pB)
    A["holm_across_two_targets"] = hA
    B["coarse_four_state"]["holm_across_two_targets"] = hB

    both = bool(
        A["supported_raw_0_05"]
        and B["coarse_four_state"]["supported_raw_0_05"]
        and hA <= .05 and hB <= .05
    )
    if both:
        verdict = "BOTH_SELECTED_ENVIRONMENTAL_SIGNATURES_TRANSPORT"
    elif A["supported_raw_0_05"]:
        verdict = "SOLAR_HETEROGENEITY_D_TRANSPORTS_BUT_BIO5_SORTING_DOES_NOT"
    elif B["coarse_four_state"]["supported_raw_0_05"]:
        verdict = "BIO5_COLOUR_SORTING_TRANSPORTS_BUT_SOLAR_HETEROGENEITY_D_DOES_NOT"
    else:
        verdict = "NEITHER_SELECTED_ENVIRONMENTAL_SIGNATURE_TRANSPORTS"

    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)
    species.to_csv(out / "third_species_targetA_metrics.csv", index=False)
    pd.DataFrame(sorting_rows).to_csv(out / "third_species_targetB_metrics.csv", index=False)

    result = {
        "schema": "fcp_selected_environmental_sorting_transport_v1",
        "date_jst": "2026-10-07",
        "status": "complete_fixed_species_disjoint_posthoc_transport",
        "confirmatory_decisions_changed": False,
        "source_commit": "7e538e5c51c05a7cc47b2fcf53eea92634c8a863",
        "third_measured_sha256": THIRD_SHA256,
        "eligible_species": int(len(species)),
        "selection_from_500_plus_500_screen": {
            "target_A": "solar-radiation heterogeneity versus species-wide D",
            "target_B": "BIO5 distance versus colour distance conditional on geography",
            "screening_Holm_support": {
                "target_A_discovery": 0.012,
                "target_A_validation": 0.003,
                "target_B_discovery": 0.050,
                "target_B_validation": 0.015,
            },
        },
        "target_A_solar_heterogeneity_vs_D": A,
        "target_B_BIO5_colour_sorting": B,
        "joint_verdict": verdict,
        "interpretation": {
            "allowed_if_B_supported": "Warm-season temperature differences repeatedly organize flower-colour turnover within species beyond geographic distance, without implying a universal direction of colour response.",
            "allowed_if_A_supported": "Species sampled across more heterogeneous solar environments show greater species-wide flower-colour-state diversity after the fixed covariate adjustment.",
            "hard_nonclaims": [
                "not untouched prospective confirmation",
                "does not measure fitness or demonstrate local adaptation",
                "does not distinguish genetic differentiation from plasticity",
                "does not prove BIO5 or solar radiation are causal agents",
                "does not establish a common direction of flower-colour response",
            ],
        },
        "technical_sensitivity_status": {
            "flower_palette_columns": int(len(flower_cols)),
            "matched_background_available": bool(have_background),
            "sensitivities_cannot_rescue_primary": True,
        },
    }
    (out / "result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

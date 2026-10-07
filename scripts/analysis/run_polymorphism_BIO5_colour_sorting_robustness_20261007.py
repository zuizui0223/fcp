#!/usr/bin/env python3
"""Robustness audit for the post hoc BIO5-associated flower-colour sorting signal."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from scipy.stats import rankdata

EARTH_RADIUS_KM = 6371.0088
MORPHS = ["white", "yellow_orange", "red_pink", "blue_purple"]
BIO9 = ["white", "yellow", "orange", "red", "pink", "magenta", "purple", "blue", "bronze"]
ALL12 = BIO9 + ["green", "brown", "black"]
THIRD_FRACTIONS = [f"flower_fraction_{x}" for x in BIO9]
PERMUTATIONS = 199
MASTER_SEED = 2026100727
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
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little", signed=False)


def as_bool(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False).astype(bool)
    return s.fillna("").astype(str).str.strip().str.casefold().isin({"true", "1", "yes", "y"})


def sample_raster(path: Path, lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    pts = [(float(x), float(y)) for x, y in zip(lon, lat)]
    out = np.full(len(pts), np.nan, dtype=float)
    good = [i for i, (x, y) in enumerate(pts) if np.isfinite(x) and np.isfinite(y)]
    if not good:
        return out
    with rasterio.open(path) as src:
        vals = list(src.sample([pts[i] for i in good]))
        for i, raw in zip(good, vals):
            v = float(raw[0])
            if src.nodata is not None and np.isclose(v, src.nodata):
                v = np.nan
            out[i] = v
    return out


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
        raise ValueError("invalid composition rows")
    mass = p.sum(axis=1)
    if np.any(mass <= 0):
        raise ValueError("zero-mass composition row")
    p = p / mass[:, None]
    a = p[:, None, :]
    b = p[None, :, :]
    m = 0.5 * (a + b)
    with np.errstate(divide="ignore", invalid="ignore"):
        ka = np.where(a > 0, a * np.log2(a / m), 0.0).sum(axis=2)
        kb = np.where(b > 0, b * np.log2(b / m), 0.0).sum(axis=2)
    return np.clip(0.5 * (ka + kb), 0.0, 1.0)


def partial_rank_vectors(geo: np.ndarray, env: np.ndarray, colour: np.ndarray) -> float:
    geo = np.asarray(geo, float)
    env = np.asarray(env, float)
    colour = np.asarray(colour, float)
    if len(geo) < 4 or np.ptp(geo) <= 1e-15 or np.ptp(env) <= 1e-15 or np.ptp(colour) <= 1e-15:
        return np.nan
    gr = rankdata(geo, method="average")
    er = rankdata(env, method="average")
    cr = rankdata(colour, method="average")
    gc = gr - gr.mean()
    g2 = float(np.dot(gc, gc))
    ec = er - er.mean()
    cc = cr - cr.mean()
    eb = 0.0 if g2 <= 1e-15 else float(np.dot(gc, ec) / g2)
    cb = 0.0 if g2 <= 1e-15 else float(np.dot(gc, cc) / g2)
    eres = ec - eb * gc
    cres = cc - cb * gc
    denom = float(np.linalg.norm(eres) * np.linalg.norm(cres))
    return np.nan if denom <= 1e-14 else float(np.dot(eres, cres) / denom)


def vertex_null(
    distance: np.ndarray,
    geo_pair: np.ndarray,
    env_pair: np.ndarray,
    *,
    cohort: str,
    taxon: int,
    representation: str,
) -> tuple[float, np.ndarray]:
    n = distance.shape[0]
    u, v = np.triu_indices(n, k=1)
    obs = partial_rank_vectors(geo_pair, env_pair, distance[u, v])
    if not np.isfinite(obs):
        return np.nan, np.full(PERMUTATIONS, np.nan)

    rng = np.random.default_rng(stable_seed("vertex", cohort, taxon, representation))
    null = np.empty(PERMUTATIONS, float)
    for i in range(PERMUTATIONS):
        p = rng.permutation(n)
        null[i] = partial_rank_vectors(geo_pair, env_pair, distance[p[u], p[v]])
    return obs, null


def same_observer_pair_test(
    observers: np.ndarray,
    colour_distance: np.ndarray,
    geo_pair_all: np.ndarray,
    env_pair_all: np.ndarray,
    u: np.ndarray,
    v: np.ndarray,
    *,
    cohort: str,
    taxon: int,
) -> tuple[float, np.ndarray, int]:
    obs_id = np.asarray(observers).astype(str)
    mask = (obs_id[u] != "") & (obs_id[u] == obs_id[v])
    if int(mask.sum()) < 5:
        return np.nan, np.full(PERMUTATIONS, np.nan), int(mask.sum())
    geo = geo_pair_all[mask]
    env = env_pair_all[mask]
    colour = colour_distance[u[mask], v[mask]]
    observed = partial_rank_vectors(geo, env, colour)
    if not np.isfinite(observed):
        return np.nan, np.full(PERMUTATIONS, np.nan), int(mask.sum())

    rng = np.random.default_rng(stable_seed("same_observer_pairs", cohort, taxon))
    null = np.empty(PERMUTATIONS, float)
    for i in range(PERMUTATIONS):
        null[i] = partial_rank_vectors(geo, env, colour[rng.permutation(len(colour))])
    return observed, null, int(mask.sum())


def aggregate_species(rows: list[dict], nulls: list[np.ndarray], value_key: str) -> dict:
    vals = np.asarray([r[value_key] for r in rows], float)
    good = np.isfinite(vals)
    arrays = [a for a, ok in zip(nulls, good) if ok and np.isfinite(a).all()]
    if not arrays or int(good.sum()) == 0:
        return {"evaluable": False, "n_species": 0}
    mat = np.vstack(arrays)
    null_mean = mat.mean(axis=0)
    observed = float(vals[good].mean())
    p = float((1 + np.count_nonzero(null_mean >= observed)) / (PERMUTATIONS + 1))
    return {
        "evaluable": True,
        "n_species": int(good.sum()),
        "observed_equal_species_mean_partial_rho": observed,
        "observed_median_species_partial_rho": float(np.median(vals[good])),
        "positive_species_fraction": float(np.mean(vals[good] > 0)),
        "matched_null_p_upper": p,
        "null_mean": float(np.mean(null_mean)),
        "null_q025": float(np.quantile(null_mean, .025)),
        "null_q975": float(np.quantile(null_mean, .975)),
        "supported_at_0_05": bool(observed > 0 and p < .05),
    }


def prepare_cohort(path: Path, cohort: str, bio5_path: Path, *, third: bool) -> pd.DataFrame:
    if file_sha256(path) != SHA[cohort]:
        raise RuntimeError(f"{cohort} input SHA256 mismatch")
    d = pd.read_csv(path, low_memory=False)
    base = {"inat_taxon_id", "species", "photo_id", "observer_id", "latitude", "longitude", "morph", "global_classifiable"}
    if third:
        required = base | set(THIRD_FRACTIONS)
    else:
        required = base | {f"palette_count_{x}" for x in ALL12}
        if cohort == "validation":
            required |= {f"background_palette_count_{x}" for x in ALL12}
    missing = sorted(required - set(d.columns))
    if missing:
        raise RuntimeError(f"{cohort}: missing required columns {missing}")

    keep = as_bool(d["global_classifiable"]) & d["morph"].astype(str).isin(MORPHS)
    d = d.loc[keep].copy()
    counts = d.groupby("inat_taxon_id").size()
    d = d.loc[d["inat_taxon_id"].isin(counts[counts >= 40].index)].copy()
    d = d.sort_values(["inat_taxon_id", "photo_id"], kind="stable").reset_index(drop=True)
    d["bio5"] = sample_raster(
        bio5_path,
        pd.to_numeric(d["longitude"], errors="coerce").to_numpy(float),
        pd.to_numeric(d["latitude"], errors="coerce").to_numpy(float),
    )
    d["cohort"] = cohort
    return d


def run_cohort(d: pd.DataFrame, cohort: str, *, third: bool) -> tuple[dict, pd.DataFrame]:
    r1_rows: list[dict] = []
    r1_nulls: list[np.ndarray] = []
    r2_rows: list[dict] = []
    r2_nulls: list[np.ndarray] = []
    r3_rows: list[dict] = []
    r3_nulls: list[np.ndarray] = []

    max_obs_per_species_observer = int(d.groupby(["inat_taxon_id", "observer_id"]).size().max())

    for taxon, g0 in d.groupby("inat_taxon_id", sort=True):
        g = g0.dropna(subset=["bio5", "latitude", "longitude"]).copy()
        if len(g) < 40:
            continue
        g = g.sort_values("photo_id", kind="stable").reset_index(drop=True)

        lat = g["latitude"].to_numpy(float)
        lon = g["longitude"].to_numpy(float)
        bio = g["bio5"].to_numpy(float)
        geo = pairwise_geo_km(lat, lon)
        u, v = np.triu_indices(len(g), k=1)
        geo_pair = geo[u, v]
        bio_pair = np.abs(bio[u] - bio[v])

        if third:
            flower9 = g[THIRD_FRACTIONS].to_numpy(float)
        else:
            flower9 = g[[f"palette_count_{x}" for x in BIO9]].to_numpy(float)
        flower9_jsd = pairwise_jsd(flower9)

        obs1, nul1 = vertex_null(
            flower9_jsd, geo_pair, bio_pair,
            cohort=cohort, taxon=int(taxon), representation="flower9",
        )
        r1_rows.append({"inat_taxon_id": int(taxon), "species": str(g["species"].iloc[0]), "r1": obs1})
        r1_nulls.append(nul1)

        obs3, nul3, same_n = same_observer_pair_test(
            g["observer_id"].fillna("").to_numpy(),
            flower9_jsd, geo_pair, bio_pair, u, v,
            cohort=cohort, taxon=int(taxon),
        )
        r3_rows.append({
            "inat_taxon_id": int(taxon),
            "species": str(g["species"].iloc[0]),
            "same_observer_pairs": same_n,
            "r3": obs3,
        })
        r3_nulls.append(nul3)

        if cohort == "validation":
            flower12 = pairwise_jsd(g[[f"palette_count_{x}" for x in ALL12]].to_numpy(float))
            back12 = pairwise_jsd(g[[f"background_palette_count_{x}" for x in ALL12]].to_numpy(float))
            differential = flower12 - back12
            obs2, nul2 = vertex_null(
                differential, geo_pair, bio_pair,
                cohort=cohort, taxon=int(taxon), representation="flower_minus_background12",
            )
            r2_rows.append({"inat_taxon_id": int(taxon), "species": str(g["species"].iloc[0]), "r2": obs2})
            r2_nulls.append(nul2)

    result = {
        "eligible_species_after_bio5": int(d.dropna(subset=["bio5"]).groupby("inat_taxon_id").filter(lambda x: len(x) >= 40)["inat_taxon_id"].nunique()),
        "maximum_photos_per_species_observer": max_obs_per_species_observer,
        "R1_continuous_nine_colour": aggregate_species(r1_rows, r1_nulls, "r1"),
        "R3_same_observer_pair": aggregate_species(r3_rows, r3_nulls, "r3"),
    }
    if cohort == "validation":
        result["R2_flower_minus_background"] = aggregate_species(r2_rows, r2_nulls, "r2")
    elif cohort == "discovery":
        result["R2_flower_minus_background"] = {
            "evaluable": False,
            "reason": "historical discovery matched-background recovery was frozen as not evaluable: 21339/21424 exact rows, 85 failures"
        }
    else:
        result["R2_flower_minus_background"] = {
            "evaluable": False,
            "reason": "third-cohort frozen measured table has no matched background palette"
        }

    details = pd.DataFrame(r1_rows).merge(pd.DataFrame(r3_rows), on=["inat_taxon_id", "species"], how="outer")
    if r2_rows:
        details = details.merge(pd.DataFrame(r2_rows), on=["inat_taxon_id", "species"], how="left")
    return result, details


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--discovery", type=Path, required=True)
    ap.add_argument("--validation", type=Path, required=True)
    ap.add_argument("--third", type=Path, required=True)
    ap.add_argument("--bio5", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()

    disc = prepare_cohort(args.discovery, "discovery", args.bio5, third=False)
    val = prepare_cohort(args.validation, "validation", args.bio5, third=False)
    third = prepare_cohort(args.third, "third", args.bio5, third=True)

    dr, dd = run_cohort(disc, "discovery", third=False)
    vr, vd = run_cohort(val, "validation", third=False)
    tr, td = run_cohort(third, "third", third=True)

    R1_all = bool(
        dr["R1_continuous_nine_colour"].get("supported_at_0_05", False)
        and vr["R1_continuous_nine_colour"].get("supported_at_0_05", False)
        and tr["R1_continuous_nine_colour"].get("supported_at_0_05", False)
    )
    R2_validation = bool(vr["R2_flower_minus_background"].get("supported_at_0_05", False))
    R3_all = bool(
        all(x["R3_same_observer_pair"].get("n_species", 0) >= 50 for x in (dr, vr, tr))
        and all(x["R3_same_observer_pair"].get("supported_at_0_05", False) for x in (dr, vr, tr))
    )

    result = {
        "schema": "fcp_BIO5_colour_sorting_robustness_v1",
        "date_jst": "2026-10-07",
        "status": "complete_posthoc_technical_robustness_audit",
        "confirmatory_decisions_changed": False,
        "permutations": PERMUTATIONS,
        "discovery": dr,
        "validation": vr,
        "third": tr,
        "cross_cohort_summary": {
            "R1_continuous_nine_colour_supported_all_three": R1_all,
            "R2_flower_minus_background_supported_validation": R2_validation,
            "R2_discovery_evaluable": False,
            "R3_same_observer_pair_supported_all_three_with_min50_species": R3_all,
        },
        "interpretation": {
            "if_R1": "BIO5-associated flower-colour turnover is not restricted to the coarse four-state classifier.",
            "if_R2": "In validation, BIO5-associated flower-colour turnover exceeds same-image background colour structure under the matched differential; discovery background inference remains unavailable by its historical fail-closed rule.",
            "if_R3": "The BIO5-colour association persists using only within-observer photo pairs, reducing observer/camera confounding.",
            "hard_nonclaims": [
                "post hoc robustness audit",
                "does not establish BIO5 causation",
                "does not measure fitness or local adaptation",
                "does not distinguish genetic differentiation from plasticity",
                "does not imply a universal direction of colour response",
            ],
        },
    }

    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)
    pd.concat([
        dd.assign(cohort="discovery"),
        vd.assign(cohort="validation"),
        td.assign(cohort="third"),
    ], ignore_index=True).to_csv(out / "species_robustness_metrics.csv", index=False)
    (out / "result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Post hoc equal-structure diagnostic for the validation D-spatial association.

Purpose
-------
Test a measurement-opportunity explanation that the frozen within-species
vertex-permutation null does not address: if every species is given the same
strength of geographic ordering while retaining its observed coordinates and
the exact multiset of observed colour vectors, does variation in species-wide
four-state diversity D alone induce a positive cross-species D-rho
association?

This is a post-outcome diagnostic. It cannot alter the frozen D-spatial
decision and it does not establish a causal spatial mechanism.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

EARTH_RADIUS_KM = 6371.0088
COLOURS = ["colour_white", "colour_yellow_orange", "colour_red_pink", "colour_blue_purple"]
MORPHS = ["white", "yellow_orange", "red_pink", "blue_purple"]
EXPECTED_SHA256 = "0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6"
EXPECTED_SPECIES = 363
EXPECTED_OBS_MEAN_RHO = 0.025482606069841617
EXPECTED_OBS_D_RHO = 0.10160084472811265
MASTER_SEED = 2026100701
DEFAULT_DESIGN = {4.0: 200, 6.0: 100, 10.0: 100}


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def as_bool(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False).astype(bool)
    return s.fillna("").astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def pairwise_geo_km(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    latr = np.deg2rad(np.asarray(lat, float))
    lonr = np.deg2rad(np.asarray(lon, float))
    c = np.cos(latr)
    xyz = np.column_stack([c * np.cos(lonr), c * np.sin(lonr), np.sin(latr)])
    dot = np.clip(xyz @ xyz.T, -1.0, 1.0)
    return np.arccos(dot) * EARTH_RADIUS_KM


def pairwise_jsd_matrix(prob: np.ndarray) -> np.ndarray:
    p = np.asarray(prob, float)
    mass = p.sum(axis=1)
    if p.ndim != 2 or np.any(~np.isfinite(p)) or np.any(p < 0) or np.any(mass <= 0):
        raise ValueError("invalid colour-composition rows")
    p = p / mass[:, None]
    a = p[:, None, :]
    b = p[None, :, :]
    m = 0.5 * (a + b)
    with np.errstate(divide="ignore", invalid="ignore"):
        ka = np.where(a > 0, a * np.log2(a / m), 0.0).sum(axis=2)
        kb = np.where(b > 0, b * np.log2(b / m), 0.0).sum(axis=2)
    return np.clip(0.5 * (ka + kb), 0.0, 1.0)


def direct_rho(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if np.ptp(y) <= 1e-15:
        return 0.0
    value = float(spearmanr(x, y).statistic)
    return value if np.isfinite(value) else 0.0


def species_D(frame: pd.DataFrame) -> float:
    counts = frame["morph"].astype(str).value_counts()
    p = np.array([counts.get(m, 0) / len(frame) for m in MORPHS], float)
    return float(1.0 - np.sum(p * p))


def spatial_axis(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    lat = np.asarray(lat, float)
    lon = np.asarray(lon, float)
    lat0 = float(np.mean(lat))
    x = EARTH_RADIUS_KM * np.cos(np.deg2rad(lat0)) * np.deg2rad(lon - np.mean(lon))
    y = EARTH_RADIUS_KM * np.deg2rad(lat - np.mean(lat))
    xy = np.column_stack([x, y])
    if np.max(np.ptp(xy, axis=0)) <= 1e-12:
        return np.zeros(len(lat), float)
    _, _, vh = np.linalg.svd(xy, full_matrices=False)
    axis = xy @ vh[0]
    sd = float(np.std(axis, ddof=0))
    return np.zeros(len(axis), float) if sd <= 1e-15 else (axis - axis.mean()) / sd


def colour_axis(jsd: np.ndarray) -> np.ndarray:
    """First classical-MDS coordinate of sqrt(JSD), used only to order vectors."""
    n = jsd.shape[0]
    if jsd.shape != (n, n):
        raise ValueError("JSD matrix must be square")
    # sqrt(JSD) is a metric; squared metric distances are therefore JSD.
    H = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * H @ jsd @ H
    eigval, eigvec = np.linalg.eigh(B)
    i = int(np.argmax(eigval))
    if eigval[i] <= 1e-14:
        return np.zeros(n, float)
    return eigvec[:, i] * np.sqrt(eigval[i])


def stable_seed(species: str, k: float, replicate: int) -> int:
    payload = f"{MASTER_SEED}|{species}|{k:g}|{replicate}".encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little", signed=False)


def impose_equal_structure(
    colours: np.ndarray,
    colour_order: np.ndarray,
    geo_pair: np.ndarray,
    spatial_score: np.ndarray,
    *,
    species: str,
    k: float,
    replicate: int,
) -> float:
    rng = np.random.default_rng(stable_seed(species, k, replicate))
    latent = spatial_score + float(k) * rng.standard_normal(len(spatial_score))

    position_order = np.argsort(latent, kind="stable")
    assigned = np.empty_like(colours)
    assigned[position_order] = colours[colour_order]

    sim_jsd = pairwise_jsd_matrix(assigned)
    u, v = np.triu_indices(len(assigned), k=1)
    return direct_rho(geo_pair, sim_jsd[u, v])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reserve", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args()

    source_sha = file_sha256(args.reserve)
    if source_sha != EXPECTED_SHA256:
        raise RuntimeError(f"validation measured-table SHA256 drift: {source_sha}")

    raw = pd.read_csv(args.reserve)
    keep = as_bool(raw["global_classifiable"]) & raw["morph"].astype(str).isin(MORPHS)
    cf = raw.loc[keep].copy()
    eligible = cf.groupby("inat_taxon_id").size()
    cf = cf.loc[cf["inat_taxon_id"].isin(eligible[eligible >= 40].index)].copy()
    if cf["inat_taxon_id"].nunique() != EXPECTED_SPECIES:
        raise RuntimeError(f"eligible species mismatch: {cf['inat_taxon_id'].nunique()}")

    species_rows = []
    cache = {}
    for taxon, g in cf.groupby("inat_taxon_id", sort=True):
        g = g.sort_values("photo_id", kind="stable").reset_index(drop=True)
        lat = g["latitude"].to_numpy(float)
        lon = g["longitude"].to_numpy(float)
        colours = g[COLOURS].to_numpy(float)
        geo = pairwise_geo_km(lat, lon)
        u, v = np.triu_indices(len(g), k=1)
        geo_pair = geo[u, v]
        jsd = pairwise_jsd_matrix(colours)
        observed_rho = direct_rho(geo_pair, jsd[u, v])
        d = species_D(g)
        species = str(g["species"].iloc[0])
        species_rows.append({
            "inat_taxon_id": int(taxon),
            "species": species,
            "n_classifiable": int(len(g)),
            "D": d,
            "observed_rho": observed_rho,
        })
        cache[int(taxon)] = {
            "species": species,
            "colours": colours,
            "colour_order": np.argsort(colour_axis(jsd), kind="stable"),
            "geo_pair": geo_pair,
            "spatial_score": spatial_axis(lat, lon),
        }

    species_df = pd.DataFrame(species_rows).sort_values("inat_taxon_id").reset_index(drop=True)
    obs_mean = float(species_df["observed_rho"].mean())
    obs_D_rho = float(spearmanr(species_df["D"], species_df["observed_rho"]).statistic)
    if abs(obs_mean - EXPECTED_OBS_MEAN_RHO) > 2e-12:
        raise RuntimeError(f"frozen validation mean rho not reproduced: {obs_mean}")
    if abs(obs_D_rho - EXPECTED_OBS_D_RHO) > 2e-12:
        raise RuntimeError(f"frozen validation D-rho not reproduced: {obs_D_rho}")

    sim_rows = []
    for k, nrep in DEFAULT_DESIGN.items():
        for rep in range(nrep):
            rhos = []
            for row in species_df.itertuples(index=False):
                item = cache[int(row.inat_taxon_id)]
                rhos.append(impose_equal_structure(
                    item["colours"],
                    item["colour_order"],
                    item["geo_pair"],
                    item["spatial_score"],
                    species=item["species"],
                    k=k,
                    replicate=rep,
                ))
            rhos = np.asarray(rhos, float)
            sim_rows.append({
                "k": k,
                "replicate": rep,
                "mean_species_rho": float(np.mean(rhos)),
                "spearman_D_rho": float(spearmanr(species_df["D"], rhos).statistic),
                "exceeds_observed_D_rho": bool(float(spearmanr(species_df["D"], rhos).statistic) >= obs_D_rho),
            })

    sims = pd.DataFrame(sim_rows)
    summaries = []
    for k, g in sims.groupby("k", sort=True):
        summaries.append({
            "k": float(k),
            "replicates": int(len(g)),
            "mean_of_mean_species_rho": float(g["mean_species_rho"].mean()),
            "min_spearman_D_rho": float(g["spearman_D_rho"].min()),
            "mean_spearman_D_rho": float(g["spearman_D_rho"].mean()),
            "max_spearman_D_rho": float(g["spearman_D_rho"].max()),
            "replicates_at_or_above_observed": int(g["exceeds_observed_D_rho"].sum()),
            "empirical_upper_exceedance_fraction": float(g["exceeds_observed_D_rho"].mean()),
        })

    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)
    species_df.to_csv(out / "validation_observed_species.csv", index=False, lineterminator="\n")
    sims.to_csv(out / "equal_structure_replicates.csv", index=False, lineterminator="\n")
    result = {
        "schema": "fcp_D_spatial_equal_structure_posthoc_v1",
        "date_jst": "2026-10-07",
        "status": "complete_posthoc_diagnostic",
        "role": "post_outcome_measurement_opportunity_diagnostic",
        "confirmatory_decisions_changed": False,
        "cohort": "species-disjoint validation",
        "source": {
            "measured_table_sha256": source_sha,
            "eligible_species": int(len(species_df)),
        },
        "observed_reproduction": {
            "mean_species_rho": obs_mean,
            "spearman_D_rho": obs_D_rho,
            "matches_frozen_values": True,
        },
        "simulation": {
            "design": "Preserve each species' observed coordinates and exact colour-vector multiset; order vectors by first classical-MDS coordinate of pairwise JSD and assign them along the first spatial principal axis plus species-specific deterministic Gaussian noise with the same k for every species.",
            "master_seed": MASTER_SEED,
            "k_interpretation": "Larger k adds more noise and therefore weaker common spatial ordering.",
            "summaries": summaries,
            "all_replicates_below_observed_D_rho": bool((~sims["exceeds_observed_D_rho"]).all()),
            "total_replicates": int(len(sims)),
            "total_replicates_at_or_above_observed": int(sims["exceeds_observed_D_rho"].sum()),
        },
        "interpretation": {
            "supported": "Under this linear-gradient equal-structure diagnostic, giving every validation species the same ordering rule does not reproduce the observed positive cross-species association between species-wide sampled colour-state diversity D and within-species spatial rho.",
            "not_supported": [
                "a confirmatory test of the D-spatial mechanism",
                "proof that high-D species are disproportionately geographically partitioned under every possible spatial process",
                "a causal explanation of geographic colour structure",
            ],
        },
        "limitations": [
            "post hoc diagnostic specified after the observed D-spatial result",
            "one spatial process family only: a linear principal-axis gradient plus Gaussian rank noise",
            "does not test patchy, threshold, multimodal or environmentally warped spatial structures",
            "simulation replicates are diagnostic rather than a preregistered Monte Carlo p-value",
        ],
    }
    (out / "result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

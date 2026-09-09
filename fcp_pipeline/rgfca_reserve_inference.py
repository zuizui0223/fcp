"""Prospectively fixed reserve replication statistics, with paired controls.

This module receives one species at a time, never selects a cohort, and never
treats photo pairs as independent replication. All permutations are vertex
permutations of complete photograph vectors.
"""
from __future__ import annotations

import hashlib

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

from fcp_pipeline.global_g3 import _upper_triangle_distances_km
from fcp_pipeline.photo_first_measurement import REFERENCE_RGB
from fcp_pipeline.rgfca_reserve_replication import COLOURS
from scripts.analysis.finalize_global_rgfca_background_control import pairwise_jsd_matrix
from scripts.analysis.run_global_rgfca_within_species_spatial_omnibus import species_observed_and_null

METRICS = ("primary", "observer_pair_exclusion", "calendar_quarter_stratification", "matched_background_differential")
PERMUTATIONS = 999


def permutation_seed(species: str, index: int, metric: str) -> int:
    if metric == "matched_background_differential":
        payload = f"202609071503|{species}|{index}"
    elif metric in METRICS[1:3]:
        payload = f"202609070902|{metric}|{species}|{index}"
    else:
        raise ValueError("unknown reserve sensitivity")
    return int.from_bytes(hashlib.sha256(payload.encode()).digest()[:8], "little")


def vertex_permutation(n: int, *, species: str, index: int, metric: str, quarters=None) -> np.ndarray:
    rng = np.random.default_rng(permutation_seed(species, index, metric))
    if metric != "calendar_quarter_stratification":
        return rng.permutation(n)
    quarters = np.asarray(quarters)
    if quarters.shape != (n,) or not np.isin(quarters, [1, 2, 3, 4]).all():
        raise ValueError("invalid fixed calendar-quarter labels")
    p = np.arange(n)
    for q in sorted(np.unique(quarters)):
        idx = np.flatnonzero(quarters == q)
        p[idx] = rng.permutation(idx)
    return p


def direct_rho(x, y) -> float:
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 3 or not np.isfinite(x).all() or np.ptp(x) <= 1e-12:
        raise ValueError("not_evaluable_pair_geometry")
    if not np.isfinite(y).all():
        raise ValueError("nonfinite colour distance")
    return 0.0 if np.ptp(y) <= 1e-15 else float(spearmanr(x, y).statistic)


def controlled_statistic(geo: np.ndarray, colour_distance: np.ndarray, *, species: str,
                         metric: str, pair_mask=None, quarters=None, permutations=PERMUTATIONS):
    n = colour_distance.shape[0]
    if colour_distance.shape != (n, n) or not np.isfinite(colour_distance).all():
        raise ValueError("invalid species pairwise colour matrix")
    u, v = np.triu_indices(n, k=1)
    if np.asarray(geo).shape != u.shape:
        raise ValueError("geographic pair vector shape mismatch")
    mask = np.ones(len(u), dtype=bool) if pair_mask is None else np.asarray(pair_mask, dtype=bool)
    if mask.shape != u.shape:
        raise ValueError("observer pair mask shape mismatch")
    x = np.asarray(geo, float)[mask]
    observed = direct_rho(x, colour_distance[u[mask], v[mask]])
    u, v = u[mask], v[mask]
    xrank = rankdata(x, method="average")
    xc = xrank - xrank.mean()
    xnorm = float(np.linalg.norm(xc))
    null = np.empty(permutations)
    checks = []
    for index in range(permutations):
        p = vertex_permutation(n, species=species, index=index, metric=metric, quarters=quarters)
        y = colour_distance[p[u], p[v]]
        # Required for the observer mask: a vertex permutation need not preserve
        # the masked distance multiset, so ranks and their norm are recomputed.
        yrank = rankdata(y, method="average")
        yc = yrank - yrank.mean()
        ynorm = float(np.linalg.norm(yc))
        null[index] = 0.0 if np.ptp(y) <= 1e-15 or ynorm <= 1e-15 else float(np.dot(xc, yc) / (xnorm * ynorm))
        if index in (0, permutations - 1):
            direct = direct_rho(x, y)
            if abs(null[index] - direct) > 2e-12:
                raise ValueError("controlled statistic differs from direct scipy")
            checks.append(abs(null[index] - direct))
    return observed, null, max(checks, default=0.0)


def evaluate_species(frame: pd.DataFrame) -> tuple[dict, np.ndarray]:
    """All four registered tests, same eligible photos, one exact taxon."""
    if frame.inat_taxon_id.nunique() != 1 or frame.species.nunique() != 1 or len(frame) < 40:
        raise ValueError("one eligible species is required")
    frame = frame.sort_values("photo_id", kind="stable")
    if frame.photo_id.duplicated().any():
        raise ValueError("duplicate photo IDs")
    species = str(frame.species.iloc[0])
    lat, lon = frame.latitude.to_numpy(float), frame.longitude.to_numpy(float)
    colours = frame[COLOURS].to_numpy(float)
    quarters = pd.to_datetime(frame.observed_on, errors="raise").dt.quarter.to_numpy()
    observers = frame.observer_id.to_numpy(float)
    if not np.isfinite(observers).all() or not np.isin(quarters, [1, 2, 3, 4]).all():
        raise ValueError("missing observer/date metadata")
    primary, primary_null, _ = species_observed_and_null(lat, lon, colours, species=species)
    geo = _upper_triangle_distances_km(lat, lon)
    jsd = pairwise_jsd_matrix(colours)
    u, v = np.triu_indices(len(frame), k=1)
    if abs(primary - direct_rho(geo, jsd[u, v])) > 2e-12:
        raise ValueError("primary direct scipy mismatch")
    flower12 = frame[[f"palette_count_{p}" for p in REFERENCE_RGB]].to_numpy(float)
    background12 = frame[[f"background_palette_count_{p}" for p in REFERENCE_RGB]].to_numpy(float)
    if (flower12 < 0).any() or (background12 < 0).any():
        raise ValueError("negative matched palette counts")
    fjsd = pairwise_jsd_matrix(flower12)
    bjsd = pairwise_jsd_matrix(background12)
    config = (
        ("observer_pair_exclusion", jsd, {"pair_mask": observers[u] != observers[v]}),
        ("calendar_quarter_stratification", jsd, {"quarters": quarters}),
        ("matched_background_differential", fjsd - bjsd, {}),
    )
    row = {"inat_taxon_id": int(frame.inat_taxon_id.iloc[0]), "species": species, "photos": len(frame),
           "rho_primary": primary, "rho_flower12_descriptive": direct_rho(geo, fjsd[u, v]),
           "rho_background12_descriptive": direct_rho(geo, bjsd[u, v])}
    nulls = [primary_null]
    errors = []
    for metric, distance, kwargs in config:
        obs, null, error = controlled_statistic(geo, distance, species=species, metric=metric, **kwargs)
        row[f"rho_{metric}"] = obs
        errors.append(error)
        nulls.append(null)
    row["maximum_direct_check_error"] = max(errors)
    row["direct_checks"] = 7
    return row, np.stack(nulls)


def summarize(species_rows: pd.DataFrame, null: np.ndarray) -> dict:
    n = len(species_rows)
    if n < 250 or null.shape != (n, len(METRICS), PERMUTATIONS) or not np.isfinite(null).all():
        raise ValueError("incomplete eligible-species null tensor")
    result = {}
    for index, metric in enumerate(METRICS):
        observed = float(species_rows[f"rho_{metric}"].mean())
        if not np.isfinite(observed):
            raise ValueError("nonfinite observed mean")
        distribution = null[:, index, :].mean(axis=0)
        p = float((1 + np.count_nonzero(distribution >= observed)) / (PERMUTATIONS + 1))
        nondegenerate = bool(np.ptp(distribution) > 1e-15)
        result[metric] = {"mean_rho": observed, "p_upper": p, "null_mean": float(distribution.mean()),
                          "null_sd": float(distribution.std(ddof=1)), "null_nondegenerate": nondegenerate,
                          "status": "evaluable" if nondegenerate else "not_evaluable_degenerate_null",
                          "positive_supported": bool(observed > 0 and p < .05 and nondegenerate)}
    rng = np.random.default_rng(202609070903)
    observed = species_rows.rho_primary.to_numpy(float)
    boots = observed[rng.integers(n, size=(4999, n))].mean(axis=1)
    return {"metrics": result,
            "directional_photo_association_replicated": result["primary"]["positive_supported"],
            "flower_specific_robust_replication": all(v["positive_supported"] for v in result.values()),
            "conditional_species_bootstrap_primary_mean_95pct": np.quantile(boots, [.025, .975]).tolist(),
            "bootstrap_ceiling": "Conditional species bootstrap, not spatially or phylogenetically independent uncertainty for all plants.",
            "cause_or_shared_boundary_inferred": False}

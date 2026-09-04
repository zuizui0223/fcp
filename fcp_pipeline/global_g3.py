"""Observed G3 prevalence and heterogeneity for the repeated global atlas.

G3 is deliberately descriptive with respect to the global distance-colour statistic:
it reports repeated species-specific Spearman rho values and between-species
heterogeneity, but it cannot rescue a null G1 and receives no confirmatory p-value.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from .global_repeated_atlas import RepeatedAtlasSchedule, build_repeated_atlas_schedule
from .global_rgfca_engine import CanonicalColourPool, canonical_colour_pool

EARTH_RADIUS_KM = 6371.0088


@dataclass(frozen=True)
class G3Result:
    schedule: RepeatedAtlasSchedule
    outer: pd.DataFrame
    species: pd.DataFrame
    tau2_fisher_z: float
    tau2_species_used: int
    median_outer_mean_rho: float
    median_outer_median_rho: float
    median_outer_positive_fraction: float


def _upper_triangle_distances_km(latitude: np.ndarray, longitude: np.ndarray) -> np.ndarray:
    lat = np.deg2rad(np.asarray(latitude, dtype=float))
    lon = np.deg2rad(np.asarray(longitude, dtype=float))
    c = np.cos(lat)
    xyz = np.column_stack([c * np.cos(lon), c * np.sin(lon), np.sin(lat)])
    dot = np.clip(xyz @ xyz.T, -1.0, 1.0)
    distance = np.arccos(dot) * EARTH_RADIUS_KM
    upper = np.triu_indices(len(lat), k=1)
    return distance[upper]


def _upper_triangle_jsd(probabilities: np.ndarray) -> np.ndarray:
    p = np.asarray(probabilities, dtype=float)
    if p.ndim != 2 or p.shape[1] != 4:
        raise ValueError("probabilities must have shape (n,4)")
    mass = p.sum(axis=1)
    if np.any(~np.isfinite(p)) or np.any(p < 0) or np.any(mass <= 0):
        raise ValueError("probabilities must be finite nonnegative rows with positive mass")
    p = p / mass[:, None]
    a = p[:, None, :]
    b = p[None, :, :]
    m = 0.5 * (a + b)
    with np.errstate(divide="ignore", invalid="ignore"):
        kl_a = np.where(a > 0, a * np.log2(a / m), 0.0).sum(axis=2)
        kl_b = np.where(b > 0, b * np.log2(b / m), 0.0).sum(axis=2)
    jsd = np.clip(0.5 * (kl_a + kl_b), 0.0, 1.0)
    upper = np.triu_indices(len(p), k=1)
    return jsd[upper]


def species_distance_colour_rho(
    latitude: Sequence[float],
    longitude: Sequence[float],
    colours: np.ndarray,
) -> float:
    """Frozen G3 rho for one fixed photo draw.

    A constant colour-distance vector has no spatial colour differentiation and is
    defined as rho=0. A constant geographic-distance vector contains no distance
    information and is returned as NaN rather than being treated as a biological zero.
    """
    lat = np.asarray(latitude, dtype=float)
    lon = np.asarray(longitude, dtype=float)
    colour = np.asarray(colours, dtype=float)
    if lat.ndim != 1 or lon.ndim != 1 or lat.shape != lon.shape or len(lat) < 3:
        raise ValueError("latitude/longitude must be equal vectors with at least three photos")
    if colour.shape != (len(lat), 4):
        raise ValueError("colours must have shape (n_photos,4)")
    distance = _upper_triangle_distances_km(lat, lon)
    jsd = _upper_triangle_jsd(colour)
    if not np.isfinite(distance).all() or not np.isfinite(jsd).all():
        raise ValueError("pairwise G3 quantities must be finite")
    if np.ptp(jsd) <= 1e-15:
        return 0.0
    if np.ptp(distance) <= 1e-12:
        return float("nan")
    rho = float(spearmanr(distance, jsd).statistic)
    return rho if np.isfinite(rho) else float("nan")


def _dersimonian_laird_tau2(
    species_rhos: dict[str, list[float]],
    *,
    variance_floor: float = 1e-6,
) -> tuple[float, int]:
    effects: list[float] = []
    variances: list[float] = []
    floor = float(variance_floor)
    if floor <= 0:
        raise ValueError("variance_floor must be positive")
    for values in species_rhos.values():
        rho = np.asarray([x for x in values if np.isfinite(x)], dtype=float)
        if len(rho) < 2:
            continue
        z = np.arctanh(np.clip(rho, -0.999999, 0.999999))
        effects.append(float(np.mean(z)))
        variances.append(max(float(np.var(z, ddof=1) / len(z)), floor))
    k = len(effects)
    if k < 2:
        return float("nan"), k
    y = np.asarray(effects, dtype=float)
    v = np.asarray(variances, dtype=float)
    w = 1.0 / v
    fixed = float(np.sum(w * y) / np.sum(w))
    q = float(np.sum(w * np.square(y - fixed)))
    c = float(np.sum(w) - np.sum(np.square(w)) / np.sum(w))
    tau2 = max(0.0, (q - (k - 1)) / c) if c > 0 else float("nan")
    return float(tau2), k


def run_g3_prevalence(
    measured_pool: pd.DataFrame,
    *,
    n_outer: int = 200,
    species_per_outer: int = 250,
    photos_per_species: int = 20,
    minimum_pool_photos_per_species: int = 40,
    species_seed: int = 2026090401,
    photo_master_seed: int = 2026090402,
    variance_floor: float = 1e-6,
) -> G3Result:
    pool: CanonicalColourPool = canonical_colour_pool(measured_pool)
    schedule = build_repeated_atlas_schedule(
        pool.photo_ids,
        pool.species,
        n_outer=int(n_outer),
        species_per_outer=int(species_per_outer),
        photos_per_species=int(photos_per_species),
        minimum_pool_photos_per_species=int(minimum_pool_photos_per_species),
        species_seed=int(species_seed),
        photo_master_seed=int(photo_master_seed),
    )
    id_to_row = {int(pid): i for i, pid in enumerate(pool.photo_ids)}
    rows: list[dict[str, object]] = []
    species_rhos: dict[str, list[float]] = {label: [] for label in schedule.species_labels}
    for outer in range(schedule.n_outer):
        for slot in range(schedule.species_per_outer):
            label = str(schedule.outer_species[outer, slot])
            ids = schedule.outer_photo_ids[outer, slot]
            idx = np.asarray([id_to_row[int(pid)] for pid in ids], dtype=np.int64)
            if np.any(pool.species[idx] != label):
                raise RuntimeError("G3 schedule assigns a photo to the wrong species")
            rho = species_distance_colour_rho(
                pool.latitude[idx],
                pool.longitude[idx],
                pool.colours[idx],
            )
            species_rhos[label].append(rho)
            rows.append({"outer": outer, "species": label, "rho": rho})
    raw = pd.DataFrame(rows)
    outer_rows: list[dict[str, object]] = []
    for outer, group in raw.groupby("outer", sort=True):
        rho = pd.to_numeric(group["rho"], errors="coerce").to_numpy(float)
        finite = rho[np.isfinite(rho)]
        outer_rows.append({
            "outer": int(outer),
            "scheduled_species": int(len(rho)),
            "evaluable_species": int(len(finite)),
            "mean_rho": float(np.mean(finite)) if len(finite) else float("nan"),
            "median_rho": float(np.median(finite)) if len(finite) else float("nan"),
            "positive_fraction": float(np.mean(finite > 0)) if len(finite) else float("nan"),
        })
    outer_frame = pd.DataFrame(outer_rows)

    species_rows: list[dict[str, object]] = []
    for label in schedule.species_labels:
        rho = np.asarray(species_rhos[label], dtype=float)
        finite = rho[np.isfinite(rho)]
        species_rows.append({
            "species": label,
            "outer_appearances": int(len(rho)),
            "evaluable_appearances": int(len(finite)),
            "mean_rho": float(np.mean(finite)) if len(finite) else float("nan"),
            "median_rho": float(np.median(finite)) if len(finite) else float("nan"),
            "q025_rho": float(np.quantile(finite, 0.025)) if len(finite) else float("nan"),
            "q975_rho": float(np.quantile(finite, 0.975)) if len(finite) else float("nan"),
        })
    species_frame = pd.DataFrame(species_rows)
    tau2, tau_k = _dersimonian_laird_tau2(species_rhos, variance_floor=float(variance_floor))
    return G3Result(
        schedule=schedule,
        outer=outer_frame,
        species=species_frame,
        tau2_fisher_z=tau2,
        tau2_species_used=tau_k,
        median_outer_mean_rho=float(np.nanmedian(outer_frame["mean_rho"])),
        median_outer_median_rho=float(np.nanmedian(outer_frame["median_rho"])),
        median_outer_positive_fraction=float(np.nanmedian(outer_frame["positive_fraction"])),
    )


__all__ = [
    "G3Result",
    "run_g3_prevalence",
    "species_distance_colour_rho",
]

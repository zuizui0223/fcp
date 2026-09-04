"""Opportunity-conditioned global recurrent-barrier field primitives.

This module contains only geometry/statistics.  It does not acquire images or
choose species.  Geometry is prepared once and reused for observed and
within-species colour permutations, so the null cannot alter where species had
an opportunity to contribute a geographic edge.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
import pandas as pd

EARTH_RADIUS_KM = 6371.0088


@dataclass(frozen=True)
class BarrierGrid:
    n_lon: int
    n_sinlat: int
    latitude: np.ndarray
    longitude: np.ndarray

    @property
    def n_cells(self) -> int:
        return int(self.n_lon * self.n_sinlat)


@dataclass(frozen=True)
class BarrierGeometry:
    grid: BarrierGrid
    species: tuple[str, ...]
    edge_species_index: np.ndarray
    edge_weighted_kernel: np.ndarray
    opportunity: np.ndarray
    distinct_species_support: np.ndarray
    kernel_km: float
    cutoff_multiplier: float


@dataclass(frozen=True)
class BarrierFieldResult:
    field: np.ndarray
    evaluable: np.ndarray
    concentration: float
    weighted_mean: float
    opportunity: np.ndarray
    distinct_species_support: np.ndarray


def equal_area_grid_centers(n_lon: int, n_sinlat: int) -> BarrierGrid:
    n_lon = int(n_lon)
    n_sinlat = int(n_sinlat)
    if n_lon < 1 or n_sinlat < 1:
        raise ValueError("grid dimensions must be positive")
    lon = -180.0 + (np.arange(n_lon, dtype=float) + 0.5) * 360.0 / n_lon
    sinlat = -1.0 + (np.arange(n_sinlat, dtype=float) + 0.5) * 2.0 / n_sinlat
    lat = np.degrees(np.arcsin(np.clip(sinlat, -1.0, 1.0)))
    lon_mesh, lat_mesh = np.meshgrid(lon, lat)
    return BarrierGrid(
        n_lon=n_lon,
        n_sinlat=n_sinlat,
        latitude=lat_mesh.ravel(),
        longitude=lon_mesh.ravel(),
    )


def _xyz(latitude: Sequence[float], longitude: Sequence[float]) -> np.ndarray:
    lat = np.deg2rad(np.asarray(latitude, dtype=float))
    lon = np.deg2rad(np.asarray(longitude, dtype=float))
    c = np.cos(lat)
    return np.column_stack([c * np.cos(lon), c * np.sin(lon), np.sin(lat)])


def spherical_midpoint(
    latitude_a: Sequence[float],
    longitude_a: Sequence[float],
    latitude_b: Sequence[float],
    longitude_b: Sequence[float],
) -> tuple[np.ndarray, np.ndarray]:
    """Great-circle midpoint using normalized Cartesian sums."""
    a = _xyz(latitude_a, longitude_a)
    b = _xyz(latitude_b, longitude_b)
    if a.shape != b.shape:
        raise ValueError("endpoint arrays must have matching shapes")
    midpoint = a + b
    norm = np.linalg.norm(midpoint, axis=1)
    if np.any(norm <= 1e-12):
        raise ValueError("antipodal endpoints do not have a unique spherical midpoint")
    midpoint = midpoint / norm[:, None]
    lat = np.degrees(np.arcsin(np.clip(midpoint[:, 2], -1.0, 1.0)))
    lon = np.degrees(np.arctan2(midpoint[:, 1], midpoint[:, 0]))
    return lat, lon


def within_species_rank_scores(
    edges: pd.DataFrame,
    *,
    score_column: str,
    species_column: str = "species",
) -> np.ndarray:
    """Mid-rank edge scores in (0, 1), with mean exactly 0.5 per species."""
    if score_column not in edges or species_column not in edges:
        raise ValueError("edges lack required score/species columns")
    frame = edges[[species_column, score_column]].copy()
    if frame[score_column].isna().any():
        raise ValueError("edge scores cannot be missing")
    rank = frame.groupby(species_column, sort=False, observed=True)[score_column].rank(method="average")
    n = frame.groupby(species_column, sort=False, observed=True)[score_column].transform("size")
    return ((rank.to_numpy(dtype=float) - 0.5) / n.to_numpy(dtype=float)).astype(float)


def _central_angle_distance_km(a_xyz: np.ndarray, b_xyz: np.ndarray) -> np.ndarray:
    dot = np.clip(a_xyz @ b_xyz.T, -1.0, 1.0)
    return np.arccos(dot) * EARTH_RADIUS_KM


def prepare_barrier_geometry(
    edges: pd.DataFrame,
    *,
    grid: BarrierGrid,
    kernel_km: float,
    cutoff_multiplier: float = 3.0,
    species_column: str = "species",
    latitude_column: str = "midpoint_latitude",
    longitude_column: str = "midpoint_longitude",
) -> BarrierGeometry:
    required = {species_column, latitude_column, longitude_column}
    missing = sorted(required - set(edges.columns))
    if missing:
        raise ValueError(f"edges lack required geometry columns: {missing}")
    if len(edges) == 0:
        raise ValueError("edges must be non-empty")
    kernel_km = float(kernel_km)
    cutoff_multiplier = float(cutoff_multiplier)
    if kernel_km <= 0 or cutoff_multiplier <= 0:
        raise ValueError("kernel parameters must be positive")

    names = edges[species_column].astype(str).to_numpy()
    species = tuple(pd.unique(names).tolist())
    species_lookup = {name: i for i, name in enumerate(species)}
    species_index = np.fromiter((species_lookup[name] for name in names), dtype=np.int64, count=len(names))
    counts = np.bincount(species_index, minlength=len(species)).astype(float)
    if np.any(counts <= 0):
        raise RuntimeError("internal species-count failure")
    edge_equal_species_weight = 1.0 / counts[species_index]

    edge_xyz = _xyz(edges[latitude_column], edges[longitude_column])
    grid_xyz = _xyz(grid.latitude, grid.longitude)
    distance = _central_angle_distance_km(edge_xyz, grid_xyz)
    cutoff_km = kernel_km * cutoff_multiplier
    kernel = np.exp(-0.5 * np.square(distance / kernel_km))
    kernel[distance > cutoff_km] = 0.0
    row_sum = kernel.sum(axis=1)
    missing_support = np.flatnonzero(row_sum <= 0)
    if len(missing_support):
        # A coarse grid can in principle have no centre within a very small
        # compact support.  Deterministically assign such an edge to its nearest
        # grid centre rather than silently deleting geographic opportunity.
        nearest = np.argmin(distance[missing_support], axis=1)
        kernel[missing_support, nearest] = 1.0
        row_sum = kernel.sum(axis=1)
    kernel = kernel / row_sum[:, None]
    weighted_kernel = kernel * edge_equal_species_weight[:, None]
    opportunity = weighted_kernel.sum(axis=0)

    support = np.zeros(grid.n_cells, dtype=np.int64)
    for species_id in range(len(species)):
        contribution = weighted_kernel[species_index == species_id].sum(axis=0)
        support += (contribution > 0).astype(np.int64)

    return BarrierGeometry(
        grid=grid,
        species=species,
        edge_species_index=species_index,
        edge_weighted_kernel=weighted_kernel,
        opportunity=opportunity,
        distinct_species_support=support,
        kernel_km=kernel_km,
        cutoff_multiplier=cutoff_multiplier,
    )


def barrier_field(
    geometry: BarrierGeometry,
    rank_scores: Sequence[float],
    *,
    minimum_distinct_species: int = 5,
) -> BarrierFieldResult:
    score = np.asarray(rank_scores, dtype=float)
    if score.shape != (geometry.edge_weighted_kernel.shape[0],):
        raise ValueError("rank_scores length does not match prepared edge geometry")
    if not np.all(np.isfinite(score)):
        raise ValueError("rank_scores must be finite")
    minimum_distinct_species = int(minimum_distinct_species)
    if minimum_distinct_species < 1:
        raise ValueError("minimum_distinct_species must be positive")
    numerator = score @ geometry.edge_weighted_kernel
    evaluable = (
        (geometry.distinct_species_support >= minimum_distinct_species)
        & np.isfinite(geometry.opportunity)
        & (geometry.opportunity > 0)
    )
    field = np.full(geometry.grid.n_cells, np.nan, dtype=float)
    field[evaluable] = numerator[evaluable] / geometry.opportunity[evaluable]
    if not np.any(evaluable):
        return BarrierFieldResult(
            field=field,
            evaluable=evaluable,
            concentration=float("nan"),
            weighted_mean=float("nan"),
            opportunity=geometry.opportunity.copy(),
            distinct_species_support=geometry.distinct_species_support.copy(),
        )
    weights = geometry.opportunity[evaluable]
    values = field[evaluable]
    weighted_mean = float(np.average(values, weights=weights))
    concentration = float(np.average(np.square(values - weighted_mean), weights=weights))
    return BarrierFieldResult(
        field=field,
        evaluable=evaluable,
        concentration=concentration,
        weighted_mean=weighted_mean,
        opportunity=geometry.opportunity.copy(),
        distinct_species_support=geometry.distinct_species_support.copy(),
    )


def permute_scores_within_species(
    score: Sequence[float],
    species_index: np.ndarray,
    *,
    rng: np.random.Generator,
) -> np.ndarray:
    score = np.asarray(score, dtype=float)
    if score.shape != species_index.shape:
        raise ValueError("score/species_index shapes differ")
    out = score.copy()
    for species_id in np.unique(species_index):
        idx = np.flatnonzero(species_index == species_id)
        out[idx] = score[rng.permutation(idx)]
    return out


def concentration_permutation_test(
    geometry: BarrierGeometry,
    rank_scores: Sequence[float],
    *,
    minimum_distinct_species: int = 5,
    permutations: int = 999,
    seed: int = 20260904,
) -> dict[str, object]:
    permutations = int(permutations)
    if permutations < 1:
        raise ValueError("permutations must be positive")
    observed = barrier_field(
        geometry,
        rank_scores,
        minimum_distinct_species=minimum_distinct_species,
    )
    if not np.isfinite(observed.concentration):
        raise ValueError("observed field is not evaluable")
    rng = np.random.default_rng(int(seed))
    null = np.empty(permutations, dtype=float)
    score = np.asarray(rank_scores, dtype=float)
    for i in range(permutations):
        permuted = permute_scores_within_species(score, geometry.edge_species_index, rng=rng)
        null[i] = barrier_field(
            geometry,
            permuted,
            minimum_distinct_species=minimum_distinct_species,
        ).concentration
    p_upper = float((1 + np.count_nonzero(null >= observed.concentration)) / (permutations + 1))
    return {
        "observed": observed,
        "null": null,
        "p_upper": p_upper,
        "null_mean": float(null.mean()),
        "null_q025": float(np.quantile(null, 0.025)),
        "null_q975": float(np.quantile(null, 0.975)),
        "permutations": permutations,
        "seed": int(seed),
    }


__all__ = [
    "BarrierFieldResult",
    "BarrierGeometry",
    "BarrierGrid",
    "barrier_field",
    "concentration_permutation_test",
    "equal_area_grid_centers",
    "prepare_barrier_geometry",
    "permute_scores_within_species",
    "spherical_midpoint",
    "within_species_rank_scores",
]

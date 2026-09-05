"""Opportunity-conditioned global recurrent-barrier field primitives.

This module contains only geometry/statistics. It does not acquire images or
choose species. Geometry is prepared once and reused for observed and null
fields, so the null cannot alter where species had an opportunity to contribute
a geographic edge.

The confirmatory RGFCA G1 null is photograph-level: complete soft colour vectors
are permuted within species and edge Jensen-Shannon divergences plus within-
species ranks are recomputed. The older edge-score permutation helper is retained
only as a legacy diagnostic and must not be used as the primary RGFCA null.
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


def edge_jensen_shannon_divergence(
    colour_vectors: Sequence[Sequence[float]],
    edge_nodes: Sequence[Sequence[int]],
) -> np.ndarray:
    """Vectorized Jensen-Shannon divergence in bits for fixed graph edges.

    Rows are complete soft colour vectors. They are normalized internally, so
    nonnegative palette masses need not sum exactly to one. The score is bounded
    to [0, 1] for base-2 logarithms.
    """
    values = np.asarray(colour_vectors, dtype=float)
    edges = np.asarray(edge_nodes)
    if values.ndim != 2 or values.shape[0] < 2 or values.shape[1] < 2:
        raise ValueError("colour_vectors must have shape (n_photos>=2, n_groups>=2)")
    if not np.all(np.isfinite(values)) or np.any(values < 0):
        raise ValueError("colour vectors must be finite and non-negative")
    mass = values.sum(axis=1)
    if np.any(mass <= 0):
        raise ValueError("every colour vector must have positive mass")
    if edges.ndim != 2 or edges.shape[1] != 2:
        raise ValueError("edge_nodes must have shape (n_edges, 2)")
    if not np.issubdtype(edges.dtype, np.integer):
        if np.any(edges != np.floor(edges)):
            raise ValueError("edge node indices must be integers")
        edges = edges.astype(np.int64)
    else:
        edges = edges.astype(np.int64, copy=False)
    if len(edges) == 0:
        return np.empty(0, dtype=float)
    if edges.min() < 0 or edges.max() >= len(values):
        raise ValueError("edge node index out of bounds")
    if np.any(edges[:, 0] == edges[:, 1]):
        raise ValueError("self edges are not allowed")

    probabilities = values / mass[:, None]
    p = probabilities[edges[:, 0]]
    q = probabilities[edges[:, 1]]
    m = 0.5 * (p + q)

    def _row_kl(a: np.ndarray, b: np.ndarray) -> np.ndarray:
        term = np.zeros_like(a, dtype=float)
        keep = a > 0
        term[keep] = a[keep] * np.log2(a[keep] / b[keep])
        return term.sum(axis=1)

    jsd = 0.5 * _row_kl(p, m) + 0.5 * _row_kl(q, m)
    return np.clip(jsd, 0.0, 1.0)


def permute_colour_vectors_within_species(
    colour_vectors: Sequence[Sequence[float]],
    photo_species: Sequence[object],
    *,
    rng: np.random.Generator,
) -> np.ndarray:
    """Permute complete photograph-level colour rows strictly within species."""
    values = np.asarray(colour_vectors, dtype=float)
    species = np.asarray(photo_species)
    if values.ndim != 2 or species.ndim != 1 or len(values) != len(species):
        raise ValueError("colour_vectors and photo_species observation counts must match")
    if not np.all(np.isfinite(values)):
        raise ValueError("colour_vectors must be finite")
    out = values.copy()
    for label in pd.unique(species):
        idx = np.flatnonzero(species == label)
        out[idx] = values[rng.permutation(idx)]
    return out


def _edge_species_from_nodes(
    photo_species: np.ndarray,
    edge_nodes: np.ndarray,
) -> np.ndarray:
    left = photo_species[edge_nodes[:, 0]]
    right = photo_species[edge_nodes[:, 1]]
    if np.any(left != right):
        raise ValueError("every barrier edge must remain within species")
    return left.astype(str)


def rank_scores_from_node_colours(
    colour_vectors: Sequence[Sequence[float]],
    photo_species: Sequence[object],
    edge_nodes: Sequence[Sequence[int]],
) -> tuple[np.ndarray, np.ndarray]:
    """Recompute edge JSD and within-species ranks from photograph colours."""
    species = np.asarray(photo_species)
    edges = np.asarray(edge_nodes, dtype=np.int64)
    if species.ndim != 1:
        raise ValueError("photo_species must be one-dimensional")
    edge_species = _edge_species_from_nodes(species, edges)
    raw = edge_jensen_shannon_divergence(colour_vectors, edges)
    frame = pd.DataFrame({"species": edge_species, "raw_colour_jsd": raw})
    rank = within_species_rank_scores(frame, score_column="raw_colour_jsd")
    return raw, rank


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
    """Legacy edge-score shuffle retained for diagnostics, not primary G1."""
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
    """Legacy edge-score null; forbidden as confirmatory RGFCA G1 primary null."""
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
        "null_unit": "legacy_edge_score",
        "primary_rgfca_allowed": False,
    }


def node_colour_concentration_permutation_test(
    geometry: BarrierGeometry,
    colour_vectors: Sequence[Sequence[float]],
    photo_species: Sequence[object],
    edge_nodes: Sequence[Sequence[int]],
    *,
    minimum_distinct_species: int = 5,
    permutations: int = 999,
    seed: int = 20260904,
) -> dict[str, object]:
    """Confirmatory single-field null using complete photo-colour vector shuffles.

    Graph geometry is fixed. For every null replicate, complete colour rows move
    only among photographs of the same species; edge JSD and within-species ranks
    are then recomputed. This preserves the dependence induced by graph edges that
    share a photograph.

    The full RGFCA runner must apply one such species-conditioned photo permutation
    to the full measured pool and reuse that assignment across its frozen 200 outer
    realizations. This helper validates the statistical primitive for one fixed
    realization.
    """
    permutations = int(permutations)
    if permutations < 1:
        raise ValueError("permutations must be positive")
    values = np.asarray(colour_vectors, dtype=float)
    species = np.asarray(photo_species)
    edges = np.asarray(edge_nodes, dtype=np.int64)
    if values.ndim != 2 or species.ndim != 1 or len(values) != len(species):
        raise ValueError("colour_vectors and photo_species counts must match")
    edge_species = _edge_species_from_nodes(species, edges)
    expected_edge_species = np.asarray([geometry.species[i] for i in geometry.edge_species_index])
    if len(edge_species) != len(expected_edge_species) or np.any(edge_species.astype(str) != expected_edge_species.astype(str)):
        raise ValueError("edge node species/order does not match prepared barrier geometry")

    _, observed_rank = rank_scores_from_node_colours(values, species, edges)
    observed = barrier_field(
        geometry,
        observed_rank,
        minimum_distinct_species=minimum_distinct_species,
    )
    if not np.isfinite(observed.concentration):
        raise ValueError("observed field is not evaluable")

    rng = np.random.default_rng(int(seed))
    null = np.empty(permutations, dtype=float)
    for b in range(permutations):
        permuted_values = permute_colour_vectors_within_species(values, species, rng=rng)
        _, rank = rank_scores_from_node_colours(permuted_values, species, edges)
        null[b] = barrier_field(
            geometry,
            rank,
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
        "null_unit": "complete_photo_colour_vector_within_species",
        "primary_rgfca_allowed": True,
    }


__all__ = [
    "BarrierFieldResult",
    "BarrierGeometry",
    "BarrierGrid",
    "barrier_field",
    "concentration_permutation_test",
    "edge_jensen_shannon_divergence",
    "equal_area_grid_centers",
    "node_colour_concentration_permutation_test",
    "permute_colour_vectors_within_species",
    "permute_scores_within_species",
    "prepare_barrier_geometry",
    "rank_scores_from_node_colours",
    "spherical_midpoint",
    "within_species_rank_scores",
]

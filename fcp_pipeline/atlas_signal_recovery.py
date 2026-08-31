"""Exact-statistic synthetic signal recovery for the frozen atlas geometry."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
from typing import Literal, Sequence

import numpy as np
from scipy import sparse

from .continuous_colour_boundaries import (
    average_rank_intensity,
    edge_colour_discontinuity,
    opportunity_weighted_concentration,
    shared_boundary_intensity,
)
from .shared_transition_surface import EdgeCellGeometry, cell_mean_intensity


@dataclass(frozen=True)
class AtlasSpeciesGeometry:
    species: str
    latitude: np.ndarray
    longitude: np.ndarray
    geometry: EdgeCellGeometry


def _species_phase(species: str) -> tuple[float, bool]:
    digest = hashlib.sha256(f"fcp-atlas-signal-v2|{species}".encode("utf-8")).digest()
    fraction = int.from_bytes(digest[:8], "big") / float(2**64)
    return fraction, bool(digest[8] & 1)


def synthetic_colour_vectors(
    geometries: Sequence[AtlasSpeciesGeometry],
    *,
    effect_size: float,
    scenario: Literal[
        "null_stationary",
        "within_species_heterogeneous_boundaries",
        "shared_geographic_boundary",
    ],
    shared_anchor_latitude: float,
    shared_anchor_longitude: float,
    rng: np.random.Generator,
) -> list[np.ndarray]:
    """Generate continuous vectors without using any observed atlas colour."""

    effect_size = float(effect_size)
    if effect_size < 0 or not np.isfinite(effect_size):
        raise ValueError("effect_size must be finite and non-negative")
    values: list[np.ndarray] = []
    for item in geometries:
        lat = np.asarray(item.latitude, dtype=float)
        lon = np.asarray(item.longitude, dtype=float)
        if scenario == "null_stationary" or effect_size == 0:
            state = np.zeros(len(lat), dtype=float)
        elif scenario == "shared_geographic_boundary":
            # A fixed great-circle-like linear separator through the highest-opportunity
            # atlas cell.  The orientation and anchor depend only on frozen geometry.
            x = np.cos(np.deg2rad(shared_anchor_latitude)) * np.sin(
                np.deg2rad(lon - shared_anchor_longitude)
            )
            y = np.sin(np.deg2rad(lat - shared_anchor_latitude))
            state = np.where(x + 0.35 * y >= 0.0, 0.5, -0.5)
        elif scenario == "within_species_heterogeneous_boundaries":
            fraction, use_latitude = _species_phase(item.species)
            coordinate = lat if use_latitude else lon
            lower, upper = np.quantile(coordinate, [0.2, 0.8])
            anchor = lower + fraction * (upper - lower)
            state = np.where(coordinate >= anchor, 0.5, -0.5)
        else:
            raise ValueError(f"unknown signal scenario: {scenario}")

        signal = effect_size * state
        noise = rng.normal(0.0, 1.0, size=(len(lat), 3))
        values.append(
            noise
            + np.column_stack(
                (
                    signal,
                    0.5 * signal,
                    -0.25 * signal,
                )
            )
        )
    return values


def atlas_shared_concentration(
    colour_vectors: Sequence[np.ndarray],
    geometries: Sequence[AtlasSpeciesGeometry],
    *,
    min_detectable_species: int,
) -> float:
    """Run the exact edge-rank, cell-mean, shared-surface concentration statistic."""

    if len(colour_vectors) != len(geometries) or not geometries:
        raise ValueError("colour vectors must match a non-empty geometry sequence")
    cell_surfaces: list[np.ndarray] = []
    detectable: list[np.ndarray] = []
    for values, item in zip(colour_vectors, geometries, strict=True):
        values = np.asarray(values, dtype=float)
        if values.shape != (len(item.latitude), 3):
            raise ValueError(f"{item.species}: synthetic vectors have the wrong shape")
        scores = edge_colour_discontinuity(values, item.geometry.retained_edges)
        intensity = average_rank_intensity(scores)
        cell_surfaces.append(cell_mean_intensity(intensity, item.geometry))
        detectable.append(item.geometry.detectable)
    shared, opportunity = shared_boundary_intensity(
        np.vstack(cell_surfaces),
        np.vstack(detectable),
        min_detectable_species=min_detectable_species,
    )
    return opportunity_weighted_concentration(
        shared,
        opportunity,
        min_opportunity=min_detectable_species,
    )


def permutation_p_value(
    colour_vectors: Sequence[np.ndarray],
    geometries: Sequence[AtlasSpeciesGeometry],
    *,
    min_detectable_species: int,
    permutations: int,
    rng: np.random.Generator,
) -> tuple[float, float]:
    """Permute complete synthetic vectors within species and rerun the full surface."""

    if permutations < 1:
        raise ValueError("permutations must be positive")
    observed = atlas_shared_concentration(
        colour_vectors,
        geometries,
        min_detectable_species=min_detectable_species,
    )
    exceed = 0
    for _ in range(permutations):
        permuted = [values[rng.permutation(len(values))] for values in colour_vectors]
        null_value = atlas_shared_concentration(
            permuted,
            geometries,
            min_detectable_species=min_detectable_species,
        )
        exceed += int(null_value >= observed)
    return observed, float((exceed + 1) / (permutations + 1))


def _batched_rank_intensity(scores: np.ndarray) -> np.ndarray:
    """Rank continuous synthetic edge scores by row, preserving rare ties exactly."""

    scores = np.asarray(scores, dtype=float)
    if scores.ndim != 2 or scores.shape[1] < 1 or not np.isfinite(scores).all():
        raise ValueError("batched edge scores must be a finite two-dimensional array")
    n_rows, n_edges = scores.shape
    if n_edges == 1:
        return np.full_like(scores, 0.5)
    order = np.argsort(scores, axis=1, kind="quicksort")
    ranked = np.empty_like(scores)
    base = np.broadcast_to(np.arange(n_edges, dtype=float) / (n_edges - 1), order.shape)
    np.put_along_axis(ranked, order, base, axis=1)

    sorted_scores = np.take_along_axis(scores, order, axis=1)
    tied_rows = np.flatnonzero(np.any(np.diff(sorted_scores, axis=1) == 0.0, axis=1))
    for row in tied_rows:
        ranked[row] = average_rank_intensity(scores[row])
    return ranked


def batched_permutation_p_value(
    colour_vectors: Sequence[np.ndarray],
    geometries: Sequence[AtlasSpeciesGeometry],
    *,
    min_detectable_species: int,
    permutations: int,
    rng: np.random.Generator,
    batch_size: int = 64,
) -> tuple[float, float]:
    """Vectorized equivalent of :func:`permutation_p_value` for formal simulations."""

    if permutations < 1 or batch_size < 1:
        raise ValueError("permutations and batch_size must be positive")
    if len(colour_vectors) != len(geometries) or not geometries:
        raise ValueError("colour vectors must match a non-empty geometry sequence")
    observed = atlas_shared_concentration(
        colour_vectors,
        geometries,
        min_detectable_species=min_detectable_species,
    )

    detectable = np.vstack([item.geometry.detectable for item in geometries])
    opportunity = detectable.sum(axis=0).astype(int)
    shared_cells = np.flatnonzero(opportunity >= min_detectable_species)
    if shared_cells.size < 1:
        raise ValueError("geometry has no shared opportunity cells")
    shared_opportunity = opportunity[shared_cells].astype(float)
    cell_index = np.full(opportunity.size, -1, dtype=int)
    cell_index[shared_cells] = np.arange(shared_cells.size)

    projections: list[sparse.csr_matrix] = []
    for item in geometries:
        geometry = item.geometry
        edge_shared = cell_index[geometry.edge_cell_id]
        valid = edge_shared >= 0
        valid &= geometry.detectable[geometry.edge_cell_id]
        edge_rows = np.flatnonzero(valid)
        columns = edge_shared[valid]
        data = 1.0 / geometry.cell_edge_count[geometry.edge_cell_id[valid]].astype(float)
        projections.append(
            sparse.csr_matrix(
                (data, (edge_rows, columns)),
                shape=(len(geometry.retained_edges), shared_cells.size),
            )
        )

    exceed = 0
    completed = 0
    while completed < permutations:
        current = min(batch_size, permutations - completed)
        numerator = np.zeros((current, shared_cells.size), dtype=float)
        for values, item, projection in zip(
            colour_vectors, geometries, projections, strict=True
        ):
            values = np.asarray(values, dtype=float)
            orders = np.empty((current, len(values)), dtype=int)
            for row in range(current):
                orders[row] = rng.permutation(len(values))
            permuted = values[orders]
            edges = item.geometry.retained_edges
            delta = permuted[:, edges[:, 0], :] - permuted[:, edges[:, 1], :]
            scores = np.sqrt(np.mean(delta * delta, axis=2))
            intensity = _batched_rank_intensity(scores)
            numerator += np.asarray(intensity @ projection)

        shared = numerator / shared_opportunity[None, :]
        weights = shared_opportunity
        means = np.sum(shared * weights[None, :], axis=1) / weights.sum()
        null = np.sum(weights[None, :] * (shared - means[:, None]) ** 2, axis=1)
        null /= weights.sum()
        exceed += int(np.count_nonzero(null >= observed))
        completed += current
    return observed, float((exceed + 1) / (permutations + 1))

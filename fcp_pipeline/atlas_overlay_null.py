"""Spatially constrained concordance null for atlas boundary overlays.

The null removes a fixed broad geographic trend and randomizes only coefficient
signs in a weighted Moran eigenbasis.  This preserves the flower surface's weighted
variance and Moran quadratic form exactly while keeping every environmental or
pollinator overlay fixed and independent of colour.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np
from scipy import linalg


@dataclass(frozen=True)
class MoranSignBasis:
    vectors: np.ndarray
    eigenvalues: np.ndarray
    sqrt_weights: np.ndarray
    weighted_design: np.ndarray


def geographic_design(latitude: np.ndarray, longitude: np.ndarray) -> np.ndarray:
    """Frozen broad-trend basis used before local boundary concordance."""

    lat = np.deg2rad(np.asarray(latitude, dtype=float))
    lon = np.deg2rad(np.asarray(longitude, dtype=float))
    if lat.ndim != 1 or lon.ndim != 1 or lat.shape != lon.shape:
        raise ValueError("latitude and longitude must be equal one-dimensional arrays")
    if not np.isfinite(lat).all() or not np.isfinite(lon).all():
        raise ValueError("coordinates must be finite")
    return np.column_stack(
        (
            np.ones(lat.size),
            np.sin(lat),
            np.sin(lat) ** 2,
            np.sin(lon),
            np.cos(lon),
            np.sin(2.0 * lon),
            np.cos(2.0 * lon),
        )
    )


def equal_area_rook_adjacency(
    cell_ids: Sequence[int] | np.ndarray,
    *,
    n_lon: int,
    n_sinlat: int,
) -> np.ndarray:
    """Rook adjacency with deterministic nearest-cell repair for subset islands."""

    ids = np.asarray(cell_ids, dtype=int)
    if ids.ndim != 1 or ids.size < 3 or len(np.unique(ids)) != ids.size:
        raise ValueError("cell_ids must contain at least three unique cells")
    if n_lon < 2 or n_sinlat < 2 or np.any((ids < 0) | (ids >= n_lon * n_sinlat)):
        raise ValueError("cell IDs or grid dimensions are invalid")
    lookup = {int(cell): index for index, cell in enumerate(ids)}
    adjacency = np.zeros((ids.size, ids.size), dtype=float)
    for index, cell in enumerate(ids):
        row, column = divmod(int(cell), n_lon)
        neighbours = [
            row * n_lon + (column - 1) % n_lon,
            row * n_lon + (column + 1) % n_lon,
        ]
        if row > 0:
            neighbours.append((row - 1) * n_lon + column)
        if row + 1 < n_sinlat:
            neighbours.append((row + 1) * n_lon + column)
        for neighbour in neighbours:
            other = lookup.get(neighbour)
            if other is not None:
                adjacency[index, other] = 1.0
                adjacency[other, index] = 1.0
    isolated = np.flatnonzero(adjacency.sum(axis=1) == 0)
    if isolated.size:
        rows, columns = np.divmod(ids, n_lon)
        sin_latitude = -1.0 + (rows + 0.5) * 2.0 / n_sinlat
        latitude = np.arcsin(np.clip(sin_latitude, -1.0, 1.0))
        longitude = -np.pi + (columns + 0.5) * 2.0 * np.pi / n_lon
        xyz = np.column_stack(
            (
                np.cos(latitude) * np.cos(longitude),
                np.cos(latitude) * np.sin(longitude),
                np.sin(latitude),
            )
        )
        similarity = xyz @ xyz.T
        np.fill_diagonal(similarity, -np.inf)
        for index in isolated:
            other = int(np.argmax(similarity[index]))
            adjacency[index, other] = 1.0
            adjacency[other, index] = 1.0
    return adjacency


def build_moran_sign_basis(
    adjacency: np.ndarray,
    weights: Sequence[float] | np.ndarray,
    design: np.ndarray,
) -> MoranSignBasis:
    """Build a weighted residual Moran basis without inspecting any overlay."""

    adjacency = np.asarray(adjacency, dtype=float)
    weights = np.asarray(weights, dtype=float)
    design = np.asarray(design, dtype=float)
    n = weights.size
    if adjacency.shape != (n, n) or not np.allclose(adjacency, adjacency.T):
        raise ValueError("adjacency must be a symmetric n x n matrix")
    if design.ndim != 2 or design.shape[0] != n or design.shape[1] < 1:
        raise ValueError("design must have n rows and at least one column")
    if (
        not np.isfinite(adjacency).all()
        or not np.isfinite(weights).all()
        or not np.isfinite(design).all()
        or np.any(weights <= 0)
    ):
        raise ValueError("basis inputs must be finite and weights positive")
    if np.any(adjacency < 0) or np.any(np.diag(adjacency) != 0):
        raise ValueError("adjacency must be non-negative without self edges")

    degree = adjacency.sum(axis=1)
    if np.any(degree <= 0):
        raise ValueError("adjacency contains an isolated cell")
    normalized = adjacency / np.sqrt(np.outer(degree, degree))
    sqrt_weights = np.sqrt(weights)
    weighted_design = sqrt_weights[:, None] * design
    q, _ = np.linalg.qr(weighted_design, mode="reduced")
    residual_space = linalg.null_space(q.T)
    if residual_space.shape[1] < 2:
        raise ValueError("too few residual dimensions for a spatial null")
    projected = residual_space.T @ normalized @ residual_space
    eigenvalues, rotation = np.linalg.eigh(projected)
    vectors = residual_space @ rotation
    return MoranSignBasis(
        vectors=vectors,
        eigenvalues=eigenvalues,
        sqrt_weights=sqrt_weights,
        weighted_design=weighted_design,
    )


def residual_coefficients(values: Sequence[float] | np.ndarray, basis: MoranSignBasis) -> np.ndarray:
    """Project one surface into the frozen weighted residual Moran basis."""

    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or values.size != basis.sqrt_weights.size:
        raise ValueError("surface must be one-dimensional and match the basis")
    if not np.isfinite(values).all():
        raise ValueError("surface values must be finite")
    weighted = basis.sqrt_weights * values
    beta, *_ = np.linalg.lstsq(basis.weighted_design, weighted, rcond=None)
    residual = weighted - basis.weighted_design @ beta
    coefficients = basis.vectors.T @ residual
    if float(np.linalg.norm(coefficients)) <= 0:
        raise ValueError("surface has no residual spatial variation")
    return coefficients


def moran_quadratic(coefficients: np.ndarray, basis: MoranSignBasis) -> float:
    coefficients = np.asarray(coefficients, dtype=float)
    if coefficients.shape != basis.eigenvalues.shape:
        raise ValueError("coefficients do not match the Moran basis")
    return float(np.dot(basis.eigenvalues, coefficients * coefficients))


def spectral_family_test(
    flower_surface: Sequence[float] | np.ndarray,
    overlays: Mapping[str, Sequence[float] | np.ndarray],
    basis: MoranSignBasis,
    *,
    randomizations: int,
    rng: np.random.Generator,
) -> dict[str, object]:
    """Test the maximum positive concordance over one frozen overlay family."""

    if not overlays or randomizations < 1:
        raise ValueError("at least one overlay and one randomization are required")
    flower = residual_coefficients(flower_surface, basis)
    flower_norm = float(np.linalg.norm(flower))
    names = tuple(sorted(overlays))
    overlay_coefficients = np.column_stack(
        [residual_coefficients(overlays[name], basis) for name in names]
    )
    overlay_norm = np.linalg.norm(overlay_coefficients, axis=0)
    observed_all = (flower @ overlay_coefficients) / (flower_norm * overlay_norm)
    observed = float(np.max(observed_all))

    null = np.empty(randomizations, dtype=float)
    batch_size = 2048
    product = flower[:, None] * overlay_coefficients
    denominators = flower_norm * overlay_norm
    completed = 0
    while completed < randomizations:
        current = min(batch_size, randomizations - completed)
        signs = rng.integers(0, 2, size=(current, flower.size), dtype=np.int8)
        signs = signs.astype(float) * 2.0 - 1.0
        correlations = (signs @ product) / denominators[None, :]
        null[completed : completed + current] = np.max(correlations, axis=1)
        completed += current
    p_value = float((np.count_nonzero(null >= observed) + 1) / (randomizations + 1))
    return {
        "overlay_names": list(names),
        "observed_by_overlay": {
            name: float(value) for name, value in zip(names, observed_all, strict=True)
        },
        "family_statistic": observed,
        "p_value": p_value,
        "null_statistics": null,
    }

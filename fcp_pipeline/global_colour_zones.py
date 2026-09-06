"""Outcome-rule-frozen extraction of persistent global flower-colour zones.

Zone geography is defined from repeated colour-field stability before any named
mountain, climate, pollinator or biogeographic layer is consulted. Components are
assigned neutral IDs (Z01, Z02, ...) and can only be ecologically annotated later.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np


@dataclass(frozen=True)
class ColourZone:
    zone_id: str
    cell_indices: tuple[int, ...]
    n_cells: int
    integrated_intensity: float
    mean_persistence: float


@dataclass(frozen=True)
class ColourZoneResult:
    persistence: np.ndarray
    consensus_field: np.ndarray
    seed_mask: np.ndarray
    evaluable_resamples: np.ndarray
    zones: tuple[ColourZone, ...]
    hotspot_quantile: float
    minimum_persistence: float
    minimum_evaluable_resamples: int


def _validate_grid_shape(n_cells: int, n_lon: int, n_sinlat: int) -> None:
    if int(n_lon) < 1 or int(n_sinlat) < 1:
        raise ValueError("grid dimensions must be positive")
    if int(n_lon) * int(n_sinlat) != int(n_cells):
        raise ValueError("grid dimensions do not match field cell count")


def persistent_hotspots(
    fields: Sequence[Sequence[float]],
    *,
    hotspot_quantile: float = 0.90,
    minimum_persistence: float = 0.60,
    minimum_evaluable_resamples: int = 100,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return persistence, consensus, seed mask and evaluable counts.

    A cell is a hotspot within one resample if it is in the upper frozen quantile
    among evaluable cells in that resample. Persistence is conditional on the cell
    being evaluable, so missing support never becomes a biological non-hotspot.
    """
    matrix = np.asarray(fields, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] < 1 or matrix.shape[1] < 1:
        raise ValueError("fields must be a non-empty resamples x cells matrix")
    hotspot_quantile = float(hotspot_quantile)
    minimum_persistence = float(minimum_persistence)
    minimum_evaluable_resamples = int(minimum_evaluable_resamples)
    if not 0 < hotspot_quantile < 1:
        raise ValueError("hotspot_quantile must lie in (0, 1)")
    if not 0 < minimum_persistence <= 1:
        raise ValueError("minimum_persistence must lie in (0, 1]")
    if minimum_evaluable_resamples < 1:
        raise ValueError("minimum_evaluable_resamples must be positive")

    evaluable = np.isfinite(matrix)
    hotspot = np.zeros(matrix.shape, dtype=bool)
    for row in range(matrix.shape[0]):
        idx = np.flatnonzero(evaluable[row])
        if len(idx) == 0:
            continue
        threshold = float(np.quantile(matrix[row, idx], hotspot_quantile))
        hotspot[row, idx] = matrix[row, idx] >= threshold

    opportunities = evaluable.sum(axis=0).astype(np.int64)
    hotspot_count = hotspot.sum(axis=0).astype(np.int64)
    persistence = np.full(matrix.shape[1], np.nan, dtype=float)
    keep = opportunities > 0
    persistence[keep] = hotspot_count[keep] / opportunities[keep]
    consensus = np.full(matrix.shape[1], np.nan, dtype=float)
    for cell in np.flatnonzero(keep):
        consensus[cell] = float(np.mean(matrix[evaluable[:, cell], cell]))
    seed = (
        (opportunities >= minimum_evaluable_resamples)
        & np.isfinite(persistence)
        & (persistence >= minimum_persistence)
    )
    return persistence, consensus, seed, opportunities


def _neighbours(cell: int, *, n_lon: int, n_sinlat: int) -> tuple[int, ...]:
    row = int(cell) // int(n_lon)
    col = int(cell) % int(n_lon)
    out: list[int] = []
    for dr in (-1, 0, 1):
        rr = row + dr
        if rr < 0 or rr >= n_sinlat:
            continue
        for dc in (-1, 0, 1):
            if dr == 0 and dc == 0:
                continue
            cc = (col + dc) % n_lon
            out.append(rr * n_lon + cc)
    return tuple(sorted(set(out)))


def connected_components(
    mask: Sequence[bool],
    *,
    n_lon: int,
    n_sinlat: int,
) -> tuple[tuple[int, ...], ...]:
    """8-neighbour components with longitude wrap across the dateline."""
    active = np.asarray(mask, dtype=bool)
    if active.ndim != 1:
        raise ValueError("mask must be one-dimensional")
    _validate_grid_shape(len(active), n_lon, n_sinlat)
    unseen = set(np.flatnonzero(active).tolist())
    components: list[tuple[int, ...]] = []
    while unseen:
        start = min(unseen)
        stack = [start]
        unseen.remove(start)
        component: list[int] = []
        while stack:
            current = stack.pop()
            component.append(current)
            for neighbour in _neighbours(current, n_lon=n_lon, n_sinlat=n_sinlat):
                if neighbour in unseen:
                    unseen.remove(neighbour)
                    stack.append(neighbour)
        components.append(tuple(sorted(component)))
    return tuple(components)


def extract_persistent_colour_zones(
    fields: Sequence[Sequence[float]],
    *,
    n_lon: int,
    n_sinlat: int,
    opportunity: Sequence[float],
    hotspot_quantile: float = 0.90,
    minimum_persistence: float = 0.60,
    minimum_evaluable_resamples: int = 100,
    minimum_zone_cells: int = 3,
) -> ColourZoneResult:
    """Extract neutral stable colour-zone IDs under the frozen primary rule."""
    matrix = np.asarray(fields, dtype=float)
    _validate_grid_shape(matrix.shape[1], n_lon, n_sinlat)
    weights = np.asarray(opportunity, dtype=float)
    if weights.shape != (matrix.shape[1],):
        raise ValueError("opportunity shape does not match field cells")
    if np.any(~np.isfinite(weights)) or np.any(weights < 0):
        raise ValueError("opportunity must be finite and non-negative")
    minimum_zone_cells = int(minimum_zone_cells)
    if minimum_zone_cells < 1:
        raise ValueError("minimum_zone_cells must be positive")

    persistence, consensus, seed, evaluable_counts = persistent_hotspots(
        matrix,
        hotspot_quantile=hotspot_quantile,
        minimum_persistence=minimum_persistence,
        minimum_evaluable_resamples=minimum_evaluable_resamples,
    )
    components = [
        component
        for component in connected_components(seed, n_lon=n_lon, n_sinlat=n_sinlat)
        if len(component) >= minimum_zone_cells
    ]
    scored: list[tuple[float, float, tuple[int, ...]]] = []
    for component in components:
        idx = np.asarray(component, dtype=int)
        finite = np.isfinite(consensus[idx])
        if not np.any(finite):
            continue
        idx = idx[finite]
        w = weights[idx]
        intensity = float(np.sum(consensus[idx] * w))
        mean_persistence = float(np.mean(persistence[idx]))
        scored.append((intensity, mean_persistence, tuple(component)))
    scored.sort(key=lambda item: (-item[0], -item[1], item[2]))
    zones = tuple(
        ColourZone(
            zone_id=f"Z{position:02d}",
            cell_indices=component,
            n_cells=len(component),
            integrated_intensity=float(intensity),
            mean_persistence=float(mean_persistence),
        )
        for position, (intensity, mean_persistence, component) in enumerate(scored, start=1)
    )
    return ColourZoneResult(
        persistence=persistence,
        consensus_field=consensus,
        seed_mask=seed,
        evaluable_resamples=evaluable_counts,
        zones=zones,
        hotspot_quantile=float(hotspot_quantile),
        minimum_persistence=float(minimum_persistence),
        minimum_evaluable_resamples=int(minimum_evaluable_resamples),
    )


__all__ = [
    "ColourZone",
    "ColourZoneResult",
    "connected_components",
    "extract_persistent_colour_zones",
    "persistent_hotspots",
]

"""Label-blind geometry and gridding for shared continuous transition surfaces.

The functions in this module keep the geometry/opportunity layer separate from flower
colour.  A species graph is built elsewhere from observation coordinates only.  Its
local edges are assigned to a fixed equal-area longitude–sin(latitude) grid by spherical
midpoint.  Detectability is then determined only by the number of retained geometry
edges supporting a species/cell combination.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

EARTH_RADIUS_KM = 6371.0088


@dataclass(frozen=True)
class EqualAreaGrid:
    """Rectangular equal-area grid in longitude and sin(latitude)."""

    n_lon: int
    n_sinlat: int

    def __post_init__(self) -> None:
        if not isinstance(self.n_lon, (int, np.integer)) or self.n_lon < 2:
            raise ValueError("n_lon must be an integer >= 2")
        if not isinstance(self.n_sinlat, (int, np.integer)) or self.n_sinlat < 2:
            raise ValueError("n_sinlat must be an integer >= 2")

    @property
    def n_cells(self) -> int:
        return int(self.n_lon * self.n_sinlat)

    @property
    def cell_area_km2(self) -> float:
        return float(4.0 * np.pi * EARTH_RADIUS_KM**2 / self.n_cells)


@dataclass(frozen=True)
class EdgeCellGeometry:
    """Fixed label-independent edge-to-cell support for one species/configuration."""

    retained_edge_indices: np.ndarray
    retained_edges: np.ndarray
    retained_edge_distance_km: np.ndarray
    edge_cell_id: np.ndarray
    cell_edge_count: np.ndarray
    detectable: np.ndarray

    def __post_init__(self) -> None:
        edge_indices = np.asarray(self.retained_edge_indices)
        edges = np.asarray(self.retained_edges)
        distance = np.asarray(self.retained_edge_distance_km, dtype=float)
        cell_id = np.asarray(self.edge_cell_id)
        counts = np.asarray(self.cell_edge_count)
        detectable = np.asarray(self.detectable, dtype=bool)
        if edge_indices.ndim != 1:
            raise ValueError("retained_edge_indices must be one-dimensional")
        if edges.ndim != 2 or edges.shape[1] != 2:
            raise ValueError("retained_edges must have shape (n_edges, 2)")
        if distance.ndim != 1 or cell_id.ndim != 1:
            raise ValueError("distance and edge_cell_id must be one-dimensional")
        if not (len(edge_indices) == len(edges) == len(distance) == len(cell_id)):
            raise ValueError("retained edge arrays must have equal length")
        if counts.ndim != 1 or detectable.ndim != 1 or counts.shape != detectable.shape:
            raise ValueError("cell_edge_count and detectable must be equal-length vectors")
        if len(edges) == 0:
            raise ValueError("at least one local edge must be retained")
        if np.any(distance < 0) or not np.isfinite(distance).all():
            raise ValueError("retained edge distances must be finite and non-negative")
        if np.any(cell_id < 0) or np.any(cell_id >= len(counts)):
            raise ValueError("edge cell id is outside the grid")
        if int(counts.sum()) != len(edges):
            raise ValueError("cell edge counts do not sum to retained edge count")


def _validate_coordinates(latitude: np.ndarray, longitude: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    lat = np.asarray(latitude, dtype=float)
    lon = np.asarray(longitude, dtype=float)
    if lat.ndim != 1 or lon.ndim != 1 or lat.shape != lon.shape:
        raise ValueError("latitude and longitude must be equal-length one-dimensional arrays")
    if len(lat) < 1 or not np.isfinite(lat).all() or not np.isfinite(lon).all():
        raise ValueError("coordinates must be non-empty and finite")
    if np.any((lat < -90.0) | (lat > 90.0)):
        raise ValueError("latitude must lie in [-90, 90]")
    if np.any((lon < -180.0) | (lon > 180.0)):
        raise ValueError("longitude must lie in [-180, 180]")
    return lat, lon


def _unit_sphere_xyz(latitude: np.ndarray, longitude: np.ndarray) -> np.ndarray:
    lat = np.deg2rad(latitude)
    lon = np.deg2rad(longitude)
    cos_lat = np.cos(lat)
    return np.column_stack((cos_lat * np.cos(lon), cos_lat * np.sin(lon), np.sin(lat)))


def equal_area_cell_ids(
    latitude: np.ndarray,
    longitude: np.ndarray,
    grid: EqualAreaGrid,
) -> np.ndarray:
    """Assign coordinates to equal-area longitude–sin(latitude) cells."""

    lat, lon = _validate_coordinates(latitude, longitude)
    lon_fraction = np.mod(lon + 180.0, 360.0) / 360.0
    col = np.floor(lon_fraction * grid.n_lon).astype(int)
    col = np.clip(col, 0, grid.n_lon - 1)

    sin_fraction = 0.5 * (np.sin(np.deg2rad(lat)) + 1.0)
    row = np.floor(sin_fraction * grid.n_sinlat).astype(int)
    row = np.clip(row, 0, grid.n_sinlat - 1)
    return row * grid.n_lon + col


def equal_area_cell_centers(grid: EqualAreaGrid) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return cell ids and latitude/longitude centres for a grid."""

    cell_id = np.arange(grid.n_cells, dtype=int)
    row = cell_id // grid.n_lon
    col = cell_id % grid.n_lon
    sin_lat = -1.0 + (row + 0.5) * (2.0 / grid.n_sinlat)
    latitude = np.rad2deg(np.arcsin(np.clip(sin_lat, -1.0, 1.0)))
    longitude = -180.0 + (col + 0.5) * (360.0 / grid.n_lon)
    return cell_id, latitude.astype(float), longitude.astype(float)


def spherical_edge_midpoints(
    latitude: np.ndarray,
    longitude: np.ndarray,
    edges: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return great-circle midpoints for non-antipodal edges."""

    lat, lon = _validate_coordinates(latitude, longitude)
    edges = np.asarray(edges)
    if edges.ndim != 2 or edges.shape[1] != 2:
        raise ValueError("edges must have shape (n_edges, 2)")
    if not np.issubdtype(edges.dtype, np.integer):
        if np.any(edges != np.floor(edges)):
            raise ValueError("edge indices must be integers")
        edges = edges.astype(int)
    else:
        edges = edges.astype(int, copy=False)
    if len(edges) == 0:
        return np.empty(0, dtype=float), np.empty(0, dtype=float)
    if edges.min() < 0 or edges.max() >= len(lat):
        raise ValueError("edge index is out of bounds")
    if np.any(edges[:, 0] == edges[:, 1]):
        raise ValueError("self edges are not allowed")

    xyz = _unit_sphere_xyz(lat, lon)
    midpoint = xyz[edges[:, 0]] + xyz[edges[:, 1]]
    norm = np.linalg.norm(midpoint, axis=1)
    if np.any(norm <= 1e-12):
        raise ValueError("antipodal edge has no unique spherical midpoint")
    midpoint = midpoint / norm[:, None]
    midpoint_lat = np.rad2deg(np.arcsin(np.clip(midpoint[:, 2], -1.0, 1.0)))
    midpoint_lon = np.rad2deg(np.arctan2(midpoint[:, 1], midpoint[:, 0]))
    return midpoint_lat, midpoint_lon


def build_edge_cell_geometry(
    latitude: np.ndarray,
    longitude: np.ndarray,
    edges: np.ndarray,
    edge_distance_km: np.ndarray,
    *,
    grid: EqualAreaGrid,
    max_edge_km: float,
    min_edges_per_cell: int,
) -> EdgeCellGeometry:
    """Freeze local edge support and detectability without using colour values."""

    lat, lon = _validate_coordinates(latitude, longitude)
    edges = np.asarray(edges, dtype=int)
    distance = np.asarray(edge_distance_km, dtype=float)
    if edges.ndim != 2 or edges.shape[1] != 2 or distance.ndim != 1 or len(edges) != len(distance):
        raise ValueError("edges and edge_distance_km must describe the same edge vector")
    max_edge_km = float(max_edge_km)
    if not np.isfinite(max_edge_km) or max_edge_km <= 0:
        raise ValueError("max_edge_km must be finite and positive")
    if not isinstance(min_edges_per_cell, (int, np.integer)) or min_edges_per_cell < 1:
        raise ValueError("min_edges_per_cell must be a positive integer")
    if np.any(distance < 0) or not np.isfinite(distance).all():
        raise ValueError("edge distances must be finite and non-negative")

    keep = distance <= max_edge_km
    retained_indices = np.flatnonzero(keep)
    retained_edges = edges[keep]
    retained_distance = distance[keep]
    if len(retained_edges) == 0:
        raise ValueError("max_edge_km removed every edge")
    midpoint_lat, midpoint_lon = spherical_edge_midpoints(lat, lon, retained_edges)
    edge_cell_id = equal_area_cell_ids(midpoint_lat, midpoint_lon, grid)
    cell_count = np.bincount(edge_cell_id, minlength=grid.n_cells).astype(int)
    detectable = cell_count >= int(min_edges_per_cell)
    return EdgeCellGeometry(
        retained_edge_indices=retained_indices,
        retained_edges=retained_edges,
        retained_edge_distance_km=retained_distance,
        edge_cell_id=edge_cell_id,
        cell_edge_count=cell_count,
        detectable=detectable,
    )


def cell_mean_intensity(
    edge_intensity: np.ndarray,
    geometry: EdgeCellGeometry,
) -> np.ndarray:
    """Average edge intensities within fixed detectable cells."""

    intensity = np.asarray(edge_intensity, dtype=float)
    if intensity.ndim != 1 or len(intensity) != len(geometry.retained_edges):
        raise ValueError("edge_intensity must match the retained geometry edges")
    if not np.isfinite(intensity).all() or np.any((intensity < 0) | (intensity > 1)):
        raise ValueError("edge intensities must be finite and lie in [0, 1]")
    numerator = np.bincount(
        geometry.edge_cell_id,
        weights=intensity,
        minlength=len(geometry.cell_edge_count),
    )
    out = np.full(len(geometry.cell_edge_count), np.nan, dtype=float)
    out[geometry.detectable] = numerator[geometry.detectable] / geometry.cell_edge_count[geometry.detectable]
    return out


def geometry_opportunity_summary(
    geometries: list[EdgeCellGeometry],
    *,
    min_detectable_species: int,
) -> dict[str, object]:
    """Summarize cross-species opportunity using only fixed geometry masks."""

    if len(geometries) == 0:
        raise ValueError("at least one species geometry is required")
    detectable = np.vstack([g.detectable for g in geometries])
    opportunity = detectable.sum(axis=0).astype(int)
    min_detectable_species = int(min_detectable_species)
    if min_detectable_species < 1:
        raise ValueError("min_detectable_species must be positive")
    shared_cells = opportunity >= min_detectable_species
    species_in_shared_cells = np.any(detectable[:, shared_cells], axis=1) if np.any(shared_cells) else np.zeros(len(geometries), dtype=bool)
    return {
        "n_cells": int(detectable.shape[1]),
        "n_cells_A_ge_1": int(np.count_nonzero(opportunity >= 1)),
        "n_cells_A_ge_2": int(np.count_nonzero(opportunity >= 2)),
        "n_cells_A_ge_3": int(np.count_nonzero(opportunity >= 3)),
        "n_cells_A_ge_4": int(np.count_nonzero(opportunity >= 4)),
        "max_A": int(opportunity.max(initial=0)),
        "species_with_any_shared_opportunity": int(np.count_nonzero(species_in_shared_cells)),
        "retained_edges_per_species": [int(len(g.retained_edges)) for g in geometries],
        "detectable_cells_per_species": [int(np.count_nonzero(g.detectable)) for g in geometries],
    }

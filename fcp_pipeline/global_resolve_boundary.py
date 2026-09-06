"""Exact pre-colour RESOLVE 2017 boundary scoring for the global RGFCA.

Scientific weights and path sampling are defined in the pre-outcome edge-mechanism
contract.  This module only implements the frozen technical mapping in
``global_rgfca_resolve_edge_mechanism_execution_v1.json``.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

EARTH_RADIUS_KM = 6371.0088


@dataclass(frozen=True)
class ResolveBoundaryIndex:
    """Spatial index plus deterministic RESOLVE classification attributes."""

    tree: object
    realm: np.ndarray
    biome_num: np.ndarray
    eco_id: np.ndarray
    source_path: str
    polygon_count: int


def _unit_xyz(latitude: float, longitude: float) -> np.ndarray:
    lat = np.deg2rad(float(latitude))
    lon = np.deg2rad(float(longitude))
    c = np.cos(lat)
    return np.asarray([c * np.cos(lon), c * np.sin(lon), np.sin(lat)], dtype=float)


def great_circle_samples(
    latitude_a: float,
    longitude_a: float,
    latitude_b: float,
    longitude_b: float,
    *,
    maximum_interval_km: float = 10.0,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Sample the shortest great-circle path at intervals no greater than the cap.

    Endpoints are always included.  Spherical linear interpolation is used except
    for coincident coordinates, where normalized linear interpolation is exactly
    equivalent for the zero-length path.  Exact antipodes have no unique shortest
    great-circle path and fail closed as missing external coverage.
    """
    maximum_interval_km = float(maximum_interval_km)
    if not np.isfinite(maximum_interval_km) or maximum_interval_km <= 0:
        raise ValueError("maximum_interval_km must be finite and positive")
    a = _unit_xyz(latitude_a, longitude_a)
    b = _unit_xyz(latitude_b, longitude_b)
    dot = float(np.clip(np.dot(a, b), -1.0, 1.0))
    angle = float(np.arccos(dot))
    distance_km = float(angle * EARTH_RADIUS_KM)
    n_intervals = max(1, int(np.ceil(distance_km / maximum_interval_km)))
    t = np.linspace(0.0, 1.0, n_intervals + 1, dtype=float)

    if angle <= 1e-14:
        xyz = np.repeat(a[None, :], n_intervals + 1, axis=0)
    else:
        sine = float(np.sin(angle))
        if abs(sine) <= 1e-12:
            # At (near-)antipodal endpoints the shortest great-circle arc is not
            # uniquely determined by endpoints alone.  The frozen analysis does not
            # permit choosing a favourable meridian, so fail this edge closed.
            return np.asarray([], dtype=float), np.asarray([], dtype=float), distance_km
        left = np.sin((1.0 - t) * angle) / sine
        right = np.sin(t * angle) / sine
        xyz = left[:, None] * a[None, :] + right[:, None] * b[None, :]
        norm = np.linalg.norm(xyz, axis=1)
        if np.any(norm <= 0) or not np.all(np.isfinite(norm)):
            return np.asarray([], dtype=float), np.asarray([], dtype=float), distance_km
        xyz = xyz / norm[:, None]

    latitude = np.rad2deg(np.arcsin(np.clip(xyz[:, 2], -1.0, 1.0)))
    longitude = np.rad2deg(np.arctan2(xyz[:, 1], xyz[:, 0]))
    return latitude.astype(float), longitude.astype(float), distance_km


def load_resolve_boundary_index(shapefile_path: str | Path) -> ResolveBoundaryIndex:
    """Load the exact verified RESOLVE shapefile with deterministic overlap order."""
    try:
        import shapefile
        from shapely.geometry import shape as shapely_shape
        from shapely.strtree import STRtree
    except ImportError as exc:  # pragma: no cover - exercised in the dedicated workflow
        raise RuntimeError("RESOLVE scoring requires pyshp and shapely") from exc

    path = Path(shapefile_path)
    if not path.exists():
        raise FileNotFoundError(path)
    reader = shapefile.Reader(str(path), encoding="latin1")
    field_names = [field[0] for field in reader.fields[1:]]
    required = {"REALM", "BIOME_NUM", "ECO_ID"}
    missing = sorted(required - set(field_names))
    if missing:
        raise RuntimeError(f"RESOLVE shapefile lacks required fields: {missing}")
    field_index = {name: field_names.index(name) for name in required}

    entries: list[tuple[int, int, str, object]] = []
    for shape_record in reader.iterShapeRecords():
        record = shape_record.record
        realm = str(record[field_index["REALM"]])
        biome_num = int(float(record[field_index["BIOME_NUM"]]))
        eco_id = int(float(record[field_index["ECO_ID"]]))
        geom = shapely_shape(shape_record.shape.__geo_interface__)
        if geom.is_empty:
            continue
        entries.append((eco_id, biome_num, realm, geom))
    if not entries:
        raise RuntimeError("RESOLVE shapefile yielded no polygons")

    # Sorting establishes the pre-frozen overlap tie-break. STRtree result indices
    # therefore have exactly the same priority order.
    entries.sort(key=lambda item: (item[0], item[1], item[2]))
    geoms = [item[3] for item in entries]
    tree = STRtree(geoms)
    return ResolveBoundaryIndex(
        tree=tree,
        realm=np.asarray([item[2] for item in entries], dtype=object),
        biome_num=np.asarray([item[1] for item in entries], dtype=np.int32),
        eco_id=np.asarray([item[0] for item in entries], dtype=np.int32),
        source_path=str(path),
        polygon_count=len(entries),
    )


def _classify_points(
    index: ResolveBoundaryIndex,
    latitude: np.ndarray,
    longitude: np.ndarray,
) -> np.ndarray:
    """Return deterministic RESOLVE polygon indices; -1 means external missingness."""
    try:
        from shapely import points
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("RESOLVE scoring requires shapely") from exc

    lat = np.asarray(latitude, dtype=float)
    lon = np.asarray(longitude, dtype=float)
    if lat.shape != lon.shape or lat.ndim != 1:
        raise ValueError("latitude and longitude must be matching 1D arrays")
    if not np.all(np.isfinite(lat)) or not np.all(np.isfinite(lon)):
        raise ValueError("classification coordinates must be finite")
    out = np.full(len(lat), -1, dtype=np.int32)
    if len(lat) == 0:
        return out

    point_geometries = points(lon, lat)
    matches = np.asarray(index.tree.query(point_geometries, predicate="covered_by"), dtype=np.int64)
    if matches.size == 0:
        return out
    if matches.ndim != 2 or matches.shape[0] != 2:
        raise RuntimeError("unexpected STRtree query result shape")
    input_index = matches[0]
    tree_index = matches[1]
    # Tree indices follow the pre-sorted ECO_ID/BIOME_NUM/REALM priority.  Sort by
    # sampled-point index and then tree index, and keep only the first match.
    order = np.lexsort((tree_index, input_index))
    ii = input_index[order]
    jj = tree_index[order]
    first = np.ones(len(ii), dtype=bool)
    first[1:] = ii[1:] != ii[:-1]
    out[ii[first]] = jj[first].astype(np.int32)
    return out


def score_resolve_boundary_pairs(
    pair_rows: Sequence[Sequence[int]],
    latitude_by_row: Sequence[float],
    longitude_by_row: Sequence[float],
    index: ResolveBoundaryIndex,
    *,
    maximum_interval_km: float = 10.0,
    maximum_points_per_batch: int = 200_000,
) -> np.ndarray:
    """Score unique unordered photo pairs under the exact frozen RESOLVE formula.

    Any sampled point outside the terrestrial RESOLVE polygons makes that edge
    externally non-evaluable (NaN). It is never converted to a zero boundary score.
    """
    pairs = np.asarray(pair_rows, dtype=np.int64)
    lat = np.asarray(latitude_by_row, dtype=float)
    lon = np.asarray(longitude_by_row, dtype=float)
    if pairs.ndim != 2 or pairs.shape[1] != 2:
        raise ValueError("pair_rows must have shape (n_pairs,2)")
    if lat.ndim != 1 or lon.shape != lat.shape:
        raise ValueError("coordinate arrays must be matching one-dimensional arrays")
    if len(pairs) and (pairs.min() < 0 or pairs.max() >= len(lat)):
        raise ValueError("pair row index outside coordinate arrays")
    maximum_points_per_batch = int(maximum_points_per_batch)
    if maximum_points_per_batch < 2:
        raise ValueError("maximum_points_per_batch must be at least two")

    scores = np.full(len(pairs), np.nan, dtype=float)
    batch_pair_indices: list[int] = []
    batch_latitude: list[np.ndarray] = []
    batch_longitude: list[np.ndarray] = []
    batch_slices: list[tuple[int, int]] = []
    batch_points = 0

    def flush() -> None:
        nonlocal batch_pair_indices, batch_latitude, batch_longitude, batch_slices, batch_points
        if not batch_pair_indices:
            return
        all_lat = np.concatenate(batch_latitude)
        all_lon = np.concatenate(batch_longitude)
        classified = _classify_points(index, all_lat, all_lon)
        for pair_index, (start, stop) in zip(batch_pair_indices, batch_slices):
            cls = classified[start:stop]
            if len(cls) < 2 or np.any(cls < 0):
                continue
            realm = index.realm[cls]
            biome = index.biome_num[cls]
            eco = index.eco_id[cls]
            value = (
                3 * int(np.count_nonzero(realm[1:] != realm[:-1]))
                + 2 * int(np.count_nonzero(biome[1:] != biome[:-1]))
                + int(np.count_nonzero(eco[1:] != eco[:-1]))
            )
            scores[pair_index] = float(value)
        batch_pair_indices = []
        batch_latitude = []
        batch_longitude = []
        batch_slices = []
        batch_points = 0

    for pair_index, (left, right) in enumerate(pairs):
        sample_lat, sample_lon, _ = great_circle_samples(
            lat[int(left)], lon[int(left)], lat[int(right)], lon[int(right)],
            maximum_interval_km=float(maximum_interval_km),
        )
        # Exact/near antipodes have no unique frozen shortest path and remain NaN.
        if len(sample_lat) == 0:
            continue
        if batch_pair_indices and batch_points + len(sample_lat) > maximum_points_per_batch:
            flush()
        start = batch_points
        stop = start + len(sample_lat)
        batch_pair_indices.append(int(pair_index))
        batch_latitude.append(sample_lat)
        batch_longitude.append(sample_lon)
        batch_slices.append((start, stop))
        batch_points = stop
        # One extremely long edge may exceed the nominal batch cap; process it alone.
        if batch_points >= maximum_points_per_batch:
            flush()
    flush()
    return scores


__all__ = [
    "EARTH_RADIUS_KM",
    "ResolveBoundaryIndex",
    "great_circle_samples",
    "load_resolve_boundary_index",
    "score_resolve_boundary_pairs",
]

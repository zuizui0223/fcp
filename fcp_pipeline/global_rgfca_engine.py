"""Scalable G1 engine for the Repeated Global Flower-Colour Atlas (RGFCA).

The engine consumes an already measured, classifiable, location-blind photo pool.
It does not acquire images, choose thresholds, or inspect external mechanisms.

Computational representation is optimized without changing the frozen statistic:
* complete photo-level four-group colour vectors are the exchangeability unit;
* one species-conditioned null assignment is reused across all outer realizations;
* pairwise colour JSD is precomputed once within each species;
* compact-support edge kernels are stored sparsely;
* all nulls assigned to one execution shard are evaluated together per outer map.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata

from .global_barrier_field import (
    BarrierGrid,
    equal_area_grid_centers,
    spherical_midpoint,
)
from .global_repeated_atlas import (
    RepeatedAtlasSchedule,
    build_repeated_atlas_schedule,
    consensus_field,
    stable_seed,
)
from .spatial_graph import spherical_knn_edges

EARTH_RADIUS_KM = 6371.0088
COLOUR_COLUMNS = ("colour_white", "colour_yellow_orange", "colour_red_pink", "colour_blue_purple")


@dataclass(frozen=True)
class CanonicalColourPool:
    frame: pd.DataFrame
    photo_ids: np.ndarray
    species: np.ndarray
    latitude: np.ndarray
    longitude: np.ndarray
    colours: np.ndarray
    species_labels: tuple[str, ...]
    species_slices: tuple[tuple[int, int], ...]


@dataclass(frozen=True)
class PairwiseJSDCache:
    flat_jsd: np.ndarray
    local_index_by_row: np.ndarray
    species_size_by_row: np.ndarray
    matrix_offset_by_row: np.ndarray


@dataclass(frozen=True)
class SparseOuterGeometry:
    grid: BarrierGrid
    edge_nodes: np.ndarray
    edge_species_index: np.ndarray
    edge_species_slices: tuple[tuple[int, int], ...]
    weighted_kernel: sparse.csr_matrix
    opportunity: np.ndarray
    distinct_species_support: np.ndarray
    evaluable: np.ndarray


@dataclass(frozen=True)
class G1ShardResult:
    schedule: RepeatedAtlasSchedule
    observed_outer_fields: np.ndarray
    observed_outer_opportunities: np.ndarray
    observed_consensus: np.ndarray
    observed_aggregate_opportunity: np.ndarray
    observed_concentration: float
    null_indices: np.ndarray
    null_consensus_fields: np.ndarray
    null_concentrations: np.ndarray
    null_unit: str


def canonical_colour_pool(frame: pd.DataFrame) -> CanonicalColourPool:
    required = {"photo_id", "species", "latitude", "longitude", *COLOUR_COLUMNS}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"measured pool lacks required columns: {missing}")
    work = frame[list(required)].copy()
    if work.isna().any().any():
        raise ValueError("canonical measured pool cannot contain missing required values")
    work["photo_id"] = pd.to_numeric(work["photo_id"], errors="raise").astype(np.int64)
    work["species"] = work["species"].astype(str).str.strip()
    if (work["species"] == "").any():
        raise ValueError("species labels must be non-empty")
    if work["photo_id"].duplicated().any():
        raise ValueError("photo_id must be globally unique")
    work["latitude"] = pd.to_numeric(work["latitude"], errors="raise").astype(float)
    work["longitude"] = pd.to_numeric(work["longitude"], errors="raise").astype(float)
    if not np.isfinite(work[["latitude", "longitude"]].to_numpy()).all():
        raise ValueError("coordinates must be finite")
    if ((work["latitude"] < -90) | (work["latitude"] > 90)).any():
        raise ValueError("latitude outside [-90,90]")
    if ((work["longitude"] < -180) | (work["longitude"] > 180)).any():
        raise ValueError("longitude outside [-180,180]")

    colours = work[list(COLOUR_COLUMNS)].apply(pd.to_numeric, errors="raise").to_numpy(float)
    if not np.isfinite(colours).all() or np.any(colours < 0):
        raise ValueError("colour masses must be finite and non-negative")
    mass = colours.sum(axis=1)
    if np.any(mass <= 0):
        raise ValueError("each classifiable photo must have positive four-group colour mass")
    colours = colours / mass[:, None]
    work.loc[:, list(COLOUR_COLUMNS)] = colours
    work = work.sort_values(["species", "photo_id"], kind="mergesort").reset_index(drop=True)

    species = work["species"].to_numpy(dtype=object)
    labels = tuple(pd.unique(species).tolist())
    slices: list[tuple[int, int]] = []
    for label in labels:
        idx = np.flatnonzero(species == label)
        if len(idx) == 0 or np.any(np.diff(idx) != 1):
            raise RuntimeError("canonical species rows are not contiguous")
        slices.append((int(idx[0]), int(idx[-1] + 1)))
    return CanonicalColourPool(
        frame=work,
        photo_ids=work["photo_id"].to_numpy(np.int64),
        species=species,
        latitude=work["latitude"].to_numpy(float),
        longitude=work["longitude"].to_numpy(float),
        colours=work[list(COLOUR_COLUMNS)].to_numpy(float),
        species_labels=labels,
        species_slices=tuple(slices),
    )


def _pairwise_jsd(probabilities: np.ndarray) -> np.ndarray:
    p = np.asarray(probabilities, dtype=float)
    if p.ndim != 2 or p.shape[1] != 4:
        raise ValueError("probabilities must have shape (n,4)")
    a = p[:, None, :]
    b = p[None, :, :]
    m = 0.5 * (a + b)
    with np.errstate(divide="ignore", invalid="ignore"):
        term_a = np.where(a > 0, a * np.log2(a / m), 0.0)
        term_b = np.where(b > 0, b * np.log2(b / m), 0.0)
    out = 0.5 * term_a.sum(axis=2) + 0.5 * term_b.sum(axis=2)
    return np.clip(out, 0.0, 1.0)


def build_pairwise_jsd_cache(pool: CanonicalColourPool) -> PairwiseJSDCache:
    n_rows = len(pool.photo_ids)
    local = np.empty(n_rows, dtype=np.int32)
    size = np.empty(n_rows, dtype=np.int32)
    offset = np.empty(n_rows, dtype=np.int64)
    pieces: list[np.ndarray] = []
    cursor = 0
    for start, stop in pool.species_slices:
        n = stop - start
        matrix = _pairwise_jsd(pool.colours[start:stop])
        pieces.append(matrix.ravel())
        local[start:stop] = np.arange(n, dtype=np.int32)
        size[start:stop] = n
        offset[start:stop] = cursor
        cursor += n * n
    return PairwiseJSDCache(
        flat_jsd=np.concatenate(pieces).astype(np.float64, copy=False),
        local_index_by_row=local,
        species_size_by_row=size,
        matrix_offset_by_row=offset,
    )


def _xyz(latitude: np.ndarray, longitude: np.ndarray) -> np.ndarray:
    lat = np.deg2rad(latitude)
    lon = np.deg2rad(longitude)
    c = np.cos(lat)
    return np.column_stack([c * np.cos(lon), c * np.sin(lon), np.sin(lat)])


def _distance_matrix_km(latitude: np.ndarray, longitude: np.ndarray, grid: BarrierGrid) -> np.ndarray:
    a = _xyz(latitude, longitude)
    b = _xyz(grid.latitude, grid.longitude)
    dot = np.clip(a @ b.T, -1.0, 1.0)
    return np.arccos(dot) * EARTH_RADIUS_KM


def prepare_sparse_outer_geometry(
    pool: CanonicalColourPool,
    photo_ids_by_species: Sequence[Sequence[int]],
    species_labels: Sequence[object],
    *,
    grid: BarrierGrid | None = None,
    k: int = 3,
    kernel_km: float = 500.0,
    cutoff_multiplier: float = 3.0,
    minimum_distinct_species: int = 5,
) -> SparseOuterGeometry:
    """Build one fixed outer geometry with sparse compact-support kernels."""
    labels = [str(x) for x in species_labels]
    draws = np.asarray(photo_ids_by_species, dtype=np.int64)
    if draws.ndim != 2 or draws.shape[0] != len(labels):
        raise ValueError("photo_ids_by_species must be species x photos")
    if draws.shape[1] <= int(k):
        raise ValueError("each outer species draw needs more photos than k")
    if len(set(labels)) != len(labels):
        raise ValueError("outer species labels must be unique")
    if grid is None:
        grid = equal_area_grid_centers(36, 18)
    kernel_km = float(kernel_km)
    cutoff_multiplier = float(cutoff_multiplier)
    minimum_distinct_species = int(minimum_distinct_species)
    if kernel_km <= 0 or cutoff_multiplier <= 0 or minimum_distinct_species < 1:
        raise ValueError("kernel/support parameters must be positive")

    id_to_row = {int(pid): i for i, pid in enumerate(pool.photo_ids)}
    label_to_index = {label: i for i, label in enumerate(pool.species_labels)}
    edge_nodes: list[np.ndarray] = []
    edge_species_index: list[np.ndarray] = []
    mid_lat: list[np.ndarray] = []
    mid_lon: list[np.ndarray] = []
    edge_slices: list[tuple[int, int]] = []
    edge_cursor = 0

    for outer_species_index, (label, ids) in enumerate(zip(labels, draws)):
        if label not in label_to_index:
            raise ValueError(f"outer species {label!r} absent from canonical pool")
        try:
            rows = np.asarray([id_to_row[int(pid)] for pid in ids], dtype=np.int64)
        except KeyError as exc:
            raise ValueError(f"outer photo_id absent from canonical pool: {exc}") from exc
        if np.any(pool.species[rows] != label):
            raise ValueError("outer schedule assigns a photo to the wrong species")
        local_edges, _ = spherical_knn_edges(pool.latitude[rows], pool.longitude[rows], k=int(k))
        global_edges = rows[local_edges]
        lat, lon = spherical_midpoint(
            pool.latitude[global_edges[:, 0]],
            pool.longitude[global_edges[:, 0]],
            pool.latitude[global_edges[:, 1]],
            pool.longitude[global_edges[:, 1]],
        )
        edge_nodes.append(global_edges)
        edge_species_index.append(np.full(len(global_edges), outer_species_index, dtype=np.int32))
        mid_lat.append(lat)
        mid_lon.append(lon)
        edge_slices.append((edge_cursor, edge_cursor + len(global_edges)))
        edge_cursor += len(global_edges)

    nodes = np.vstack(edge_nodes)
    species_index = np.concatenate(edge_species_index)
    latitude = np.concatenate(mid_lat)
    longitude = np.concatenate(mid_lon)
    distance = _distance_matrix_km(latitude, longitude, grid)
    cutoff = kernel_km * cutoff_multiplier
    kernel = np.exp(-0.5 * np.square(distance / kernel_km))
    kernel[distance > cutoff] = 0.0
    row_sum = kernel.sum(axis=1)
    missing = np.flatnonzero(row_sum <= 0)
    if len(missing):
        nearest = np.argmin(distance[missing], axis=1)
        kernel[missing, nearest] = 1.0
        row_sum = kernel.sum(axis=1)
    kernel = kernel / row_sum[:, None]

    species_edge_counts = np.bincount(species_index, minlength=len(labels)).astype(float)
    edge_equal_species_weight = 1.0 / species_edge_counts[species_index]
    kernel *= edge_equal_species_weight[:, None]
    weighted = sparse.csr_matrix(kernel)
    opportunity = np.asarray(weighted.sum(axis=0)).ravel()

    support = np.zeros(grid.n_cells, dtype=np.int32)
    for start, stop in edge_slices:
        contribution = np.asarray(weighted[start:stop].sum(axis=0)).ravel()
        support += (contribution > 0).astype(np.int32)
    evaluable = (support >= minimum_distinct_species) & (opportunity > 0)
    return SparseOuterGeometry(
        grid=grid,
        edge_nodes=nodes,
        edge_species_index=species_index,
        edge_species_slices=tuple(edge_slices),
        weighted_kernel=weighted,
        opportunity=opportunity,
        distinct_species_support=support,
        evaluable=evaluable,
    )


def _raw_jsd_from_source_rows(
    cache: PairwiseJSDCache,
    edge_nodes: np.ndarray,
    source_rows: np.ndarray,
) -> np.ndarray:
    left_target = edge_nodes[:, 0]
    right_target = edge_nodes[:, 1]
    left_source = source_rows[left_target]
    right_source = source_rows[right_target]
    # The null mapping is species conditioned, so the target species block is also
    # the correct source JSD matrix block.
    n = cache.species_size_by_row[left_target].astype(np.int64)
    offset = cache.matrix_offset_by_row[left_target]
    left_local = cache.local_index_by_row[left_source].astype(np.int64)
    right_local = cache.local_index_by_row[right_source].astype(np.int64)
    index = offset + left_local * n + right_local
    return cache.flat_jsd[index]


def _raw_jsd_matrix_from_source_rows(
    cache: PairwiseJSDCache,
    edge_nodes: np.ndarray,
    source_rows_matrix: np.ndarray,
) -> np.ndarray:
    source = np.asarray(source_rows_matrix, dtype=np.int64)
    if source.ndim != 2:
        raise ValueError("source_rows_matrix must be nulls x pool_rows")
    left_target = edge_nodes[:, 0]
    right_target = edge_nodes[:, 1]
    left_source = source[:, left_target]
    right_source = source[:, right_target]
    n = cache.species_size_by_row[left_target].astype(np.int64)[None, :]
    offset = cache.matrix_offset_by_row[left_target][None, :]
    left_local = cache.local_index_by_row[left_source].astype(np.int64)
    right_local = cache.local_index_by_row[right_source].astype(np.int64)
    index = offset + left_local * n + right_local
    return cache.flat_jsd[index]


def _rank_edges(raw: np.ndarray, slices: Sequence[tuple[int, int]]) -> np.ndarray:
    scores = np.asarray(raw, dtype=float)
    if scores.ndim != 1:
        raise ValueError("raw edge scores must be one-dimensional")
    out = np.empty_like(scores)
    for start, stop in slices:
        n = stop - start
        ranks = rankdata(scores[start:stop], method="average")
        out[start:stop] = (ranks - 0.5) / n
    return out


def _rank_edges_matrix(raw: np.ndarray, slices: Sequence[tuple[int, int]]) -> np.ndarray:
    scores = np.asarray(raw, dtype=float)
    if scores.ndim != 2:
        raise ValueError("raw edge score matrix must be two-dimensional")
    out = np.empty_like(scores)
    for start, stop in slices:
        n = stop - start
        ranks = rankdata(scores[:, start:stop], method="average", axis=1)
        out[:, start:stop] = (ranks - 0.5) / n
    return out


def _field_from_rank(geometry: SparseOuterGeometry, rank_scores: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    numerator = np.asarray(geometry.weighted_kernel.T.dot(np.asarray(rank_scores, dtype=float))).ravel()
    field = np.full(geometry.grid.n_cells, np.nan, dtype=float)
    field[geometry.evaluable] = numerator[geometry.evaluable] / geometry.opportunity[geometry.evaluable]
    return field, numerator


def _field_matrix_from_rank(geometry: SparseOuterGeometry, rank_scores: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    scores = np.asarray(rank_scores, dtype=float)
    numerator = np.asarray(geometry.weighted_kernel.T.dot(scores.T).T)
    field = np.full((scores.shape[0], geometry.grid.n_cells), np.nan, dtype=float)
    field[:, geometry.evaluable] = numerator[:, geometry.evaluable] / geometry.opportunity[geometry.evaluable]
    return field, numerator


def null_source_row_matrix(
    pool: CanonicalColourPool,
    null_indices: Sequence[int],
    *,
    master_seed: int = 2026090403,
) -> np.ndarray:
    """Build deterministic species-conditioned source-row maps for null IDs."""
    indices = np.asarray(null_indices, dtype=np.int64)
    if indices.ndim != 1 or len(indices) == 0 or np.any(indices < 0):
        raise ValueError("null_indices must be a non-empty vector of nonnegative integers")
    if len(np.unique(indices)) != len(indices):
        raise ValueError("null_indices must be unique within one shard")
    source = np.tile(np.arange(len(pool.photo_ids), dtype=np.int64), (len(indices), 1))
    for null_row, null_index in enumerate(indices):
        for label, (start, stop) in zip(pool.species_labels, pool.species_slices):
            rng = np.random.default_rng(stable_seed(int(master_seed), int(null_index), label))
            rows = np.arange(start, stop, dtype=np.int64)
            source[null_row, start:stop] = rows[rng.permutation(len(rows))]
    return source


def _weighted_concentration(field: np.ndarray, opportunity: np.ndarray) -> float:
    keep = np.isfinite(field) & np.isfinite(opportunity) & (opportunity > 0)
    if not np.any(keep):
        return float("nan")
    x = field[keep]
    w = opportunity[keep]
    mean = float(np.average(x, weights=w))
    return float(np.average(np.square(x - mean), weights=w))


def run_g1_shard(
    measured_pool: pd.DataFrame,
    *,
    null_indices: Sequence[int],
    n_outer: int = 200,
    species_per_outer: int = 250,
    photos_per_species: int = 20,
    minimum_pool_photos_per_species: int = 40,
    k: int = 3,
    n_lon: int = 36,
    n_sinlat: int = 18,
    kernel_km: float = 500.0,
    cutoff_multiplier: float = 3.0,
    minimum_distinct_species: int = 5,
    species_seed: int = 2026090401,
    photo_master_seed: int = 2026090402,
    null_master_seed: int = 2026090403,
) -> G1ShardResult:
    """Run observed G1 plus a deterministic subset of the frozen 999 nulls.

    This function is suitable for deterministic null sharding. Each shard repeats
    the observed outer program and must yield an identical observed digest; a final
    reducer can then require null indices 0..998 exactly once before computing p.
    """
    pool = canonical_colour_pool(measured_pool)
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
    cache = build_pairwise_jsd_cache(pool)
    null_indices_arr = np.asarray(null_indices, dtype=np.int64)
    source_null = null_source_row_matrix(pool, null_indices_arr, master_seed=int(null_master_seed))
    identity_source = np.arange(len(pool.photo_ids), dtype=np.int64)
    grid = equal_area_grid_centers(int(n_lon), int(n_sinlat))

    observed_fields = np.full((int(n_outer), grid.n_cells), np.nan, dtype=float)
    observed_opportunities = np.zeros((int(n_outer), grid.n_cells), dtype=float)
    null_aggregate_numerator = np.zeros((len(null_indices_arr), grid.n_cells), dtype=float)
    aggregate_opportunity = np.zeros(grid.n_cells, dtype=float)

    for outer in range(int(n_outer)):
        geometry = prepare_sparse_outer_geometry(
            pool,
            schedule.outer_photo_ids[outer],
            schedule.outer_species[outer],
            grid=grid,
            k=int(k),
            kernel_km=float(kernel_km),
            cutoff_multiplier=float(cutoff_multiplier),
            minimum_distinct_species=int(minimum_distinct_species),
        )
        raw_observed = _raw_jsd_from_source_rows(cache, geometry.edge_nodes, identity_source)
        rank_observed = _rank_edges(raw_observed, geometry.edge_species_slices)
        field_observed, _ = _field_from_rank(geometry, rank_observed)
        observed_fields[outer] = field_observed
        observed_opportunities[outer, geometry.evaluable] = geometry.opportunity[geometry.evaluable]
        aggregate_opportunity[geometry.evaluable] += geometry.opportunity[geometry.evaluable]

        raw_null = _raw_jsd_matrix_from_source_rows(cache, geometry.edge_nodes, source_null)
        rank_null = _rank_edges_matrix(raw_null, geometry.edge_species_slices)
        _, numerator_null = _field_matrix_from_rank(geometry, rank_null)
        null_aggregate_numerator[:, geometry.evaluable] += numerator_null[:, geometry.evaluable]

    observed = consensus_field(observed_fields, observed_opportunities)
    if not np.allclose(observed.aggregate_opportunity, aggregate_opportunity):
        raise RuntimeError("observed consensus aggregate opportunity mismatch")

    null_consensus = np.full((len(null_indices_arr), grid.n_cells), np.nan, dtype=float)
    keep = aggregate_opportunity > 0
    null_consensus[:, keep] = null_aggregate_numerator[:, keep] / aggregate_opportunity[keep]
    null_concentration = np.asarray(
        [_weighted_concentration(row, aggregate_opportunity) for row in null_consensus],
        dtype=float,
    )
    return G1ShardResult(
        schedule=schedule,
        observed_outer_fields=observed_fields,
        observed_outer_opportunities=observed_opportunities,
        observed_consensus=observed.field,
        observed_aggregate_opportunity=observed.aggregate_opportunity,
        observed_concentration=float(observed.concentration),
        null_indices=null_indices_arr,
        null_consensus_fields=null_consensus,
        null_concentrations=null_concentration,
        null_unit="complete_photo_colour_vector_within_species_reused_across_outer_realizations",
    )


__all__ = [
    "COLOUR_COLUMNS",
    "CanonicalColourPool",
    "G1ShardResult",
    "PairwiseJSDCache",
    "SparseOuterGeometry",
    "build_pairwise_jsd_cache",
    "canonical_colour_pool",
    "null_source_row_matrix",
    "prepare_sparse_outer_geometry",
    "run_g1_shard",
]

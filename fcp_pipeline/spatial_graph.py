"""Label-blind spatial graph helpers for Chapter 1 colour-structure tests.

Graphs are constructed from observation geometry and species identity only. Colour values
are never used to choose neighbours or truncate edges. Complete colour vectors can then
be permuted strictly within species while the graph remains fixed.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.spatial import cKDTree

from .continuous_colour_boundaries import (
    edge_colour_discontinuity,
    species_conditioned_vector_permutation,
    weighted_graph_discontinuity,
)

EARTH_RADIUS_KM = 6371.0088


@dataclass(frozen=True)
class SpeciesGraph:
    """Undirected within-species graph built without colour information."""

    edges: np.ndarray
    edge_species: np.ndarray
    edge_distance_km: np.ndarray

    def __post_init__(self) -> None:
        edges = np.asarray(self.edges)
        edge_species = np.asarray(self.edge_species)
        distance = np.asarray(self.edge_distance_km, dtype=float)
        if edges.ndim != 2 or edges.shape[1] != 2:
            raise ValueError("edges must have shape (n_edges, 2)")
        if edge_species.ndim != 1 or distance.ndim != 1:
            raise ValueError("edge_species and edge_distance_km must be one-dimensional")
        if not (len(edges) == len(edge_species) == len(distance)):
            raise ValueError("graph arrays must have equal edge counts")
        if len(edges) == 0:
            raise ValueError("graph must contain at least one edge")
        if np.any(distance < 0) or not np.isfinite(distance).all():
            raise ValueError("edge distances must be finite and non-negative")


def _validate_lat_lon(latitude: np.ndarray, longitude: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    lat = np.asarray(latitude, dtype=float)
    lon = np.asarray(longitude, dtype=float)
    if lat.ndim != 1 or lon.ndim != 1 or lat.shape != lon.shape:
        raise ValueError("latitude and longitude must be equal-length one-dimensional arrays")
    if len(lat) < 2:
        raise ValueError("at least two observations are required")
    if not np.isfinite(lat).all() or not np.isfinite(lon).all():
        raise ValueError("latitude and longitude must be finite")
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


def _great_circle_edge_distances(xyz: np.ndarray, edges: np.ndarray) -> np.ndarray:
    dot = np.sum(xyz[edges[:, 0]] * xyz[edges[:, 1]], axis=1)
    angle = np.arccos(np.clip(dot, -1.0, 1.0))
    return EARTH_RADIUS_KM * angle


def spherical_knn_edges(
    latitude: np.ndarray,
    longitude: np.ndarray,
    *,
    k: int,
    max_edge_km: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Build a symmetrized spherical k-nearest-neighbour graph.

    Great-circle neighbour order is obtained from Euclidean chord distance on the unit
    sphere, which is monotonic in great-circle distance. The graph depends only on
    coordinates. An optional maximum edge length is applied after kNN construction and
    is therefore also label-blind.
    """

    lat, lon = _validate_lat_lon(latitude, longitude)
    n = len(lat)
    if not isinstance(k, (int, np.integer)) or k < 1 or k >= n:
        raise ValueError("k must be an integer in [1, n_observations-1]")
    if max_edge_km is not None:
        max_edge_km = float(max_edge_km)
        if not np.isfinite(max_edge_km) or max_edge_km <= 0:
            raise ValueError("max_edge_km must be finite and positive when provided")

    xyz = _unit_sphere_xyz(lat, lon)
    tree = cKDTree(xyz)
    # k+1 is enough whether self is included or tied duplicates reorder the zero-distance
    # results: after removing self there are at least k candidates.
    _, neighbours = tree.query(xyz, k=k + 1)
    if neighbours.ndim == 1:
        neighbours = neighbours[:, None]

    undirected: set[tuple[int, int]] = set()
    for i in range(n):
        candidates = [int(j) for j in np.atleast_1d(neighbours[i]) if int(j) != i]
        if len(candidates) < k:
            # This can occur only under an extreme duplicate-coordinate tie ordering.
            # Query all observations rather than silently reducing the requested degree.
            _, all_neighbours = tree.query(xyz[i], k=n)
            candidates = [int(j) for j in np.atleast_1d(all_neighbours) if int(j) != i]
        for j in candidates[:k]:
            a, b = (i, j) if i < j else (j, i)
            undirected.add((a, b))

    edges = np.asarray(sorted(undirected), dtype=int)
    if edges.size == 0:
        raise ValueError("kNN construction produced no edges")
    distance = _great_circle_edge_distances(xyz, edges)
    if max_edge_km is not None:
        keep = distance <= max_edge_km
        edges = edges[keep]
        distance = distance[keep]
        if len(edges) == 0:
            raise ValueError("max_edge_km removed every graph edge")
    return edges, distance


def species_conditioned_knn_graph(
    latitude: np.ndarray,
    longitude: np.ndarray,
    species: np.ndarray,
    *,
    k: int,
    max_edge_km: float | None = None,
) -> SpeciesGraph:
    """Build one colour-blind kNN graph per species and combine global row indices.

    Every species is required to contribute at least one retained edge. This prevents a
    species from silently disappearing from the inferential null after edge truncation.
    """

    lat, lon = _validate_lat_lon(latitude, longitude)
    sp = np.asarray(species)
    if sp.ndim != 1 or len(sp) != len(lat):
        raise ValueError("species must be one-dimensional and match the coordinate rows")

    all_edges: list[np.ndarray] = []
    all_species: list[np.ndarray] = []
    all_distance: list[np.ndarray] = []
    for label in np.unique(sp):
        idx = np.flatnonzero(sp == label)
        if len(idx) <= k:
            raise ValueError(f"species {label!r} has {len(idx)} rows, insufficient for k={k}")
        local_edges, local_distance = spherical_knn_edges(
            lat[idx],
            lon[idx],
            k=k,
            max_edge_km=max_edge_km,
        )
        global_edges = idx[local_edges]
        if len(global_edges) == 0:
            raise ValueError(f"species {label!r} contributed no retained graph edges")
        all_edges.append(global_edges)
        all_species.append(np.full(len(global_edges), label, dtype=sp.dtype))
        all_distance.append(local_distance)

    return SpeciesGraph(
        edges=np.vstack(all_edges),
        edge_species=np.concatenate(all_species),
        edge_distance_km=np.concatenate(all_distance),
    )


def equal_species_graph_discontinuity(
    edge_scores: np.ndarray,
    edge_species: np.ndarray,
    *,
    expected_species: np.ndarray,
    weights: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Return species-specific Q_i values and their equal-species global mean.

    Edge counts are normalized within each species first. The final global statistic is
    the arithmetic mean across species, so a species with many retained edges does not
    receive more global weight than a species with few retained edges.
    """

    scores = np.asarray(edge_scores, dtype=float)
    edge_sp = np.asarray(edge_species)
    expected = np.asarray(expected_species)
    if scores.ndim != 1 or edge_sp.ndim != 1 or len(scores) != len(edge_sp):
        raise ValueError("edge_scores and edge_species must be equal-length vectors")
    if expected.ndim != 1 or len(expected) == 0:
        raise ValueError("expected_species must be a non-empty one-dimensional array")
    w = None if weights is None else np.asarray(weights, dtype=float)
    if w is not None and (w.ndim != 1 or w.shape != scores.shape):
        raise ValueError("weights must match edge_scores")

    expected_unique = np.unique(expected)
    edge_unique = np.unique(edge_sp)
    if not np.array_equal(edge_unique, expected_unique):
        raise ValueError("edge graph species do not exactly match expected_species; no species may disappear")

    q = []
    for label in expected_unique:
        idx = edge_sp == label
        q.append(
            weighted_graph_discontinuity(
                scores[idx],
                weights=None if w is None else w[idx],
            )
        )
    q_arr = np.asarray(q, dtype=float)
    return expected_unique, q_arr, float(q_arr.mean())


def species_conditioned_graph_permutation_null(
    standardized_values: np.ndarray,
    species: np.ndarray,
    graph: SpeciesGraph,
    *,
    n_permutations: int,
    rng: np.random.Generator,
    weights: np.ndarray | None = None,
) -> dict[str, np.ndarray | float]:
    """Compute the Stage-A graph null with complete within-species vector shuffles.

    The graph is kept fixed because it is constructed from species and observation
    geometry only. Colour vectors are permuted as indivisible rows within species for
    every replicate. Species-specific Q_i values are calculated first and the global
    statistic gives every species exactly one vote.
    """

    values = np.asarray(standardized_values, dtype=float)
    sp = np.asarray(species)
    if values.ndim != 2 or sp.ndim != 1 or values.shape[0] != len(sp):
        raise ValueError("standardized_values and species observation counts must match")
    if not isinstance(n_permutations, (int, np.integer)) or n_permutations < 1:
        raise ValueError("n_permutations must be a positive integer")
    if graph.edges.min() < 0 or graph.edges.max() >= len(sp):
        raise ValueError("graph edge index is out of bounds for the supplied observations")
    if np.any(sp[graph.edges[:, 0]] != graph.edge_species) or np.any(
        sp[graph.edges[:, 1]] != graph.edge_species
    ):
        raise ValueError("graph contains an edge that crosses species")

    expected_species = np.unique(sp)
    observed_edge_scores = edge_colour_discontinuity(values, graph.edges)
    species_order, observed_q, observed_global = equal_species_graph_discontinuity(
        observed_edge_scores,
        graph.edge_species,
        expected_species=expected_species,
        weights=weights,
    )

    null_species = np.empty((n_permutations, len(species_order)), dtype=float)
    null_global = np.empty(n_permutations, dtype=float)
    for b in range(n_permutations):
        permuted = species_conditioned_vector_permutation(values, sp, rng=rng)
        edge_scores = edge_colour_discontinuity(permuted, graph.edges)
        order_b, q_b, global_b = equal_species_graph_discontinuity(
            edge_scores,
            graph.edge_species,
            expected_species=expected_species,
            weights=weights,
        )
        if not np.array_equal(order_b, species_order):
            raise RuntimeError("species order changed across permutation replicates")
        null_species[b] = q_b
        null_global[b] = global_b

    return {
        "species": species_order,
        "observed_species_q": observed_q,
        "observed_global_equal_species_mean": observed_global,
        "null_species_q": null_species,
        "null_global_equal_species_mean": null_global,
    }

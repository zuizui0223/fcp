"""Continuous within-species flower-colour boundary helpers.

This module generalizes the discrete graph mismatch statistic without changing the
species-conditioned null.  A complete colour vector moves as one indivisible row under
within-species permutation, locations and geometry remain fixed, and boundary
opportunity is still defined independently of colour.

For a one-dimensional binary colour state encoded as 0/1, the RMS edge distance is
exactly ``1(y_a != y_b)``.  The weighted graph discontinuity statistic therefore reduces
exactly to the original discrete mismatch statistic in that special case.
"""

from __future__ import annotations

import numpy as np


def species_conditioned_vector_permutation(
    values: np.ndarray,
    species: np.ndarray,
    *,
    rng: np.random.Generator,
) -> np.ndarray:
    """Permute complete colour vectors strictly within species.

    The entire feature row is moved together.  Component values from different
    observations are never recombined, and no vector can cross a species boundary.
    Observation row order—and therefore spatial locations—remains fixed.
    """

    values = np.asarray(values)
    species = np.asarray(species)
    if values.ndim != 2:
        raise ValueError("values must be a two-dimensional (n_observations, n_features) array")
    if species.ndim != 1:
        raise ValueError("species must be one-dimensional")
    if values.shape[0] != species.shape[0]:
        raise ValueError("values and species must have equal observation counts")
    if not np.isfinite(values).all():
        raise ValueError("values must be finite before permutation")

    out = values.copy()
    for sp in np.unique(species):
        idx = np.flatnonzero(species == sp)
        order = rng.permutation(idx.size)
        out[idx] = values[idx][order]
    return out


def edge_colour_discontinuity(
    standardized_values: np.ndarray,
    edges: np.ndarray,
) -> np.ndarray:
    """Return RMS colour-vector distance for each geometry-defined edge.

    ``standardized_values`` must already use a species-specific transform frozen from
    calibration data.  The RMS normalization by ``sqrt(p)`` makes the score comparable
    across feature dimensions without changing the binary one-feature special case.

    The edge list is assumed to have been constructed from spatial geometry only.
    This function never creates or edits the graph.
    """

    values = np.asarray(standardized_values, dtype=float)
    edges = np.asarray(edges)
    if values.ndim != 2 or values.shape[1] < 1:
        raise ValueError("standardized_values must be a non-empty two-dimensional array")
    if not np.isfinite(values).all():
        raise ValueError("standardized_values must be finite")
    if edges.ndim != 2 or edges.shape[1] != 2:
        raise ValueError("edges must have shape (n_edges, 2)")
    if not np.issubdtype(edges.dtype, np.integer):
        if np.any(edges != np.floor(edges)):
            raise ValueError("edge indices must be integers")
        edges = edges.astype(int)
    else:
        edges = edges.astype(int, copy=False)
    if edges.shape[0] == 0:
        return np.empty(0, dtype=float)
    if edges.min() < 0 or edges.max() >= values.shape[0]:
        raise ValueError("edge index is out of bounds")
    if np.any(edges[:, 0] == edges[:, 1]):
        raise ValueError("self edges are not allowed")

    delta = values[edges[:, 0]] - values[edges[:, 1]]
    return np.sqrt(np.mean(delta * delta, axis=1))


def weighted_graph_discontinuity(
    edge_scores: np.ndarray,
    *,
    weights: np.ndarray | None = None,
) -> float:
    """Calculate the species-level graph discontinuity statistic Q_i.

    Lower values than the within-species permutation null indicate that nearby
    observations are more colour-similar than expected from random labelling.  With a
    scalar 0/1 state this is exactly the weighted discrete edge-mismatch fraction.
    """

    scores = np.asarray(edge_scores, dtype=float)
    if scores.ndim != 1 or scores.size == 0:
        raise ValueError("edge_scores must be a non-empty one-dimensional array")
    if not np.isfinite(scores).all() or np.any(scores < 0):
        raise ValueError("edge_scores must be finite and non-negative")

    if weights is None:
        return float(scores.mean())
    w = np.asarray(weights, dtype=float)
    if w.ndim != 1 or w.shape != scores.shape:
        raise ValueError("weights must have the same one-dimensional shape as edge_scores")
    if not np.isfinite(w).all() or np.any(w < 0):
        raise ValueError("weights must be finite and non-negative")
    total = float(w.sum())
    if total <= 0:
        raise ValueError("weights must have positive total weight")
    return float(np.dot(w, scores) / total)


def average_rank_intensity(edge_scores: np.ndarray) -> np.ndarray:
    """Map edge discontinuities to within-species [0, 1] average-rank intensity.

    Equal edge scores receive equal average ranks.  This avoids arbitrary label-based
    tie breaking.  Rank intensity is useful for cross-species shared-boundary maps
    because it removes species-specific colour scale while retaining relative spatial
    transition strength.  A single-edge graph receives intensity 0.5.
    """

    scores = np.asarray(edge_scores, dtype=float)
    if scores.ndim != 1 or scores.size == 0:
        raise ValueError("edge_scores must be a non-empty one-dimensional array")
    if not np.isfinite(scores).all():
        raise ValueError("edge_scores must be finite")
    n = scores.size
    if n == 1:
        return np.array([0.5], dtype=float)

    order = np.argsort(scores, kind="mergesort")
    sorted_scores = scores[order]
    ranks = np.empty(n, dtype=float)
    start = 0
    while start < n:
        end = start + 1
        while end < n and sorted_scores[end] == sorted_scores[start]:
            end += 1
        # Zero-based average rank, then map 0..n-1 to 0..1.
        average_rank = 0.5 * (start + end - 1)
        ranks[order[start:end]] = average_rank / (n - 1)
        start = end
    return ranks


def shared_boundary_intensity(
    intensity: np.ndarray,
    detectable: np.ndarray,
    *,
    min_detectable_species: int = 1,
) -> tuple[np.ndarray, np.ndarray]:
    """Average continuous boundary intensity over label-blind species opportunities.

    ``intensity`` has shape ``(n_species, n_cells)`` and is expected to be in [0, 1].
    ``detectable`` has the same shape and must be derived from sampling geometry/support
    without colour.  Cells below the opportunity threshold return NaN rather than zero.
    """

    intensity = np.asarray(intensity, dtype=float)
    detectable = np.asarray(detectable, dtype=bool)
    if intensity.ndim != 2 or detectable.ndim != 2 or intensity.shape != detectable.shape:
        raise ValueError("intensity and detectable must be identically shaped two-dimensional arrays")
    if min_detectable_species < 1:
        raise ValueError("min_detectable_species must be >= 1")
    if np.any(~np.isfinite(intensity[detectable])):
        raise ValueError("intensity must be finite wherever a species is detectable")
    if np.any((intensity[detectable] < 0) | (intensity[detectable] > 1)):
        raise ValueError("detectable intensity values must lie in [0, 1]")

    A = detectable.sum(axis=0).astype(int)
    numerator = np.where(detectable, intensity, 0.0).sum(axis=0)
    shared = np.full(A.shape, np.nan, dtype=float)
    valid = A >= min_detectable_species
    shared[valid] = numerator[valid] / A[valid]
    return shared, A


def opportunity_weighted_concentration(
    shared_intensity: np.ndarray,
    opportunity: np.ndarray,
    *,
    min_opportunity: int = 1,
) -> float:
    """Global concentration statistic for co-located continuous boundaries.

    The statistic is the opportunity-weighted variance of the shared intensity across
    evaluable cells.  Large values mean transition intensity is geographically
    concentrated rather than spatially even.  Geometry-induced heterogeneity is handled
    by comparing the observed value with the complete within-species permutation
    pipeline under the same opportunity surface.
    """

    shared = np.asarray(shared_intensity, dtype=float)
    A = np.asarray(opportunity)
    if shared.ndim != 1 or A.ndim != 1 or shared.shape != A.shape:
        raise ValueError("shared_intensity and opportunity must be equal-length vectors")
    if min_opportunity < 1:
        raise ValueError("min_opportunity must be >= 1")
    valid = (A >= min_opportunity) & np.isfinite(shared)
    if not np.any(valid):
        raise ValueError("no evaluable cells for global concentration")
    w = A[valid].astype(float)
    x = shared[valid]
    mean = float(np.dot(w, x) / w.sum())
    return float(np.dot(w, (x - mean) ** 2) / w.sum())

"""Pre-image shared-transition scan for the terminal FCP image-first atlas.

The scan operates on within-species ranks of transition-cell scores.  Geometry and
opportunity are frozen before colour.  Each species contributes at most one high-rank
indicator to a cell, so globally common taxa cannot dominate through edge count.

The real-colour null remains the complete within-species colour-vector permutation
specified by the atlas contract.  This module also exposes a conditional high-rank
placement null used *before pixels* to qualify the cross-species scan layer on the exact
terminal opportunity geometry.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np


EARTH_RADIUS_KM = 6371.0088


@dataclass(frozen=True)
class CoexceedanceReference:
    detectable: np.ndarray
    valid_cell_ids: np.ndarray
    valid_mask: np.ndarray
    high_counts: np.ndarray
    high_probabilities: np.ndarray
    expected_count: np.ndarray
    standard_deviation: np.ndarray

    def __post_init__(self) -> None:
        detectable = np.asarray(self.detectable, dtype=bool)
        if detectable.ndim != 2 or detectable.shape[0] < 1 or detectable.shape[1] < 1:
            raise ValueError("detectable must be a non-empty species x cell matrix")
        if np.asarray(self.valid_cell_ids).ndim != 1:
            raise ValueError("valid_cell_ids must be one-dimensional")
        if np.asarray(self.valid_mask).shape != (detectable.shape[1],):
            raise ValueError("valid_mask does not match the cell dimension")
        if np.asarray(self.high_counts).shape != (detectable.shape[0],):
            raise ValueError("high_counts does not match the species dimension")
        if np.asarray(self.high_probabilities).shape != (detectable.shape[0],):
            raise ValueError("high_probabilities does not match the species dimension")
        n_valid = int(np.count_nonzero(self.valid_mask))
        if np.asarray(self.expected_count).shape != (n_valid,):
            raise ValueError("expected_count does not match valid cells")
        if np.asarray(self.standard_deviation).shape != (n_valid,):
            raise ValueError("standard_deviation does not match valid cells")
        if n_valid == 0:
            raise ValueError("no cells satisfy the frozen cross-species opportunity gate")
        if np.any(np.asarray(self.standard_deviation) <= 0):
            raise ValueError("valid cells must have positive null standard deviation")


def build_detectability_matrix(
    species_detectable_cells: Sequence[Sequence[int]], *, n_cells: int
) -> np.ndarray:
    """Construct a frozen species x cell opportunity matrix from cell-id ledgers."""
    n_cells = int(n_cells)
    if n_cells < 1 or len(species_detectable_cells) < 1:
        raise ValueError("species_detectable_cells and n_cells must be non-empty")
    matrix = np.zeros((len(species_detectable_cells), n_cells), dtype=bool)
    for i, cells in enumerate(species_detectable_cells):
        ids = np.asarray(list(cells), dtype=int)
        if ids.ndim != 1 or ids.size == 0:
            raise ValueError(f"species index {i} has no detectable cells")
        if np.any(ids < 0) or np.any(ids >= n_cells):
            raise ValueError(f"species index {i} contains an out-of-grid cell id")
        if len(np.unique(ids)) != len(ids):
            raise ValueError(f"species index {i} repeats a detectable cell id")
        matrix[i, ids] = True
    return matrix


def build_coexceedance_reference(
    detectable: np.ndarray,
    *,
    high_transition_quantile: float = 0.9,
    min_detectable_species: int = 4,
) -> CoexceedanceReference:
    """Freeze the conditional expectation for top-tail transition co-exceedance.

    For species ``i`` with ``d_i`` detectable cells, exactly
    ``ceil((1-q) d_i)`` cells are labelled high transition.  Under the conditional
    placement null those labels are uniformly distributed among the species' frozen
    detectable cells.  The cellwise expectation and variance therefore account for
    unequal opportunity without weighting by raw observation or edge count.
    """
    D = np.asarray(detectable, dtype=bool)
    if D.ndim != 2 or D.shape[0] < 1 or D.shape[1] < 1:
        raise ValueError("detectable must be a non-empty species x cell matrix")
    q = float(high_transition_quantile)
    if not (0.5 < q < 1.0):
        raise ValueError("high_transition_quantile must lie in (0.5, 1)")
    min_detectable_species = int(min_detectable_species)
    if min_detectable_species < 2:
        raise ValueError("min_detectable_species must be >= 2")

    detectable_counts = D.sum(axis=1).astype(int)
    if np.any(detectable_counts < 1):
        raise ValueError("every species must retain at least one detectable cell")
    tail_fraction = 1.0 - q
    high_counts = np.maximum(1, np.ceil(tail_fraction * detectable_counts).astype(int))
    high_probabilities = high_counts / detectable_counts

    opportunity = D.sum(axis=0).astype(int)
    valid_mask = opportunity >= min_detectable_species
    valid_ids = np.flatnonzero(valid_mask)
    if len(valid_ids) == 0:
        raise ValueError("no cells satisfy min_detectable_species")
    D_valid = D[:, valid_mask]
    expected = (D_valid * high_probabilities[:, None]).sum(axis=0)
    variance = (
        D_valid
        * (high_probabilities * (1.0 - high_probabilities))[:, None]
    ).sum(axis=0)
    if np.any(variance <= 0):
        raise ValueError("frozen valid cells have zero conditional variance")

    return CoexceedanceReference(
        detectable=D,
        valid_cell_ids=valid_ids.astype(int),
        valid_mask=valid_mask,
        high_counts=high_counts,
        high_probabilities=high_probabilities,
        expected_count=expected.astype(float),
        standard_deviation=np.sqrt(variance).astype(float),
    )


def high_transition_mask_from_scores(
    transition_scores: np.ndarray,
    reference: CoexceedanceReference,
) -> np.ndarray:
    """Convert species-cell transition scores to the prospectively frozen top tail."""
    scores = np.asarray(transition_scores, dtype=float)
    D = np.asarray(reference.detectable, dtype=bool)
    if scores.shape != D.shape:
        raise ValueError("transition_scores must match the frozen detectability matrix")
    if np.any(~np.isfinite(scores[D])):
        raise ValueError("transition scores must be finite wherever a species is detectable")

    high = np.zeros_like(D, dtype=bool)
    for i in range(D.shape[0]):
        cells = np.flatnonzero(D[i])
        count = int(reference.high_counts[i])
        values = scores[i, cells]
        # Stable two-key ordering makes ties deterministic without using species labels.
        order = np.lexsort((cells, values))
        chosen = cells[order[-count:]]
        high[i, chosen] = True
    return high


def coexceedance_z(
    high_transition: np.ndarray,
    reference: CoexceedanceReference,
) -> np.ndarray:
    high = np.asarray(high_transition, dtype=bool)
    D = np.asarray(reference.detectable, dtype=bool)
    if high.shape != D.shape:
        raise ValueError("high_transition must match the frozen detectability matrix")
    if np.any(high & ~D):
        raise ValueError("high-transition flags cannot occur outside detectable cells")
    if not np.array_equal(high.sum(axis=1), reference.high_counts):
        raise ValueError("each species must contribute exactly its frozen high-tail count")
    K = high[:, reference.valid_mask].sum(axis=0).astype(float)
    return (K - reference.expected_count) / reference.standard_deviation


def coexceedance_scan_statistic(
    high_transition: np.ndarray,
    reference: CoexceedanceReference,
) -> float:
    """Maximum standardized excess co-exceedance across prospectively valid cells."""
    return float(np.max(coexceedance_z(high_transition, reference)))


def conditional_rank_scan_null(
    reference: CoexceedanceReference,
    *,
    n_permutations: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Conditional high-rank placement null for pre-image scan qualification.

    The high-tail count for every species is held fixed and placed uniformly among that
    species' frozen detectable cells.  This null is intentionally colour-free and is
    used only to qualify the cross-species scan layer before terminal pixels are opened.
    """
    B = int(n_permutations)
    if B < 1:
        raise ValueError("n_permutations must be positive")
    D = np.asarray(reference.detectable, dtype=bool)
    valid_ids = reference.valid_cell_ids
    full_to_valid = np.full(D.shape[1], -1, dtype=int)
    full_to_valid[valid_ids] = np.arange(len(valid_ids), dtype=int)
    counts = np.zeros((B, len(valid_ids)), dtype=np.int16)

    for i in range(D.shape[0]):
        cells = np.flatnonzero(D[i])
        m = int(reference.high_counts[i])
        random_scores = rng.random((B, len(cells)))
        chosen_local = np.argpartition(random_scores, m - 1, axis=1)[:, :m]
        chosen_cells = cells[chosen_local]
        mapped = full_to_valid[chosen_cells]
        keep = mapped >= 0
        row_ids = np.broadcast_to(np.arange(B)[:, None], mapped.shape)[keep]
        col_ids = mapped[keep]
        np.add.at(counts, (row_ids, col_ids), 1)

    z = (
        counts.astype(float) - reference.expected_count[None, :]
    ) / reference.standard_deviation[None, :]
    return np.max(z, axis=1)


def monte_carlo_p(observed: float, null_distribution: np.ndarray) -> float:
    null = np.asarray(null_distribution, dtype=float)
    observed = float(observed)
    if null.ndim != 1 or null.size < 1 or not np.isfinite(null).all():
        raise ValueError("null_distribution must be a finite non-empty vector")
    if not np.isfinite(observed):
        raise ValueError("observed statistic must be finite")
    return float((1 + np.count_nonzero(null >= observed)) / (len(null) + 1))


def equal_area_cell_xyz(*, n_lon: int, n_sinlat: int) -> np.ndarray:
    """Unit-sphere coordinates for equal-area longitude-sin(latitude) cell centres."""
    n_lon = int(n_lon)
    n_sinlat = int(n_sinlat)
    if n_lon < 2 or n_sinlat < 2:
        raise ValueError("grid dimensions must be >= 2")
    cell = np.arange(n_lon * n_sinlat, dtype=int)
    row = cell // n_lon
    col = cell % n_lon
    sin_lat = -1.0 + (row + 0.5) * (2.0 / n_sinlat)
    lat = np.arcsin(np.clip(sin_lat, -1.0, 1.0))
    lon = -np.pi + (col + 0.5) * (2.0 * np.pi / n_lon)
    cos_lat = np.cos(lat)
    return np.column_stack((cos_lat * np.cos(lon), cos_lat * np.sin(lon), np.sin(lat)))


def _random_unit_normal(rng: np.random.Generator) -> np.ndarray:
    vector = rng.normal(size=3)
    norm = float(np.linalg.norm(vector))
    while norm <= 1e-12:
        vector = rng.normal(size=3)
        norm = float(np.linalg.norm(vector))
    return vector / norm


def _boundary_proximity(
    xyz: np.ndarray, normal: np.ndarray, *, sigma_radians: float
) -> np.ndarray:
    angular_distance = np.arcsin(
        np.clip(np.abs(np.asarray(xyz, dtype=float) @ np.asarray(normal, dtype=float)), 0.0, 1.0)
    )
    return np.exp(-0.5 * (angular_distance / float(sigma_radians)) ** 2)


def simulate_transition_scores(
    reference: CoexceedanceReference,
    xyz: np.ndarray,
    *,
    scenario: str,
    effect_size: float,
    boundary_sigma_radians: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """Simulate colour-free cell transition scores on exact frozen opportunity geometry."""
    if scenario not in {
        "null_stationary",
        "within_species_heterogeneous_nonshared_boundaries",
        "shared_geographic_boundary",
    }:
        raise ValueError(f"unknown signal-recovery scenario: {scenario}")
    effect = float(effect_size)
    sigma = float(boundary_sigma_radians)
    if effect < 0 or not np.isfinite(effect):
        raise ValueError("effect_size must be finite and non-negative")
    if sigma <= 0 or not np.isfinite(sigma):
        raise ValueError("boundary_sigma_radians must be finite and positive")
    xyz = np.asarray(xyz, dtype=float)
    D = np.asarray(reference.detectable, dtype=bool)
    if xyz.shape != (D.shape[1], 3):
        raise ValueError("xyz must contain one unit-sphere centre per grid cell")

    scores = np.full(D.shape, np.nan, dtype=float)
    common = _random_unit_normal(rng) if scenario == "shared_geographic_boundary" else None
    for i in range(D.shape[0]):
        cells = np.flatnonzero(D[i])
        noise = rng.normal(size=len(cells))
        if scenario == "null_stationary":
            values = noise
        else:
            normal = common if common is not None else _random_unit_normal(rng)
            values = noise + effect * _boundary_proximity(
                xyz[cells], normal, sigma_radians=sigma
            )
        scores[i, cells] = values
    return scores


def signal_recovery_rates(
    reference: CoexceedanceReference,
    xyz: np.ndarray,
    *,
    n_repetitions: int,
    n_permutations: int,
    alpha: float,
    effect_sizes: Sequence[float],
    boundary_sigma_radians: float,
    seed: int,
) -> dict[str, object]:
    """Run the frozen scan-layer signal-recovery qualification.

    The randomization null is generated once because, conditional on the frozen
    opportunity matrix and each species' fixed top-tail count, it is independent of the
    synthetic transition-score realization.
    """
    repetitions = int(n_repetitions)
    if repetitions < 1:
        raise ValueError("n_repetitions must be positive")
    alpha = float(alpha)
    if not (0.0 < alpha < 1.0):
        raise ValueError("alpha must lie in (0, 1)")
    effects = tuple(float(value) for value in effect_sizes)
    if effects != (0.0, 0.5, 1.0, 2.0):
        raise ValueError("effect_sizes must remain frozen at 0, 0.5, 1, 2")

    rng = np.random.default_rng(int(seed))
    null_scan = conditional_rank_scan_null(
        reference, n_permutations=int(n_permutations), rng=rng
    )

    def evaluate(scenario: str, effect: float) -> dict[str, object]:
        statistics = np.empty(repetitions, dtype=float)
        p_values = np.empty(repetitions, dtype=float)
        for r in range(repetitions):
            scores = simulate_transition_scores(
                reference,
                xyz,
                scenario=scenario,
                effect_size=effect,
                boundary_sigma_radians=boundary_sigma_radians,
                rng=rng,
            )
            high = high_transition_mask_from_scores(scores, reference)
            statistic = coexceedance_scan_statistic(high, reference)
            statistics[r] = statistic
            p_values[r] = monte_carlo_p(statistic, null_scan)
        return {
            "scenario": scenario,
            "effect_size": effect,
            "rejections": int(np.count_nonzero(p_values <= alpha)),
            "rate": float(np.mean(p_values <= alpha)),
            "median_p": float(np.median(p_values)),
            "median_statistic": float(np.median(statistics)),
            "p_values": p_values,
            "statistics": statistics,
        }

    null_result = evaluate("null_stationary", 0.0)
    heterogeneous = evaluate(
        "within_species_heterogeneous_nonshared_boundaries", 2.0
    )
    shared = [evaluate("shared_geographic_boundary", effect) for effect in effects[1:]]
    return {
        "null_distribution": null_scan,
        "null_result": null_result,
        "heterogeneous_result": heterogeneous,
        "shared_results": shared,
    }

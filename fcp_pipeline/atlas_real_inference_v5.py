"""Exact real-colour statistics for the terminal FCP atlas v5 cascade.

This module begins only after location-blind measurement and the sealed coordinate join.
It never selects species, observations, graph edges, cells, scales or cohorts from colour.
The same standardized three-component Lab rows are reused across reached branches and
complete vectors move together strictly within species under every biological null.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Sequence

import numpy as np

from .atlas_colour_inference import robust_standardize_lab
from .atlas_inference_cascade import classify_shared_transition, classify_spatial_organization
from .atlas_shared_transition_v5 import (
    build_coexceedance_reference,
    coexceedance_scan_statistic,
    high_transition_mask_from_scores,
    monte_carlo_p,
)
from .continuous_colour_boundaries import average_rank_intensity, edge_colour_discontinuity
from .shared_transition_surface import (
    EdgeCellGeometry,
    EqualAreaGrid,
    build_edge_cell_geometry,
    cell_mean_intensity,
)
from .spatial_graph import spherical_knn_edges


PROTOCOL = "jbi-atlas-real-colour-inference-amendment-v5-v1"
COHORTS = tuple(f"C{i:02d}" for i in range(1, 9))


@dataclass(frozen=True)
class FrozenSpeciesColourState:
    species_id: str
    cohort_id: str
    standardized_lab: np.ndarray
    geometry: EdgeCellGeometry

    def __post_init__(self) -> None:
        values = np.asarray(self.standardized_lab, dtype=float)
        if not self.species_id or self.cohort_id not in COHORTS:
            raise ValueError("species state requires a frozen species and C01-C08 cohort")
        if values.ndim != 2 or values.shape[1] != 3 or values.shape[0] < 2:
            raise ValueError("standardized Lab must have shape (n>=2, 3)")
        if not np.isfinite(values).all():
            raise ValueError("standardized Lab must be finite")
        if int(np.max(self.geometry.retained_edges)) >= values.shape[0]:
            raise ValueError("geometry edge index exceeds the species colour denominator")


def validate_real_inference_amendment(contract: Mapping[str, Any]) -> None:
    if contract.get("protocol") != PROTOCOL:
        raise ValueError("unexpected real-colour inference amendment")
    if (
        contract.get("status") != "prospectively_frozen_before_any_terminal_scaleout_candidate_pixel"
        or contract.get("pixel_status_at_freeze") != "not_revealed"
    ):
        raise ValueError("real-colour inference was not frozen before terminal pixels")
    denominator = contract.get("terminal_denominator", {})
    if (
        denominator.get("species") != 200
        or denominator.get("cohorts") != 8
        or denominator.get("species_per_cohort") != 25
        or denominator.get("frozen_observations") != 60000
        or denominator.get("minimum_evaluable_species_per_cohort") != 20
        or denominator.get("minimum_evaluable_species_total") != 160
        or denominator.get("all_eight_cohorts_required") is not True
        or denominator.get("replacement_or_resampling_after_colour") is not False
    ):
        raise ValueError("terminal real-colour denominator changed")
    colour = contract.get("frozen_colour_field", {})
    if (
        colour.get("primary_components") != ["flower_L_mean", "flower_a_mean", "flower_b_mean"]
        or colour.get("standardize_once_before_any_branch") is not True
        or colour.get("same_standardized_vector_table_reused_by_all_reached_branches") is not True
        or colour.get("complete_three_component_vector_moves_together_under_null") is not True
        or colour.get("permutation_strictly_within_species") is not True
        or colour.get("coordinates_graphs_opportunity_and_cohort_membership_fixed_under_null") is not True
        or colour.get("downstream_remeasurement_or_restandardization_forbidden") is not True
    ):
        raise ValueError("real-colour field firewall changed")
    geometry = contract.get("geometry", {})
    if (
        geometry.get("knn_k") != 5
        or geometry.get("primary_scale_km") != 100
        or geometry.get("sensitivity_scales_km") != [250, 500]
        or geometry.get("equal_area_grids", {}).get("100") != [320, 160]
        or geometry.get("minimum_edges_per_species_cell") != 2
        or geometry.get("minimum_retained_edges_per_species_configuration") != 100
        or geometry.get("minimum_detectable_cells_per_species_configuration") != 10
        or geometry.get("geometry_is_coordinate_only") is not True
        or geometry.get("colour_must_not_select_edges_cells_or_scales") is not True
    ):
        raise ValueError("real-colour primary geometry changed")
    spatial = contract.get("species_conditioned_spatial_organization", {})
    if (
        spatial.get("randomizations") != 9999
        or spatial.get("seed") != 20260901
        or spatial.get("alpha") != 0.05
        or spatial.get("minimum_directionally_positive_cohorts") != 6
        or spatial.get("alternative") != "greater"
    ):
        raise ValueError("spatial-organization inference rule changed")
    shared = contract.get("shared_transition", {})
    if (
        shared.get("high_transition_quantile") != 0.9
        or shared.get("minimum_detectable_species_per_tested_cell") != 4
        or shared.get("randomizations") != 9999
        or shared.get("seed") != 20260902
        or shared.get("alpha") != 0.05
        or shared.get("minimum_directionally_positive_cohorts") != 6
        or shared.get("preimage_conditional_rank_placement_null_must_not_be_used_for_biological_p_values") is not True
    ):
        raise ValueError("shared-transition real-colour inference rule changed")
    environment = contract.get("environmental_concordance", {})
    if (
        environment.get("randomizations") != 9999
        or environment.get("seed") != 20260831
        or environment.get("primary_scale_km") != 100
        or environment.get("primary_season") != "all_dates"
        or environment.get("primary_families_in_fixed_order")
        != ["macroclimate_boundary", "land_cover_boundary", "ecoregion_boundary"]
        or environment.get("terrain_boundary") != "not_evaluable_precolour_and_not_substituted"
    ):
        raise ValueError("environmental v5 inference rule changed")
    if contract.get("cascade", {}).get("not_evaluable_never_advances_confirmatory") is not True:
        raise ValueError("not_evaluable firewall changed")


def prepare_species_colour_state(
    *,
    species_id: str,
    cohort_id: str,
    latitude: Sequence[float] | np.ndarray,
    longitude: Sequence[float] | np.ndarray,
    lab_values: np.ndarray,
    grid: EqualAreaGrid,
    knn_k: int = 5,
    max_edge_km: float = 100.0,
    minimum_edges_per_cell: int = 2,
    minimum_retained_edges: int = 100,
    minimum_detectable_cells: int = 10,
) -> FrozenSpeciesColourState:
    """Freeze one species transform and coordinate-only primary geometry once."""
    latitude = np.asarray(latitude, dtype=float)
    longitude = np.asarray(longitude, dtype=float)
    values, _ = robust_standardize_lab(np.asarray(lab_values, dtype=float))
    if values.shape[0] != latitude.size or longitude.shape != latitude.shape:
        raise ValueError("coordinate and Lab denominators differ")
    edges, distances = spherical_knn_edges(latitude, longitude, k=knn_k)
    geometry = build_edge_cell_geometry(
        latitude,
        longitude,
        edges,
        distances,
        grid=grid,
        max_edge_km=max_edge_km,
        min_edges_per_cell=minimum_edges_per_cell,
    )
    if len(geometry.retained_edges) < minimum_retained_edges:
        raise ValueError("species has too few retained primary edges")
    if int(np.count_nonzero(geometry.detectable)) < minimum_detectable_cells:
        raise ValueError("species has too few detectable primary cells")
    return FrozenSpeciesColourState(
        species_id=str(species_id),
        cohort_id=str(cohort_id),
        standardized_lab=values,
        geometry=geometry,
    )


def _validate_states(states: Sequence[FrozenSpeciesColourState], *, require_terminal: bool) -> tuple[FrozenSpeciesColourState, ...]:
    values = tuple(states)
    if not values:
        raise ValueError("at least one frozen species state is required")
    species = [item.species_id for item in values]
    if len(set(species)) != len(species):
        raise ValueError("species IDs must be unique across frozen states")
    counts = {cohort: 0 for cohort in COHORTS}
    for item in values:
        counts[item.cohort_id] += 1
    if require_terminal:
        if len(values) < 160 or any(counts[cohort] < 20 for cohort in COHORTS):
            raise ValueError("terminal inference requires >=160 species and >=20 in every cohort")
    elif any(counts[cohort] == 0 for cohort in COHORTS):
        raise ValueError("all eight cohorts must be represented")
    n_cells = {len(item.geometry.detectable) for item in values}
    if len(n_cells) != 1:
        raise ValueError("species states use different grids")
    return values


def _spatial_species_stat(lab: np.ndarray, geometry: EdgeCellGeometry) -> float:
    return -float(np.mean(edge_colour_discontinuity(lab, geometry.retained_edges)))


def _cohort_means(species_stats: np.ndarray, states: Sequence[FrozenSpeciesColourState]) -> np.ndarray:
    return np.asarray(
        [
            np.mean([species_stats[i] for i, item in enumerate(states) if item.cohort_id == cohort])
            for cohort in COHORTS
        ],
        dtype=float,
    )


def run_spatial_organization_test(
    states: Sequence[FrozenSpeciesColourState],
    *,
    inference_v5: Mapping[str, Any],
    randomizations: int = 9999,
    seed: int = 20260901,
    require_terminal: bool = True,
) -> dict[str, Any]:
    """Run the frozen equal-species primary spatial-organization permutation test."""
    states = _validate_states(states, require_terminal=require_terminal)
    B = int(randomizations)
    if B < 1:
        raise ValueError("randomizations must be positive")
    observed_species = np.asarray(
        [_spatial_species_stat(item.standardized_lab, item.geometry) for item in states],
        dtype=float,
    )
    observed_pooled = float(np.mean(observed_species))
    observed_cohorts = _cohort_means(observed_species, states)
    null_pooled = np.empty(B, dtype=float)
    null_cohorts = np.empty((B, 8), dtype=float)
    rng = np.random.default_rng(int(seed))
    for b in range(B):
        permuted_stats = np.empty(len(states), dtype=float)
        for i, item in enumerate(states):
            order = rng.permutation(item.standardized_lab.shape[0])
            permuted_stats[i] = _spatial_species_stat(item.standardized_lab[order], item.geometry)
        null_pooled[b] = float(np.mean(permuted_stats))
        null_cohorts[b] = _cohort_means(permuted_stats, states)
    p_value = monte_carlo_p(observed_pooled, null_pooled)
    null_medians = np.median(null_cohorts, axis=0)
    directions = observed_cohorts - null_medians
    outcome = classify_spatial_organization(
        pooled_p_value=p_value,
        cohort_directions=directions.tolist(),
        contract=inference_v5,
    )
    return {
        "branch": "species_conditioned_spatial_organization",
        "outcome": outcome,
        "randomizations": B,
        "seed": int(seed),
        "evaluable_species": len(states),
        "observed_pooled_statistic": observed_pooled,
        "pooled_p_value": p_value,
        "observed_cohort_statistics": dict(zip(COHORTS, observed_cohorts.tolist(), strict=True)),
        "null_cohort_medians": dict(zip(COHORTS, null_medians.tolist(), strict=True)),
        "cohort_directions": dict(zip(COHORTS, directions.tolist(), strict=True)),
        "positive_cohorts": int(np.count_nonzero(directions > 0)),
        "null_pooled_median": float(np.median(null_pooled)),
        "null_pooled_q95": float(np.quantile(null_pooled, 0.95)),
    }


def _transition_score_matrix(states: Sequence[FrozenSpeciesColourState], labs: Sequence[np.ndarray]) -> np.ndarray:
    n_cells = len(states[0].geometry.detectable)
    scores = np.full((len(states), n_cells), np.nan, dtype=float)
    for i, (item, lab) in enumerate(zip(states, labs, strict=True)):
        edge_scores = edge_colour_discontinuity(lab, item.geometry.retained_edges)
        edge_ranks = average_rank_intensity(edge_scores)
        scores[i] = cell_mean_intensity(edge_ranks, item.geometry)
    return scores


def _detectability(states: Sequence[FrozenSpeciesColourState]) -> np.ndarray:
    return np.vstack([item.geometry.detectable for item in states]).astype(bool)


def run_shared_transition_test(
    states: Sequence[FrozenSpeciesColourState],
    *,
    inference_v5: Mapping[str, Any],
    qualification_passed: bool,
    randomizations: int = 9999,
    seed: int = 20260902,
    high_transition_quantile: float = 0.9,
    min_detectable_species: int = 4,
    require_terminal: bool = True,
) -> dict[str, Any]:
    """Run the qualified co-exceedance scan using the real Lab-vector permutation null."""
    states = _validate_states(states, require_terminal=require_terminal)
    if qualification_passed is not True:
        return {
            "branch": "shared_transition",
            "outcome": "not_evaluable",
            "reason": "preimage_qualification_not_passed",
        }
    B = int(randomizations)
    if B < 1:
        raise ValueError("randomizations must be positive")
    D = _detectability(states)
    pooled_reference = build_coexceedance_reference(
        D,
        high_transition_quantile=high_transition_quantile,
        min_detectable_species=min_detectable_species,
    )
    cohort_indices = {
        cohort: np.asarray([i for i, item in enumerate(states) if item.cohort_id == cohort], dtype=int)
        for cohort in COHORTS
    }
    cohort_references = {
        cohort: build_coexceedance_reference(
            D[index],
            high_transition_quantile=high_transition_quantile,
            min_detectable_species=min_detectable_species,
        )
        for cohort, index in cohort_indices.items()
    }

    observed_scores = _transition_score_matrix(states, [item.standardized_lab for item in states])
    observed_high = high_transition_mask_from_scores(observed_scores, pooled_reference)
    observed_pooled = coexceedance_scan_statistic(observed_high, pooled_reference)
    observed_cohorts = np.asarray(
        [
            coexceedance_scan_statistic(
                observed_high[cohort_indices[cohort]], cohort_references[cohort]
            )
            for cohort in COHORTS
        ],
        dtype=float,
    )

    null_pooled = np.empty(B, dtype=float)
    null_cohorts = np.empty((B, 8), dtype=float)
    rng = np.random.default_rng(int(seed))
    for b in range(B):
        labs = []
        for item in states:
            order = rng.permutation(item.standardized_lab.shape[0])
            labs.append(item.standardized_lab[order])
        scores = _transition_score_matrix(states, labs)
        high = high_transition_mask_from_scores(scores, pooled_reference)
        null_pooled[b] = coexceedance_scan_statistic(high, pooled_reference)
        for j, cohort in enumerate(COHORTS):
            idx = cohort_indices[cohort]
            null_cohorts[b, j] = coexceedance_scan_statistic(
                high[idx], cohort_references[cohort]
            )

    p_value = monte_carlo_p(observed_pooled, null_pooled)
    null_medians = np.median(null_cohorts, axis=0)
    directions = observed_cohorts - null_medians
    outcome = classify_shared_transition(
        pooled_p_value=p_value,
        cohort_directions=directions.tolist(),
        qualification_passed=True,
        contract=inference_v5,
    )
    return {
        "branch": "shared_transition",
        "outcome": outcome,
        "qualification_passed": True,
        "real_colour_null": "complete_within_species_lab_vector_permutation",
        "randomizations": B,
        "seed": int(seed),
        "evaluable_species": len(states),
        "observed_pooled_statistic": float(observed_pooled),
        "pooled_p_value": p_value,
        "observed_cohort_statistics": dict(zip(COHORTS, observed_cohorts.tolist(), strict=True)),
        "null_cohort_medians": dict(zip(COHORTS, null_medians.tolist(), strict=True)),
        "cohort_directions": dict(zip(COHORTS, directions.tolist(), strict=True)),
        "positive_cohorts": int(np.count_nonzero(directions > 0)),
        "null_pooled_median": float(np.median(null_pooled)),
        "null_pooled_q95": float(np.quantile(null_pooled, 0.95)),
        "valid_pooled_cells": int(len(pooled_reference.valid_cell_ids)),
    }


__all__ = [
    "COHORTS",
    "FrozenSpeciesColourState",
    "prepare_species_colour_state",
    "run_shared_transition_test",
    "run_spatial_organization_test",
    "validate_real_inference_amendment",
]

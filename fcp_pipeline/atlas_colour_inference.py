"""Frozen species-conditioned colour surfaces and joint environmental inference."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Sequence

import numpy as np

from .atlas_overlay_null import (
    MoranSignBasis,
    build_moran_sign_basis,
    equal_area_rook_adjacency,
    geographic_design,
    residual_coefficients,
)
from .continuous_colour_boundaries import (
    average_rank_intensity,
    edge_colour_discontinuity,
)
from .shared_transition_surface import (
    EdgeCellGeometry,
    EqualAreaGrid,
    build_edge_cell_geometry,
    cell_mean_intensity,
)
from .spatial_graph import spherical_knn_edges


COLOUR_INFERENCE_PROTOCOL = "jbi-atlas-colour-surface-and-environment-inference-v1"
SEASON_CONFIGURATIONS = (
    "all_dates",
    "same_calendar_month_edges",
    "same_local_solar_quarter_edges",
)


@dataclass(frozen=True)
class SpeciesTransitionSurface:
    status: str
    surface: np.ndarray
    detectable: np.ndarray
    retained_edges: int
    detectable_cells: int
    nonconstant_components: int
    geometry: EdgeCellGeometry | None


@dataclass(frozen=True)
class SpectralCohortTest:
    name: str
    basis: MoranSignBasis
    flower_coefficients: np.ndarray
    overlay_names: tuple[str, ...]
    overlay_coefficients: np.ndarray
    cells: int


def validate_colour_inference_contract(contract: Mapping[str, Any]) -> None:
    if contract.get("protocol") != COLOUR_INFERENCE_PROTOCOL:
        raise ValueError("unexpected atlas colour-inference protocol")
    if contract.get("status") != (
        "prospectively_frozen_before_any_scaleout_candidate_pixel_or_colour_environment_join"
    ):
        raise ValueError("colour inference was not frozen before scale-out pixels")
    if any(value is not False for value in contract.get("outcome_firewall", {}).values()):
        raise ValueError("colour-inference outcome firewall is open")
    vectors = contract.get("colour_vectors", {})
    if vectors.get("primary") != ["flower_L_mean", "flower_a_mean", "flower_b_mean"]:
        raise ValueError("primary atlas Lab vector changed")
    if vectors.get("background_diagnostic") != [
        "background_L_mean",
        "background_a_mean",
        "background_b_mean",
    ]:
        raise ValueError("background diagnostic vector changed")
    transition = contract.get("transition_surface", {})
    if (
        transition.get("knn_k") != 5
        or transition.get("scales_km") != [100, 250, 500]
        or transition.get("minimum_edges_per_species_cell") != 2
        or transition.get("minimum_retained_edges_per_species_configuration") != 100
        or transition.get("minimum_detectable_cells_per_species_configuration") != 10
        or transition.get("minimum_evaluable_species_per_primary_cohort") != 20
    ):
        raise ValueError("species-transition rules changed")
    if tuple(contract.get("season_configurations", {}).get("sensitivities", ())) != (
        "same_calendar_month_edges",
        "same_local_solar_quarter_edges",
    ):
        raise ValueError("season sensitivity rules changed")
    environment = contract.get("environmental_families", {})
    if environment.get("primary") != [
        "macroclimate_boundary",
        "land_cover_boundary",
        "ecoregion_boundary",
    ] or environment.get("minimum_test_cells") != 20:
        raise ValueError("environmental test family changed")
    inference = contract.get("joint_inference", {})
    if (
        inference.get("randomizations") != 9999
        or inference.get("seed") != 20260831
        or inference.get("alpha") != 0.05
    ):
        raise ValueError("joint environmental inference changed")


def robust_standardize_lab(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Median/IQR-standardize Lab, retaining zero-IQR components as exact zero."""

    values = np.asarray(values, dtype=float)
    if values.ndim != 2 or values.shape[1] != 3 or values.shape[0] < 2:
        raise ValueError("Lab values must have shape (n>=2, 3)")
    if not np.isfinite(values).all():
        raise ValueError("Lab values must be finite")
    median = np.median(values, axis=0)
    q25, q75 = np.quantile(values, [0.25, 0.75], axis=0)
    iqr = q75 - q25
    variable = iqr > 0
    if not np.any(variable):
        raise ValueError("all Lab components have zero IQR")
    standardized = np.zeros_like(values, dtype=float)
    standardized[:, variable] = (values[:, variable] - median[variable]) / iqr[variable]
    return standardized, variable


def _season_edge_subset(
    edges: np.ndarray,
    distances: np.ndarray,
    season_labels: np.ndarray | None,
) -> tuple[np.ndarray, np.ndarray]:
    if season_labels is None:
        return edges, distances
    labels = np.asarray(season_labels)
    if labels.ndim != 1 or len(labels) <= int(edges.max(initial=-1)):
        raise ValueError("season labels must match observations")
    if any(str(value) == "" for value in labels):
        raise ValueError("season labels cannot be empty")
    keep = labels[edges[:, 0]] == labels[edges[:, 1]]
    return edges[keep], distances[keep]


def build_species_transition_surface(
    latitude: Sequence[float] | np.ndarray,
    longitude: Sequence[float] | np.ndarray,
    lab_values: np.ndarray,
    *,
    grid: EqualAreaGrid,
    scale_km: int,
    season_labels: Sequence[object] | np.ndarray | None = None,
    knn_k: int = 5,
    minimum_edges_per_cell: int = 2,
    minimum_retained_edges: int = 100,
    minimum_detectable_cells: int = 10,
) -> SpeciesTransitionSurface:
    """Build one species surface without pooling colours across taxa."""

    latitude = np.asarray(latitude, dtype=float)
    longitude = np.asarray(longitude, dtype=float)
    empty = np.full(grid.n_cells, np.nan, dtype=float)
    detectable_empty = np.zeros(grid.n_cells, dtype=bool)
    try:
        standardized, variable = robust_standardize_lab(lab_values)
        if standardized.shape[0] != latitude.size or longitude.shape != latitude.shape:
            raise ValueError("colour and coordinate denominators differ")
        edges, distances = spherical_knn_edges(latitude, longitude, k=knn_k)
        if season_labels is not None and len(season_labels) != latitude.size:
            raise ValueError("season labels must match observations")
        edges, distances = _season_edge_subset(
            edges,
            distances,
            None if season_labels is None else np.asarray(season_labels),
        )
        geometry = build_edge_cell_geometry(
            latitude,
            longitude,
            edges,
            distances,
            grid=grid,
            max_edge_km=float(scale_km),
            min_edges_per_cell=minimum_edges_per_cell,
        )
    except ValueError:
        return SpeciesTransitionSurface(
            status="not_evaluable",
            surface=empty,
            detectable=detectable_empty,
            retained_edges=0,
            detectable_cells=0,
            nonconstant_components=0,
            geometry=None,
        )
    retained_edges = len(geometry.retained_edges)
    detectable_cells = int(np.count_nonzero(geometry.detectable))
    if retained_edges < minimum_retained_edges or detectable_cells < minimum_detectable_cells:
        return SpeciesTransitionSurface(
            status="not_evaluable",
            surface=empty,
            detectable=geometry.detectable,
            retained_edges=retained_edges,
            detectable_cells=detectable_cells,
            nonconstant_components=int(np.count_nonzero(variable)),
            geometry=geometry,
        )
    edge_scores = edge_colour_discontinuity(
        standardized, geometry.retained_edges
    )
    surface = cell_mean_intensity(average_rank_intensity(edge_scores), geometry)
    return SpeciesTransitionSurface(
        status="evaluable",
        surface=surface,
        detectable=geometry.detectable,
        retained_edges=retained_edges,
        detectable_cells=detectable_cells,
        nonconstant_components=int(np.count_nonzero(variable)),
        geometry=geometry,
    )


def equal_species_cohort_surface(
    surfaces: Sequence[SpeciesTransitionSurface],
    *,
    minimum_species: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Average evaluable within-species ranks with one vote per detectable species."""

    evaluable = [item for item in surfaces if item.status == "evaluable"]
    if len(evaluable) < minimum_species:
        raise ValueError("cohort has too few evaluable species")
    cell_count = {len(item.surface) for item in evaluable}
    if len(cell_count) != 1:
        raise ValueError("species surfaces use different grids")
    detectable = np.vstack([item.detectable for item in evaluable])
    values = np.vstack([item.surface for item in evaluable])
    opportunity = detectable.sum(axis=0).astype(int)
    numerator = np.where(detectable, values, 0.0).sum(axis=0)
    cohort = np.full(opportunity.shape, np.nan, dtype=float)
    present = opportunity > 0
    cohort[present] = numerator[present] / opportunity[present]
    return cohort, opportunity


def prepare_spectral_cohort_test(
    name: str,
    flower_surface: Sequence[float] | np.ndarray,
    opportunity: Sequence[float] | np.ndarray,
    cell_ids: Sequence[int] | np.ndarray,
    latitude: Sequence[float] | np.ndarray,
    longitude: Sequence[float] | np.ndarray,
    overlays: Mapping[str, Sequence[float] | np.ndarray],
    *,
    n_lon: int,
    n_sinlat: int,
    minimum_cells: int = 20,
) -> SpectralCohortTest:
    """Prepare one fixed cohort/family basis on exact finite common cells."""

    flower = np.asarray(flower_surface, dtype=float)
    weight = np.asarray(opportunity, dtype=float)
    cells = np.asarray(cell_ids, dtype=int)
    lat = np.asarray(latitude, dtype=float)
    lon = np.asarray(longitude, dtype=float)
    if not overlays:
        raise ValueError("at least one environmental overlay is required")
    arrays = {key: np.asarray(value, dtype=float) for key, value in overlays.items()}
    expected = flower.shape
    if (
        flower.ndim != 1
        or any(array.shape != expected for array in (weight, cells, lat, lon))
        or any(array.shape != expected for array in arrays.values())
    ):
        raise ValueError("spectral test inputs must be equal one-dimensional vectors")
    finite = np.isfinite(flower) & np.isfinite(weight) & (weight > 0)
    finite &= np.isfinite(lat) & np.isfinite(lon)
    for array in arrays.values():
        finite &= np.isfinite(array)
    if int(np.count_nonzero(finite)) < minimum_cells:
        raise ValueError("too few common finite cells for environmental inference")
    cells = cells[finite]
    if len(np.unique(cells)) != len(cells):
        raise ValueError("environmental test cell IDs must be unique")
    adjacency = equal_area_rook_adjacency(
        cells, n_lon=n_lon, n_sinlat=n_sinlat
    )
    basis = build_moran_sign_basis(
        adjacency,
        weight[finite],
        geographic_design(lat[finite], lon[finite]),
    )
    names = tuple(sorted(arrays))
    return SpectralCohortTest(
        name=name,
        basis=basis,
        flower_coefficients=residual_coefficients(flower[finite], basis),
        overlay_names=names,
        overlay_coefficients=np.column_stack(
            [residual_coefficients(arrays[key][finite], basis) for key in names]
        ),
        cells=len(cells),
    )


def _observed_family_statistic(test: SpectralCohortTest) -> tuple[float, dict[str, float]]:
    flower = test.flower_coefficients
    overlays = test.overlay_coefficients
    denominator = np.linalg.norm(flower) * np.linalg.norm(overlays, axis=0)
    correlations = (flower @ overlays) / denominator
    return float(np.max(correlations)), {
        name: float(value)
        for name, value in zip(test.overlay_names, correlations, strict=True)
    }


def _sign_null_family(
    test: SpectralCohortTest,
    *,
    randomizations: int,
    rng: np.random.Generator,
) -> np.ndarray:
    flower = test.flower_coefficients
    overlays = test.overlay_coefficients
    product = flower[:, None] * overlays
    denominator = np.linalg.norm(flower) * np.linalg.norm(overlays, axis=0)
    output = np.empty(randomizations, dtype=float)
    completed = 0
    while completed < randomizations:
        current = min(2048, randomizations - completed)
        signs = rng.integers(0, 2, size=(current, len(flower)), dtype=np.int8)
        signs = signs.astype(float) * 2.0 - 1.0
        correlations = (signs @ product) / denominator[None, :]
        output[completed : completed + current] = np.max(correlations, axis=1)
        completed += current
    return output


def joint_equal_cohort_spectral_test(
    groups: Mapping[str, Sequence[SpectralCohortTest]],
    *,
    randomizations: int,
    rng: np.random.Generator,
    expected_cohorts: int = 8,
) -> dict[str, Any]:
    """Apply one maximum null to all cohort and equal-cohort family statistics."""

    if not groups or randomizations < 1 or expected_cohorts < 1:
        raise ValueError("joint inference needs groups, cohorts and randomizations")
    global_max = np.full(randomizations, -np.inf, dtype=float)
    working: dict[str, dict[str, Any]] = {}
    for group_name in sorted(groups):
        tests = tuple(groups[group_name])
        if len(tests) != expected_cohorts or len({test.name for test in tests}) != len(tests):
            raise ValueError(f"{group_name}: expected unique complete cohort tests")
        cohort_observed = []
        cohort_null = []
        cohort_details = []
        for test in tests:
            observed, by_overlay = _observed_family_statistic(test)
            null = _sign_null_family(
                test, randomizations=randomizations, rng=rng
            )
            global_max = np.maximum(global_max, null)
            cohort_observed.append(observed)
            cohort_null.append(null)
            cohort_details.append(
                {
                    "cohort_test": test.name,
                    "cells": test.cells,
                    "overlay_names": list(test.overlay_names),
                    "observed_by_overlay": by_overlay,
                    "family_statistic": observed,
                }
            )
        aggregate_observed = float(np.mean(cohort_observed))
        aggregate_null = np.mean(np.vstack(cohort_null), axis=0)
        global_max = np.maximum(global_max, aggregate_null)
        working[group_name] = {
            "aggregate_observed": aggregate_observed,
            "aggregate_null": aggregate_null,
            "cohorts": cohort_details,
        }
    results: dict[str, Any] = {}
    for group_name, item in working.items():
        observed = item["aggregate_observed"]
        own_null = item.pop("aggregate_null")
        item["unadjusted_p"] = float(
            (np.count_nonzero(own_null >= observed) + 1) / (randomizations + 1)
        )
        item["familywise_adjusted_p"] = float(
            (np.count_nonzero(global_max >= observed) + 1) / (randomizations + 1)
        )
        results[group_name] = item
    return {
        "randomizations": randomizations,
        "groups": results,
        "joint_null_minimum": float(np.min(global_max)),
        "joint_null_median": float(np.median(global_max)),
        "joint_null_maximum": float(np.max(global_max)),
    }

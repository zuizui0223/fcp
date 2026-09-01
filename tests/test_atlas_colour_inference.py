from __future__ import annotations

import copy
import json
from pathlib import Path

import numpy as np
import pytest

from fcp_pipeline.atlas_colour_inference import (
    SpectralCohortTest,
    build_species_transition_surface,
    equal_species_cohort_surface,
    joint_equal_cohort_spectral_test,
    prepare_spectral_cohort_test,
    robust_standardize_lab,
    validate_colour_inference_contract,
)
from fcp_pipeline.shared_transition_surface import EqualAreaGrid
from scripts.data.run_jbi_atlas_environmental_inference import species_free_map_row


CONTRACT = json.loads(
    Path("docs/supporting/jbi_atlas_colour_surface_contract_v1.json").read_text(
        encoding="utf-8"
    )
)


def test_colour_surface_contract_is_prospective_and_pinned() -> None:
    validate_colour_inference_contract(CONTRACT)
    changed = copy.deepcopy(CONTRACT)
    changed["joint_inference"]["randomizations"] = 999
    with pytest.raises(ValueError, match="joint environmental inference"):
        validate_colour_inference_contract(changed)


def test_robust_lab_keeps_zero_iqr_component_as_zero() -> None:
    values = np.array(
        [[10.0, 1.0, 2.0], [20.0, 1.0, 4.0], [30.0, 1.0, 8.0], [40.0, 1.0, 10.0]]
    )
    standardized, variable = robust_standardize_lab(values)
    assert variable.tolist() == [True, False, True]
    assert np.all(standardized[:, 1] == 0.0)
    with pytest.raises(ValueError, match="all Lab components"):
        robust_standardize_lab(np.ones((5, 3)))


def test_public_atlas_map_withholds_species_and_cohort_labels() -> None:
    public = species_free_map_row(
        {
            "measurement_id": "FCPM-1",
            "species": "Must stay hidden",
            "cohort_id": "C01",
            "latitude": 35.0,
            "longitude": 135.0,
        },
        [50.0, 10.0, -5.0],
    )
    assert set(public) == {
        "measurement_id",
        "latitude",
        "longitude",
        "flower_L_mean",
        "flower_a_mean",
        "flower_b_mean",
    }


def test_species_transition_filters_season_edges_before_cell_ranks() -> None:
    rng = np.random.default_rng(12)
    latitude = rng.uniform(-2.0, 2.0, 180)
    longitude = rng.uniform(-2.0, 2.0, 180)
    values = rng.normal(size=(180, 3))
    grid = EqualAreaGrid(32, 16)
    all_dates = build_species_transition_surface(
        latitude,
        longitude,
        values,
        grid=grid,
        scale_km=500,
        minimum_retained_edges=1,
        minimum_detectable_cells=1,
    )
    same_month = build_species_transition_surface(
        latitude,
        longitude,
        values,
        grid=grid,
        scale_km=500,
        season_labels=np.arange(180) % 12,
        minimum_retained_edges=1,
        minimum_detectable_cells=1,
    )
    assert all_dates.status == "evaluable"
    assert same_month.retained_edges < all_dates.retained_edges


def test_equal_species_surface_uses_opportunity_not_missing_as_zero() -> None:
    grid = EqualAreaGrid(2, 2)
    first = build_species_transition_surface(
        np.array([0.0, 0.1, 0.2, 0.3]),
        np.array([0.0, 0.1, 0.2, 0.3]),
        np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]], dtype=float),
        grid=grid,
        scale_km=500,
        knn_k=2,
        minimum_edges_per_cell=1,
        minimum_retained_edges=1,
        minimum_detectable_cells=1,
    )
    cohort, opportunity = equal_species_cohort_surface([first], minimum_species=1)
    assert np.all(np.isnan(cohort[opportunity == 0]))
    assert np.all(np.isfinite(cohort[opportunity > 0]))


def make_spectral_test(name: str, seed: int) -> SpectralCohortTest:
    rng = np.random.default_rng(seed)
    n_lon, n_sinlat = 8, 4
    cells = np.arange(n_lon * n_sinlat)
    row, column = np.divmod(cells, n_lon)
    latitude = np.rad2deg(np.arcsin(-1.0 + (row + 0.5) * 2.0 / n_sinlat))
    longitude = -180.0 + (column + 0.5) * 360.0 / n_lon
    flower = rng.normal(size=len(cells))
    return prepare_spectral_cohort_test(
        name,
        flower,
        np.ones(len(cells)),
        cells,
        latitude,
        longitude,
        {"aligned": flower, "noise": rng.normal(size=len(cells))},
        n_lon=n_lon,
        n_sinlat=n_sinlat,
    )


def test_joint_null_protects_equal_cohort_family_and_detects_alignment() -> None:
    tests = [make_spectral_test(f"C{index:02d}", index) for index in range(1, 9)]
    result = joint_equal_cohort_spectral_test(
        {"flower|100|all_dates|macroclimate": tests},
        randomizations=199,
        rng=np.random.default_rng(99),
    )
    group = result["groups"]["flower|100|all_dates|macroclimate"]
    assert group["aggregate_observed"] == pytest.approx(1.0)
    assert group["familywise_adjusted_p"] <= 0.01
    assert len(group["cohorts"]) == 8

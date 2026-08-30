from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import numpy as np
import pytest

from fcp_pipeline.shared_transition_surface import EqualAreaGrid, build_edge_cell_geometry

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "analysis" / "run_jbi_ch1_stage_b_shared_transition_concentration.py"
SPEC = spec_from_file_location("stage_b_runner", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
STAGE_B = module_from_spec(SPEC)
SPEC.loader.exec_module(STAGE_B)


class ColourAccessTrap(dict):
    def __getitem__(self, key):
        if key == "values":
            raise AssertionError("geometry selection accessed observed colour values")
        return super().__getitem__(key)


def _geometry_only_species_inputs():
    latitude = np.array([0.0, 0.0, 0.1, 0.1])
    longitude = np.array([0.0, 0.1, 0.0, 0.1])
    edges = np.array([[0, 1], [0, 2], [2, 3]])
    distance = np.array([10.0, 20.0, 30.0])
    return {
        species: ColourAccessTrap(
            {
                "values": object(),
                "latitude": latitude,
                "longitude": longitude,
                "base_edges": edges,
                "base_edge_distance_km": distance,
            }
        )
        for species in STAGE_B.EXPECTED_SPECIES
    }


def test_geometry_primary_selection_never_reads_observed_colour_values():
    contract = {
        "geometry_only_primary_selection": {
            "minimum_retained_edges_per_species_cell": 1,
            "candidate_max_edge_km_in_priority_order": [100],
            "candidate_equal_area_grids_in_priority_order_within_edge_cap": [
                {"n_lon": 8, "n_sinlat": 4}
            ],
            "priority_rule": "first passing configuration",
            "passing_criteria": {
                "minimum_retained_edges_per_species": 3,
                "minimum_cells_A_ge_2": 1,
                "minimum_cells_A_ge_3": 1,
                "minimum_species_with_any_A_ge_2_opportunity": 6,
            },
        },
        "shared_surface": {"minimum_detectable_species": 2},
    }

    audit, cache, selected = STAGE_B.freeze_geometry_candidates(
        _geometry_only_species_inputs(),
        contract,
    )

    assert audit["selection_used_colour_values"] is False
    assert selected == "cap_100km_grid_8x4"
    assert audit["candidate_configurations"][0]["passes_geometry_only_criteria"] is True
    assert selected in cache


def test_monte_carlo_upper_summary_uses_plus_one_correction():
    summary = STAGE_B.mc_upper_summary(4.0, np.array([1.0, 2.0, 3.0, 5.0]))
    assert summary["p_upper_tail"] == pytest.approx((1 + 1) / (4 + 1))
    assert summary["concentration_excess"] == pytest.approx(1.25)
    assert summary["null_quantiles"]["p50"] == pytest.approx(2.5)


def test_complete_surface_keeps_fixed_label_independent_opportunity():
    grid = EqualAreaGrid(n_lon=8, n_sinlat=4)
    latitude = np.array([0.0, 0.0, 0.1, 0.1])
    longitude = np.array([0.0, 0.1, 0.0, 0.1])
    edges = np.array([[0, 1], [0, 2], [2, 3]])
    distance = np.array([10.0, 20.0, 30.0])
    geometries = [
        build_edge_cell_geometry(
            latitude,
            longitude,
            edges,
            distance,
            grid=grid,
            max_edge_km=100.0,
            min_edges_per_cell=1,
        )
        for _ in STAGE_B.EXPECTED_SPECIES
    ]
    species_inputs = {
        species: {
            "values": np.array([[0.0], [0.1], [0.9], [1.0]]),
        }
        for species in STAGE_B.EXPECTED_SPECIES
    }

    observed = STAGE_B.compute_surface(
        species_inputs,
        geometries,
        permuted_values=None,
        min_detectable_species=2,
    )
    permuted_vectors = [
        species_inputs[species]["values"][[3, 2, 1, 0]]
        for species in STAGE_B.EXPECTED_SPECIES
    ]
    permuted = STAGE_B.compute_surface(
        species_inputs,
        geometries,
        permuted_values=permuted_vectors,
        min_detectable_species=2,
    )

    np.testing.assert_array_equal(observed[1], permuted[1])
    np.testing.assert_array_equal(observed[3], permuted[3])
    assert np.count_nonzero(np.isfinite(observed[2])) >= 1

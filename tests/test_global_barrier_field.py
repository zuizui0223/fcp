from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from fcp_pipeline.global_barrier_field import (
    barrier_field,
    concentration_permutation_test,
    equal_area_grid_centers,
    prepare_barrier_geometry,
    spherical_midpoint,
    within_species_rank_scores,
)


def _repeated_positions_fixture(n_species: int = 40) -> pd.DataFrame:
    # Every species has the same eight broad geographic opportunities.  This is
    # deliberately simple: the null must retain those opportunities and change
    # only which within-species colour discontinuity rank lands at each place.
    longitude = np.array([-150.0, -100.0, -50.0, -5.0, 5.0, 50.0, 100.0, 150.0])
    latitude = np.zeros_like(longitude)
    rows = []
    for species_id in range(n_species):
        for edge_id, (lat, lon) in enumerate(zip(latitude, longitude)):
            rows.append(
                {
                    "species": f"species_{species_id:03d}",
                    "edge_id": edge_id,
                    "midpoint_latitude": lat,
                    "midpoint_longitude": lon,
                }
            )
    return pd.DataFrame(rows)


def test_spherical_midpoint_handles_dateline_without_planar_artifact():
    lat, lon = spherical_midpoint([0.0], [170.0], [0.0], [-170.0])
    assert lat[0] == pytest.approx(0.0, abs=1e-8)
    assert abs(abs(lon[0]) - 180.0) < 1e-8


def test_within_species_midrank_has_mean_half_per_species():
    edges = pd.DataFrame(
        {
            "species": ["a", "a", "a", "b", "b"],
            "raw": [1.0, 3.0, 2.0, 5.0, 5.0],
        }
    )
    rank = within_species_rank_scores(edges, score_column="raw")
    for species in ["a", "b"]:
        idx = edges.index[edges["species"] == species].to_numpy()
        assert rank[idx].mean() == pytest.approx(0.5)
    assert np.all((rank > 0.0) & (rank < 1.0))


def test_opportunity_is_equal_species_normalized_even_when_edge_counts_differ():
    edges = pd.DataFrame(
        {
            "species": ["a", "a", "b"],
            "midpoint_latitude": [0.0, 0.0, 0.0],
            "midpoint_longitude": [-10.0, 10.0, 0.0],
        }
    )
    geometry = prepare_barrier_geometry(
        edges,
        grid=equal_area_grid_centers(36, 18),
        kernel_km=500.0,
    )
    # Kernel rows each sum to one; equal-species edge weights therefore make
    # total opportunity exactly one per species irrespective of edge count.
    assert geometry.opportunity.sum() == pytest.approx(2.0)


def test_planted_recurrent_barrier_is_detected_against_species_conditioned_null():
    edges = _repeated_positions_fixture(n_species=40)
    # Two central longitudes are the recurrent colour-transition zone for every
    # species; other geographic opportunities carry weaker discontinuities.
    raw = np.tile(np.array([1.0, 2.0, 3.0, 7.0, 8.0, 4.0, 5.0, 6.0]), 40)
    score_edges = edges.copy()
    score_edges["raw_colour_divergence"] = raw
    rank = within_species_rank_scores(score_edges, score_column="raw_colour_divergence")
    geometry = prepare_barrier_geometry(
        score_edges,
        grid=equal_area_grid_centers(36, 18),
        kernel_km=500.0,
        cutoff_multiplier=3.0,
    )
    result = concentration_permutation_test(
        geometry,
        rank,
        minimum_distinct_species=5,
        permutations=99,
        seed=20260904,
    )
    assert result["observed"].concentration > result["null_mean"]
    assert result["p_upper"] <= 0.02


def test_balanced_no_barrier_fixture_is_not_declared_supported():
    edges = _repeated_positions_fixture(n_species=40)
    # Rotate the complete within-species rank set across species.  Because 40 is
    # divisible by eight, every geographic opportunity receives every rank the
    # same number of times, so the observed consensus is deliberately flat.
    base = np.arange(1.0, 9.0)
    raw = []
    for species_id in range(40):
        raw.extend(np.roll(base, species_id % 8))
    score_edges = edges.copy()
    score_edges["raw_colour_divergence"] = np.asarray(raw)
    rank = within_species_rank_scores(score_edges, score_column="raw_colour_divergence")
    geometry = prepare_barrier_geometry(
        score_edges,
        grid=equal_area_grid_centers(36, 18),
        kernel_km=500.0,
        cutoff_multiplier=3.0,
    )
    observed = barrier_field(geometry, rank, minimum_distinct_species=5)
    result = concentration_permutation_test(
        geometry,
        rank,
        minimum_distinct_species=5,
        permutations=99,
        seed=20260904,
    )
    assert observed.concentration >= 0.0
    assert result["p_upper"] > 0.05


def test_cells_below_distinct_species_support_are_not_biological_zeroes():
    edges = pd.DataFrame(
        {
            "species": ["a", "b", "c", "d"],
            "midpoint_latitude": [0.0, 0.0, 0.0, 0.0],
            "midpoint_longitude": [0.0, 0.0, 0.0, 0.0],
        }
    )
    geometry = prepare_barrier_geometry(
        edges,
        grid=equal_area_grid_centers(36, 18),
        kernel_km=500.0,
    )
    field = barrier_field(geometry, [0.2, 0.4, 0.6, 0.8], minimum_distinct_species=5)
    assert not field.evaluable.any()
    assert np.isnan(field.field).all()
    assert np.isnan(field.concentration)

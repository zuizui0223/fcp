import numpy as np
import pytest

from fcp_pipeline.continuous_colour_boundaries import (
    average_rank_intensity,
    edge_colour_discontinuity,
    opportunity_weighted_concentration,
    shared_boundary_intensity,
    species_conditioned_vector_permutation,
    weighted_graph_discontinuity,
)
from fcp_pipeline.spatial_boundaries import shared_boundary_strength


def test_vector_permutation_preserves_complete_rows_within_species_only():
    values = np.array(
        [
            [1.0, 10.0],
            [2.0, 20.0],
            [3.0, 30.0],
            [101.0, 110.0],
            [102.0, 120.0],
            [103.0, 130.0],
        ]
    )
    species = np.array(["a", "a", "a", "b", "b", "b"], dtype=object)

    permuted = species_conditioned_vector_permutation(
        values,
        species,
        rng=np.random.default_rng(17),
    )

    assert permuted.shape == values.shape
    for sp in np.unique(species):
        idx = np.flatnonzero(species == sp)
        before = sorted(map(tuple, values[idx].tolist()))
        after = sorted(map(tuple, permuted[idx].tolist()))
        assert after == before

    # The two columns remain coupled as original observation rows; feature values are
    # never independently shuffled or recombined.
    original_rows = {tuple(row) for row in values.tolist()}
    assert all(tuple(row) in original_rows for row in permuted.tolist())


def test_binary_scalar_edge_distance_reduces_exactly_to_discrete_mismatch():
    labels = np.array([0.0, 1.0, 1.0, 0.0])
    edges = np.array([[0, 1], [1, 2], [2, 3], [0, 3]])

    continuous = edge_colour_discontinuity(labels[:, None], edges)
    discrete = (labels[edges[:, 0]] != labels[edges[:, 1]]).astype(float)

    np.testing.assert_array_equal(continuous, discrete)

    weights = np.array([1.0, 2.0, 3.0, 4.0])
    q_continuous = weighted_graph_discontinuity(continuous, weights=weights)
    q_discrete = float(np.dot(weights, discrete) / weights.sum())
    assert q_continuous == pytest.approx(q_discrete)


def test_multidimensional_edge_distance_uses_rms_feature_difference():
    values = np.array([[0.0, 0.0], [3.0, 4.0], [3.0, 0.0]])
    edges = np.array([[0, 1], [1, 2]])

    scores = edge_colour_discontinuity(values, edges)

    np.testing.assert_allclose(
        scores,
        np.array([np.sqrt((3.0**2 + 4.0**2) / 2.0), np.sqrt((0.0**2 + 4.0**2) / 2.0)]),
    )


def test_average_rank_intensity_uses_average_ranks_for_ties():
    scores = np.array([1.0, 3.0, 3.0, 5.0])

    intensity = average_rank_intensity(scores)

    np.testing.assert_allclose(intensity, np.array([0.0, 0.5, 0.5, 1.0]))
    np.testing.assert_allclose(average_rank_intensity(np.array([7.0])), np.array([0.5]))


def test_shared_continuous_intensity_reduces_to_discrete_shared_strength_for_binary_input():
    boundary = np.array(
        [
            [True, False, True, False],
            [True, True, False, False],
            [False, True, True, True],
        ]
    )
    detectable = np.array(
        [
            [True, True, False, False],
            [True, False, True, False],
            [True, True, True, False],
        ]
    )

    continuous, A_continuous = shared_boundary_intensity(
        boundary.astype(float),
        detectable,
        min_detectable_species=2,
    )
    discrete, A_discrete = shared_boundary_strength(
        boundary,
        detectable,
        min_detectable_species=2,
    )

    np.testing.assert_array_equal(A_continuous, A_discrete)
    np.testing.assert_allclose(continuous, discrete, equal_nan=True)


def test_shared_intensity_ignores_undetectable_values_and_preserves_not_evaluable_cells():
    intensity = np.array(
        [
            [0.2, np.nan, 1.0],
            [0.8, 0.4, np.nan],
            [0.6, 0.9, np.nan],
        ]
    )
    detectable = np.array(
        [
            [True, False, True],
            [True, True, False],
            [False, True, False],
        ]
    )

    shared, A = shared_boundary_intensity(
        intensity,
        detectable,
        min_detectable_species=2,
    )

    np.testing.assert_array_equal(A, np.array([2, 2, 1]))
    np.testing.assert_allclose(shared[:2], np.array([0.5, 0.65]))
    assert np.isnan(shared[2])


def test_opportunity_weighted_concentration_is_zero_for_spatially_even_intensity():
    shared = np.array([0.4, 0.4, 0.4, np.nan])
    opportunity = np.array([1, 2, 7, 0])

    assert opportunity_weighted_concentration(shared, opportunity) == pytest.approx(0.0)


def test_opportunity_weighted_concentration_matches_weighted_variance():
    shared = np.array([0.0, 1.0, 0.5])
    opportunity = np.array([1, 3, 2])

    expected_mean = (1 * 0.0 + 3 * 1.0 + 2 * 0.5) / 6
    expected = (
        1 * (0.0 - expected_mean) ** 2
        + 3 * (1.0 - expected_mean) ** 2
        + 2 * (0.5 - expected_mean) ** 2
    ) / 6

    assert opportunity_weighted_concentration(shared, opportunity) == pytest.approx(expected)


def test_invalid_edges_and_non_evaluable_global_surface_fail_explicitly():
    values = np.array([[0.0], [1.0]])
    with pytest.raises(ValueError, match="self edges"):
        edge_colour_discontinuity(values, np.array([[0, 0]]))

    with pytest.raises(ValueError, match="no evaluable cells"):
        opportunity_weighted_concentration(np.array([np.nan]), np.array([0]))

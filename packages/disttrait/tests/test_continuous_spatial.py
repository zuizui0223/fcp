import numpy as np
import pytest

from disttrait import (
    continuous_spatial_permutation_null,
    continuous_spatial_rho,
    pairwise_absolute_difference,
    spatial_permutation_null_from_distance_matrix,
    spatial_rho_from_distance_matrix,
)


def test_continuous_spatial_rho_uses_absolute_trait_difference() -> None:
    lat = [0.0, 0.0, 0.0, 0.0, 0.0]
    lon = [0.0, 1.0, 2.0, 3.0, 4.0]
    values = [0.0, 1.0, 2.0, 3.0, 4.0]
    assert continuous_spatial_rho(lat, lon, values) == pytest.approx(1.0, abs=1e-14)
    pairwise = pairwise_absolute_difference(values)
    assert len(pairwise) == 10


def test_distance_matrix_api_matches_continuous_wrapper() -> None:
    lat = [35.0, 35.2, 35.5, 36.0, 36.5]
    lon = [135.0, 135.1, 135.4, 135.8, 136.2]
    values = np.array([1.0, 2.5, 2.0, 5.0, 7.0])
    matrix = np.abs(values[:, None] - values[None, :])
    observed_generic = spatial_rho_from_distance_matrix(lat, lon, matrix)
    observed_wrapper = continuous_spatial_rho(lat, lon, values)
    assert observed_generic == pytest.approx(observed_wrapper, abs=1e-14)


def test_continuous_permutation_null_is_deterministic_and_matches_matrix_api() -> None:
    lat = [35.0, 35.2, 35.5, 36.0, 36.5]
    lon = [135.0, 135.1, 135.4, 135.8, 136.2]
    values = np.array([1.0, 2.5, 2.0, 5.0, 7.0])
    matrix = np.abs(values[:, None] - values[None, :])
    obs1, null1 = continuous_spatial_permutation_null(
        lat, lon, values, n_permutations=19, seed=19, key="continuous"
    )
    obs2, null2 = spatial_permutation_null_from_distance_matrix(
        lat, lon, matrix, n_permutations=19, seed=19, key="continuous"
    )
    assert obs1 == pytest.approx(obs2, abs=1e-14)
    np.testing.assert_allclose(null1, null2, rtol=0, atol=1e-14)
    assert len(null1) == 19

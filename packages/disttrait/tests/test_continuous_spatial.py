import numpy as np
import pytest

from disttrait import (
    absolute_pairwise,
    continuous_spatial_permutation_null,
    continuous_spatial_rho,
)


def test_absolute_pairwise_scalar_trait():
    observed = absolute_pairwise([1.0, 3.0, 6.0])
    np.testing.assert_allclose(observed, [2.0, 5.0, 3.0], rtol=0, atol=0)


def test_continuous_spatial_rho_detects_monotone_spatial_change():
    lat = [0.0] * 5
    lon = [0.0, 1.0, 2.0, 3.0, 4.0]
    values = [1.0, 1.5, 2.5, 4.0, 6.0]
    rho = continuous_spatial_rho(lat, lon, values)
    assert rho > 0.8


def test_continuous_spatial_null_is_deterministic():
    lat = [0.0] * 5
    lon = [0.0, 1.0, 2.0, 3.0, 4.0]
    values = [1.0, 1.5, 2.5, 4.0, 6.0]
    obs1, null1 = continuous_spatial_permutation_null(
        lat, lon, values, n_permutations=19, seed=20260919, key="fixture"
    )
    obs2, null2 = continuous_spatial_permutation_null(
        lat, lon, values, n_permutations=19, seed=20260919, key="fixture"
    )
    assert obs1 == pytest.approx(obs2, abs=0, rel=0)
    np.testing.assert_array_equal(null1, null2)
    assert len(null1) == 19

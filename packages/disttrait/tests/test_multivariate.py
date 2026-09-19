import numpy as np

from disttrait import (
    continuous_spatial_rho,
    euclidean_pairwise,
    multivariate_spatial_permutation_null,
    multivariate_spatial_rho,
)


def test_euclidean_pairwise_known_geometry():
    values = np.array([
        [0.0, 0.0],
        [3.0, 4.0],
        [3.0, 0.0],
    ])
    observed = euclidean_pairwise(values)
    np.testing.assert_allclose(observed, [5.0, 3.0, 4.0], rtol=0, atol=1e-14)


def test_one_column_multivariate_rho_matches_scalar_continuous_rho():
    lat = [0.0, 0.0, 0.0, 0.0, 0.0]
    lon = [0.0, 1.0, 2.0, 3.0, 4.0]
    values = np.array([0.0, 0.1, 0.4, 0.9, 1.6])
    scalar = continuous_spatial_rho(lat, lon, values)
    multi = multivariate_spatial_rho(lat, lon, values[:, None])
    assert scalar == multi


def test_multivariate_permutation_null_is_deterministic():
    lat = [0.0] * 5
    lon = [0.0, 1.0, 2.0, 3.0, 4.0]
    values = np.array([
        [0.0, 0.0],
        [0.2, -0.1],
        [0.5, 0.3],
        [0.9, 0.7],
        [1.4, 1.2],
    ])
    obs1, null1 = multivariate_spatial_permutation_null(
        lat,
        lon,
        values,
        n_permutations=19,
        seed=17,
        key="fixture",
    )
    obs2, null2 = multivariate_spatial_permutation_null(
        lat,
        lon,
        values,
        n_permutations=19,
        seed=17,
        key="fixture",
    )
    assert obs1 == obs2
    np.testing.assert_array_equal(null1, null2)


def test_complete_rows_are_the_permutation_unit():
    lat = [0.0] * 4
    lon = [0.0, 1.0, 2.0, 3.0]
    # Perfect within-row relationship. A complete-row permutation retains the
    # geometry of the four 2-D points rather than creating impossible cross-
    # dimension combinations.
    values = np.array([
        [0.0, 0.0],
        [1.0, 2.0],
        [2.0, 4.0],
        [3.0, 6.0],
    ])
    _, null = multivariate_spatial_permutation_null(
        lat,
        lon,
        values,
        n_permutations=7,
        seed=21,
        key="rows",
    )
    assert np.isfinite(null).all()
    assert len(null) == 7

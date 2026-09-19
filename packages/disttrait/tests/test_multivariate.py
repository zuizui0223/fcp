import hashlib

import numpy as np
from scipy.stats import spearmanr

from disttrait import (
    continuous_spatial_rho,
    euclidean_pairwise,
    great_circle_pairwise_km,
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


def test_optimized_multivariate_null_matches_direct_reranking():
    lat = [0.0] * 6
    lon = [0.0, 0.7, 1.9, 3.1, 4.8, 6.0]
    values = np.array([
        [0.0, 0.0],
        [0.2, -0.1],
        [0.5, 0.4],
        [0.9, 0.2],
        [1.4, 1.1],
        [1.8, 0.7],
    ])
    seed = 31
    key = "equivalence"
    n_permutations = 13

    observed, optimized = multivariate_spatial_permutation_null(
        lat,
        lon,
        values,
        n_permutations=n_permutations,
        seed=seed,
        key=key,
    )

    geo = great_circle_pairwise_km(lat, lon)
    u, v = np.triu_indices(len(values), k=1)
    direct_observed = spearmanr(
        geo,
        np.linalg.norm(values[u] - values[v], axis=1),
    ).statistic
    direct = []
    for idx in range(n_permutations):
        payload = f"{seed}|{key}|{idx}".encode()
        permutation_seed = int.from_bytes(
            hashlib.sha256(payload).digest()[:8],
            "little",
        )
        p = np.random.default_rng(permutation_seed).permutation(len(values))
        direct.append(
            spearmanr(
                geo,
                np.linalg.norm(values[p[u]] - values[p[v]], axis=1),
            ).statistic
        )

    assert observed == direct_observed
    np.testing.assert_allclose(optimized, direct, rtol=0, atol=2e-15)

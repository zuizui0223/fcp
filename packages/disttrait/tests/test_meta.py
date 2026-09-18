import numpy as np
import pytest

from disttrait import (
    random_effects_from_groups,
    random_effects_slope_permutation_test,
    random_effects_slope_summary,
    species_slope_estimate,
)


def test_species_slope_estimate_recovers_signed_ols_slope():
    estimate = species_slope_estimate(
        [-1.0, 0.0, 1.0, 2.0],
        [0.0, 1.0, 2.0, 3.0],
    )
    assert estimate.slope == pytest.approx(1.0, abs=1e-14)
    assert estimate.n == 4
    assert estimate.variance > 0


def test_random_effects_summary_detects_shared_mean_without_heterogeneity():
    result = random_effects_slope_summary(
        [0.8, 1.0, 1.2, 1.0],
        [0.04, 0.04, 0.04, 0.04],
    )
    assert result.fixed_mean == pytest.approx(1.0, abs=1e-14)
    assert result.random_mean == pytest.approx(1.0, abs=1e-14)
    assert result.tau2 == 0.0
    assert result.p_heterogeneity > 0.5
    assert result.p_mean_two_sided < 1e-10
    assert result.p_omnibus < 1e-9


def test_random_effects_summary_detects_opposing_slope_heterogeneity():
    result = random_effects_slope_summary(
        [1.0, -1.0, 1.0, -1.0],
        [0.01, 0.01, 0.01, 0.01],
    )
    assert result.random_mean == pytest.approx(0.0, abs=1e-14)
    assert result.p_mean_two_sided == pytest.approx(1.0, abs=1e-14)
    assert result.tau2 > 1.0
    assert result.p_heterogeneity < 1e-20
    assert result.p_omnibus < 1e-19


def test_random_effects_from_groups_preserves_signed_species_slopes():
    x = np.linspace(-1.0, 1.0, 12)
    groups = [
        (x, 0.8 * x + np.linspace(-0.1, 0.1, 12)),
        (x, -0.8 * x + np.linspace(0.1, -0.1, 12)),
        (x, 1.1 * x + np.linspace(-0.05, 0.05, 12)),
        (x, -1.1 * x + np.linspace(0.05, -0.05, 12)),
    ]
    result, estimates = random_effects_from_groups(groups)
    signs = np.sign([x.slope for x in estimates]).tolist()
    assert signs == [1.0, -1.0, 1.0, -1.0]
    assert result.p_heterogeneity < 0.05


def test_permutation_calibrated_meta_detects_shared_signed_effect_reproducibly():
    x = np.linspace(-1.0, 1.0, 14)
    noise = 0.08 * np.sin(np.arange(14, dtype=float))
    groups = [(x, 0.8 * x + noise + 0.02 * i) for i in range(8)]
    a = random_effects_slope_permutation_test(
        groups, n_permutations=99, seed=17, key="shared"
    )
    b = random_effects_slope_permutation_test(
        groups, n_permutations=99, seed=17, key="shared"
    )
    assert a.p_mean_permutation <= 0.02
    assert a.p_omnibus_permutation <= 0.04
    np.testing.assert_array_equal(a.null_random_mean, b.null_random_mean)
    np.testing.assert_array_equal(a.null_q, b.null_q)


def test_permutation_calibrated_meta_detects_opposing_direction_heterogeneity():
    x = np.linspace(-1.0, 1.0, 14)
    noise = 0.08 * np.cos(np.arange(14, dtype=float))
    groups = []
    for i in range(10):
        sign = 1.0 if i % 2 == 0 else -1.0
        groups.append((x, sign * 0.9 * x + noise))
    result = random_effects_slope_permutation_test(
        groups, n_permutations=99, seed=23, key="opposing"
    )
    assert result.p_mean_permutation > 0.2
    assert result.p_heterogeneity_permutation <= 0.02
    assert result.p_omnibus_permutation <= 0.04
    assert result.observed.tau2 > 0

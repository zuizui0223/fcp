import numpy as np

from disttrait import (
    alignment_statistic,
    categorical_diversity,
    gini_simpson,
    hellinger_rows,
    one_vs_rest_contrast,
    permute_rows_within_strata,
    two_mode_axis,
)


def test_gini_simpson_accepts_counts_or_probabilities():
    assert gini_simpson([10, 0]) == 0.0
    assert gini_simpson([1, 1]) == 0.5
    assert gini_simpson([0.5, 0.5]) == 0.5
    assert categorical_diversity(["a", "a", "b", "b"]) == 0.5


def test_hellinger_rows_are_normalized_before_transform():
    x = hellinger_rows(np.array([[2.0, 2.0], [9.0, 1.0]]))
    assert np.allclose(x[0] ** 2, [0.5, 0.5])
    assert np.allclose(x[1] ** 2, [0.9, 0.1])


def test_two_mode_axis_recovers_one_vs_rest_geometry():
    x = np.array([
        [0.95, 0.03, 0.02],
        [0.90, 0.05, 0.05],
        [0.92, 0.04, 0.04],
        [0.05, 0.50, 0.45],
        [0.02, 0.48, 0.50],
        [0.04, 0.46, 0.50],
    ])
    result = two_mode_axis(x)
    q = one_vs_rest_contrast(3, focal_index=0)
    assert result.minor_fraction == 0.5
    assert alignment_statistic(result.unit_axis[None, :], q) > 0.98


def test_stratified_row_permutation_preserves_each_stratum_multiset():
    values = np.array([[1], [2], [3], [4], [5], [6]])
    strata = ["a", "a", "a", "b", "b", "b"]
    out = permute_rows_within_strata(values, strata, rng=np.random.default_rng(2))
    assert sorted(out[:3, 0].tolist()) == [1, 2, 3]
    assert sorted(out[3:, 0].tolist()) == [4, 5, 6]

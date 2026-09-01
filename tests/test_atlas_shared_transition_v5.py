from __future__ import annotations

import numpy as np
import pytest

from fcp_pipeline.atlas_shared_transition_v5 import (
    build_coexceedance_reference,
    build_detectability_matrix,
    coexceedance_scan_statistic,
    coexceedance_z,
    conditional_rank_scan_null,
    high_transition_mask_from_scores,
    monte_carlo_p,
    signal_recovery_rates,
)


def toy_detectability() -> np.ndarray:
    return build_detectability_matrix(
        [
            [0, 1, 2, 3, 4],
            [0, 1, 2, 3, 5],
            [0, 1, 2, 4, 5],
            [0, 1, 3, 4, 5],
            [0, 2, 3, 4, 5],
            [1, 2, 3, 4, 5],
        ],
        n_cells=6,
    )


def test_reference_matches_conditional_expectation() -> None:
    D = toy_detectability()
    ref = build_coexceedance_reference(
        D, high_transition_quantile=0.8, min_detectable_species=4
    )
    assert np.array_equal(ref.high_counts, np.ones(6, dtype=int))
    assert np.allclose(ref.high_probabilities, np.full(6, 0.2))
    assert len(ref.valid_cell_ids) == 6
    expected = D.sum(axis=0) * 0.2
    variance = D.sum(axis=0) * 0.2 * 0.8
    assert np.allclose(ref.expected_count, expected)
    assert np.allclose(ref.standard_deviation, np.sqrt(variance))


def test_top_tail_is_exact_and_deterministic_under_ties() -> None:
    D = toy_detectability()
    ref = build_coexceedance_reference(
        D, high_transition_quantile=0.8, min_detectable_species=4
    )
    scores = np.full(D.shape, np.nan)
    scores[D] = 1.0
    high1 = high_transition_mask_from_scores(scores, ref)
    high2 = high_transition_mask_from_scores(scores, ref)
    assert np.array_equal(high1, high2)
    assert np.array_equal(high1.sum(axis=1), ref.high_counts)
    assert not np.any(high1 & ~D)


def test_coexceedance_rejects_flags_outside_opportunity() -> None:
    D = toy_detectability()
    ref = build_coexceedance_reference(
        D, high_transition_quantile=0.8, min_detectable_species=4
    )
    scores = np.full(D.shape, np.nan)
    scores[D] = np.arange(D.sum(), dtype=float)
    high = high_transition_mask_from_scores(scores, ref)
    bad = high.copy()
    row, cell = np.argwhere(~D)[0]
    bad[row, cell] = True
    with pytest.raises(ValueError, match="outside detectable"):
        coexceedance_z(bad, ref)


def test_handcrafted_shared_cell_has_positive_scan_excess() -> None:
    D = toy_detectability()
    ref = build_coexceedance_reference(
        D, high_transition_quantile=0.8, min_detectable_species=4
    )
    scores = np.full(D.shape, np.nan)
    scores[D] = 0.0
    for i in range(D.shape[0]):
        shared = np.flatnonzero(D[i])[0]
        scores[i, shared] = 100.0
    high = high_transition_mask_from_scores(scores, ref)
    assert coexceedance_scan_statistic(high, ref) > 1.0


def test_conditional_null_is_reproducible_and_finite() -> None:
    ref = build_coexceedance_reference(
        toy_detectability(), high_transition_quantile=0.8, min_detectable_species=4
    )
    a = conditional_rank_scan_null(
        ref, n_permutations=99, rng=np.random.default_rng(123)
    )
    b = conditional_rank_scan_null(
        ref, n_permutations=99, rng=np.random.default_rng(123)
    )
    assert a.shape == (99,)
    assert np.isfinite(a).all()
    assert np.array_equal(a, b)


def test_monte_carlo_p_uses_plus_one_correction() -> None:
    null = np.array([1.0, 2.0, 3.0, 4.0])
    assert monte_carlo_p(5.0, null) == pytest.approx(1 / 5)
    assert monte_carlo_p(3.0, null) == pytest.approx(3 / 5)


def test_signal_recovery_protocol_is_seed_reproducible() -> None:
    D = np.ones((12, 64), dtype=bool)
    ref = build_coexceedance_reference(
        D, high_transition_quantile=0.9, min_detectable_species=4
    )
    from fcp_pipeline.atlas_shared_transition_v5 import equal_area_cell_xyz

    xyz = equal_area_cell_xyz(n_lon=16, n_sinlat=4)
    kwargs = dict(
        n_repetitions=3,
        n_permutations=49,
        alpha=0.05,
        effect_sizes=(0.0, 0.5, 1.0, 2.0),
        boundary_sigma_radians=0.1,
        seed=20260901,
    )
    first = signal_recovery_rates(ref, xyz, **kwargs)
    second = signal_recovery_rates(ref, xyz, **kwargs)
    assert np.array_equal(first["null_distribution"], second["null_distribution"])
    assert first["null_result"]["rate"] == second["null_result"]["rate"]
    assert [x["rate"] for x in first["shared_results"]] == [
        x["rate"] for x in second["shared_results"]
    ]

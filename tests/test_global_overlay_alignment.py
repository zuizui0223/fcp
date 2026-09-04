from __future__ import annotations

import numpy as np
import pytest

from fcp_pipeline.global_overlay_alignment import (
    evaluate_primary_overlays,
    holm_adjust,
    overlay_alignment_permutation_test,
    weighted_spearman,
)


def test_weighted_spearman_recovers_monotone_alignment():
    x = np.arange(1.0, 11.0)
    y = x**3
    w = np.linspace(1.0, 2.0, len(x))
    assert weighted_spearman(x, y, w) == pytest.approx(1.0)
    assert weighted_spearman(x, -y, w) == pytest.approx(-1.0)


def test_planted_external_barrier_alignment_beats_colour_permutation_null():
    rng = np.random.default_rng(20260904)
    n_cells = 120
    predictor = np.linspace(-1.0, 1.0, n_cells)
    observed = predictor + rng.normal(0.0, 0.05, n_cells)
    # Null fields preserve the same cell support but no systematic alignment to
    # the fixed external surface.
    null = np.vstack([rng.permutation(observed) for _ in range(199)])
    result = overlay_alignment_permutation_test(
        predictor_name="planted_barrier",
        observed_field=observed,
        null_fields=null,
        predictor_surface=predictor,
        opportunity=np.ones(n_cells),
        minimum_cells=50,
    )
    assert result.status == "evaluated"
    assert result.observed_rho > 0.95
    assert result.p_upper <= 0.01


def test_anti_alignment_is_not_supported_by_upper_tail_test():
    rng = np.random.default_rng(7)
    n_cells = 100
    predictor = np.arange(n_cells, dtype=float)
    observed = -predictor
    null = np.vstack([rng.permutation(observed) for _ in range(99)])
    result = overlay_alignment_permutation_test(
        predictor_name="wrong_direction",
        observed_field=observed,
        null_fields=null,
        predictor_surface=predictor,
        opportunity=np.ones(n_cells),
        minimum_cells=50,
    )
    assert result.observed_rho < -0.99
    assert result.p_upper > 0.95


def test_insufficient_external_coverage_is_not_recast_as_zero_alignment():
    observed = np.linspace(0.0, 1.0, 100)
    predictor = np.full(100, np.nan)
    predictor[:20] = np.linspace(0.0, 1.0, 20)
    null = np.tile(observed, (9, 1))
    result = overlay_alignment_permutation_test(
        predictor_name="sparse",
        observed_field=observed,
        null_fields=null,
        predictor_surface=predictor,
        opportunity=np.ones(100),
        minimum_cells=50,
    )
    assert result.status == "not_evaluable_external_surface_coverage"
    assert np.isnan(result.observed_rho)
    assert result.n_cells == 20


def test_holm_adjustment_is_monotone_and_family_wise():
    adjusted = holm_adjust({"a": 0.01, "b": 0.03, "c": 0.20, "d": 0.80})
    assert adjusted["a"] == pytest.approx(0.04)
    assert adjusted["b"] == pytest.approx(0.09)
    assert adjusted["c"] == pytest.approx(0.40)
    assert adjusted["d"] == pytest.approx(0.80)


def test_primary_overlays_cannot_run_before_colour_field_and_stability_gates():
    args = dict(
        observed_field=np.arange(60, dtype=float),
        null_fields=np.tile(np.arange(60, dtype=float), (9, 1)),
        predictor_surfaces={"x": np.arange(60, dtype=float)},
        opportunity=np.ones(60),
        minimum_cells=50,
    )
    assert evaluate_primary_overlays(g1_supported=False, g2_stable=True, **args)["status"] == "not_run_g1_hierarchical_gate"
    assert evaluate_primary_overlays(g1_supported=True, g2_stable=False, **args)["status"] == "not_run_g2_stability_gate"


def test_fixed_primary_family_uses_holm_and_keeps_sparse_predictor_not_evaluable():
    rng = np.random.default_rng(13)
    n_cells = 100
    observed = np.linspace(-1.0, 1.0, n_cells)
    null = np.vstack([rng.permutation(observed) for _ in range(199)])
    sparse = np.full(n_cells, np.nan)
    sparse[:10] = observed[:10]
    payload = evaluate_primary_overlays(
        g1_supported=True,
        g2_stable=True,
        observed_field=observed,
        null_fields=null,
        predictor_surfaces={
            "strong": observed,
            "noise": rng.normal(size=n_cells),
            "sparse": sparse,
        },
        opportunity=np.ones(n_cells),
        minimum_cells=50,
        alpha=0.05,
    )
    assert payload["status"] == "evaluated"
    assert payload["results"]["strong"]["supported"] is True
    assert payload["results"]["sparse"]["status"] == "not_evaluable_external_surface_coverage"
    assert np.isnan(payload["results"]["sparse"]["p_holm"])

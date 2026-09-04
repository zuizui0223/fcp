from __future__ import annotations

import numpy as np
import pytest

from fcp_pipeline.global_pair_context_profiles import build_pair_context_profiles


def _reference():
    x = np.linspace(0.0, 1.0, 200)
    return {
        "pollinator_turnover": x,
        "climate_turnover": x,
        "edaphic_turnover": x,
        "terrain_turnover": x,
        "marine_gap": x,
        "major_river_crossing": x,
        "biogeographic_boundary_crossing": x,
        "mountain_boundary_crossing": x,
    }


def test_sympatric_pair_gets_low_turnover_similarity_labels():
    reference = _reference()
    scores = {name: [0.10] for name in reference}
    payload = build_pair_context_profiles(
        pair_ids=["p1"],
        geography=["sympatric"],
        pair_scores=scores,
        reference_scores=reference,
        minimum_reference=100,
    )
    profile = payload["profiles"]["p1"]
    assert profile.geography == "sympatric"
    assert "pollinator-similar" in profile.labels
    assert "climate-similar" in profile.labels
    assert "edaphic-similar" in profile.labels
    assert "terrain-similar" in profile.labels
    assert all(not label.endswith("-separated") for label in profile.labels)


def test_allopatric_pair_gets_high_separation_labels():
    reference = _reference()
    scores = {name: [0.90] for name in reference}
    payload = build_pair_context_profiles(
        pair_ids=["p2"],
        geography=["allopatric"],
        pair_scores=scores,
        reference_scores=reference,
        minimum_reference=100,
    )
    profile = payload["profiles"]["p2"]
    assert "climate-separated" in profile.labels
    assert "edaphic-separated" in profile.labels
    assert "marine-gap-separated" in profile.labels
    assert "biogeographic-boundary-separated" in profile.labels
    assert all(not label.endswith("-similar") for label in profile.labels)


def test_midrange_pair_can_remain_unlabelled():
    reference = _reference()
    scores = {name: [0.50] for name in reference}
    payload = build_pair_context_profiles(
        pair_ids=["p3", "p4"],
        geography=["sympatric", "allopatric"],
        pair_scores={name: [0.50, 0.50] for name in reference},
        reference_scores=reference,
        minimum_reference=100,
    )
    assert payload["profiles"]["p3"].labels == ()
    assert payload["profiles"]["p4"].labels == ()


def test_sparse_external_reference_is_not_evaluable_not_zero():
    reference = _reference()
    reference["climate_turnover"] = np.arange(20.0)
    payload = build_pair_context_profiles(
        pair_ids=["p1"],
        geography=["allopatric"],
        pair_scores={name: [1.0] for name in reference},
        reference_scores=reference,
        minimum_reference=100,
    )
    assert payload["thresholds"]["climate_turnover"]["status"] == "not_evaluable_reference_coverage"
    assert "climate-separated" not in payload["profiles"]["p1"].labels


def test_constant_reference_cannot_generate_a_favourable_label():
    reference = _reference()
    reference["marine_gap"] = np.ones(200)
    payload = build_pair_context_profiles(
        pair_ids=["p1"],
        geography=["allopatric"],
        pair_scores={name: [1.0] for name in reference},
        reference_scores=reference,
        minimum_reference=100,
    )
    assert payload["thresholds"]["marine_gap"]["status"] == "not_evaluable_constant_reference"
    assert "marine-gap-separated" not in payload["profiles"]["p1"].labels


def test_invalid_geography_is_rejected():
    with pytest.raises(ValueError, match="sympatric or allopatric"):
        build_pair_context_profiles(
            pair_ids=["p1"],
            geography=["parapatric"],
            pair_scores={},
            reference_scores={},
        )

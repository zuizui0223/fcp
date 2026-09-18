import copy
import json
from pathlib import Path

import numpy as np
import pytest

from fcp_pipeline.p500_white_measurement_control import (
    RECEIPT_PATH,
    classify_measurement_control,
    digital_highlight_metrics,
    inspect_candidate_metadata,
    response_blind_high_clip_mask,
    validate_preopening_firewall,
)


ROOT = Path(__file__).resolve().parents[1]
PREFLIGHT = ROOT / "results/polymorphism_h2_p500_candidate_metadata_20260913/result.json"
RECEIPT = ROOT / RECEIPT_PATH
PROTOCOL = ROOT / "docs/POLYMORPHISM_H2_P500_WHITE_MEASUREMENT_CONTROL_PROTOCOL_20260913.md"


def _load(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def _synthetic_preopening():
    # Only for unit testing the accepted schema; never save as real provenance.
    return {"gate_results": {"overall_pass": True},
            "blindness": dict(_load(RECEIPT)["required_blind_flags"])}


def test_real_acquisition_is_not_a_preopening_execution_receipt():
    assert PROTOCOL.exists()
    with pytest.raises(RuntimeError, match="metadata gate did not pass"):
        validate_preopening_firewall(_load(PREFLIGHT), _load(RECEIPT))
    report = inspect_candidate_metadata(_load(PREFLIGHT))
    assert report["candidate_metadata_supported"] is True
    assert report["execution_chronology_verified"] is False
    assert report["opening_authorized"] is False


def test_synthetic_complete_schema_with_unchanged_receipt():
    validate_preopening_firewall(_synthetic_preopening(), _load(RECEIPT))


def test_firewall_fails_closed_if_any_required_outcome_was_opened():
    preflight = _synthetic_preopening()
    receipt = _load(RECEIPT)
    tampered = copy.deepcopy(preflight)
    tampered["blindness"]["image_bytes_opened"] = True
    with pytest.raises(RuntimeError, match="already opened"):
        validate_preopening_firewall(tampered, receipt)


def test_firewall_rejects_estimand_mutation():
    preflight = _load(PREFLIGHT)
    receipt = _load(RECEIPT)
    tampered = copy.deepcopy(receipt)
    tampered["h2_estimand_mutated"] = True
    with pytest.raises(RuntimeError, match="estimand mutation"):
        validate_preopening_firewall(preflight, tampered)


def test_digital_highlight_metrics_are_non_circular_rgb_diagnostics():
    rgb = np.array(
        [
            [255, 255, 255],
            [250, 120, 20],
            [249, 249, 249],
            [10, 20, 30],
        ],
        dtype=np.uint8,
    )
    out = digital_highlight_metrics(rgb)
    assert out["bit_depth"] == 8
    assert out["clip_fraction"] == pytest.approx(0.25)
    assert out["near_clip_fraction"] == pytest.approx(0.50)
    assert 0.0 <= out["luminance_q99"] <= 1.0


def test_response_blind_high_clip_set_uses_only_technical_values():
    values = np.array([0.0] * 18 + [0.02, 0.20])
    mask, threshold = response_blind_high_clip_mask(values)
    assert threshold >= 0.01
    assert mask.dtype == bool
    assert mask[-1]


def test_clear_requires_equivalence_ci_and_h2_stability():
    decision = classify_measurement_control(
        or_ci_low=0.90,
        or_ci_high=1.10,
        primary_h2_support=True,
        sensitivity_h2_support=True,
        primary_h2_species_n=400,
        sensitivity_h2_species_n=390,
    )
    assert decision == "CLEAR"


def test_large_coupling_or_support_loss_is_flagged():
    assert (
        classify_measurement_control(
            or_ci_low=1.30,
            or_ci_high=1.60,
            primary_h2_support=True,
            sensitivity_h2_support=True,
            primary_h2_species_n=400,
            sensitivity_h2_species_n=390,
        )
        == "FLAGGED"
    )
    assert (
        classify_measurement_control(
            or_ci_low=0.90,
            or_ci_high=1.10,
            primary_h2_support=True,
            sensitivity_h2_support=False,
            primary_h2_species_n=400,
            sensitivity_h2_species_n=390,
        )
        == "FLAGGED"
    )


def test_overlap_or_insufficient_retention_is_indeterminate():
    assert (
        classify_measurement_control(
            or_ci_low=0.75,
            or_ci_high=1.10,
            primary_h2_support=True,
            sensitivity_h2_support=True,
            primary_h2_species_n=400,
            sensitivity_h2_species_n=390,
        )
        == "INDETERMINATE"
    )


@pytest.mark.parametrize("field,value", [
    ("primary_technical_predictor", "white_fraction"),
    ("equivalence_or_interval", [0.01, 100]),
    ("minimum_h2_species_retention", 0),
    ("current_h2_result", {"opened": True}),
    ("current_measurement_control_result", "CLEAR"),
    ("candidate_species_n", 499),
    ("preflight_overall_pass", 1),
    ("technical_metrics", []),
    ("decision_states", ["CLEAR"]),
])
def test_reject_frozen_receipt_mutations(field, value):
    receipt = _load(RECEIPT)
    receipt[field] = value
    with pytest.raises(RuntimeError, match="field mismatch"):
        validate_preopening_firewall(_synthetic_preopening(), receipt)


@pytest.mark.parametrize("field,value", [
    ("primary_h2_species_n", 0), ("primary_h2_species_n", True),
    ("primary_h2_species_n", 400.5), ("sensitivity_h2_species_n", 401),
    ("sensitivity_h2_species_n", -1), ("sensitivity_h2_species_n", float("nan")),
    ("primary_h2_support", "false"), ("sensitivity_h2_support", 1),
    ("or_ci_low", "0.9"), ("or_ci_high", float("inf")),
    ("or_ci_low", True),
])
def test_invalid_decision_inputs_never_clear(field, value):
    args = dict(or_ci_low=.9, or_ci_high=1.1, primary_h2_support=True,
                sensitivity_h2_support=True, primary_h2_species_n=400,
                sensitivity_h2_species_n=390)
    args[field] = value
    assert classify_measurement_control(**args) == "INDETERMINATE"


@pytest.mark.parametrize("bit_depth", [True, 8.5, 0, -1, 33])
def test_invalid_decode_depth(bit_depth):
    with pytest.raises(ValueError, match="bit_depth"):
        digital_highlight_metrics(np.array([[1, 2, 3]]), bit_depth=bit_depth)


def test_fractional_rgb_is_not_an_integer_decode():
    with pytest.raises(ValueError, match="integral"):
        digital_highlight_metrics(np.array([[1.5, 2, 3]]))


def test_16_bit_highlight_boundaries():
    result = digital_highlight_metrics(np.array([[65535, 0, 0], [64250, 0, 0],
                                                [64249, 0, 0]], dtype=np.uint16), bit_depth=16)
    assert result["clip_fraction"] == pytest.approx(1/3)
    assert result["near_clip_fraction"] == pytest.approx(2/3)


@pytest.mark.parametrize("flag", ["image_pixels_opened", "flower_colour_opened", "palette_opened", "D_opened", "H2_W_opened"])
def test_actual_metadata_open_flags_fail_closed(flag):
    candidate = _load(PREFLIGHT)
    candidate["outcome_firewall"][flag] = True
    report = inspect_candidate_metadata(candidate)
    assert report["candidate_metadata_supported"] is False
    assert report["opening_authorized"] is False


def test_incomplete_metadata_fails_closed():
    assert inspect_candidate_metadata({})["candidate_metadata_supported"] is False
    assert (
        classify_measurement_control(
            or_ci_low=0.90,
            or_ci_high=1.10,
            primary_h2_support=True,
            sensitivity_h2_support=True,
            primary_h2_species_n=400,
            sensitivity_h2_species_n=350,
        )
        == "INDETERMINATE"
    )

import copy
import json
from pathlib import Path

import numpy as np
import pytest

from fcp_pipeline.p500_white_measurement_control import (
    RECEIPT_PATH,
    classify_measurement_control,
    digital_highlight_metrics,
    response_blind_high_clip_mask,
    validate_preopening_firewall,
)


ROOT = Path(__file__).resolve().parents[1]
PREFLIGHT = ROOT / "results/polymorphism_h2_p500_candidate_metadata_20260913/result.json"
RECEIPT = ROOT / RECEIPT_PATH
PROTOCOL = ROOT / "docs/POLYMORPHISM_H2_P500_WHITE_MEASUREMENT_CONTROL_PROTOCOL_20260913.md"


def _load(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def test_frozen_preopening_receipt_is_valid_and_protocol_exists():
    assert PROTOCOL.exists()
    validate_preopening_firewall(_load(PREFLIGHT), _load(RECEIPT))


def test_firewall_fails_closed_if_any_required_outcome_was_opened():
    preflight = _load(PREFLIGHT)
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

"""Prospective P500 white-state measurement-control firewall.

This module contains only pre-specified technical diagnostics and gate logic.
It does not fetch images, classify flower colour, estimate D, or run H2.
"""

from __future__ import annotations

from collections.abc import Mapping
import json
from numbers import Integral, Real
from typing import Any

import numpy as np


FREEZE_ID = "polymorphism-h2-p500-white-measurement-control-20260913"
BASE_SHA = "c55fb858579074d2e3a4ce1a71b5a2d266e11011"
PROTOCOL_PATH = "docs/POLYMORPHISM_H2_P500_WHITE_MEASUREMENT_CONTROL_PROTOCOL_20260913.md"
PREFLIGHT_PATH = "results/polymorphism_h2_p500_candidate_metadata_20260913/result.json"
RECEIPT_PATH = "results/polymorphism_h2_p500_white_measurement_control_freeze_20260913/receipt.json"

REQUIRED_BLIND_FLAGS = (
    "image_results_requested",
    "filtered_specimens_opened",
    "image_urls_opened",
    "image_bytes_opened",
    "flower_colour_opened",
    "palette_opened",
    "D_opened",
    "H1_opened",
    "H2_W_opened",
)

PRIMARY_PREDICTOR = "near_clip_fraction"
FORBIDDEN_TECHNICAL_PREDICTORS = frozenset(
    {
        "hsv_saturation",
        "hsl_saturation",
        "lab_chroma",
        "palette_distance",
        "white_fraction",
    }
)
EQUIVALENCE_OR_LOW = 0.80
EQUIVALENCE_OR_HIGH = 1.25
MIN_H2_SPECIES_RETENTION = 0.90


def validate_preopening_firewall(
    preflight: Mapping[str, Any], receipt: Mapping[str, Any]
) -> None:
    """Fail closed unless the frozen pre-opening contract is intact."""

    if not isinstance(preflight, Mapping) or not isinstance(receipt, Mapping):
        raise RuntimeError("preflight and receipt must be mappings")
    if receipt.get("freeze_id") != FREEZE_ID:
        raise RuntimeError("P500 measurement-control freeze_id mismatch")
    if receipt.get("base_sha") != BASE_SHA:
        raise RuntimeError("P500 measurement-control base SHA mismatch")
    if receipt.get("protocol_path") != PROTOCOL_PATH:
        raise RuntimeError("P500 measurement-control protocol path mismatch")
    if receipt.get("preflight_result_path") != PREFLIGHT_PATH:
        raise RuntimeError("P500 measurement-control preflight path mismatch")
    if receipt.get("h2_estimand_mutated") is not False:
        raise RuntimeError("H2 estimand mutation is forbidden")
    if receipt.get("image_or_colour_outcomes_opened_at_freeze") is not False:
        raise RuntimeError("freeze must precede P500 image/colour opening")

    expected = {
        "status": "preopening_protocol_frozen",
        "base_branch": "analysis/polymorphism-42111-h1-h2-gates-20260912",
        "freeze_branch": "analysis/p500-white-measurement-control-freeze-20260913",
        "preflight_overall_pass": True,
        "candidate_species_n": 500,
        "eligible_species_n": 499,
        "request_error_species_n": 0,
        "primary_technical_predictor": PRIMARY_PREDICTOR,
        "technical_metrics": ["clip_fraction", "near_clip_fraction", "luminance_q99"],
        "decision_states": ["CLEAR", "FLAGGED", "INDETERMINATE"],
        "equivalence_or_interval": [EQUIVALENCE_OR_LOW, EQUIVALENCE_OR_HIGH],
        "minimum_h2_species_retention": MIN_H2_SPECIES_RETENTION,
        "current_measurement_control_result": None,
        "current_h2_result": None,
        "allowed_next_step": "response_blind_technical_image_materialization_after_freeze_verification",
    }
    for key, value in expected.items():
        # JSON identity rejects false/0 and true/1 substitutions as well.
        if key not in receipt or json.dumps(receipt[key], sort_keys=True) != json.dumps(value, sort_keys=True):
            raise RuntimeError(f"frozen receipt field mismatch: {key}")

    gates = preflight.get("gate_results", {})
    if not isinstance(gates, Mapping) or gates.get("overall_pass") is not True:
        raise RuntimeError("P500 candidate metadata gate did not pass")

    blindness = preflight.get("blindness", {})
    if not isinstance(blindness, Mapping):
        raise RuntimeError("preflight blindness must be a mapping")
    missing = [key for key in REQUIRED_BLIND_FLAGS if key not in blindness]
    if missing:
        raise RuntimeError(f"missing preflight blindness flags: {missing}")
    opened = [key for key in REQUIRED_BLIND_FLAGS if blindness.get(key) is not False]
    if opened:
        raise RuntimeError(f"P500 outcomes were already opened: {opened}")

    frozen_flags = receipt.get("required_blind_flags", {})
    if not isinstance(frozen_flags, Mapping) or set(frozen_flags) != set(REQUIRED_BLIND_FLAGS):
        raise RuntimeError("receipt blindness flag set mismatch")
    if any(frozen_flags[key] is not False for key in REQUIRED_BLIND_FLAGS):
        raise RuntimeError("receipt does not freeze all required flags as false")

    forbidden = set(receipt.get("forbidden_technical_predictors", ()))
    if forbidden != set(FORBIDDEN_TECHNICAL_PREDICTORS):
        raise RuntimeError("forbidden technical predictor set mismatch")


def inspect_candidate_metadata(candidate: Mapping[str, Any]) -> dict[str, Any]:
    """Inspect the actual acquisition schema, never fabricate opening evidence.

    This is separate from validate_preopening_firewall. Acquisition provenance
    and a static freeze declaration are not a verified execution chronology.
    """
    try:
        passed = (
            candidate["analysis"] == "polymorphism_h2_p500_candidate_metadata_acquisition"
            and candidate["status"] == "complete_fresh_metadata_draw_before_pixels"
            and candidate["decision"]["verdict"] == "P500_CANDIDATE_METADATA_GATE_PASS"
            and candidate["decision"]["no_species_replacement"] is True
            and candidate["decision"]["no_target_relaxation"] is True
            and type(candidate["P500"]["species"]) is int
            and candidate["P500"]["species"] == 500
            and type(candidate["capacity"]["full100_species"]) is int
            and candidate["capacity"]["full100_species"] == 499
            and type(candidate["transport"]["request_errors"]) is int
            and candidate["transport"]["request_errors"] == 0
            and all(candidate["outcome_firewall"][key] is False for key in (
                "image_pixels_opened", "flower_colour_opened", "palette_opened",
                "D_opened", "H2_W_opened", "morph_opened", "H3_predictors_used"))
        )
    except (KeyError, TypeError):
        passed = False
    return {
        "candidate_metadata_supported": bool(passed),
        "execution_chronology_verified": False,
        "opening_authorized": False,
        "reason": "candidate_metadata_is_not_a_verified_preopening_execution_record",
    }


def _srgb_to_linear(x: np.ndarray) -> np.ndarray:
    return np.where(x <= 0.04045, x / 12.92, ((x + 0.055) / 1.055) ** 2.4)


def digital_highlight_metrics(
    rgb_pixels: np.ndarray, *, bit_depth: int = 8
) -> dict[str, float | int]:
    """Return non-circular digital highlight diagnostics for flower-mask pixels."""

    pixels = np.asarray(rgb_pixels)
    if pixels.ndim != 2 or pixels.shape[1] != 3:
        raise ValueError("rgb_pixels must have shape (n_pixels, 3)")
    if len(pixels) == 0:
        raise ValueError("rgb_pixels must contain at least one pixel")
    if isinstance(bit_depth, (bool, np.bool_)) or not isinstance(bit_depth, Integral) or not 1 <= bit_depth <= 32:
        raise ValueError("bit_depth must be an integer in [1, 32]")
    bit_depth = int(bit_depth)
    ceiling = float((1 << bit_depth) - 1)
    x = pixels.astype(np.float64, copy=False)
    if not np.isfinite(x).all() or np.any((x < 0.0) | (x > ceiling)):
        raise ValueError("rgb values must be finite and within the decode range")
    if np.any(x != np.floor(x)):
        raise ValueError("decoded integer RGB samples must be integral")

    channel_max = x.max(axis=1)
    clip_fraction = float(np.mean(channel_max == ceiling))
    near_clip_threshold = ceiling * (250.0 / 255.0)
    near_clip_fraction = float(np.mean(channel_max >= near_clip_threshold))

    normalized = x / ceiling
    linear = _srgb_to_linear(normalized)
    luminance = (
        0.2126729 * linear[:, 0]
        + 0.7151522 * linear[:, 1]
        + 0.0721750 * linear[:, 2]
    )
    return {
        "bit_depth": bit_depth,
        "clip_fraction": clip_fraction,
        "near_clip_fraction": near_clip_fraction,
        "luminance_q99": float(np.quantile(luminance, 0.99)),
    }


def response_blind_high_clip_mask(near_clip_fraction: np.ndarray) -> tuple[np.ndarray, float]:
    """Freeze the technical-only high-clipping exclusion set."""

    values = np.asarray(near_clip_fraction, dtype=float)
    if values.ndim != 1 or len(values) == 0:
        raise ValueError("near_clip_fraction must be a non-empty 1D array")
    if not np.isfinite(values).all() or np.any((values < 0.0) | (values > 1.0)):
        raise ValueError("near_clip_fraction must lie in [0, 1]")
    threshold = float(max(0.01, np.quantile(values, 0.95)))
    return values > threshold, threshold


def classify_measurement_control(
    *,
    or_ci_low: float | None,
    or_ci_high: float | None,
    primary_h2_support: bool | None,
    sensitivity_h2_support: bool | None,
    primary_h2_species_n: int,
    sensitivity_h2_species_n: int,
) -> str:
    """Apply the frozen CLEAR/FLAGGED/INDETERMINATE rule."""

    if (
        not isinstance(primary_h2_support, (bool, np.bool_))
        or not isinstance(sensitivity_h2_support, (bool, np.bool_))
        or any(isinstance(n, (bool, np.bool_)) or not isinstance(n, Integral)
               for n in (primary_h2_species_n, sensitivity_h2_species_n))
    ):
        return "INDETERMINATE"
    if (
        or_ci_low is None
        or or_ci_high is None
        or primary_h2_support is None
        or sensitivity_h2_support is None
        or primary_h2_species_n <= 0
        or sensitivity_h2_species_n < 0
        or sensitivity_h2_species_n > primary_h2_species_n
    ):
        return "INDETERMINATE"
    if any(isinstance(x, (bool, np.bool_)) or not isinstance(x, Real)
           for x in (or_ci_low, or_ci_high)):
        return "INDETERMINATE"
    low = float(or_ci_low)
    high = float(or_ci_high)
    if not np.isfinite([low, high]).all() or low <= 0.0 or high < low:
        return "INDETERMINATE"

    retention = sensitivity_h2_species_n / primary_h2_species_n
    if retention < MIN_H2_SPECIES_RETENTION:
        return "INDETERMINATE"

    h2_changed_from_support = bool(primary_h2_support) and not bool(sensitivity_h2_support)
    ci_entirely_outside = high < EQUIVALENCE_OR_LOW or low > EQUIVALENCE_OR_HIGH
    if h2_changed_from_support or ci_entirely_outside:
        return "FLAGGED"

    ci_inside = low >= EQUIVALENCE_OR_LOW and high <= EQUIVALENCE_OR_HIGH
    h2_unchanged = bool(primary_h2_support) == bool(sensitivity_h2_support)
    if ci_inside and h2_unchanged:
        return "CLEAR"
    return "INDETERMINATE"


__all__ = [
    "BASE_SHA",
    "EQUIVALENCE_OR_HIGH",
    "EQUIVALENCE_OR_LOW",
    "FORBIDDEN_TECHNICAL_PREDICTORS",
    "FREEZE_ID",
    "MIN_H2_SPECIES_RETENTION",
    "PREFLIGHT_PATH",
    "PRIMARY_PREDICTOR",
    "PROTOCOL_PATH",
    "RECEIPT_PATH",
    "REQUIRED_BLIND_FLAGS",
    "classify_measurement_control",
    "digital_highlight_metrics",
    "inspect_candidate_metadata",
    "response_blind_high_clip_mask",
    "validate_preopening_firewall",
]

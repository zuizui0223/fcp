"""Recompute and validate immutable ROI-v4 JRC gate evidence."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

from .flower_roi_v4 import summarize_composite_gate, validate_roi_v4_contract


INTEGER_FIELDS = (
    "detector_predictions",
    "retained_instances",
    "flip_detector_predictions",
    "flip_retained_instances",
    "true_positive",
    "false_positive",
    "false_negative",
    "mask_pixels",
    "mask_pixels_inside_reference_box_union",
    "background_pixels",
    "flip_background_pixels",
    "source_annotation_boxes",
    "reference_boxes",
    "source_not_evaluable_boxes",
    "small_reference_boxes",
    "small_hit_boxes",
    "medium_reference_boxes",
    "medium_hit_boxes",
    "large_reference_boxes",
    "large_hit_boxes",
)
FLOAT_FIELDS = (
    "image_mask_pixels_inside_reference_box_union",
    "horizontal_flip_mask_iou",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def parse_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    folded = str(value).strip().casefold()
    if folded == "true":
        return True
    if folded == "false":
        return False
    raise ValueError(f"expected a boolean, got {value!r}")


def normalize_gate_row(row: Mapping[str, Any]) -> dict[str, Any]:
    normalized = dict(row)
    for field in INTEGER_FIELDS:
        normalized[field] = int(row[field])
    for field in FLOAT_FIELDS:
        normalized[field] = float(row[field])
    normalized["estimator_admitted"] = parse_bool(row["estimator_admitted"])
    return normalized


def _same(observed: Any, expected: Any) -> bool:
    if isinstance(expected, float):
        return math.isclose(float(observed), expected, rel_tol=1e-12, abs_tol=1e-12)
    return observed == expected


def validate_gate_artifacts(
    rows_path: Path,
    result_path: Path,
    contract: Mapping[str, Any],
    *,
    phase: str,
    trained_weight_sha256: str,
) -> dict[str, Any]:
    """Recompute a complete gate and enforce its firewall fields."""

    validate_roi_v4_contract(contract)
    rows = read_csv(rows_path)
    expected_images = 400 if phase == "development" else 100
    identities = {
        (row.get("image_id", ""), row.get("file_name", ""), row.get("image_sha256", ""))
        for row in rows
    }
    if (
        len(rows) != expected_images
        or len(identities) != expected_images
        or any(not all(identity) for identity in identities)
    ):
        raise RuntimeError(f"ROI v4 {phase} denominator or identity changed")
    recomputed = summarize_composite_gate(
        [normalize_gate_row(row) for row in rows], contract, phase=phase
    )
    result = json.loads(result_path.read_text(encoding="utf-8"))
    if (
        result.get("protocol") != contract["protocol"]
        or result.get("phase") != phase
        or result.get("trained_weight_sha256") != trained_weight_sha256
        or result.get("rows_sha256") != sha256(rows_path)
        or result.get("source_annotation_sha256")
        != contract["jrc_source"][f"{'train' if phase == 'development' else 'test'}_annotation_sha256"]
        or result.get("scaleout_candidate_pixels_opened") is not False
    ):
        raise RuntimeError(f"ROI v4 {phase} result provenance changed")
    expected_test_opened = phase == "locked_test"
    if result.get("jrc_test_images_decoded_or_scored") is not expected_test_opened:
        raise RuntimeError("JRC locked-test firewall changed")
    for section in ("metrics", "checks"):
        if set(result.get(section, {})) != set(recomputed[section]):
            raise RuntimeError(f"ROI v4 {phase} {section} fields changed")
        for key, expected in recomputed[section].items():
            if not _same(result[section][key], expected):
                raise RuntimeError(f"ROI v4 {phase} changed: {section}.{key}")
    for field in (
        "status",
        "jrc_locked_test_permitted",
        "scaleout_candidate_pixels_permitted",
    ):
        if result.get(field) != recomputed[field]:
            raise RuntimeError(f"ROI v4 {phase} decision changed: {field}")
    if phase == "development" and result["scaleout_candidate_pixels_permitted"] is not False:
        raise RuntimeError("development cannot authorize scale-out pixels")
    return result

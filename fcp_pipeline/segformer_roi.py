"""Prospectively frozen prompt-free flower-candidate ROI qualification."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
import math
from typing import Any

import numpy as np
from skimage.color import rgb2lab


PROTOCOL = "jbi-atlas-roi-estimator-v3"
BOX_EDGE_AMENDMENT_PROTOCOL = "jbi-atlas-roi-v3-jrc-box-edge-amendment-v1"
BOX_EDGE_AMENDMENT_V2_PROTOCOL = "jbi-atlas-roi-v3-jrc-box-edge-amendment-v2"
CANVAS_SIZE = 512


def validate_roi_v3_contract(contract: Mapping[str, Any]) -> None:
    if contract.get("protocol") != PROTOCOL:
        raise ValueError("unexpected ROI estimator protocol")
    firewall = contract.get("outcome_firewall", {})
    if (
        firewall.get("scaleout_candidate_pixels_opened") is not False
        or firewall.get("jrc_estimator_predictions_run") is not False
        or firewall.get("oxford102_estimator_predictions_run") is not False
        or firewall.get("jrc_test_images_decoded_or_scored") is not False
    ):
        raise ValueError("ROI v3 was not frozen before protected outcomes")
    estimator = contract.get("estimator", {})
    expected_estimator = {
        "model_id": "nvidia/segformer-b0-finetuned-ade-512-512",
        "revision": "489d5cd81a0b59fab9b7ea758d3548ebe99677da",
        "model_safetensors_bytes": 15_036_944,
        "model_safetensors_sha256": (
            "6ae39addd01de6b1b8bde2cf677d43a5cd733424b8d186de3f95d1c51fee23f9"
        ),
        "flower_label": 66,
        "plant_background_control_label": 17,
        "plant_fallback_permitted": False,
        "weights_redistributed_in_repository": False,
    }
    for key, value in expected_estimator.items():
        if estimator.get(key) != value:
            raise ValueError(f"ROI estimator identity changed: {key}")
    rules = contract.get("preprocessing_and_roi", {})
    expected_rules = {
        "minimum_flower_pixels": 100,
        "minimum_plant_background_control_pixels": 100,
        "horizontal_flip_mask_iou_minimum": 0.5,
        "horizontal_flip_colour_delta_e_maximum": 5.0,
    }
    for key, value in expected_rules.items():
        if rules.get(key) != value:
            raise ValueError(f"ROI post-processing rule changed: {key}")
    jrc = contract.get("jrc_field_gate", {})
    if (
        jrc.get("train_images") != 400
        or jrc.get("test_images") != 100
        or jrc.get("train_boxes") != 6992
        or jrc.get("test_boxes") != 2524
    ):
        raise ValueError("JRC source denominator changed")
    acceptance = jrc.get("acceptance", {})
    expected_acceptance = {
        "minimum_images": 100,
        "minimum_admitted_fraction": 0.8,
        "minimum_pooled_predicted_pixel_precision_inside_box_union": 0.7,
        "minimum_median_image_predicted_pixel_precision_inside_box_union": 0.7,
        "minimum_pooled_object_recall": 0.35,
        "minimum_medium_object_recall": 0.35,
        "minimum_large_object_recall": 0.5,
        "small_object_recall": "reported_without_gate",
    }
    if acceptance != expected_acceptance:
        raise ValueError("JRC acceptance criteria changed")
    proxy = contract.get("oxford102_proxy", {})
    if (
        proxy.get("locked_images") != 2040
        or proxy.get("locked_selection_salt")
        != "fcp-atlas-roi-v3-oxford102-balanced-proxy"
        or proxy.get("drop_exact_duplicates_against_full_oxford17_ids")
        != [3448, 3456, 4657, 4691, 6241]
    ):
        raise ValueError("Oxford-102 proxy selection changed")


def validate_jrc_box_edge_amendment(amendment: Mapping[str, Any]) -> None:
    if amendment.get("protocol") != BOX_EDGE_AMENDMENT_PROTOCOL:
        raise ValueError("unexpected JRC box-edge amendment")
    evidence = amendment.get("stop_evidence", {})
    audit = amendment.get("annotation_only_audit", {})
    if (
        evidence.get("locked_jrc_test_images_decoded_or_scored") is not False
        or evidence.get("scaleout_candidate_pixels_opened") is not False
        or audit.get("train_boxes") != 6992
        or audit.get("train_boxes_crossing_image_edge") != 206
        or audit.get("test_boxes") != 2524
        or audit.get("test_boxes_crossing_image_edge") != 0
    ):
        raise ValueError("JRC box-edge stop evidence changed")
    if not str(amendment.get("frozen_correction", "")).startswith(
        "intersect every official COCO box with the closed image extent"
    ):
        raise ValueError("JRC box clipping rule changed")


def validate_jrc_box_edge_amendment_v2(amendment: Mapping[str, Any]) -> None:
    if amendment.get("protocol") != BOX_EDGE_AMENDMENT_V2_PROTOCOL:
        raise ValueError("unexpected JRC box-edge v2 amendment")
    evidence = amendment.get("stop_evidence", {})
    audit = amendment.get("annotation_only_audit", {})
    zero = audit.get("zero_area_train_annotation", {})
    if (
        evidence.get("locked_jrc_test_images_decoded_or_scored") is not False
        or evidence.get("scaleout_candidate_pixels_opened") is not False
        or audit.get("train_boxes") != 6992
        or audit.get("train_boxes_with_positive_clipped_area") != 6991
        or audit.get("train_boxes_with_zero_clipped_area") != 1
        or audit.get("test_boxes") != 2524
        or audit.get("test_boxes_with_zero_clipped_area") != 0
        or zero.get("annotation_id") != 5998
        or zero.get("image_id") != 350
        or zero.get("bbox") != [1110.0, -1.0, 63.0, 1.0]
    ):
        raise ValueError("JRC zero-area annotation evidence changed")
    if not str(amendment.get("frozen_correction", "")).startswith(
        "retain the original source-box denominator"
    ):
        raise ValueError("JRC zero-area annotation rule changed")


def class_masks_from_labels(
    labels: np.ndarray, *, flower_label: int = 66, plant_label: int = 17
) -> tuple[np.ndarray, np.ndarray]:
    labels = np.asarray(labels)
    if labels.shape != (CANVAS_SIZE, CANVAS_SIZE) or labels.ndim != 2:
        raise ValueError("SegFormer labels must be a 512 by 512 canvas")
    if not np.issubdtype(labels.dtype, np.integer):
        raise ValueError("SegFormer labels must be integer class IDs")
    return labels == int(flower_label), labels == int(plant_label)


def _mask_mean_lab(rgb: np.ndarray, mask: np.ndarray) -> np.ndarray | None:
    rgb = np.asarray(rgb)
    mask = np.asarray(mask, dtype=bool)
    if rgb.shape != (CANVAS_SIZE, CANVAS_SIZE, 3) or mask.shape != rgb.shape[:2]:
        raise ValueError("RGB and mask must share the frozen 512 canvas")
    if not np.any(mask):
        return None
    lab = rgb2lab(rgb.astype(np.float32) / 255.0)
    return np.mean(lab[mask], axis=0)


def evaluate_flip_stable_admission(
    rgb: np.ndarray,
    flower_mask: np.ndarray,
    plant_mask: np.ndarray,
    flipped_flower_mask_unflipped: np.ndarray,
    contract: Mapping[str, Any],
) -> dict[str, Any]:
    """Apply the frozen size, control and horizontal-flip admission rules."""

    validate_roi_v3_contract(contract)
    flower = np.asarray(flower_mask, dtype=bool)
    plant = np.asarray(plant_mask, dtype=bool)
    flipped = np.asarray(flipped_flower_mask_unflipped, dtype=bool)
    if flower.shape != (CANVAS_SIZE, CANVAS_SIZE) or plant.shape != flower.shape:
        raise ValueError("ROI masks must use the frozen 512 canvas")
    if flipped.shape != flower.shape:
        raise ValueError("flipped ROI mask has the wrong shape")
    intersection = int(np.count_nonzero(flower & flipped))
    union = int(np.count_nonzero(flower | flipped))
    flip_iou = intersection / union if union else 0.0
    original_lab = _mask_mean_lab(rgb, flower)
    flipped_lab = _mask_mean_lab(rgb, flipped)
    colour_delta_e = (
        float(np.linalg.norm(original_lab - flipped_lab))
        if original_lab is not None and flipped_lab is not None
        else math.inf
    )
    rules = contract["preprocessing_and_roi"]
    failures = []
    if int(flower.sum()) < int(rules["minimum_flower_pixels"]):
        failures.append("insufficient_flower_pixels")
    if int(plant.sum()) < int(rules["minimum_plant_background_control_pixels"]):
        failures.append("insufficient_plant_background_control_pixels")
    if flip_iou < float(rules["horizontal_flip_mask_iou_minimum"]):
        failures.append("horizontal_flip_mask_instability")
    if colour_delta_e > float(rules["horizontal_flip_colour_delta_e_maximum"]):
        failures.append("horizontal_flip_colour_instability")
    return {
        "estimator_admitted": not failures,
        "failure_reasons": ";".join(failures),
        "flower_pixels": int(flower.sum()),
        "plant_background_control_pixels": int(plant.sum()),
        "horizontal_flip_mask_iou": float(flip_iou),
        "horizontal_flip_colour_delta_e": float(colour_delta_e),
    }


def _scaled_box(
    bbox: Sequence[float], *, source_width: int, source_height: int
) -> tuple[int, int, int, int, float] | None:
    if len(bbox) != 4 or source_width < 1 or source_height < 1:
        raise ValueError("invalid COCO box or source dimensions")
    x, y, width, height = (float(value) for value in bbox)
    if not all(math.isfinite(value) for value in (x, y, width, height)):
        raise ValueError("COCO box contains a non-finite value")
    if width <= 0 or height <= 0:
        raise ValueError("COCO box has non-positive area")
    clipped_x0 = max(0.0, min(float(source_width), x))
    clipped_y0 = max(0.0, min(float(source_height), y))
    clipped_x1 = max(0.0, min(float(source_width), x + width))
    clipped_y1 = max(0.0, min(float(source_height), y + height))
    if clipped_x1 <= clipped_x0 or clipped_y1 <= clipped_y0:
        return None
    x0 = max(0, min(CANVAS_SIZE, math.floor(clipped_x0 * CANVAS_SIZE / source_width)))
    y0 = max(0, min(CANVAS_SIZE, math.floor(clipped_y0 * CANVAS_SIZE / source_height)))
    x1 = max(0, min(CANVAS_SIZE, math.ceil(clipped_x1 * CANVAS_SIZE / source_width)))
    y1 = max(0, min(CANVAS_SIZE, math.ceil(clipped_y1 * CANVAS_SIZE / source_height)))
    if x1 <= x0 or y1 <= y0:
        raise ValueError("COCO box vanishes on the frozen canvas")
    area = ((clipped_x1 - clipped_x0) * CANVAS_SIZE / source_width) * (
        (clipped_y1 - clipped_y0) * CANVAS_SIZE / source_height
    )
    return x0, y0, x1, y1, float(area)


def score_jrc_boxes(
    flower_mask: np.ndarray,
    boxes: Sequence[Sequence[float]],
    *,
    source_width: int,
    source_height: int,
) -> dict[str, Any]:
    """Score class-66 pixels against every manually drawn JRC flower box."""

    flower = np.asarray(flower_mask, dtype=bool)
    if flower.shape != (CANVAS_SIZE, CANVAS_SIZE) or not boxes:
        raise ValueError("JRC scoring requires a 512 mask and at least one box")
    union = np.zeros_like(flower)
    hits = {"small": 0, "medium": 0, "large": 0}
    totals = {"small": 0, "medium": 0, "large": 0}
    source_not_evaluable = 0
    for bbox in boxes:
        scaled = _scaled_box(
            bbox, source_width=source_width, source_height=source_height
        )
        if scaled is None:
            source_not_evaluable += 1
            continue
        x0, y0, x1, y1, area = scaled
        union[y0:y1, x0:x1] = True
        size_bin = "small" if area < 1024 else "medium" if area < 9216 else "large"
        totals[size_bin] += 1
        intersection = int(np.count_nonzero(flower[y0:y1, x0:x1]))
        if intersection >= max(9, math.ceil(0.05 * area)):
            hits[size_bin] += 1
    predicted = int(flower.sum())
    inside = int(np.count_nonzero(flower & union))
    return {
        "predicted_flower_pixels": predicted,
        "predicted_flower_pixels_inside_box_union": inside,
        "image_predicted_pixel_precision_inside_box_union": (
            inside / predicted if predicted else 0.0
        ),
        "reference_boxes": sum(totals.values()),
        "source_annotation_boxes": len(boxes),
        "source_not_evaluable_boxes": source_not_evaluable,
        "hit_boxes": sum(hits.values()),
        **{f"{name}_reference_boxes": totals[name] for name in totals},
        **{f"{name}_hit_boxes": hits[name] for name in hits},
    }


def summarize_jrc_gate(
    rows: Sequence[Mapping[str, Any]],
    contract: Mapping[str, Any],
    *,
    phase: str,
) -> dict[str, Any]:
    """Apply the same predeclared gates to development or locked JRC rows."""

    validate_roi_v3_contract(contract)
    if phase not in {"development", "locked_test"} or not rows:
        raise ValueError("JRC phase or rows are invalid")
    acceptance = contract["jrc_field_gate"]["acceptance"]
    admitted = sum(bool(row["estimator_admitted"]) for row in rows)
    predicted = sum(int(row["predicted_flower_pixels"]) for row in rows)
    inside = sum(int(row["predicted_flower_pixels_inside_box_union"]) for row in rows)
    reference = sum(int(row["reference_boxes"]) for row in rows)
    hit = sum(int(row["hit_boxes"]) for row in rows)
    recalls = {}
    for size in ("small", "medium", "large"):
        denominator = sum(int(row[f"{size}_reference_boxes"]) for row in rows)
        numerator = sum(int(row[f"{size}_hit_boxes"]) for row in rows)
        recalls[size] = numerator / denominator if denominator else math.nan
    metrics = {
        "images": len(rows),
        "admitted_images": admitted,
        "admitted_fraction": admitted / len(rows),
        "pooled_predicted_pixel_precision_inside_box_union": (
            inside / predicted if predicted else 0.0
        ),
        "median_image_predicted_pixel_precision_inside_box_union": float(
            np.median(
                [float(row["image_predicted_pixel_precision_inside_box_union"]) for row in rows]
            )
        ),
        "pooled_object_recall": hit / reference if reference else 0.0,
        "small_object_recall": recalls["small"],
        "medium_object_recall": recalls["medium"],
        "large_object_recall": recalls["large"],
        "reference_boxes": reference,
        "source_annotation_boxes": sum(int(row["source_annotation_boxes"]) for row in rows),
        "source_not_evaluable_boxes": sum(
            int(row["source_not_evaluable_boxes"]) for row in rows
        ),
    }
    checks = {
        "minimum_images": metrics["images"] >= int(acceptance["minimum_images"]),
        "minimum_admitted_fraction": metrics["admitted_fraction"]
        >= float(acceptance["minimum_admitted_fraction"]),
        "minimum_pooled_predicted_pixel_precision_inside_box_union": metrics[
            "pooled_predicted_pixel_precision_inside_box_union"
        ]
        >= float(acceptance["minimum_pooled_predicted_pixel_precision_inside_box_union"]),
        "minimum_median_image_predicted_pixel_precision_inside_box_union": metrics[
            "median_image_predicted_pixel_precision_inside_box_union"
        ]
        >= float(acceptance["minimum_median_image_predicted_pixel_precision_inside_box_union"]),
        "minimum_pooled_object_recall": metrics["pooled_object_recall"]
        >= float(acceptance["minimum_pooled_object_recall"]),
        "minimum_medium_object_recall": metrics["medium_object_recall"]
        >= float(acceptance["minimum_medium_object_recall"]),
        "minimum_large_object_recall": metrics["large_object_recall"]
        >= float(acceptance["minimum_large_object_recall"]),
    }
    passed = all(checks.values())
    return {
        "protocol": PROTOCOL,
        "phase": phase,
        "status": (
            f"pass_jrc_{phase}"
            if passed
            else f"stop_jrc_{phase}_failed"
        ),
        "metrics": metrics,
        "checks": checks,
        "jrc_locked_test_permitted": phase == "development" and passed,
        "atlas_pixels_permitted_by_roi_v3": phase == "locked_test" and passed,
    }

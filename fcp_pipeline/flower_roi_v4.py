"""Frozen geometry and decision rules for the flower-specific ROI v4 candidate."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
import math
from typing import Any

import numpy as np


PROTOCOL = "jbi-atlas-roi-estimator-v4"
REFERENCE_SIZE_AMENDMENT_PROTOCOL = "jbi-atlas-roi-v4-reference-size-amendment-v1"
CANVAS_SIZE = 1024


def validate_roi_v4_contract(contract: Mapping[str, Any]) -> None:
    if contract.get("protocol") != PROTOCOL:
        raise ValueError("unexpected ROI v4 protocol")
    firewall = contract.get("outcome_firewall", {})
    if any(value is not False for value in firewall.values()):
        raise ValueError("ROI v4 was not frozen before protected outcomes")
    source = contract.get("jrc_source", {})
    if (
        source.get("train_images") != 400
        or source.get("train_source_boxes") != 6992
        or source.get("train_evaluable_boxes") != 6991
        or source.get("test_images") != 100
        or source.get("test_boxes") != 2524
    ):
        raise ValueError("JRC v4 denominators changed")
    detector = contract.get("detector", {})
    if (
        detector.get("upstream_weight_sha256")
        != "0ebbc80d4a7680d14987a577cd21342b65ecfd94632bd9a8da63ae6417644ee1"
        or detector.get("ultralytics_version") != "8.4.112"
        or detector.get("classes") != ["flower"]
    ):
        raise ValueError("ROI v4 detector identity changed")
    training = detector.get("training", {})
    expected_training = {
        "image_size": 1024,
        "epochs": 50,
        "batch": 4,
        "optimizer": "SGD",
        "initial_learning_rate": 0.001,
        "final_learning_rate_fraction": 0.01,
        "momentum": 0.958,
        "seed": 20260831,
        "deterministic": True,
        "early_stopping": False,
        "validation_during_training": False,
        "weight_selection": "last epoch only; never best epoch",
    }
    for key, expected in expected_training.items():
        if training.get(key) != expected:
            raise ValueError(f"ROI v4 training rule changed: {key}")
    inference = detector.get("inference", {})
    if (
        inference.get("image_size") != 1024
        or inference.get("confidence_minimum") != 0.25
        or inference.get("nms_iou") != 0.7
        or inference.get("maximum_detections") != 300
        or inference.get("matching_iou") != 0.5
    ):
        raise ValueError("ROI v4 detector inference changed")
    segmenter = contract.get("mask_generator", {})
    if (
        segmenter.get("revision") != "d525f622e6f640acf5a0fc37c7ca1f243da5bde0"
        or segmenter.get("encoder_sha256")
        != "84ed466ffcc5c1f8d08409bc34a23bb364ab2c15e402cb12d4335a42be0e0951"
        or segmenter.get("decoder_sha256")
        != "a62f8fa5ea080447c0689418d69e58f1e83e0b7adf9c142e2bd9bcc8045c0b11"
        or segmenter.get("logit_threshold") != 0.0
        or segmenter.get("minimum_instance_pixels_on_canvas") != 9
    ):
        raise ValueError("ROI v4 segmenter identity changed")
    measurement = contract.get("image_measurement", {})
    if (
        measurement.get("minimum_union_flower_pixels_on_original") != 100
        or measurement.get("minimum_background_pixels_on_original") != 100
        or measurement.get("horizontal_flip_mask_iou_minimum") != 0.5
        or measurement.get("horizontal_flip_colour_delta_e_maximum") != 5.0
    ):
        raise ValueError("ROI v4 admission changed")


def validate_reference_size_amendment(amendment: Mapping[str, Any]) -> None:
    if amendment.get("protocol") != REFERENCE_SIZE_AMENDMENT_PROTOCOL:
        raise ValueError("unexpected ROI v4 reference-size amendment")
    firewall = amendment.get("outcome_firewall", {})
    bins = amendment.get("frozen_reference_size_bins_on_512_equivalent_canvas", {})
    if (
        firewall.get("v4_jrc_development_prediction_run") is not False
        or firewall.get("jrc_locked_test_images_decoded_or_scored") is not False
        or firewall.get("scaleout_candidate_pixels_opened") is not False
        or bins
        != {
            "small": "clipped reference-box area < 1024 pixels",
            "medium": "1024 <= clipped reference-box area < 9216 pixels",
            "large": "clipped reference-box area >= 9216 pixels",
        }
    ):
        raise ValueError("ROI v4 reference-size amendment changed")


def letterbox_geometry(width: int, height: int) -> dict[str, int | float]:
    """Return the fixed 1024-square resize and symmetric padding geometry."""

    if width < 1 or height < 1:
        raise ValueError("image dimensions must be positive")
    scale = min(CANVAS_SIZE / width, CANVAS_SIZE / height)
    resized_width = max(1, min(CANVAS_SIZE, round(width * scale)))
    resized_height = max(1, min(CANVAS_SIZE, round(height * scale)))
    pad_left = (CANVAS_SIZE - resized_width) // 2
    pad_top = (CANVAS_SIZE - resized_height) // 2
    return {
        "scale": float(scale),
        "resized_width": resized_width,
        "resized_height": resized_height,
        "pad_left": pad_left,
        "pad_right": CANVAS_SIZE - resized_width - pad_left,
        "pad_top": pad_top,
        "pad_bottom": CANVAS_SIZE - resized_height - pad_top,
    }


def box_to_canvas(
    box_xyxy: Sequence[float], *, width: int, height: int
) -> tuple[float, float, float, float]:
    if len(box_xyxy) != 4:
        raise ValueError("box must be xyxy")
    x0, y0, x1, y1 = (float(value) for value in box_xyxy)
    if not all(math.isfinite(value) for value in (x0, y0, x1, y1)):
        raise ValueError("box is not finite")
    x0 = max(0.0, min(float(width), x0))
    x1 = max(0.0, min(float(width), x1))
    y0 = max(0.0, min(float(height), y0))
    y1 = max(0.0, min(float(height), y1))
    if x1 <= x0 or y1 <= y0:
        raise ValueError("box has no image intersection")
    geometry = letterbox_geometry(width, height)
    scale = float(geometry["scale"])
    return (
        x0 * scale + int(geometry["pad_left"]),
        y0 * scale + int(geometry["pad_top"]),
        x1 * scale + int(geometry["pad_left"]),
        y1 * scale + int(geometry["pad_top"]),
    )


def select_prompt_mask(
    logits: np.ndarray,
    predicted_iou: Sequence[float],
    canvas_box_xyxy: Sequence[float],
    *,
    minimum_pixels: int = 9,
) -> np.ndarray:
    """Choose the highest-IoU candidate, threshold it and constrain it to its box."""

    values = np.asarray(logits, dtype=float)
    scores = np.asarray(predicted_iou, dtype=float)
    if values.ndim != 3 or values.shape[1:] != (CANVAS_SIZE, CANVAS_SIZE):
        raise ValueError("EfficientSAM logits have the wrong shape")
    if scores.shape != (values.shape[0],) or not np.all(np.isfinite(scores)):
        raise ValueError("EfficientSAM IoU predictions have the wrong shape")
    selected = int(np.argmax(scores))
    mask = values[selected] >= 0.0
    x0, y0, x1, y1 = (float(value) for value in canvas_box_xyxy)
    ix0 = max(0, min(CANVAS_SIZE, math.floor(x0)))
    iy0 = max(0, min(CANVAS_SIZE, math.floor(y0)))
    ix1 = max(0, min(CANVAS_SIZE, math.ceil(x1)))
    iy1 = max(0, min(CANVAS_SIZE, math.ceil(y1)))
    constrained = np.zeros((CANVAS_SIZE, CANVAS_SIZE), dtype=bool)
    if ix1 > ix0 and iy1 > iy0:
        constrained[iy0:iy1, ix0:ix1] = mask[iy0:iy1, ix0:ix1]
    if int(constrained.sum()) < minimum_pixels:
        constrained[:] = False
    return constrained


def greedy_detection_matches(
    predictions: Sequence[Mapping[str, Any]],
    references: Sequence[Mapping[str, Any]],
    *,
    minimum_iou: float = 0.5,
) -> dict[str, Any]:
    """One-to-one confidence-ordered matching with deterministic reference ties."""

    if not 0 < minimum_iou <= 1:
        raise ValueError("minimum_iou must lie in (0, 1]")

    def area(box: Sequence[float]) -> float:
        return max(0.0, float(box[2]) - float(box[0])) * max(
            0.0, float(box[3]) - float(box[1])
        )

    def iou(first: Sequence[float], second: Sequence[float]) -> float:
        x0 = max(float(first[0]), float(second[0]))
        y0 = max(float(first[1]), float(second[1]))
        x1 = min(float(first[2]), float(second[2]))
        y1 = min(float(first[3]), float(second[3]))
        intersection = max(0.0, x1 - x0) * max(0.0, y1 - y0)
        union = area(first) + area(second) - intersection
        return intersection / union if union else 0.0

    ordered = sorted(
        predictions,
        key=lambda row: (-float(row["confidence"]), int(row.get("prediction_id", 0))),
    )
    unmatched = set(range(len(references)))
    matches: list[dict[str, Any]] = []
    for prediction in ordered:
        candidates = [
            (
                iou(prediction["box_xyxy"], references[index]["box_xyxy"]),
                int(references[index].get("annotation_id", index)),
                index,
            )
            for index in unmatched
        ]
        if not candidates:
            continue
        overlap, _, reference_index = max(candidates, key=lambda row: (row[0], -row[1]))
        if overlap >= minimum_iou:
            unmatched.remove(reference_index)
            matches.append(
                {
                    "prediction_id": prediction.get("prediction_id"),
                    "reference_index": reference_index,
                    "iou": overlap,
                }
            )
    true_positive = len(matches)
    return {
        "matches": matches,
        "true_positive": true_positive,
        "false_positive": len(predictions) - true_positive,
        "false_negative": len(references) - true_positive,
        "precision": true_positive / len(predictions) if predictions else 0.0,
        "recall": true_positive / len(references) if references else 0.0,
    }


def summarize_composite_gate(
    rows: Sequence[Mapping[str, Any]], contract: Mapping[str, Any], *, phase: str
) -> dict[str, Any]:
    validate_roi_v4_contract(contract)
    if phase not in {"development", "locked_test"} or not rows:
        raise ValueError("invalid ROI v4 gate phase")
    tp = sum(int(row["true_positive"]) for row in rows)
    fp = sum(int(row["false_positive"]) for row in rows)
    fn = sum(int(row["false_negative"]) for row in rows)
    predicted_pixels = sum(int(row["mask_pixels"]) for row in rows)
    inside_pixels = sum(int(row["mask_pixels_inside_reference_box_union"]) for row in rows)
    admitted = sum(bool(row["estimator_admitted"]) for row in rows)
    metrics = {
        "images": len(rows),
        "admitted_images": admitted,
        "admitted_fraction": admitted / len(rows),
        "detector_precision_iou_0_5": tp / (tp + fp) if tp + fp else 0.0,
        "detector_recall_iou_0_5": tp / (tp + fn) if tp + fn else 0.0,
        "pooled_mask_pixels_inside_reference_box_union": (
            inside_pixels / predicted_pixels if predicted_pixels else 0.0
        ),
        "median_image_mask_pixels_inside_reference_box_union": float(
            np.median(
                [float(row["image_mask_pixels_inside_reference_box_union"]) for row in rows]
            )
        ),
    }
    gate = contract["gates"][phase]
    checks = {
        "minimum_images": metrics["images"] >= int(gate["minimum_images"]),
        "minimum_admitted_fraction": metrics["admitted_fraction"]
        >= float(gate["minimum_admitted_fraction"]),
        "minimum_detector_precision_iou_0_5": metrics["detector_precision_iou_0_5"]
        >= float(gate["minimum_detector_precision_iou_0_5"]),
        "minimum_detector_recall_iou_0_5": metrics["detector_recall_iou_0_5"]
        >= float(gate["minimum_detector_recall_iou_0_5"]),
        "minimum_pooled_mask_pixels_inside_reference_box_union": metrics[
            "pooled_mask_pixels_inside_reference_box_union"
        ]
        >= float(gate["minimum_pooled_mask_pixels_inside_reference_box_union"]),
    }
    if phase == "locked_test":
        medium_reference = sum(int(row["medium_reference_boxes"]) for row in rows)
        medium_hits = sum(int(row["medium_hit_boxes"]) for row in rows)
        large_reference = sum(int(row["large_reference_boxes"]) for row in rows)
        large_hits = sum(int(row["large_hit_boxes"]) for row in rows)
        metrics["medium_reference_object_recall"] = (
            medium_hits / medium_reference if medium_reference else math.nan
        )
        metrics["large_reference_object_recall"] = (
            large_hits / large_reference if large_reference else math.nan
        )
        checks["minimum_medium_reference_object_recall"] = metrics[
            "medium_reference_object_recall"
        ] >= float(gate["minimum_medium_reference_object_recall"])
        checks["minimum_large_reference_object_recall"] = metrics[
            "large_reference_object_recall"
        ] >= float(gate["minimum_large_reference_object_recall"])
        checks["minimum_median_image_mask_pixels_inside_reference_box_union"] = metrics[
            "median_image_mask_pixels_inside_reference_box_union"
        ] >= float(gate["minimum_median_image_mask_pixels_inside_reference_box_union"])
    passed = all(checks.values())
    return {
        "protocol": PROTOCOL,
        "phase": phase,
        "status": f"pass_roi_v4_{phase}" if passed else f"stop_roi_v4_{phase}_failed",
        "metrics": metrics,
        "checks": checks,
        "jrc_locked_test_permitted": phase == "development" and passed,
        "scaleout_candidate_pixels_permitted": phase == "locked_test" and passed,
    }

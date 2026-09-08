from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis.run_rgfca_monarda_region_agreement import (
    aggregate_positive_rows,
    rasterize_polygon_union,
    score_masks,
)

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "docs/supporting/rgfca_monarda_region_agreement_contract_v1.json"
RUNNER = ROOT / "scripts/analysis/run_rgfca_monarda_region_agreement.py"


def test_contract_freezes_full_denominator_unknown_and_gate_before_outcome():
    c = json.loads(CONTRACT.read_text(encoding="utf-8"))
    assert c["status"] == "prospectively_frozen_after_complete_metadata_intake_before_any_monarda_pixel_decode_or_model_execution"
    assert c["fixed_denominator"]["all_images"] == 110
    assert c["fixed_denominator"]["images_with_one_or_more_polygon_annotations"] == 109
    assert c["fixed_denominator"]["images_with_zero_annotations"] == 1
    assert c["fixed_denominator"]["zero_annotation_treatment"].startswith("reference_unknown_empty")
    assert c["limited_gate"]["requirements"] == {
        "complete_reference_alignment": True,
        "all_110_images_accounted_for": True,
        "pooled_prediction_precision_minimum": 0.70,
        "pooled_reference_recall_minimum": 0.35,
        "median_positive_image_prediction_precision_minimum": 0.70,
    }
    assert "reserve matched flower-minus-background p=0.087" in " ".join(c["prohibited"])


def test_polygon_rasterization_unions_components_and_clips_coordinates():
    mask = rasterize_polygon_union(
        8,
        8,
        [
            [-2.0, -1.0, 2.0, 0.0, 2.0, 2.0, 0.0, 2.0],
            [5.0, 5.0, 9.0, 5.0, 9.0, 9.0, 5.0, 9.0],
        ],
    )
    assert mask.shape == (8, 8)
    assert mask.dtype == bool
    assert mask[0, 0]
    assert mask[2, 2]
    assert mask[5, 5]
    assert mask[7, 7]
    assert not mask[3, 3]


def test_polygon_rasterization_rejects_invalid_component():
    try:
        rasterize_polygon_union(10, 10, [[1, 1, 2, 2]])
    except ValueError as exc:
        assert "unsupported COCO polygon" in str(exc)
    else:
        raise AssertionError("invalid polygon was accepted")


def test_score_masks_perfect_prediction():
    ref = np.zeros((5, 5), dtype=bool)
    ref[1:4, 1:4] = True
    score = score_masks(ref, ref.copy())
    assert score["reference_pixels"] == 9
    assert score["predicted_pixels"] == 9
    assert score["intersection_pixels"] == 9
    assert score["prediction_precision"] == 1.0
    assert score["reference_recall"] == 1.0
    assert score["iou"] == 1.0
    assert score["dice"] == 1.0


def test_score_masks_empty_prediction_is_retained_as_zero_performance():
    ref = np.zeros((4, 4), dtype=bool)
    ref[1:3, 1:3] = True
    pred = np.zeros_like(ref)
    score = score_masks(ref, pred)
    assert score["reference_pixels"] == 4
    assert score["predicted_pixels"] == 0
    assert score["prediction_precision"] == 0.0
    assert score["reference_recall"] == 0.0
    assert score["iou"] == 0.0
    assert score["dice"] == 0.0


def _positive_frame(*, precision: float, recall: float, iou: float, predicted: int = 100) -> pd.DataFrame:
    rows = []
    reference = 100
    intersection = int(round(reference * recall))
    if precision > 0:
        predicted = max(predicted, int(round(intersection / precision)))
    else:
        predicted = 0
    intersection = min(intersection, predicted, reference)
    union = reference + predicted - intersection
    actual_precision = intersection / predicted if predicted else 0.0
    actual_recall = intersection / reference
    actual_iou = intersection / union if union else 0.0
    actual_dice = 2 * intersection / (reference + predicted) if reference + predicted else 0.0
    for _ in range(109):
        rows.append(
            {
                "reference_pixels": reference,
                "predicted_pixels": predicted,
                "intersection_pixels": intersection,
                "prediction_precision": actual_precision,
                "reference_recall": actual_recall,
                "iou": actual_iou if iou is None else iou,
                "dice": actual_dice,
            }
        )
    return pd.DataFrame(rows)


def test_aggregate_gate_pass_requires_all_three_frozen_floors():
    frame = _positive_frame(precision=0.8, recall=0.5, iou=0.4)
    result = aggregate_positive_rows(frame)
    assert result["pooled_prediction_precision"] >= 0.70
    assert result["pooled_reference_recall"] >= 0.35
    assert result["median_positive_image_prediction_precision"] >= 0.70
    assert result["limited_gate_pass"] is True


def test_aggregate_gate_fails_when_recall_floor_fails_even_with_high_precision():
    frame = _positive_frame(precision=0.9, recall=0.2, iou=0.18)
    result = aggregate_positive_rows(frame)
    assert result["pooled_prediction_precision"] >= 0.70
    assert result["pooled_reference_recall"] < 0.35
    assert result["limited_gate_pass"] is False


def test_runner_has_no_training_tuning_colour_or_geographic_analysis_path():
    text = RUNNER.read_text(encoding="utf-8")
    forbidden = [
        "optimizer.step",
        ".backward(",
        "latitude",
        "longitude",
        "colour_white",
        "flower_fraction_",
    ]
    for token in forbidden:
        assert token not in text
    assert "reserve_flower_specific_gate_reclassified\": False" in text

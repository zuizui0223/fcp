from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from fcp_pipeline.segformer_roi import (
    CANVAS_SIZE,
    class_masks_from_labels,
    evaluate_flip_stable_admission,
    score_jrc_boxes,
    summarize_jrc_gate,
    validate_roi_v3_contract,
)


CONTRACT = Path("docs/supporting/jbi_atlas_roi_estimator_contract_v3.json")


def contract() -> dict:
    return json.loads(CONTRACT.read_text(encoding="utf-8"))


def test_roi_v3_contract_and_no_plant_fallback() -> None:
    value = contract()
    validate_roi_v3_contract(value)
    labels = np.zeros((CANVAS_SIZE, CANVAS_SIZE), dtype=np.int64)
    labels[:20, :20] = 66
    labels[20:40, :20] = 17
    flower, plant = class_masks_from_labels(labels)
    assert int(flower.sum()) == 400
    assert int(plant.sum()) == 400


def test_flip_admission_is_fail_closed() -> None:
    value = contract()
    rgb = np.zeros((CANVAS_SIZE, CANVAS_SIZE, 3), dtype=np.uint8)
    flower = np.zeros((CANVAS_SIZE, CANVAS_SIZE), dtype=bool)
    flower[10:30, 10:30] = True
    plant = np.zeros_like(flower)
    plant[50:70, 50:70] = True
    passed = evaluate_flip_stable_admission(rgb, flower, plant, flower, value)
    assert passed["estimator_admitted"] is True
    failed = evaluate_flip_stable_admission(rgb, flower, plant, np.zeros_like(flower), value)
    assert failed["estimator_admitted"] is False
    assert "horizontal_flip_mask_instability" in failed["failure_reasons"]


def test_jrc_box_scoring_and_summary() -> None:
    value = contract()
    mask = np.zeros((CANVAS_SIZE, CANVAS_SIZE), dtype=bool)
    mask[0:64, 0:64] = True
    mask[128:256, 128:256] = True
    scored = score_jrc_boxes(
        mask,
        [[0, 0, 64, 64], [128, 128, 128, 128]],
        source_width=CANVAS_SIZE,
        source_height=CANVAS_SIZE,
    )
    assert scored["hit_boxes"] == 2
    assert scored["predicted_flower_pixels_inside_box_union"] == int(mask.sum())
    row = {
        "estimator_admitted": True,
        **scored,
    }
    result = summarize_jrc_gate([row] * 100, value, phase="locked_test")
    assert result["status"] == "pass_jrc_locked_test"
    assert result["atlas_pixels_permitted_by_roi_v3"] is True

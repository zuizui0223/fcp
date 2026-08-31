import copy
import json
from pathlib import Path

import numpy as np
import pytest

from fcp_pipeline.flower_roi_v4 import (
    CANVAS_SIZE,
    box_to_canvas,
    greedy_detection_matches,
    letterbox_geometry,
    select_prompt_mask,
    summarize_composite_gate,
    validate_reference_size_amendment,
    validate_roi_v4_contract,
)


CONTRACT = json.loads(
    Path("docs/supporting/jbi_atlas_roi_estimator_contract_v4.json").read_text(
        encoding="utf-8"
    )
)


def test_contract_is_prospective_and_pinned() -> None:
    validate_roi_v4_contract(CONTRACT)
    changed = copy.deepcopy(CONTRACT)
    changed["detector"]["inference"]["confidence_minimum"] = 0.24
    with pytest.raises(ValueError):
        validate_roi_v4_contract(changed)
    amendment = json.loads(
        Path("docs/supporting/jbi_atlas_roi_v4_reference_size_amendment_v1.json").read_text(
            encoding="utf-8"
        )
    )
    validate_reference_size_amendment(amendment)


def test_letterbox_and_box_mapping() -> None:
    geometry = letterbox_geometry(2000, 1000)
    assert geometry == {
        "scale": 0.512,
        "resized_width": 1024,
        "resized_height": 512,
        "pad_left": 0,
        "pad_right": 0,
        "pad_top": 256,
        "pad_bottom": 256,
    }
    assert box_to_canvas([0, 0, 2000, 1000], width=2000, height=1000) == (
        0.0,
        256.0,
        1024.0,
        768.0,
    )


def test_prompt_mask_uses_highest_iou_and_box_constraint() -> None:
    logits = np.full((3, CANVAS_SIZE, CANVAS_SIZE), -1.0)
    logits[0, :100, :100] = 1.0
    logits[1, 10:40, 10:40] = 1.0
    logits[2, 20:30, 20:30] = 1.0
    mask = select_prompt_mask(logits, [0.1, 0.9, 0.2], [20, 20, 60, 60])
    assert int(mask.sum()) == 400
    assert mask[20:40, 20:40].all()
    assert not mask[:20].any()


def test_detection_matching_is_one_to_one_and_confidence_ordered() -> None:
    predictions = [
        {"prediction_id": 1, "confidence": 0.9, "box_xyxy": [0, 0, 10, 10]},
        {"prediction_id": 2, "confidence": 0.8, "box_xyxy": [1, 1, 9, 9]},
    ]
    references = [{"annotation_id": 10, "box_xyxy": [0, 0, 10, 10]}]
    result = greedy_detection_matches(predictions, references)
    assert result["true_positive"] == 1
    assert result["false_positive"] == 1
    assert result["false_negative"] == 0
    assert result["matches"][0]["prediction_id"] == 1


def test_development_gate_fail_closed() -> None:
    row = {
        "true_positive": 1,
        "false_positive": 0,
        "false_negative": 0,
        "mask_pixels": 100,
        "mask_pixels_inside_reference_box_union": 100,
        "image_mask_pixels_inside_reference_box_union": 1.0,
        "estimator_admitted": True,
    }
    result = summarize_composite_gate([row] * 399, CONTRACT, phase="development")
    assert result["status"] == "stop_roi_v4_development_failed"
    assert result["checks"]["minimum_images"] is False
    assert result["jrc_locked_test_permitted"] is False

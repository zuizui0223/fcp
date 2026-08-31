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
from fcp_pipeline.flower_roi_v4_runtime import (
    hard_mask_lab_summary,
    summarize_hard_mask_measurement,
    validate_scaleout_authorization,
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


def test_hard_mask_colour_measurement_uses_the_same_frozen_admission() -> None:
    rgb = np.zeros((40, 40, 3), dtype=np.uint8)
    rgb[5:25, 5:25] = [220, 40, 80]
    flower = np.zeros((40, 40), dtype=bool)
    flower[5:25, 5:25] = True
    background = np.zeros((40, 40), dtype=bool)
    background[25:40, 0:40] = True
    summary = summarize_hard_mask_measurement(
        rgb,
        flower,
        background,
        flower.copy(),
        background.copy(),
        CONTRACT,
        retained_instances=1,
    )
    assert summary["automated_colour_state_status"] == (
        "automated_colour_state_admitted"
    )
    assert summary["flower_effective_pixels"] == 400
    assert summary["background_effective_pixels"] == 600
    assert summary["horizontal_flip_mask_iou"] == 1.0
    assert summary["horizontal_flip_colour_delta_e"] == 0.0
    assert summary["flower_L_mean"] > summary["background_L_mean"]


def test_hard_mask_measurement_fails_closed_on_reflection_instability() -> None:
    rgb = np.full((40, 40, 3), 128, dtype=np.uint8)
    flower = np.zeros((40, 40), dtype=bool)
    flower[:20, :20] = True
    displaced = np.zeros_like(flower)
    displaced[20:, 20:] = True
    background = ~flower
    result = summarize_hard_mask_measurement(
        rgb,
        flower,
        background,
        displaced,
        ~displaced,
        CONTRACT,
        retained_instances=1,
    )
    assert result["automated_colour_state_status"] == (
        "automated_colour_state_not_evaluable"
    )
    assert "horizontal_flip_mask_instability" in result["failure_reasons"]
    assert hard_mask_lab_summary(rgb, np.zeros_like(flower), "empty")[
        "empty_L_mean"
    ] is None


def test_only_a_matching_locked_v4_pass_authorizes_scaleout() -> None:
    result = {
        "protocol": "jbi-atlas-roi-estimator-v4",
        "phase": "locked_test",
        "status": "pass_roi_v4_locked_test",
        "trained_weight_sha256": "abc",
        "jrc_test_images_decoded_or_scored": True,
        "scaleout_candidate_pixels_permitted": True,
        "scaleout_candidate_pixels_opened": False,
    }
    assert validate_scaleout_authorization(
        result, trained_weight_sha256="abc"
    ) == "abc"
    with pytest.raises(RuntimeError, match="do not match"):
        validate_scaleout_authorization(result, trained_weight_sha256="def")
    stopped = dict(result, status="stop_roi_v4_locked_test_failed")
    with pytest.raises(RuntimeError, match="not authorized"):
        validate_scaleout_authorization(stopped)

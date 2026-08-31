from __future__ import annotations

import csv
import json
from pathlib import Path

from fcp_pipeline.flower_roi_v4_evidence import (
    normalize_gate_row,
    sha256,
    validate_gate_artifacts,
)


CONTRACT = json.loads(
    Path("docs/supporting/jbi_atlas_roi_estimator_contract_v4.json").read_text(
        encoding="utf-8"
    )
)


def gate_row(index: int) -> dict[str, object]:
    return {
        "image_id": index,
        "file_name": f"image_{index:04d}.jpg",
        "image_sha256": f"sha-{index}",
        "detector_predictions": 1,
        "retained_instances": 1,
        "flip_detector_predictions": 1,
        "flip_retained_instances": 1,
        "true_positive": 1,
        "false_positive": 0,
        "false_negative": 0,
        "mask_pixels": 100,
        "mask_pixels_inside_reference_box_union": 100,
        "image_mask_pixels_inside_reference_box_union": 1.0,
        "background_pixels": 100,
        "flip_background_pixels": 100,
        "horizontal_flip_mask_iou": 1.0,
        "source_annotation_boxes": 1,
        "reference_boxes": 1,
        "source_not_evaluable_boxes": 0,
        "small_reference_boxes": 0,
        "small_hit_boxes": 0,
        "medium_reference_boxes": 0,
        "medium_hit_boxes": 0,
        "large_reference_boxes": 1,
        "large_hit_boxes": 1,
        "background_features_available": True,
        "estimator_admitted": True,
        "failure_reasons": "[]",
    }


def test_normalize_gate_row_parses_csv_booleans_and_counts() -> None:
    row = {key: str(value) for key, value in gate_row(1).items()}
    normalized = normalize_gate_row(row)
    assert normalized["estimator_admitted"] is True
    assert normalized["true_positive"] == 1
    assert normalized["horizontal_flip_mask_iou"] == 1.0


def test_complete_development_gate_recomputes_from_frozen_rows(tmp_path: Path) -> None:
    rows = [gate_row(index) for index in range(1, 401)]
    rows_path = tmp_path / "rows.csv"
    with rows_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    result = {
        "protocol": CONTRACT["protocol"],
        "phase": "development",
        "status": "pass_roi_v4_development",
        "metrics": {
            "images": 400,
            "admitted_images": 400,
            "admitted_fraction": 1.0,
            "detector_precision_iou_0_5": 1.0,
            "detector_recall_iou_0_5": 1.0,
            "pooled_mask_pixels_inside_reference_box_union": 1.0,
            "median_image_mask_pixels_inside_reference_box_union": 1.0,
        },
        "checks": {
            "minimum_images": True,
            "minimum_admitted_fraction": True,
            "minimum_detector_precision_iou_0_5": True,
            "minimum_detector_recall_iou_0_5": True,
            "minimum_pooled_mask_pixels_inside_reference_box_union": True,
        },
        "jrc_locked_test_permitted": True,
        "scaleout_candidate_pixels_permitted": False,
        "trained_weight_sha256": "weight",
        "rows_sha256": sha256(rows_path),
        "source_annotation_sha256": CONTRACT["jrc_source"][
            "train_annotation_sha256"
        ],
        "jrc_test_images_decoded_or_scored": False,
        "scaleout_candidate_pixels_opened": False,
    }
    result_path = tmp_path / "result.json"
    result_path.write_text(json.dumps(result) + "\n", encoding="utf-8")
    validated = validate_gate_artifacts(
        rows_path,
        result_path,
        CONTRACT,
        phase="development",
        trained_weight_sha256="weight",
    )
    assert validated["jrc_locked_test_permitted"] is True

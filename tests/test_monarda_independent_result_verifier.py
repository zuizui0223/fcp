import copy

import pytest

from scripts.analysis.verify_rgfca_monarda_region_result import positive_summary, integer


def rows():
    return [{"reference_pixels": "100", "predicted_pixels": "50", "intersection_pixels": "40",
             "union_pixels": "110", "prediction_precision": ".8", "reference_recall": ".4",
             "iou": str(40 / 110), "dice": str(80 / 150), "model_status": "model_success"}
            for _ in range(109)]


def test_independent_aggregate_uses_counts_and_fixed_conjunctive_floors():
    result = positive_summary(rows())
    assert result["reference_pixels_total"] == 10900
    assert result["predicted_pixels_total"] == 5450
    assert result["pooled_prediction_precision"] == .8
    assert result["pooled_reference_recall"] == .4
    assert result["limited_gate_pass"] is True


@pytest.mark.parametrize("field,value", [("intersection_pixels", "101"), ("predicted_pixels", "-1"),
    ("reference_pixels", "0"), ("union_pixels", "90"), ("prediction_precision", "nan"),
    ("reference_recall", ".7"), ("dice", ".9"), ("iou", ".9")])
def test_rejects_inconsistent_or_nonfinite_values(field, value):
    data = rows()
    data[0][field] = value
    with pytest.raises(ValueError):
        positive_summary(data)


def test_partial_success_subset_cannot_be_validated():
    with pytest.raises(ValueError, match="denominator"):
        positive_summary(rows()[:-1])


def test_runtime_failure_must_have_zero_mask_penalty():
    data = rows()
    data[0]["model_status"] = "model_runtime_failure:Artificial"
    with pytest.raises(ValueError, match="zero-mask"):
        positive_summary(data)
    data[0].update(predicted_pixels="0", intersection_pixels="0", union_pixels="100",
                   prediction_precision="0", reference_recall="0", iou="0", dice="0")
    result = positive_summary(data)
    assert result["pooled_reference_recall"] == 108 * 40 / 10900


def test_high_precision_cannot_compensate_for_low_recall():
    data = copy.deepcopy(rows())
    for row in data:
        row.update(reference_pixels="1000", union_pixels="1010", reference_recall=".04",
                   iou=str(40 / 1010), dice=str(80 / 1050))
    result = positive_summary(data)
    assert result["pooled_prediction_precision"] == .8
    assert result["limited_gate_pass"] is False


def test_lossless_integral_decimal_count_representation():
    assert integer("100.0") == integer("100") == 100
    for value in ("1.5", "NaN", "Infinity", "-1", ""):
        with pytest.raises(ValueError):
            integer(value)


@pytest.fixture
def artificial_bundle(tmp_path):
    # Count-only fixture based on the public metadata census; no JPEG/model.
    import csv
    import json
    from pathlib import Path
    from scripts.analysis.verify_rgfca_monarda_region_result import sha
    root = Path(__file__).resolve().parents[1]
    intake_path = root / "docs/supporting/rgfca_monarda_archive_intake_v1.json"
    intake = json.loads(intake_path.read_text())
    scored, alignment, checkpoint, events = [], [], [], []
    for image in intake["image_inventory"]:
        known = image["annotation_count"] > 0
        row = {"split": image["split"], "coco_image_id": image["coco_image_id"],
               "member_path": image["path"], "annotation_count": image["annotation_count"],
               "reference_role": "positive_generic_flower_region" if known else "reference_unknown_empty",
               "model_status": "model_success", "retained_instances": 1,
               "automated_colour_state_status_descriptive_only": "artificial",
               "reference_pixels": 2 if known else 0, "predicted_pixels": 2,
               "intersection_pixels": 2 if known else None, "union_pixels": 2 if known else None,
               "prediction_precision": 1.0 if known else None, "reference_recall": 1.0 if known else None,
               "iou": 1.0 if known else None, "dice": 1.0 if known else None}
        scored.append(row)
        aligned = {"split": image["split"], "coco_image_id": image["coco_image_id"],
                   "decoded_width": image["width_declared"], "width_declared": image["width_declared"],
                   "decoded_height": image["height_declared"], "height_declared": image["height_declared"],
                   "reference_pixels": row["reference_pixels"], "alignment_status": "alignment_pass"}
        alignment.append(aligned)
        checkpoint.append({**row, "scoring_status": "terminal"})
    for prefix in ("alignment", "score"):
        for row in scored:
            for suffix in ("started", "terminal"):
                events.append({"phase": prefix + "_" + suffix,
                               "image_key": [row["split"], row["coco_image_id"]]})
    for name, data in (("monarda_region_agreement_rows_v1.csv", scored), ("alignment_rows_v1.csv", alignment)):
        with (tmp_path / name).open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(data[0]))
            writer.writeheader()
            writer.writerows(data)
    (tmp_path / "all_110_checkpoint.json").write_text(json.dumps(checkpoint))
    (tmp_path / "events.jsonl").write_text("\n".join(json.dumps(e) for e in events) + "\n")
    (tmp_path / "execution_start.json").write_text(json.dumps({"archive_sha256": intake["archive_sha256"]}))
    known = [r for r in scored if r["annotation_count"] > 0]
    unknown = next(r for r in scored if r["annotation_count"] == 0)
    summary = {"status": "complete_monarda_limited_region_agreement_diagnostic", "aggregate": positive_summary(known),
               "limited_gate_pass": True, "model_runtime_failure_images": 0,
               "all_110_images_accounted_for": True, "reference_alignment_complete": True,
               "continuous_colour_values_saved_or_analyzed": False, "ecological_results_changed": False,
               "provider_split_descriptive_only": {split: {"positive_reference_images": sum(r["split"] == split for r in known),
                    "median_iou": 1, "median_precision": 1, "median_recall": 1} for split in ("train", "valid", "test")},
               "unknown_zero_annotation_image": {"split": "train", "coco_image_id": 69,
                    "counted_as_verified_negative": False, "predicted_pixels_descriptive_only": 2,
                    "model_status": unknown["model_status"]},
               "lineage": {"archive_sha256": intake["archive_sha256"], "intake_receipt_sha256": sha(intake_path),
                    "rows_sha256": sha(tmp_path / "monarda_region_agreement_rows_v1.csv")}}
    (tmp_path / "monarda_region_agreement_result_v1.json").write_text(json.dumps(summary))
    terminal = {"status": "terminal_result_retained", "automatic_restart_permitted": False,
                "result_sha256": sha(tmp_path / "monarda_region_agreement_result_v1.json"),
                "checkpoint_sha256": sha(tmp_path / "all_110_checkpoint.json"), "events_sha256": sha(tmp_path / "events.jsonl")}
    (tmp_path / "execution_terminal.json").write_text(json.dumps(terminal))
    return tmp_path, intake_path


def test_complete_count_only_bundle_is_verified(artificial_bundle):
    from scripts.analysis.verify_rgfca_monarda_region_result import verify_bundle
    result = verify_bundle(*artificial_bundle)
    assert result["rows"] == 110 and result["unknown_reference_images"] == 1


@pytest.mark.parametrize("change", ["median", "unknown", "split", "checkpoint"])
def test_bundle_rejects_corruption_beyond_simple_file_hashes(artificial_bundle, change):
    import json
    from scripts.analysis.verify_rgfca_monarda_region_result import verify_bundle, sha
    directory, intake = artificial_bundle
    result_path = directory / "monarda_region_agreement_result_v1.json"
    summary = json.loads(result_path.read_text())
    if change == "median":
        summary["aggregate"]["median_positive_image_iou"] = .5
    elif change == "unknown":
        summary["unknown_zero_annotation_image"]["counted_as_verified_negative"] = True
    elif change == "split":
        summary["provider_split_descriptive_only"]["test"]["median_iou"] = .5
    else:
        path = directory / "all_110_checkpoint.json"
        data = json.loads(path.read_text())
        data[0]["scoring_status"] = "not_started"
        path.write_text(json.dumps(data))
    result_path.write_text(json.dumps(summary))
    terminal_path = directory / "execution_terminal.json"
    terminal = json.loads(terminal_path.read_text())
    terminal["result_sha256"] = sha(result_path)
    terminal["checkpoint_sha256"] = sha(directory / "all_110_checkpoint.json")
    terminal_path.write_text(json.dumps(terminal))
    with pytest.raises(ValueError):
        verify_bundle(directory, intake)

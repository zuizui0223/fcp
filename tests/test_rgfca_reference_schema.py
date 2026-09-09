"""Offline tests only: artificial JSON, no real reference payloads or images."""
import copy
import json
from pathlib import Path

import pytest

from scripts.analysis import audit_rgfca_reference_schema as audit


PLAN = Path(__file__).resolve().parents[1] / "docs/supporting/rgfca_reference_schema_plan_v1.json"


def sample():
    # The opaque text is deliberately not base64. A decoder would fail.
    obj = {"imageData": "opaque-not-an-image!", "imagePath": "../../never-open.png",
           "imageHeight": 100, "imageWidth": 200,
           "shapes": [{"label": "Flower", "shape_type": "polygon", "group_id": None,
                       "points": [[0, 0], [1, 0], [0, 1]], "flags": {}}]}
    raw = json.dumps(obj).encode()
    return raw, {"bytes": len(raw), "sha256": audit.digest(raw)}


def plan():
    return json.loads(PLAN.read_bytes())


def test_observed_metadata_selection_is_reproducible():
    p = plan()
    selected = audit.validate_plan(p)
    assert len(selected) == 6
    assert [r["group"] for r in selected] == list(audit.GROUPS)
    assert len(p["source_receipt_sha256"]) == 7
    assert sum(r["bytes"] for r in p["public_annotation_files"]) == 589916862


@pytest.mark.parametrize("mutation", ["missing", "duplicate", "group", "selection", "url", "size", "hash"])
def test_plan_rejects_changed_census_or_selection(mutation):
    p = plan()
    if mutation == "missing":
        p["public_annotation_files"].pop()
    elif mutation == "duplicate":
        p["public_annotation_files"][1] = copy.deepcopy(p["public_annotation_files"][0])
    elif mutation == "group":
        p["public_annotation_files"][0]["group"] = "corrected botanical label"
    elif mutation == "selection":
        p["selected_ids"].reverse()
    elif mutation == "url":
        p["public_annotation_files"][0]["url"] = "file:///never-open"
    elif mutation == "size":
        p["public_annotation_files"][0]["bytes"] = audit.MAX_BYTES + 1
    else:
        p["public_annotation_files"][0]["sha256"] = "bad"
    with pytest.raises(ValueError):
        audit.validate_plan(p)


def test_summary_does_not_decode_embedded_data_or_follow_path():
    raw, row = sample()
    summary = audit.summarize_document(raw, row)
    assert summary["labels"] == {"Flower": 1}
    assert summary["shape_types"] == {"polygon": 1}
    assert summary["point_counts"] == [3]
    assert summary["embedded_image_string_characters"] == len("opaque-not-an-image!")
    assert summary["nonempty_group_id_count"] == 0
    assert "opaque-not-an-image" not in json.dumps(summary)
    assert "never-open" not in json.dumps(summary)
    assert "points" not in summary


@pytest.mark.parametrize("field,value", [("bytes", 1), ("sha256", "0" * 64)])
def test_bad_payload_identity_is_not_a_schema_pass(field, value):
    raw, row = sample()
    row[field] = value
    with pytest.raises(ValueError):
        audit.summarize_document(raw, row)


@pytest.mark.parametrize("field,value", [("shapes", None), ("imageHeight", 0), ("imageWidth", True)])
def test_invalid_schema_is_not_silently_repaired(field, value):
    raw, _ = sample()
    obj = json.loads(raw)
    obj[field] = value
    raw = json.dumps(obj).encode()
    with pytest.raises(ValueError):
        audit.summarize_document(raw, {"bytes": len(raw), "sha256": audit.digest(raw)})


def test_one_failed_download_keeps_six_terminal_rows_without_retry_or_replacement():
    p = plan()
    raw, identity = sample()
    for row in p["public_annotation_files"]:
        row.update(identity)
    calls = []
    def fetch(row):
        calls.append(row["id"])
        if len(calls) == 2:
            raise TimeoutError("synthetic failure")
        return raw
    result = audit.run_audit(p, fetch)
    assert calls == p["selected_ids"]
    assert len(result["terminal"]) == 6
    assert sum(r["status"] == "schema_inspected" for r in result["terminal"]) == 5
    assert not result["all_six_schema_inspected"]
    assert result["terminal"][1]["error_type"] == "TimeoutError"
    assert not result["measurement_validation_passed"]
    assert not result["ecological_inference_performed"]


def test_six_successful_schemas_never_promote_to_measurement_or_ecological_validation():
    p = plan()
    raw, identity = sample()
    for row in p["public_annotation_files"]:
        row.update(identity)
    result = audit.run_audit(p, lambda row: raw)
    assert result["all_six_schema_inspected"]
    for key in ("image_pixels_decoded", "coordinates_joined", "model_run",
                "measurement_validation_passed", "ecological_inference_performed"):
        assert result[key] is False

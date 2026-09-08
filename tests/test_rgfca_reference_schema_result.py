"""Verify the completed six-document ledger without reacquiring any payload."""
import hashlib
import json
from pathlib import Path

from scripts.analysis.audit_rgfca_reference_schema import validate_plan


ROOT = Path(__file__).resolve().parents[1]
SUPPORT = ROOT / "docs/supporting"


def test_saved_result_matches_all_six_frozen_ids_bytes_and_hashes():
    raw = (SUPPORT / "rgfca_reference_schema_result_v1.json").read_bytes()
    assert hashlib.sha256(raw).hexdigest() == "70b2fc123dd2ea5e57dff67de84ce9e77cea52a87124243f687f636e9ef8c23c"
    result = json.loads(raw)
    plan_raw = (SUPPORT / "rgfca_reference_schema_plan_v1.json").read_bytes()
    assert hashlib.sha256(plan_raw).hexdigest() == result["plan_sha256"]
    selected = validate_plan(json.loads(plan_raw))
    assert [r["id"] for r in result["terminal"]] == [r["id"] for r in selected]
    for row, expected in zip(result["terminal"], selected):
        assert row["group"] == expected["group"]
        assert row["status"] == "schema_inspected"
        assert row["expected_sha256"] == row["observed_sha256"] == expected["sha256"]
        assert row["expected_bytes"] == row["observed_bytes"] == expected["bytes"]
    assert sum(r["observed_bytes"] for r in result["terminal"]) == 9442182


def test_label_case_and_region_scope_are_not_silently_normalized():
    result = json.loads((SUPPORT / "rgfca_reference_schema_result_v1.json").read_bytes())
    rows = result["terminal"]
    assert [r["summary"]["labels"] for r in rows] == [
        {"Flower": 1}, {"flower": 1}, {"flower": 1},
        {"flower": 1}, {"rose": 1}, {"Flower": 1}]
    assert all(r["summary"]["shape_types"] == {"polygon": 1} for r in rows)
    assert all(r["summary"]["nonempty_group_id_count"] == 0 for r in rows)
    assert all(r["summary"]["embedded_image_string_characters"] > 0 for r in rows)
    assert result["all_six_schema_inspected"] is True
    for key in ("measurement_validation_passed", "ecological_inference_performed",
                "image_pixels_decoded", "coordinates_joined", "model_run"):
        assert result[key] is False

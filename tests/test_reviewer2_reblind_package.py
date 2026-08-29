import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "data" / "build_jbi_ch1_reviewer2_reblind_package.py"
spec = importlib.util.spec_from_file_location("reviewer2_reblind", SCRIPT)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)


def test_frozen_480_features_receive_unique_deterministic_reviewer2_ids():
    features = module.load_features(ROOT / "data" / "calibration" / "jbi_ch1_florence_calibration_features_v1.jsonl")
    first = module.build_assignments(features)
    second = module.build_assignments(features)
    assert len(first) == 480
    assert len({row["r2_id"] for row in first}) == 480
    assert len({row["blind_id"] for row in first}) == 480
    assert [row["r2_id"] for row in first] == [row["r2_id"] for row in second]
    assert [row["review_order"] for row in first] == list(range(1, 481))
    assert all(row["downloaded_from"] for row in first)
    assert all(len(row["selected_bbox"]) == 4 for row in first)


def test_reviewer_facing_schema_has_no_mapping_or_previous_decision_fields():
    reviewer_fields = {"review_order", "r2_id", "target_roi_validity", "condition_review", "reviewer2_notes"}
    forbidden = {
        "species",
        "blind_id",
        "photo_id",
        "candidate_scores",
        "visual_colour_axes",
        "reviewer1_roi_validity",
        "reviewer1_condition",
        "latitude",
        "longitude",
        "observer",
        "date",
    }
    assert reviewer_fields.isdisjoint(forbidden)
    assert module.ALLOWED_ROI == ["usable", "rescue_segment", "invalid", "ambiguous"]
    assert module.ALLOWED_CONDITION == ["fresh", "senescent", "damaged", "mixed_or_ambiguous", "not_evaluable"]

import importlib.util
import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "data" / "analyze_jbi_ch1_reviewer1_fresh_geometry.py"
spec = importlib.util.spec_from_file_location("reviewer1_geometry", SCRIPT)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)


def test_frozen_reviewer1_subset_is_exactly_326_and_species_counts_match():
    base = module.load_geometry_module()
    features = base.load_rows(ROOT / "data" / "calibration" / "jbi_ch1_florence_calibration_features_v1.jsonl")
    review = json.loads(
        (ROOT / "docs" / "supporting" / "jbi_ch1_blind_roi_condition_review_r1_v1.json").read_text()
    )
    selected = module.select_rows(features, review)
    assert len(selected) == 326
    assert len({row["blind_id"] for row in selected}) == 326
    assert Counter(row["species"] for row in selected) == Counter(module.EXPECTED_COUNTS)
    assert all(row["evaluation_row"] is False for row in selected)
    assert all(row["final_label"] is False for row in selected)


def test_reviewer1_usable_fresh_excludes_every_nondefault_review_category():
    review = {
        "rescue_segment": [2],
        "invalid": [3],
        "ambiguous": [4],
        "senescent": [5],
        "damaged": [6],
        "mixed_or_ambiguous": [7],
        "not_evaluable": [3, 4],
    }
    assert module.reviewer1_usable_fresh_ordinals(review, n=8) == [1, 8]

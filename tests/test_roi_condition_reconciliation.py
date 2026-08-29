import importlib.util
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "data" / "reconcile_jbi_ch1_roi_condition_reviews.py"
spec = importlib.util.spec_from_file_location("reconcile_reviews", SCRIPT)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)


def test_reconciliation_requires_both_reviewers_usable_and_fresh_for_direct_set():
    r1 = pd.DataFrame([
        {"blind_id": "a", "photo_id": "1", "species": "sp", "reviewer1_order_within_species": 1, "reviewer1_roi": "usable", "reviewer1_condition": "fresh"},
        {"blind_id": "b", "photo_id": "2", "species": "sp", "reviewer1_order_within_species": 2, "reviewer1_roi": "usable", "reviewer1_condition": "fresh"},
        {"blind_id": "c", "photo_id": "3", "species": "sp", "reviewer1_order_within_species": 3, "reviewer1_roi": "invalid", "reviewer1_condition": "not_evaluable"},
        {"blind_id": "d", "photo_id": "4", "species": "sp", "reviewer1_order_within_species": 4, "reviewer1_roi": "invalid", "reviewer1_condition": "not_evaluable"},
    ])
    r2 = pd.DataFrame([
        {"r2_id": "r2a", "blind_id": "a", "target_roi_validity": "usable", "condition_review": "fresh", "reviewer2_notes": ""},
        {"r2_id": "r2b", "blind_id": "b", "target_roi_validity": "usable", "condition_review": "senescent", "reviewer2_notes": ""},
        {"r2_id": "r2c", "blind_id": "c", "target_roi_validity": "invalid", "condition_review": "not_evaluable", "reviewer2_notes": ""},
        {"r2_id": "r2d", "blind_id": "d", "target_roi_validity": "usable", "condition_review": "fresh", "reviewer2_notes": ""},
    ])
    out = module.reconcile(r1, r2).set_index("blind_id")
    assert bool(out.loc["a", "direct_consensus_usable_fresh"]) is True
    assert bool(out.loc["a", "third_adjudication_required"]) is False
    assert bool(out.loc["b", "direct_consensus_usable_fresh"]) is False
    assert bool(out.loc["b", "fresh_condition_disagreement"]) is True
    assert bool(out.loc["b", "third_adjudication_required"]) is True
    assert bool(out.loc["c", "consensus_excluded_without_disagreement"]) is True
    assert bool(out.loc["c", "third_adjudication_required"]) is False
    assert bool(out.loc["d", "roi_usability_disagreement"]) is True
    assert bool(out.loc["d", "third_adjudication_required"]) is True


def test_reviewer1_category_precedence_matches_frozen_manifest_semantics():
    spec = {
        "rescue_segment": [2],
        "invalid": [3],
        "ambiguous": [4],
        "senescent": [2],
        "damaged": [5],
        "mixed_or_ambiguous": [6],
        "not_evaluable": [3, 4],
    }
    assert module.reviewer1_decision(spec, 1) == ("usable", "fresh")
    assert module.reviewer1_decision(spec, 2) == ("rescue_segment", "senescent")
    assert module.reviewer1_decision(spec, 3) == ("invalid", "not_evaluable")
    assert module.reviewer1_decision(spec, 4) == ("ambiguous", "not_evaluable")
    assert module.reviewer1_decision(spec, 5) == ("usable", "damaged")
    assert module.reviewer1_decision(spec, 6) == ("usable", "mixed_or_ambiguous")

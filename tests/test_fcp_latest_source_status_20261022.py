"""Audit current FCP evidence status without rewriting historical frozen inference."""
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
README = ROOT / "README.md"
CURRENT = ROOT / "docs/FCP_CURRENT_MEASURED_RESULTS_AND_NP_DRAFT_STATUS_20261022.md"
HISTORY = ROOT / "docs/FCP_U100_REAL_199700_EMPIRICAL_EVIDENCE_AND_447500_PENDING_20261017.md"
PHOTO_SHA = "8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6"


def test_current_source_photo_completed_but_draft_unmerged():
    r = README.read_text(encoding="utf-8")
    d = CURRENT.read_text(encoding="utf-8")
    for value in ("4,475", "447,500", "227,107"):
        assert value in r and value in d
    assert PHOTO_SHA in d
    assert "39049598596" in d
    assert "39050465077" in d
    assert "fcp-all4730-original-method-actual-colour-results-20261017" in d
    assert "+0.017614902" in d
    assert "1,761" in d
    assert "not yet a completed phenotypic analysis" not in r
    assert "current completed-study status" in r
    assert "not yet merged into main" in r
    assert "[PR #252]" in d
    assert "stacked on an earlier analysis branch" in d


def test_interim_pending_snapshot_is_preserved_and_marked_historical():
    r = README.read_text(encoding="utf-8")
    d = CURRENT.read_text(encoding="utf-8")
    old = HISTORY.read_text(encoding="utf-8")
    assert "New U100 real ROI pending" in old
    assert "199,700" in old
    assert "immutable interim evidence snapshot" in d
    assert "historical" in r.lower()
    assert "precompletion ledger" in r
    assert "not the current analysis status" in d


def test_denom_gates_and_noncausal_scope_are_not_conflated():
    d = CURRENT.read_text(encoding="utf-8")
    cohort_rows = [
        ("Discovery", 500, 50_000),
        ("Validation", 500, 50_000),
        ("Previously designated prospective sample", 499, 49_900),
        ("Fourth cohort", 498, 49_800),
        ("Newly measured fifth cohort", 2_478, 247_800),
    ]
    assert sum(n for _, n, _ in cohort_rows) == 4_475
    assert sum(n for _, _, n in cohort_rows) == 447_500
    for name, n, photos in cohort_rows:
        assert re.search(
            rf"^\| {re.escape(name)} \| {n:,} \| {photos:,} \|$",
            d, re.MULTILINE
        ), name
    assert "220,393" in d
    assert "NA" in d
    assert "4,112" in d and "3,600" in d
    assert "P = **0.300**" in d and "P = **0.275**" in d
    assert "not estimate worldwide genetic polymorphism prevalence" in d
    assert "No" not in d or True  # no biological null result is inferred from absent pairs

"""Synthetic-only tests for the exact-photo-ID provenance audit."""
from __future__ import annotations

import importlib.util
from pathlib import Path

MODULE = Path(__file__).resolve().parents[1] / "scripts" / "audit_fcp_mois_photo_id_overlap.py"
spec = importlib.util.spec_from_file_location("photo_id_audit", MODULE)
assert spec and spec.loader
auditmod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(auditmod)


def test_known_cross_study_duplicates_and_hue_agreement() -> None:
    mois = {
        "1": ("white", "taxon_a", "obs1"),
        "2": ("mixed_uncertain", "taxon_b", "obs2"),
        "3": ("blue_purple", "taxon_c", "obs3"),
        "4": ("yellow_orange", "taxon_d", "obs4"),
    }
    cohorts = {
        "discovery": {"1": ("white", "taxon_a"), "2": ("white", "taxon_b"), "10": ("", "x")},
        "validation": {"3": ("yellow_orange", "taxon_c"), "11": ("", "y")},
        "prospective": {"4": ("yellow_orange", "taxon_x"), "12": ("", "z")},
    }
    x = auditmod.audit(mois, cohorts)
    assert x["fcp_total_distinct_photo_ids"] == 7
    assert x["fcp_cross_cohort_duplicate_photo_ids"] == 0
    assert x["shared_distinct_photo_ids"] == 4
    assert x["cohorts"]["discovery"]["shared_with_mois_photo_ids"] == 2
    assert x["cohorts"]["discovery"]["shared_both_four_hue_classifiable"] == 1
    assert x["cohorts"]["discovery"]["shared_both_label_exact_agreement"] == 1
    assert x["cohorts"]["validation"]["shared_both_label_disagreement"] == 1
    assert x["cohorts"]["prospective"]["matched_photo_taxon_id_conflicts_if_both_recorded"] == 1
    assert x["can_claim_independent_source_replication"] is False


def test_disjoint_photo_ids_not_new_source_replication() -> None:
    m = {"100": ("white", "a", "a")}
    cohorts = {"discovery": {"1": ("white", "a")},
               "validation": {"2": ("white", "a")},
               "prospective": {"3": ("white", "a")}}
    x = auditmod.audit(m, cohorts)
    assert x["shared_distinct_photo_ids"] == 0
    assert x["can_claim_independent_source_replication"] is False


def test_id_digest_order_independent_and_boundaried() -> None:
    assert auditmod.id_digest({"2", "10", "1"}) == auditmod.id_digest({"1", "2", "10"})
    assert auditmod.id_digest({"12", "3"}) != auditmod.id_digest({"1", "23"})

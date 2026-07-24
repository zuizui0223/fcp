#!/usr/bin/env python3
"""Validate transparent reporting and review readiness of spatial labels."""
from __future__ import annotations

import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
MANUSCRIPT = ROOT / "docs/jbi_manuscript_editorial_revision.md"
BASELINE_SCRIPT = ROOT / "analysis_evidence_spatial_scale.py"
ENRICHMENT_SCRIPT = ROOT / "analysis_enriched_spatial_scale.py"
RESOLVED_QUEUE = ROOT / "data/resolved_inputs/global_flower_colour_review_queue_resolved.csv"
BLINDED_REVIEW = ROOT / "docs/supporting/jbi_table_s18_blinded_classification_review.csv"
RULE_KEY = ROOT / "docs/supporting/jbi_table_s19_rule_classification_key.csv"
REVIEW_PROTOCOL = ROOT / "docs/jbi_classification_review_protocol.md"
RULE_AUDIT = ROOT / "docs/jbi_classification_rule_audit.md"


def fail(message: str) -> None:
    raise SystemExit(message)


def require(text: str, phrase: str, label: str) -> None:
    if phrase not in text:
        fail(f"{label} missing transparency phrase: {phrase}")


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), list(reader)


def main() -> None:
    required_paths = (
        MANUSCRIPT,
        BASELINE_SCRIPT,
        ENRICHMENT_SCRIPT,
        RESOLVED_QUEUE,
        BLINDED_REVIEW,
        RULE_KEY,
        REVIEW_PROTOCOL,
        RULE_AUDIT,
    )
    for path in required_paths:
        if not path.is_file() or path.stat().st_size == 0:
            fail(f"Missing or empty classification-transparency source: {path}")

    manuscript = MANUSCRIPT.read_text(encoding="utf-8")
    baseline = BASELINE_SCRIPT.read_text(encoding="utf-8")
    enrichment = ENRICHMENT_SCRIPT.read_text(encoding="utf-8")
    protocol = REVIEW_PROTOCOL.read_text(encoding="utf-8")
    audit = RULE_AUDIT.read_text(encoding="utf-8")

    for phrase in (
        "source-traceable, rule-derived literature classifications",
        "rule-based screening of retained source text",
        "screening labels rather than final biological annotations",
        "review_status = unreviewed",
        "no field documenting completed independent human screening",
        "direct colour signal, no artificial signal and an evidence score of at least 20",
        "### Baseline-unambiguous association",
    ):
        require(manuscript, phrase, "manuscript")

    if re.search(r"\baudited\b", manuscript, flags=re.IGNORECASE):
        fail("Manuscript still uses `audited`, which may imply completed human review")

    for phrase in (
        "Automated labels are screening labels, not final biological annotations",
        "manual_review_required",
        "['within_population','among_population','mixed']",
        "polymorphic populations?",
    ):
        require(baseline, phrase, "baseline classifier")

    compact_enrichment = enrichment.replace(" ", "")
    for phrase in (
        "direct_colour_signal==1",
        "artificial_signal==0",
        "score>=20",
        "high_confidence_enrichment",
        "manual_review_required",
        "polymorphicpopulations?",
    ):
        require(compact_enrichment, phrase.replace(" ", ""), "enrichment classifier")

    queue_fields, queue_rows = read_csv(RESOLVED_QUEUE)
    if "review_status" not in queue_fields:
        fail("Resolved queue has no review_status field")
    if not queue_rows:
        fail("Resolved queue is empty")
    statuses = {str(row.get("review_status") or "").strip() for row in queue_rows}
    if statuses != {"unreviewed"}:
        fail(f"Resolved queue review_status values changed: {sorted(statuses)}")

    forbidden_reviewer_fields = {
        "reviewed_by",
        "reviewer",
        "second_reviewer",
        "adjudicator",
        "agreement",
        "agreement_statistic",
    }
    unexpected = sorted(set(queue_fields) & forbidden_reviewer_fields)
    if unexpected:
        fail(f"Reviewer metadata fields appeared and manuscript must be reassessed: {unexpected}")

    blind_fields, blind_rows = read_csv(BLINDED_REVIEW)
    key_fields, key_rows = read_csv(RULE_KEY)
    if len(blind_rows) != 34 or len(key_rows) != 34:
        fail(f"Review packet must contain 34 rows: S18={len(blind_rows)}, S19={len(key_rows)}")
    if {row["review_id"] for row in blind_rows} != {row["review_id"] for row in key_rows}:
        fail("S18 and S19 review IDs do not match")
    if "frozen_rule_label" in blind_fields or "queue_recomputed_label" in blind_fields:
        fail("Blinded S18 exposes rule labels")
    if any(
        str(row.get(field) or "").strip()
        for row in blind_rows
        for field in ("reviewer_label", "reviewer_name_or_initials", "review_date", "review_notes")
    ):
        fail("Blinded S18 has been partly filled; freeze completed reviews outside the template before replacing it")

    required_key_fields = {
        "frozen_rule_label",
        "queue_recomputed_label",
        "queue_rule_comparison",
        "source_match_status",
        "queue_review_status",
    }
    if not required_key_fields.issubset(key_fields):
        fail(f"S19 missing fields: {sorted(required_key_fields - set(key_fields))}")
    comparisons = {row["queue_rule_comparison"] for row in key_rows}
    if comparisons != {"matches_frozen_label"}:
        failures = [row["canonical_name"] for row in key_rows if row["queue_rule_comparison"] != "matches_frozen_label"]
        fail(f"Current queue rules do not reproduce frozen labels: {failures}")
    source_matches = {row["source_match_status"] for row in key_rows}
    if not source_matches.issubset({"matches_queue_best_doi", "matches_queue_best_openalex"}):
        failures = [row["canonical_name"] for row in key_rows if not row["source_match_status"].startswith("matches_queue_best")]
        fail(f"Classification source differs from queue best source: {failures}")
    if {row["queue_review_status"] for row in key_rows} != {"unreviewed"}:
        fail("S19 does not preserve the unreviewed queue status")

    for phrase in (
        "Keep S19 unavailable until all first-pass labels are frozen",
        "one reviewer verifies all 34 species",
        "source-traceable, rule-derived classifications",
    ):
        require(protocol, phrase, "review protocol")

    for phrase in (
        "plural `polymorphic populations`",
        "Number of frozen baseline classifications changed: **0**",
        "code-parity and reproducibility correction",
    ):
        require(audit, phrase, "rule audit")

    print(
        {
            "status": "pass",
            "resolved_queue_rows": len(queue_rows),
            "review_status_values": sorted(statuses),
            "documented_completed_human_review": False,
            "baseline_labels": "automated screening labels",
            "manual_review_queue_classes": ["mixed", "unclear"],
            "enrichment_rule": "direct colour signal; no artificial signal; score >= 20",
            "blinded_review_rows": len(blind_rows),
            "rule_key_rows": len(key_rows),
            "rule_reproduction": "34/34",
            "classification_source_matches": "34/34",
            "review_template_completed": False,
        }
    )


if __name__ == "__main__":
    main()

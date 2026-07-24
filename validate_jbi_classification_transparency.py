#!/usr/bin/env python3
"""Validate transparent reporting of the literature-derived spatial labels."""
from __future__ import annotations

import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
MANUSCRIPT = ROOT / "docs/jbi_manuscript_editorial_revision.md"
BASELINE_SCRIPT = ROOT / "analysis_evidence_spatial_scale.py"
ENRICHMENT_SCRIPT = ROOT / "analysis_enriched_spatial_scale.py"
RESOLVED_QUEUE = ROOT / "data/resolved_inputs/global_flower_colour_review_queue_resolved.csv"


def fail(message: str) -> None:
    raise SystemExit(message)


def require(text: str, phrase: str, label: str) -> None:
    if phrase not in text:
        fail(f"{label} missing transparency phrase: {phrase}")


def main() -> None:
    for path in (MANUSCRIPT, BASELINE_SCRIPT, ENRICHMENT_SCRIPT, RESOLVED_QUEUE):
        if not path.is_file() or path.stat().st_size == 0:
            fail(f"Missing or empty classification-transparency source: {path}")

    manuscript = MANUSCRIPT.read_text(encoding="utf-8")
    baseline = BASELINE_SCRIPT.read_text(encoding="utf-8")
    enrichment = ENRICHMENT_SCRIPT.read_text(encoding="utf-8")

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
    ):
        require(baseline, phrase, "baseline classifier")

    for phrase in (
        "direct_colour_signal==1",
        "artificial_signal==0",
        "score>=20",
        "high_confidence_enrichment",
        "manual_review_required",
    ):
        require(enrichment.replace(" ", ""), phrase.replace(" ", ""), "enrichment classifier")

    with RESOLVED_QUEUE.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = set(reader.fieldnames or [])
        if "review_status" not in fields:
            fail("Resolved queue has no review_status field")
        rows = list(reader)

    if not rows:
        fail("Resolved queue is empty")
    statuses = {str(row.get("review_status") or "").strip() for row in rows}
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
    unexpected = sorted(fields & forbidden_reviewer_fields)
    if unexpected:
        fail(f"Reviewer metadata fields appeared and manuscript must be reassessed: {unexpected}")

    print(
        {
            "status": "pass",
            "resolved_queue_rows": len(rows),
            "review_status_values": sorted(statuses),
            "documented_completed_human_review": False,
            "baseline_labels": "automated screening labels",
            "manual_review_queue_classes": ["mixed", "unclear"],
            "enrichment_rule": "direct colour signal; no artificial signal; score >= 20",
        }
    )


if __name__ == "__main__":
    main()

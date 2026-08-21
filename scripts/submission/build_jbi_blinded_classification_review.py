#!/usr/bin/env python3
"""Build the blinded 34-species classification-review sheet and separate rule key.

This utility never changes frozen classifications. It exposes source provenance and
evidence excerpts for human review while keeping the current rule label separate.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

from fcp_pipeline.evidence import (
    ALLOWED_SPATIAL_LABELS,
    clean_excerpt,
    rule_label,
    source_match,
)

BLINDED_FIELDS = [
    "review_id", "canonical_name", "family", "classification_source_id",
    "queue_best_doi", "queue_best_openalex_id", "queue_best_title",
    "evidence_excerpt", "source_match_status", "reviewer_label",
    "reviewer_name_or_initials", "review_date", "review_notes",
]
KEY_FIELDS = [
    "review_id", "canonical_name", "frozen_rule_label", "queue_recomputed_label",
    "queue_rule_comparison", "within_signal_in_queue_text", "geographic_signal_in_queue_text",
    "classification_source_id", "queue_best_doi", "queue_best_openalex_id",
    "source_match_status", "queue_review_status", "queue_review_reason",
]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", default="docs/supporting/frozen_classification_manifest.csv")
    ap.add_argument("--resolved-queue", default="data/resolved_inputs/global_flower_colour_review_queue_resolved.csv")
    ap.add_argument("--blinded-out", default="docs/supporting/blinded_classification_review.csv")
    ap.add_argument("--key-out", default="docs/supporting/rule_classification_key.csv")
    args = ap.parse_args()

    manifest = read_csv(Path(args.manifest))
    queue = read_csv(Path(args.resolved_queue))
    queue_by_name = {row["canonical_name"].strip(): row for row in queue}
    if len(manifest) != 34 or len({row["canonical_name"] for row in manifest}) != 34:
        raise SystemExit("Frozen classification manifest must contain 34 unique species")

    blinded: list[dict[str, str]] = []
    key: list[dict[str, str]] = []
    for index, row in enumerate(sorted(manifest, key=lambda x: x["canonical_name"]), start=1):
        name = row["canonical_name"].strip()
        q = queue_by_name.get(name)
        if q is None:
            raise SystemExit(f"Baseline species absent from resolved queue: {name}")
        frozen = row["spatial_scale"].strip()
        if frozen not in ALLOWED_SPATIAL_LABELS:
            raise SystemExit(f"Unexpected frozen label for {name}: {frozen}")

        title = q.get("best_title", "")
        evidence = q.get("best_match_evidence", "")
        reason = q.get("review_reason", "")
        recomputed, within, geographic = rule_label(" ".join((title, evidence, reason)))
        review_id = f"JBI-{index:03d}"
        match = source_match(row.get("source_id", ""), q.get("best_doi", ""), q.get("best_openalex_id", ""))

        blinded.append({
            "review_id": review_id,
            "canonical_name": name,
            "family": row.get("family", ""),
            "classification_source_id": row.get("source_id", ""),
            "queue_best_doi": q.get("best_doi", ""),
            "queue_best_openalex_id": q.get("best_openalex_id", ""),
            "queue_best_title": title,
            "evidence_excerpt": clean_excerpt(evidence),
            "source_match_status": match,
            "reviewer_label": "",
            "reviewer_name_or_initials": "",
            "review_date": "",
            "review_notes": "",
        })
        key.append({
            "review_id": review_id,
            "canonical_name": name,
            "frozen_rule_label": frozen,
            "queue_recomputed_label": recomputed,
            "queue_rule_comparison": "matches_frozen_label" if recomputed == frozen else "differs_from_frozen_label",
            "within_signal_in_queue_text": str(within),
            "geographic_signal_in_queue_text": str(geographic),
            "classification_source_id": row.get("source_id", ""),
            "queue_best_doi": q.get("best_doi", ""),
            "queue_best_openalex_id": q.get("best_openalex_id", ""),
            "source_match_status": match,
            "queue_review_status": q.get("review_status", ""),
            "queue_review_reason": reason,
        })

    write_csv(Path(args.blinded_out), BLINDED_FIELDS, blinded)
    write_csv(Path(args.key_out), KEY_FIELDS, key)
    print({
        "status": "pass",
        "rows": len(blinded),
        "source_match_counts": {x: sum(r["source_match_status"] == x for r in blinded) for x in sorted({r["source_match_status"] for r in blinded})},
        "queue_rule_comparison_counts": {x: sum(r["queue_rule_comparison"] == x for r in key) for x in sorted({r["queue_rule_comparison"] for r in key})},
    })


if __name__ == "__main__":
    main()

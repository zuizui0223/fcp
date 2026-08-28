#!/usr/bin/env python3
"""Build a compact adjudication queue from merged duplicate record-screen reviews."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

FIELDS = [
    "record_relevance",
    "natural_intraspecific_variation",
    "floral_display_colour",
    "full_text_required",
    "focal_taxon_text",
]


def clean(value):
    return " ".join(str(value or "").split()).strip()


def norm(value):
    return clean(value).lower()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--reviewed", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--summary-out", required=True)
    args = p.parse_args()

    with Path(args.reviewed).open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise SystemExit("No reviewed rows")

    out_rows = []
    by_field = {f: 0 for f in FIELDS}
    for row in rows:
        disagree = []
        for field in FIELDS:
            a = norm(row.get(f"reviewer_1_{field}"))
            b = norm(row.get(f"reviewer_2_{field}"))
            if a and b and a != b:
                disagree.append(field)
                by_field[field] += 1
        if not disagree:
            continue
        out_rows.append({
            "record_review_id": row.get("record_review_id", ""),
            "batch_id": row.get("batch_id", ""),
            "source_id": row.get("source_id", ""),
            "title": row.get("title", ""),
            "abstract": row.get("abstract", ""),
            "language": row.get("language", ""),
            "disagreement_fields": ";".join(disagree),
            "reviewer_1_record_relevance": row.get("reviewer_1_record_relevance", ""),
            "reviewer_2_record_relevance": row.get("reviewer_2_record_relevance", ""),
            "reviewer_1_natural_intraspecific_variation": row.get("reviewer_1_natural_intraspecific_variation", ""),
            "reviewer_2_natural_intraspecific_variation": row.get("reviewer_2_natural_intraspecific_variation", ""),
            "reviewer_1_floral_display_colour": row.get("reviewer_1_floral_display_colour", ""),
            "reviewer_2_floral_display_colour": row.get("reviewer_2_floral_display_colour", ""),
            "reviewer_1_focal_taxon_text": row.get("reviewer_1_focal_taxon_text", ""),
            "reviewer_2_focal_taxon_text": row.get("reviewer_2_focal_taxon_text", ""),
            "reviewer_1_full_text_required": row.get("reviewer_1_full_text_required", ""),
            "reviewer_2_full_text_required": row.get("reviewer_2_full_text_required", ""),
            "adjudicated_record_relevance": "",
            "adjudicated_natural_intraspecific_variation": "",
            "adjudicated_floral_display_colour": "",
            "adjudicated_focal_taxon_text": "",
            "adjudicated_full_text_required": "",
            "adjudication_notes": "",
        })

    fields = list(out_rows[0].keys()) if out_rows else [
        "record_review_id", "batch_id", "source_id", "title", "abstract", "language",
        "disagreement_fields", "reviewer_1_record_relevance", "reviewer_2_record_relevance",
        "reviewer_1_natural_intraspecific_variation", "reviewer_2_natural_intraspecific_variation",
        "reviewer_1_floral_display_colour", "reviewer_2_floral_display_colour",
        "reviewer_1_focal_taxon_text", "reviewer_2_focal_taxon_text",
        "reviewer_1_full_text_required", "reviewer_2_full_text_required",
        "adjudicated_record_relevance", "adjudicated_natural_intraspecific_variation",
        "adjudicated_floral_display_colour", "adjudicated_focal_taxon_text",
        "adjudicated_full_text_required", "adjudication_notes"
    ]
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with Path(args.out).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        if out_rows:
            writer.writerows(out_rows)

    summary = {
        "status": "complete",
        "reviewed_records": len(rows),
        "records_requiring_adjudication": len(out_rows),
        "disagreements_by_field": by_field,
        "semantic_guard": "Queue creation only; disagreements remain unresolved until explicit adjudication.",
    }
    Path(args.summary_out).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

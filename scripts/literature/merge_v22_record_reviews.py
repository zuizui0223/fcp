#!/usr/bin/env python3
"""Merge independently completed Reviewer 1/2 v2.2 record-screening sheets.

The merger refuses silent metadata drift, duplicate IDs, or reviewer-column cross-contamination.
It does not adjudicate disagreements.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

R1_PREFIX = "reviewer_1_"
R2_PREFIX = "reviewer_2_"
CORE_SUFFIXES = [
    "record_relevance",
    "natural_intraspecific_variation",
    "floral_display_colour",
    "full_text_required",
]


def clean(value):
    return " ".join(str(value or "").split()).strip()


def read_rows(path: str):
    with Path(path).open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise SystemExit(f"No rows: {path}")
    ids = [clean(r.get("record_review_id")) for r in rows]
    if any(not x for x in ids):
        raise SystemExit(f"Missing record_review_id: {path}")
    if len(ids) != len(set(ids)):
        raise SystemExit(f"Duplicate record_review_id: {path}")
    return rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--reviewer1", required=True)
    p.add_argument("--reviewer2", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--summary-out", required=True)
    args = p.parse_args()

    r1 = read_rows(args.reviewer1)
    r2 = read_rows(args.reviewer2)
    m1 = {r["record_review_id"]: r for r in r1}
    m2 = {r["record_review_id"]: r for r in r2}
    if set(m1) != set(m2):
        only1 = sorted(set(m1) - set(m2))[:20]
        only2 = sorted(set(m2) - set(m1))[:20]
        raise SystemExit(f"Reviewer ID sets differ: only_r1={only1}, only_r2={only2}")

    all_fields = list(dict.fromkeys(list(r1[0].keys()) + list(r2[0].keys())))
    immutable = [
        f for f in all_fields
        if not f.startswith(R1_PREFIX)
        and not f.startswith(R2_PREFIX)
        and not f.startswith("adjudicated_")
        and f not in {"adjudication_notes", "review_status"}
    ]

    metadata_mismatches = []
    cross_contamination = []
    merged = []
    for rid in [r["record_review_id"] for r in r1]:
        a, b = m1[rid], m2[rid]
        for field in immutable:
            if clean(a.get(field)) != clean(b.get(field)):
                metadata_mismatches.append({"record_review_id": rid, "field": field})
        for field in all_fields:
            if field.startswith(R2_PREFIX) and clean(a.get(field)):
                cross_contamination.append({"record_review_id": rid, "sheet": "reviewer1", "field": field})
            if field.startswith(R1_PREFIX) and clean(b.get(field)):
                cross_contamination.append({"record_review_id": rid, "sheet": "reviewer2", "field": field})

        row = {field: a.get(field, "") for field in all_fields}
        for field in all_fields:
            if field.startswith(R2_PREFIX):
                row[field] = b.get(field, "")
        row["review_status"] = "double_review_complete" if all(
            clean(row.get(R1_PREFIX + s)) and clean(row.get(R2_PREFIX + s))
            for s in CORE_SUFFIXES
        ) else "double_review_incomplete"
        merged.append(row)

    if metadata_mismatches:
        raise SystemExit(f"Immutable metadata changed in {len(metadata_mismatches)} cells; first={metadata_mismatches[:10]}")
    if cross_contamination:
        raise SystemExit(f"Reviewer-column cross-contamination detected; first={cross_contamination[:10]}")

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with Path(args.out).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=all_fields)
        writer.writeheader()
        writer.writerows(merged)

    completeness = {}
    for suffix in CORE_SUFFIXES:
        completeness[suffix] = {
            "reviewer_1_nonblank": sum(bool(clean(r.get(R1_PREFIX + suffix))) for r in merged),
            "reviewer_2_nonblank": sum(bool(clean(r.get(R2_PREFIX + suffix))) for r in merged),
            "double_coded": sum(
                bool(clean(r.get(R1_PREFIX + suffix))) and bool(clean(r.get(R2_PREFIX + suffix)))
                for r in merged
            ),
        }
    summary = {
        "status": "complete",
        "records": len(merged),
        "fully_double_reviewed_records": sum(r["review_status"] == "double_review_complete" for r in merged),
        "completeness": completeness,
        "metadata_mismatches": 0,
        "cross_contamination": 0,
        "semantic_guard": "Merge only; no disagreement is resolved and no biological inclusion is inferred.",
    }
    Path(args.summary_out).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

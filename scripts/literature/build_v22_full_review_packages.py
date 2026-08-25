#!/usr/bin/env python3
"""Split the canonical 12,064-record blind sheet into reviewer-specific batch files."""
from __future__ import annotations

import argparse
import csv
import json
import shutil
from collections import Counter
from pathlib import Path

R1_PREFIX = "reviewer_1_"
R2_PREFIX = "reviewer_2_"


def read_rows(path):
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write(path: Path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--blind", required=True)
    p.add_argument("--codebook", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--summary-out", required=True)
    args = p.parse_args()

    rows = read_rows(args.blind)
    if len(rows) != 12064:
        raise SystemExit(f"Expected 12064 rows, found {len(rows)}")
    ids = [r.get("record_review_id", "") for r in rows]
    if len(ids) != len(set(ids)):
        raise SystemExit("Duplicate record_review_id in canonical blind sheet")

    batches = Counter(r.get("batch_id", "") for r in rows)
    expected = [f"B{i:02d}" for i in range(1, 14)]
    if sorted(batches) != expected:
        raise SystemExit(f"Expected batches {expected}, found {sorted(batches)}")
    if any(batches[b] > 1000 for b in expected):
        raise SystemExit(f"Batch over 1000 records: {dict(batches)}")
    if sum(batches.values()) != 12064:
        raise SystemExit("Batch counts do not sum to 12064")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    base_fields = list(rows[0].keys())
    r1_fields = [
        f for f in base_fields
        if not f.startswith(R2_PREFIX) and not f.startswith("adjudicated_") and f != "adjudication_notes"
    ]
    r2_fields = [
        f for f in base_fields
        if not f.startswith(R1_PREFIX) and not f.startswith("adjudicated_") and f != "adjudication_notes"
    ]

    for batch in expected:
        subset = [r for r in rows if r.get("batch_id") == batch]
        # Reviewer-specific sheets deliberately omit the other reviewer's columns.
        write(outdir / "Reviewer1" / f"FCP_{batch}_Reviewer1.csv", subset, r1_fields)
        write(outdir / "Reviewer2" / f"FCP_{batch}_Reviewer2.csv", subset, r2_fields)

    shutil.copyfile(args.codebook, outdir / "JBI_V22_RECORD_SCREENING_CODEBOOK.md")
    readme = (
        "# FCP JBI v2.2 full record-screening packages\n\n"
        "- Reviewer 1 uses only files in `Reviewer1/`.\n"
        "- Reviewer 2 uses only files in `Reviewer2/`.\n"
        "- Reviewers work independently and must not exchange decisions before a batch is locked.\n"
        "- Do not provide reviewers with coordinator keys, query memberships, historical labels, automated taxon hints, climatic results, or the other reviewer's files.\n"
        "- Complete Wave 0 calibration and pass the prespecified gate before starting B01.\n"
        "- Return completed files without changing bibliographic columns or record_review_id.\n"
    )
    (outdir / "README_REVIEWERS.md").write_text(readme, encoding="utf-8")

    summary = {
        "status": "complete",
        "records": len(rows),
        "review_batches": 13,
        "batch_counts": dict(sorted(batches.items())),
        "reviewer1_files": 13,
        "reviewer2_files": 13,
        "other_reviewer_columns_present": False,
        "adjudication_columns_present": False,
        "gate_before_B01": "Wave 0 calibration must pass the prespecified gate",
        "semantic_guard": "Packaging changes workload only; it does not determine record inclusion.",
    }
    Path(args.summary_out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.summary_out).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

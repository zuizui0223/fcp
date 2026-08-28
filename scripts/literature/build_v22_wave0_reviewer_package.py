#!/usr/bin/env python3
"""Build reviewer-only Wave 0 packages without exposing coordinator-only columns."""
from __future__ import annotations

import argparse
import csv
import json
import shutil
from pathlib import Path


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def project_write(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--blind", required=True)
    parser.add_argument("--codebook", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--summary-out", required=True)
    args = parser.parse_args()

    rows = read_rows(Path(args.blind))
    if len(rows) != 384:
        raise SystemExit(f"Expected 384 Wave 0 rows, found {len(rows)}")
    if len({row.get('record_review_id', '') for row in rows}) != 384:
        raise SystemExit("Wave 0 record_review_id values are not unique")

    all_fields = list(rows[0].keys())
    forbidden_common = {
        "calibration_stratum",
        "query_ids",
        "representative_query_id",
        "historical_34_exact_source_match",
        "detected_binomial_strings",
        "screen_priority",
        "automated_within_signal",
        "automated_among_signal",
    }
    leaked = forbidden_common.intersection(all_fields)
    if leaked:
        raise SystemExit(f"Coordinator-only fields leaked into Wave 0 blind input: {sorted(leaked)}")

    reviewer1_fields = [
        field for field in all_fields
        if not field.startswith("reviewer_2_")
        and not field.startswith("adjudicated_")
        and field != "adjudication_notes"
    ]
    reviewer2_fields = [
        field for field in all_fields
        if not field.startswith("reviewer_1_")
        and not field.startswith("adjudicated_")
        and field != "adjudication_notes"
    ]

    if any(field.startswith("reviewer_2_") for field in reviewer1_fields):
        raise SystemExit("Reviewer 2 fields remain in Reviewer 1 package")
    if any(field.startswith("reviewer_1_") for field in reviewer2_fields):
        raise SystemExit("Reviewer 1 fields remain in Reviewer 2 package")
    if any(field.startswith("adjudicated_") for field in reviewer1_fields + reviewer2_fields):
        raise SystemExit("Adjudication fields remain in reviewer package")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    project_write(outdir / "FCP_Wave0_Reviewer1.csv", rows, reviewer1_fields)
    project_write(outdir / "FCP_Wave0_Reviewer2.csv", rows, reviewer2_fields)
    shutil.copyfile(args.codebook, outdir / "JBI_V22_RECORD_SCREENING_CODEBOOK.md")

    readme = """# FCP JBI v2.2 — Wave 0 independent calibration

This directory is reviewer-facing only. No coordinator key is included.

- Reviewer 1 uses only `FCP_Wave0_Reviewer1.csv` and fills `reviewer_1_*` fields.
- Reviewer 2 uses only `FCP_Wave0_Reviewer2.csv` and fills `reviewer_2_*` fields.
- Reviewers must work independently and must not compare decisions before both files are locked.
- Do not change `record_review_id` or bibliographic columns.
- Use `uncertain` and/or `full_text_required=yes` when title/abstract evidence is insufficient.
- Return both completed CSV files to the coordinator.
- B01–B13 must not start until Wave 0 passes the prespecified gate in the codebook.
"""
    (outdir / "README_REVIEWERS.md").write_text(readme, encoding="utf-8")

    summary = {
        "status": "complete",
        "records_per_reviewer": 384,
        "reviewer1_files": ["FCP_Wave0_Reviewer1.csv"],
        "reviewer2_files": ["FCP_Wave0_Reviewer2.csv"],
        "coordinator_key_included": False,
        "other_reviewer_columns_present": False,
        "adjudication_columns_present": False,
        "semantic_guard": "Packaging changes visibility only; it does not adjudicate or classify records.",
    }
    Path(args.summary_out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.summary_out).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

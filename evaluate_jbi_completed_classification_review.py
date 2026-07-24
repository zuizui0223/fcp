#!/usr/bin/env python3
"""Evaluate completed blinded classification-review sheets.

The tool compares one or more completed copies of S18 with each other and with the
separate S19 key. It reports agreement and discrepancies but never changes the
frozen classification manifest or analysis inputs.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from itertools import combinations
from pathlib import Path

ALLOWED = {"within_population", "among_population", "mixed", "unclear"}


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), list(reader)


def cohen_kappa(a: list[str], b: list[str]) -> tuple[float, float, float]:
    if len(a) != len(b) or not a:
        raise ValueError("Kappa inputs must be non-empty and equal length")
    n = len(a)
    observed = sum(x == y for x, y in zip(a, b)) / n
    ca = Counter(a)
    cb = Counter(b)
    expected = sum((ca[label] / n) * (cb[label] / n) for label in ALLOWED)
    kappa = (observed - expected) / (1 - expected) if expected < 1 else 1.0
    return observed, expected, kappa


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--review", action="append", required=True, help="Completed copy of S18; repeat for independent reviewers")
    parser.add_argument("--key", default="docs/supporting/jbi_table_s19_rule_classification_key.csv")
    parser.add_argument("--outdir", default="analysis_outputs/jbi_classification_review")
    args = parser.parse_args()

    key_path = Path(args.key)
    key_fields, key_rows = read_csv(key_path)
    if len(key_rows) != 34:
        raise SystemExit(f"Expected 34 key rows; found {len(key_rows)}")
    required_key = {"review_id", "canonical_name", "frozen_rule_label"}
    if not required_key.issubset(key_fields):
        raise SystemExit(f"Key missing fields: {sorted(required_key - set(key_fields))}")
    key = {row["review_id"]: row for row in key_rows}

    review_sets: list[dict[str, object]] = []
    for path_text in args.review:
        path = Path(path_text)
        fields, rows = read_csv(path)
        required = {"review_id", "canonical_name", "reviewer_label", "reviewer_name_or_initials", "review_date", "review_notes"}
        if not required.issubset(fields):
            raise SystemExit(f"Review {path} missing fields: {sorted(required - set(fields))}")
        if len(rows) != 34:
            raise SystemExit(f"Review {path} must contain 34 rows; found {len(rows)}")
        ids = [row["review_id"] for row in rows]
        if len(set(ids)) != 34 or set(ids) != set(key):
            raise SystemExit(f"Review IDs in {path} do not match S19")
        missing_labels = [row["review_id"] for row in rows if row["reviewer_label"].strip() not in ALLOWED]
        if missing_labels:
            raise SystemExit(f"Review {path} has blank or invalid labels: {missing_labels}")
        missing_identity = [row["review_id"] for row in rows if not row["reviewer_name_or_initials"].strip()]
        missing_date = [row["review_id"] for row in rows if not row["review_date"].strip()]
        if missing_identity or missing_date:
            raise SystemExit(f"Review {path} is missing reviewer identity/date: identity={missing_identity}, date={missing_date}")
        by_id = {row["review_id"]: row for row in rows}
        reviewer_names = sorted({row["reviewer_name_or_initials"].strip() for row in rows})
        if len(reviewer_names) != 1:
            raise SystemExit(f"Each review file must represent one reviewer; found {reviewer_names} in {path}")
        review_sets.append({"path": str(path), "reviewer": reviewer_names[0], "rows": by_id})

    ordered_ids = sorted(key)
    comparison_rows: list[dict[str, str]] = []
    for review_id in ordered_ids:
        labels = [str(review["rows"][review_id]["reviewer_label"]).strip() for review in review_sets]
        all_agree = len(set(labels)) == 1
        consensus = labels[0] if all_agree else "unresolved_disagreement"
        row = {
            "review_id": review_id,
            "canonical_name": key[review_id]["canonical_name"],
            "frozen_rule_label": key[review_id]["frozen_rule_label"],
            "reviewer_labels": ";".join(labels),
            "reviewers_agree": str(all_agree).lower(),
            "consensus_label": consensus,
            "consensus_matches_frozen": str(all_agree and consensus == key[review_id]["frozen_rule_label"]).lower(),
            "requires_adjudication": str(not all_agree).lower(),
        }
        comparison_rows.append(row)

    pairwise: list[dict[str, object]] = []
    for left, right in combinations(review_sets, 2):
        labels_left = [left["rows"][rid]["reviewer_label"].strip() for rid in ordered_ids]
        labels_right = [right["rows"][rid]["reviewer_label"].strip() for rid in ordered_ids]
        observed, expected, kappa = cohen_kappa(labels_left, labels_right)
        pairwise.append(
            {
                "reviewer_1": left["reviewer"],
                "reviewer_2": right["reviewer"],
                "raw_agreement": observed,
                "expected_agreement": expected,
                "cohen_kappa": kappa,
                "n": len(ordered_ids),
            }
        )

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    comparison_path = outdir / "jbi_classification_review_comparison.csv"
    with comparison_path.open("w", newline="", encoding="utf-8") as handle:
        fields = list(comparison_rows[0])
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(comparison_rows)

    disagreements = [row for row in comparison_rows if row["requires_adjudication"] == "true"]
    rule_disagreements = [
        row for row in comparison_rows
        if row["reviewers_agree"] == "true" and row["consensus_matches_frozen"] == "false"
    ]
    summary = {
        "status": "complete",
        "review_files": [review["path"] for review in review_sets],
        "reviewers": [review["reviewer"] for review in review_sets],
        "n_species": len(ordered_ids),
        "n_reviewers": len(review_sets),
        "pairwise_agreement": pairwise,
        "n_interreviewer_disagreements": len(disagreements),
        "interreviewer_disagreement_ids": [row["review_id"] for row in disagreements],
        "n_consensus_disagreements_with_frozen_rule": len(rule_disagreements),
        "consensus_disagreement_ids": [row["review_id"] for row in rule_disagreements],
        "all_species_resolved": not disagreements,
        "automatic_analysis_update_performed": False,
        "next_action": "Adjudicate disagreements, document accepted decisions, update the correction log, then rerun analyses."
        if disagreements or rule_disagreements
        else "Record verification in the correction/audit record before describing labels as human verified.",
    }
    (outdir / "jbi_classification_review_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()

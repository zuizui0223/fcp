#!/usr/bin/env python3
"""Build a blinded human-review sheet and a separate rule-label key.

The blinded sheet does not expose the current rule-derived spatial label. Human
reviewers should fill it before opening the key. This script does not overwrite
classifications or analysis inputs.
"""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

from analysis_evidence_spatial_scale import GEOGRAPHIC, WITHIN

BLINDED_FIELDS = [
    "review_id",
    "canonical_name",
    "family",
    "classification_source_id",
    "queue_best_doi",
    "queue_best_openalex_id",
    "queue_best_title",
    "evidence_excerpt",
    "source_match_status",
    "reviewer_label",
    "reviewer_name_or_initials",
    "review_date",
    "review_notes",
]

KEY_FIELDS = [
    "review_id",
    "canonical_name",
    "rule_label",
    "within_signal",
    "geographic_signal",
    "classification_source_id",
    "queue_best_doi",
    "queue_best_openalex_id",
    "source_match_status",
    "queue_review_status",
    "queue_review_reason",
]

ALLOWED_LABELS = {"within_population", "among_population", "mixed", "unclear"}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def normalize_doi(value: str) -> str:
    x = str(value or "").strip().lower()
    x = re.sub(r"^https?://(?:dx\.)?doi\.org/", "", x)
    return x


def clean_excerpt(value: str, limit: int = 1800) -> str:
    x = re.sub(r"\s+", " ", str(value or "")).strip()
    return x[:limit]


def source_match(source_id: str, doi: str, openalex_id: str) -> str:
    source = str(source_id or "").strip()
    if not source:
        return "classification_source_missing"
    if normalize_doi(source) and normalize_doi(source) == normalize_doi(doi):
        return "matches_queue_best_doi"
    if source.rstrip("/").lower() == str(openalex_id or "").strip().rstrip("/").lower():
        return "matches_queue_best_openalex"
    return "classification_source_differs_from_queue_best"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest",
        default="docs/supporting/jbi_table_s6_frozen_classification_manifest.csv",
    )
    parser.add_argument(
        "--resolved-queue",
        default="data/resolved_inputs/global_flower_colour_review_queue_resolved.csv",
    )
    parser.add_argument(
        "--blinded-out",
        default="docs/supporting/jbi_table_s18_blinded_classification_review.csv",
    )
    parser.add_argument(
        "--key-out",
        default="docs/supporting/jbi_table_s19_rule_classification_key.csv",
    )
    args = parser.parse_args()

    manifest = read_csv(Path(args.manifest))
    queue_rows = read_csv(Path(args.resolved_queue))
    queue_by_name = {row["canonical_name"].strip(): row for row in queue_rows}

    if len(manifest) != 34:
        raise SystemExit(f"Expected 34 baseline rows; found {len(manifest)}")
    if len({row["canonical_name"] for row in manifest}) != 34:
        raise SystemExit("Baseline manifest contains duplicate species")

    blinded: list[dict[str, str]] = []
    key: list[dict[str, str]] = []

    for index, row in enumerate(sorted(manifest, key=lambda x: x["canonical_name"]), start=1):
        name = row["canonical_name"].strip()
        q = queue_by_name.get(name)
        if q is None:
            raise SystemExit(f"Baseline species absent from resolved queue: {name}")
        label = row["spatial_scale"].strip()
        if label not in ALLOWED_LABELS:
            raise SystemExit(f"Unexpected baseline label for {name}: {label}")

        title = q.get("best_title", "")
        evidence = q.get("best_match_evidence", "")
        reason = q.get("review_reason", "")
        rule_text = " ".join((title, evidence, reason))
        within = int(bool(WITHIN.search(rule_text)))
        geographic = int(bool(GEOGRAPHIC.search(rule_text)))
        recomputed = (
            "within_population" if within and not geographic else
            "among_population" if geographic and not within else
            "mixed" if within and geographic else
            "unclear"
        )
        if recomputed != label:
            raise SystemExit(
                f"Frozen label does not match current rule for {name}: {label} != {recomputed}"
            )

        review_id = f"JBI-{index:03d}"
        match = source_match(
            row.get("source_id", ""),
            q.get("best_doi", ""),
            q.get("best_openalex_id", ""),
        )
        blinded.append(
            {
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
            }
        )
        key.append(
            {
                "review_id": review_id,
                "canonical_name": name,
                "rule_label": label,
                "within_signal": str(within),
                "geographic_signal": str(geographic),
                "classification_source_id": row.get("source_id", ""),
                "queue_best_doi": q.get("best_doi", ""),
                "queue_best_openalex_id": q.get("best_openalex_id", ""),
                "source_match_status": match,
                "queue_review_status": q.get("review_status", ""),
                "queue_review_reason": q.get("review_reason", ""),
            }
        )

    for output, fields, rows in (
        (Path(args.blinded_out), BLINDED_FIELDS, blinded),
        (Path(args.key_out), KEY_FIELDS, key),
    ):
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)

    print(
        {
            "status": "pass",
            "rows": len(blinded),
            "blinded_output": str(args.blinded_out),
            "key_output": str(args.key_out),
            "source_match_counts": {
                status: sum(row["source_match_status"] == status for row in blinded)
                for status in sorted({row["source_match_status"] for row in blinded})
            },
        }
    )


if __name__ == "__main__":
    main()

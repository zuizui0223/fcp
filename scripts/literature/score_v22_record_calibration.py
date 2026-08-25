#!/usr/bin/env python3
"""Score duplicate blind v2.2 record-screening calibration decisions."""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

FIELDS = [
    ("record_relevance", "reviewer_1_record_relevance", "reviewer_2_record_relevance"),
    (
        "natural_intraspecific_variation",
        "reviewer_1_natural_intraspecific_variation",
        "reviewer_2_natural_intraspecific_variation",
    ),
    ("floral_display_colour", "reviewer_1_floral_display_colour", "reviewer_2_floral_display_colour"),
    ("full_text_required", "reviewer_1_full_text_required", "reviewer_2_full_text_required"),
]


def clean(value):
    return " ".join(str(value or "").split()).strip().lower()


def kappa(values_a, values_b):
    pairs = [(clean(a), clean(b)) for a, b in zip(values_a, values_b) if clean(a) and clean(b)]
    n = len(pairs)
    if not n:
        return {"n_double_coded": 0, "raw_agreement": None, "cohen_kappa": None, "labels": []}
    observed = sum(a == b for a, b in pairs) / n
    counts_a = Counter(a for a, _ in pairs)
    counts_b = Counter(b for _, b in pairs)
    labels = set(counts_a) | set(counts_b)
    expected = sum((counts_a[label] / n) * (counts_b[label] / n) for label in labels)
    coefficient = None if abs(1 - expected) < 1e-12 else (observed - expected) / (1 - expected)
    return {
        "n_double_coded": n,
        "raw_agreement": observed,
        "cohen_kappa": coefficient,
        "labels": sorted(labels),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reviewed", required=True)
    parser.add_argument("--summary-out", required=True)
    parser.add_argument("--disagreements-out", required=True)
    args = parser.parse_args()

    with Path(args.reviewed).open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    metrics = {}
    disagreements = []
    for name, field_1, field_2 in FIELDS:
        metrics[name] = kappa(
            [row.get(field_1, "") for row in rows],
            [row.get(field_2, "") for row in rows],
        )
        for row in rows:
            value_1 = clean(row.get(field_1, ""))
            value_2 = clean(row.get(field_2, ""))
            if value_1 and value_2 and value_1 != value_2:
                disagreements.append({
                    "record_review_id": row.get("record_review_id", ""),
                    "field": name,
                    "reviewer_1": value_1,
                    "reviewer_2": value_2,
                    "title": row.get("title", ""),
                })

    taxon_pairs = []
    for row in rows:
        value_1 = clean(row.get("reviewer_1_focal_taxon_text", ""))
        value_2 = clean(row.get("reviewer_2_focal_taxon_text", ""))
        if value_1 and value_2:
            taxon_pairs.append((value_1, value_2))
    taxon_agreement = (
        None if not taxon_pairs else sum(a == b for a, b in taxon_pairs) / len(taxon_pairs)
    )

    summary = {
        "status": "complete",
        "records": len(rows),
        "agreement": metrics,
        "focal_taxon_text": {
            "n_double_coded": len(taxon_pairs),
            "normalized_exact_agreement": taxon_agreement,
        },
        "n_disagreements": len(disagreements),
        "semantic_guard": (
            "Agreement statistics assess duplicate coding consistency only. They do not adjudicate "
            "record inclusion or spatial state."
        ),
    }
    Path(args.summary_out).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    fields = ["record_review_id", "field", "reviewer_1", "reviewer_2", "title"]
    with Path(args.disagreements_out).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(disagreements)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

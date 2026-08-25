#!/usr/bin/env python3
"""Prepare source-level blinded review materials for the upstream FCP re-audit.

The human-review unit is species × source. Reviewer-facing output contains only
bibliographic/source context and empty review fields. Automated screening signals,
priority labels and other navigation metadata are written to a separate coordinator-only
key and must not be shown to reviewers before independent coding is complete.

No climatic variables or ecological model outputs are read here.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
from collections import defaultdict
from pathlib import Path
from typing import Any

PRIORITY_ORDER = {
    "P1_high_natural_itv": 0,
    "P1_high_population_itv": 1,
    "P2_possible_itv": 2,
    "P3_likely_exclusion_review": 3,
    "P4_flower_colour_context_only": 4,
    "P5_low_relevance": 5,
}

FORBIDDEN_REVIEW_PREFIXES = ("automated_", "historical_", "systematic_screen_")
FORBIDDEN_REVIEW_COLUMNS = {"screen_priority", "spatial_scale", "automated_evidence_state"}


def clean(value: Any) -> str:
    return " ".join(str(value or "").split())


def truthy(value: Any) -> bool:
    return clean(value).lower() in {"1", "true", "yes"}


def integer(value: Any) -> int:
    try:
        return int(float(value or 0))
    except (TypeError, ValueError):
        return 0


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def source_key(row: dict[str, str]) -> str:
    return clean(row.get("doi") or row.get("record_id") or row.get("title"))


def review_fields() -> dict[str, str]:
    return {
        "review_status": "unreviewed",
        "reviewer_1_initials": "",
        "reviewer_1_date": "",
        "reviewer_1_display_colour_relevant": "",
        "reviewer_1_natural_eligibility": "",
        "reviewer_1_variation_form": "",
        "reviewer_1_local_coexistence": "",
        "reviewer_1_geographic_structure": "",
        "reviewer_1_local_evidence": "",
        "reviewer_1_geographic_evidence": "",
        "reviewer_1_notes": "",
        "reviewer_2_initials": "",
        "reviewer_2_date": "",
        "reviewer_2_display_colour_relevant": "",
        "reviewer_2_natural_eligibility": "",
        "reviewer_2_variation_form": "",
        "reviewer_2_local_coexistence": "",
        "reviewer_2_geographic_structure": "",
        "reviewer_2_local_evidence": "",
        "reviewer_2_geographic_evidence": "",
        "reviewer_2_notes": "",
        "adjudicator_initials": "",
        "adjudication_date": "",
        "adjudicated_display_colour_relevant": "",
        "adjudicated_natural_eligibility": "",
        "adjudicated_variation_form": "",
        "adjudicated_local_coexistence": "",
        "adjudicated_geographic_structure": "",
        "adjudicated_local_evidence": "",
        "adjudicated_geographic_evidence": "",
        "adjudication_notes": "",
    }


def merge_group(accepted_name: str, family: str, rows: list[dict[str, str]]) -> tuple[dict[str, Any], dict[str, Any]]:
    best = min(
        rows,
        key=lambda row: (
            PRIORITY_ORDER.get(clean(row.get("screen_priority")), 99),
            -len(clean(row.get("evidence_excerpt"))),
        ),
    )
    excerpts = sorted(
        {clean(row.get("evidence_excerpt")) for row in rows if clean(row.get("evidence_excerpt"))},
        key=len,
        reverse=True,
    )
    source = source_key(best)
    review = {
        "canonical_name": accepted_name,
        "family": family,
        "source_id": source,
        "doi": clean(best.get("doi")),
        "record_id": clean(best.get("record_id")),
        "title": clean(best.get("title")),
        "navigation_excerpt": excerpts[0][:2000] if excerpts else "",
        **review_fields(),
    }
    key = {
        "canonical_name": accepted_name,
        "family": family,
        "source_id": source,
        "input_names": ";".join(sorted({clean(row.get("input_name")) for row in rows if clean(row.get("input_name"))})),
        "screen_priority": clean(best.get("screen_priority")),
        "automated_within_signal": max(integer(row.get("within_signal")) for row in rows),
        "automated_among_signal": max(integer(row.get("among_signal")) for row in rows),
        "automated_natural_signal": max(integer(row.get("natural_signal")) for row in rows),
        "automated_cultivated_signal": max(integer(row.get("cultivated_signal")) for row in rows),
        "automated_induced_signal": max(integer(row.get("induced_signal")) for row in rows),
        "automated_ontogenetic_signal": max(integer(row.get("ontogenetic_signal")) for row in rows),
        "automated_non_display_floral_signal": max(integer(row.get("non_display_floral_signal")) for row in rows),
        "automated_variation_form": ";".join(sorted({clean(row.get("provisional_variation_form")) for row in rows if clean(row.get("provisional_variation_form"))})),
        "n_merged_candidate_links": len(rows),
    }
    return review, key


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise SystemExit(f"No rows to write: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--record-links", required=True)
    parser.add_argument("--resolution-audit", required=True)
    parser.add_argument("--review-out", required=True)
    parser.add_argument("--key-out", required=True)
    args = parser.parse_args()

    links = read_csv(Path(args.record_links))
    audit = read_csv(Path(args.resolution_audit))
    resolution = {
        clean(row.get("input_name")): row
        for row in audit
        if truthy(row.get("accepted")) and clean(row.get("accepted_name"))
    }

    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    family_by_species: dict[str, str] = {}
    for row in links:
        input_name = clean(row.get("input_name"))
        resolved = resolution.get(input_name)
        if not resolved:
            continue
        accepted_name = clean(resolved.get("accepted_name"))
        family = clean(resolved.get("family"))
        key = source_key(row)
        if not key:
            continue
        family_by_species[accepted_name] = family
        grouped[(accepted_name, key)].append(row)

    pairs: list[tuple[dict[str, Any], dict[str, Any]]] = []
    for accepted_name, source in sorted(grouped):
        review, key = merge_group(
            accepted_name,
            family_by_species.get(accepted_name, ""),
            grouped[(accepted_name, source)],
        )
        stable = hashlib.sha256(f"{accepted_name}|{source}".encode("utf-8")).hexdigest()[:12]
        review_id = f"JSRC-{stable}"
        review["source_review_id"] = review_id
        key["source_review_id"] = review_id
        pairs.append((review, key))

    # Stable pseudo-random order prevents grouping by automated state or family while
    # preserving exactly the same review IDs across reruns.
    pairs.sort(key=lambda pair: pair[0]["source_review_id"])
    review_rows = [pair[0] for pair in pairs]
    key_rows = [pair[1] for pair in pairs]

    review_columns = set(review_rows[0])
    leaked = sorted(
        column for column in review_columns
        if column in FORBIDDEN_REVIEW_COLUMNS or column.startswith(FORBIDDEN_REVIEW_PREFIXES)
    )
    if leaked:
        raise SystemExit(f"Blinding failure: reviewer-facing columns leak navigation labels: {leaked}")

    write_csv(Path(args.review_out), review_rows)
    write_csv(Path(args.key_out), key_rows)

    mixed_navigation_species = {
        row["canonical_name"]
        for row in key_rows
        if row["automated_within_signal"] and row["automated_among_signal"]
    }
    species_with_within = {row["canonical_name"] for row in key_rows if row["automated_within_signal"]}
    species_with_among = {row["canonical_name"] for row in key_rows if row["automated_among_signal"]}
    print({
        "status": "complete",
        "source_review_rows": len(review_rows),
        "species": len({row['canonical_name'] for row in review_rows}),
        "species_with_automated_within_signal_hidden": len(species_with_within),
        "species_with_automated_among_signal_hidden": len(species_with_among),
        "species_with_same_source_mixed_signal_hidden": len(mixed_navigation_species),
        "reviewer_facing_automated_columns": 0,
        "semantic_guard": "Automated flags are coordinator-only and absent from the reviewer-facing sheet.",
    })


if __name__ == "__main__":
    main()

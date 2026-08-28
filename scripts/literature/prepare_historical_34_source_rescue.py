#!/usr/bin/env python3
"""Prepare direct blinded source-review rows for the historical 34-species freeze.

This rescue lane ignores archived systematic-search priority. Each historical species is
looked up by its frozen classification source_id. Reviewer-facing output hides the old
spatial label, automated signals, search priority and recovery status. Those fields are
written only to a coordinator key.

No historical label is validated or changed by this script.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path
from typing import Any


def clean(value: Any) -> str:
    return " ".join(str(value or "").split()).strip()


def normalize_doi(value: Any) -> str:
    text = clean(value).lower()
    text = re.sub(r"^https?://(?:dx\.)?doi\.org/", "", text)
    text = re.sub(r"^doi:\s*", "", text)
    return text.rstrip(".")


def normalize_openalex(value: Any) -> str:
    text = clean(value).rstrip("/").lower()
    if "openalex.org/" in text:
        return text.rsplit("/", 1)[-1]
    if re.fullmatch(r"w\d+", text):
        return text
    return ""


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def source_key_variants(row: dict[str, str]) -> set[str]:
    keys: set[str] = set()
    doi = normalize_doi(row.get("doi"))
    if doi:
        keys.add("doi:" + doi)
    oa = normalize_openalex(row.get("record_id"))
    if oa:
        keys.add("oa:" + oa)
    return keys


def historical_key(source_id: str) -> str:
    doi = normalize_doi(source_id)
    if doi.startswith("10."):
        return "doi:" + doi
    oa = normalize_openalex(source_id)
    if oa:
        return "oa:" + oa
    return "raw:" + clean(source_id).lower()


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
    parser.add_argument("--manifest", default="docs/supporting/frozen_classification_manifest.csv")
    parser.add_argument("--systematic-queue", required=True)
    parser.add_argument("--review-out", required=True)
    parser.add_argument("--key-out", required=True)
    parser.add_argument("--summary-out", required=True)
    args = parser.parse_args()

    manifest = read_csv(Path(args.manifest))
    systematic = read_csv(Path(args.systematic_queue))
    if len(manifest) != 34 or len({row.get('canonical_name') for row in manifest}) != 34:
        raise SystemExit("Historical manifest must contain exactly 34 unique species")

    by_key: dict[str, list[dict[str, str]]] = {}
    for row in systematic:
        for key in source_key_variants(row):
            by_key.setdefault(key, []).append(row)

    review_rows: list[dict[str, Any]] = []
    key_rows: list[dict[str, Any]] = []
    exact_matches = 0
    missing_sources: list[str] = []
    non_p1_exact = 0

    priority_rank = {
        "P1_high_natural_itv": 0,
        "P1_high_population_itv": 1,
        "P2_possible_itv": 2,
        "P3_likely_exclusion_review": 3,
        "P4_flower_colour_context_only": 4,
        "P5_low_relevance": 5,
    }

    for row in sorted(manifest, key=lambda x: clean(x.get("canonical_name"))):
        species = clean(row.get("canonical_name"))
        source_id = clean(row.get("source_id"))
        key = historical_key(source_id)
        matches = by_key.get(key, [])
        if matches:
            exact_matches += 1
            best = min(
                matches,
                key=lambda x: (
                    priority_rank.get(clean(x.get("screen_priority")), 99),
                    -len(clean(x.get("abstract"))),
                ),
            )
            priority = clean(best.get("screen_priority"))
            if not priority.startswith("P1_"):
                non_p1_exact += 1
            recovery = "exact_source_recovered"
        else:
            best = {}
            priority = ""
            recovery = "source_not_in_archived_systematic_corpus"
            missing_sources.append(species)

        stable = hashlib.sha256(f"{species}|{source_id}".encode("utf-8")).hexdigest()[:12]
        review_id = f"JHIST-{stable}"
        reviewer_source_id = clean(best.get("doi") or best.get("record_id") or source_id)
        review_rows.append({
            "historical_review_id": review_id,
            "canonical_name": species,
            "family": clean(row.get("family")),
            "source_id": reviewer_source_id,
            "doi": clean(best.get("doi")),
            "record_id": clean(best.get("record_id")),
            "title": clean(best.get("title")),
            "navigation_excerpt": clean(best.get("abstract"))[:2400],
            **review_fields(),
        })
        key_rows.append({
            "historical_review_id": review_id,
            "canonical_name": species,
            "family": clean(row.get("family")),
            "historical_spatial_scale": clean(row.get("spatial_scale")),
            "historical_source_id": source_id,
            "historical_source_key": key,
            "systematic_source_recovery": recovery,
            "systematic_screen_priority": priority,
            "systematic_record_id": clean(best.get("record_id")),
            "systematic_doi": clean(best.get("doi")),
            "systematic_title": clean(best.get("title")),
            "automated_within_signal": clean(best.get("within_signal")),
            "automated_among_signal": clean(best.get("among_signal")),
            "automated_natural_signal": clean(best.get("natural_signal")),
            "automated_cultivated_signal": clean(best.get("cultivated_signal")),
            "automated_induced_signal": clean(best.get("induced_signal")),
            "automated_ontogenetic_signal": clean(best.get("ontogenetic_signal")),
            "audit_guard": "Coordinator-only key; historical labels and automated signals must not be shown during independent review.",
        })

    review_rows.sort(key=lambda row: row["historical_review_id"])
    key_rows.sort(key=lambda row: row["historical_review_id"])

    forbidden = {
        "historical_spatial_scale", "historical_source_id", "systematic_screen_priority",
        "automated_within_signal", "automated_among_signal", "automated_natural_signal",
    }
    leaked = sorted(forbidden.intersection(review_rows[0].keys()))
    if leaked:
        raise SystemExit(f"Historical blinding failure: {leaked}")

    write_csv(Path(args.review_out), review_rows)
    write_csv(Path(args.key_out), key_rows)

    summary = {
        "status": "complete",
        "historical_species": 34,
        "exact_historical_sources_recovered_in_systematic_corpus": exact_matches,
        "exact_historical_sources_outside_p1": non_p1_exact,
        "historical_sources_not_in_archived_systematic_corpus": len(missing_sources),
        "missing_source_species": missing_sources,
        "reviewer_facing_historical_label_columns": 0,
        "semantic_guard": "Historical source rescue is priority-independent; old labels and automated signals are coordinator-only.",
    }
    Path(args.summary_out).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)

    if len(review_rows) != 34 or len(key_rows) != 34:
        raise SystemExit("Historical source review materials did not retain all 34 species")


if __name__ == "__main__":
    main()

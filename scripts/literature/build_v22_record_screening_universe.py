#!/usr/bin/env python3
"""Build the canonical record-level blind screening universe from the v2.2 search.

Every deduplicated v2.2 record is retained. No record is excluded because a binomial
cannot be extracted automatically. Search-query membership, automated navigation hints,
and historical benchmark identity are coordinator-only and are never written to the
reviewer-facing sheet.

This is deliberately upstream of taxon validation and spatial-state coding.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path
from typing import Any

BINOMIAL_IN_TEXT = re.compile(r"\b([A-Z][a-z-]{2,})\s+([a-z][a-z-]{2,})\b")


def clean(value: Any) -> str:
    return " ".join(str(value or "").split()).strip()


def normalize_doi(value: Any) -> str:
    text = clean(value).lower()
    text = re.sub(r"^https?://(?:dx\.)?doi\.org/", "", text)
    text = re.sub(r"^doi:\s*", "", text)
    return text.rstrip(".")


def normalize_openalex(value: Any) -> str:
    text = clean(value).rstrip("/")
    if "openalex.org/" in text.lower():
        return text.rsplit("/", 1)[-1].upper()
    if re.fullmatch(r"W\d+", text, re.I):
        return text.upper()
    return ""


def source_key(row: dict[str, str]) -> str:
    doi = normalize_doi(row.get("doi"))
    if doi:
        return "doi:" + doi
    oa = normalize_openalex(row.get("record_id"))
    if oa:
        return "openalex:" + oa
    dedup = clean(row.get("dedup_key"))
    if dedup:
        return "dedup:" + dedup
    return "title:" + clean(row.get("title")).lower()


def review_id(row: dict[str, str]) -> str:
    digest = hashlib.sha256(source_key(row).encode("utf-8")).hexdigest()[:16]
    return "V22R-" + digest.upper()


def detected_binomials(row: dict[str, str]) -> list[str]:
    text = f"{row.get('title', '')} {row.get('abstract', '')}"
    seen: set[str] = set()
    out: list[str] = []
    for genus, epithet in BINOMIAL_IN_TEXT.findall(text):
        name = f"{genus} {epithet}"
        if name not in seen:
            seen.add(name)
            out.append(name)
    return out[:25]


def read_historical_keys(path: Path) -> set[str]:
    keys: set[str] = set()
    if not path.exists():
        return keys
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            sid = clean(row.get("source_id"))
            doi = normalize_doi(sid)
            if doi.startswith("10."):
                keys.add("doi:" + doi)
                continue
            oa = normalize_openalex(sid)
            if oa:
                keys.add("openalex:" + oa)
    return keys


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--works", required=True)
    parser.add_argument("--historical-manifest", default="docs/supporting/frozen_classification_manifest.csv")
    parser.add_argument("--blind-out", required=True)
    parser.add_argument("--key-out", required=True)
    parser.add_argument("--summary-out", required=True)
    parser.add_argument("--batch-size", type=int, default=1000)
    args = parser.parse_args()

    with Path(args.works).open(newline="", encoding="utf-8") as handle:
        works = list(csv.DictReader(handle))
    if not works:
        raise SystemExit("No v2.2 works supplied")

    historical_keys = read_historical_keys(Path(args.historical_manifest))
    ids = [review_id(row) for row in works]
    if len(set(ids)) != len(ids):
        raise SystemExit("Record review IDs are not unique")

    # Hash order makes batch allocation deterministic but unrelated to query membership,
    # publication year, historical status, or any biological navigation signal.
    order = sorted(range(len(works)), key=lambda i: hashlib.sha256(ids[i].encode("ascii")).hexdigest())
    batch_by_index: dict[int, str] = {}
    size = max(1, args.batch_size)
    for position, index in enumerate(order):
        batch_by_index[index] = f"B{position // size + 1:02d}"

    blind_rows: list[dict[str, Any]] = []
    key_rows: list[dict[str, Any]] = []
    historical_hits = 0
    records_with_detected_binomial = 0

    for index, row in enumerate(works):
        rid = ids[index]
        skey = source_key(row)
        bins = detected_binomials(row)
        if bins:
            records_with_detected_binomial += 1
        hist = int(skey in historical_keys)
        historical_hits += hist

        blind_rows.append({
            "record_review_id": rid,
            "batch_id": batch_by_index[index],
            "source_id": clean(row.get("doi") or row.get("record_id")),
            "title": clean(row.get("title")),
            "abstract": clean(row.get("abstract")),
            "year": clean(row.get("year")),
            "journal": clean(row.get("journal")),
            "work_type": clean(row.get("work_type")),
            "language": clean(row.get("language")),
            "url": clean(row.get("url")),
            "review_status": "unreviewed",
            "reviewer_1_record_relevance": "",
            "reviewer_1_natural_intraspecific_variation": "",
            "reviewer_1_floral_display_colour": "",
            "reviewer_1_focal_taxon_text": "",
            "reviewer_1_full_text_required": "",
            "reviewer_1_exclusion_reason": "",
            "reviewer_1_notes": "",
            "reviewer_2_record_relevance": "",
            "reviewer_2_natural_intraspecific_variation": "",
            "reviewer_2_floral_display_colour": "",
            "reviewer_2_focal_taxon_text": "",
            "reviewer_2_full_text_required": "",
            "reviewer_2_exclusion_reason": "",
            "reviewer_2_notes": "",
            "adjudicated_record_relevance": "",
            "adjudicated_natural_intraspecific_variation": "",
            "adjudicated_floral_display_colour": "",
            "adjudicated_focal_taxon_text": "",
            "adjudicated_full_text_required": "",
            "adjudicated_exclusion_reason": "",
            "adjudication_notes": "",
        })

        key_rows.append({
            "record_review_id": rid,
            "source_key": skey,
            "dedup_key": clean(row.get("dedup_key")),
            "record_id": clean(row.get("record_id")),
            "doi": normalize_doi(row.get("doi")),
            "query_ids": clean(row.get("query_ids")),
            "representative_query_id": clean(row.get("query_id")),
            "cited_by_count": clean(row.get("cited_by_count")),
            "detected_binomial_strings": ";".join(bins),
            "historical_34_exact_source_match": hist,
        })

    blind_fields = list(blind_rows[0].keys())
    key_fields = list(key_rows[0].keys())
    write_csv(Path(args.blind_out), blind_rows, blind_fields)
    write_csv(Path(args.key_out), key_rows, key_fields)

    summary = {
        "status": "complete",
        "input_v22_deduplicated_works": len(works),
        "blind_record_screening_rows": len(blind_rows),
        "coordinator_key_rows": len(key_rows),
        "records_with_detected_binomial_navigation_hint_hidden": records_with_detected_binomial,
        "records_without_detected_binomial_retained_for_review": len(works) - records_with_detected_binomial,
        "historical_34_exact_source_rows_present": historical_hits,
        "review_batches": len(set(row["batch_id"] for row in blind_rows)),
        "batch_size_target": size,
        "reviewer_facing_query_membership_columns": 0,
        "reviewer_facing_historical_status_columns": 0,
        "reviewer_facing_automated_taxon_hint_columns": 0,
        "semantic_guard": (
            "All v2.2 records are retained for blind record screening. Absence of an automatically detected "
            "binomial cannot exclude a record. Taxon validation and spatial evidence coding occur only after "
            "record-level screening/adjudication."
        ),
    }
    Path(args.summary_out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.summary_out).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)

    if len(blind_rows) != 12064:
        raise SystemExit(f"Expected 12064 v2.2 records, found {len(blind_rows)}")
    if historical_hits != 34:
        raise SystemExit(f"Expected 34 historical exact source rows, found {historical_hits}")


if __name__ == "__main__":
    main()

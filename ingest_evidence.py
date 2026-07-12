#!/usr/bin/env python3
"""Merge reviewed evidence records into the angiosperm census.

Safety rule: absence of a literature hit never becomes monomorphism. Only rows
with an explicit reviewed evidence status can overwrite unknown fields.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, List

ALLOWED_OUTCOMES = {
    "unknown", "monomorphic_confirmed", "geographic_mosaic", "mixed",
    "local_coexistence", "excluded",
}
ALLOWED_EVIDENCE = {"unreviewed", "screening", "probable", "verified", "excluded"}


def read_rows(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def key(row: Dict[str, str]) -> str:
    if row.get("gbif_key"):
        return f"gbif:{row['gbif_key'].strip()}"
    name = row.get("canonical_name") or row.get("scientific_name") or row.get("species") or ""
    return "name:" + " ".join(name.strip().lower().split())


def validate_evidence(row: Dict[str, str]) -> None:
    outcome = (row.get("outcome_class") or "unknown").strip()
    status = (row.get("evidence_status") or "unreviewed").strip()
    if outcome not in ALLOWED_OUTCOMES:
        raise ValueError(f"invalid outcome_class: {outcome}")
    if status not in ALLOWED_EVIDENCE:
        raise ValueError(f"invalid evidence_status: {status}")
    if outcome == "monomorphic_confirmed" and status not in {"verified", "probable"}:
        raise ValueError("monomorphic_confirmed requires reviewed positive evidence, not absence of hits")
    if outcome in {"geographic_mosaic", "mixed", "local_coexistence"} and status == "unreviewed":
        raise ValueError("polymorphism outcomes cannot be unreviewed")


def merge(census: List[Dict[str, str]], evidence: List[Dict[str, str]]) -> List[Dict[str, str]]:
    index = {key(row): row for row in evidence}
    out: List[Dict[str, str]] = []
    editable = {
        "pollination_mode", "display_opportunity", "outcome_class",
        "evidence_status", "sampling_effort", "notes",
    }
    for row in census:
        merged = dict(row)
        ev = index.get(key(row))
        if ev:
            validate_evidence(ev)
            for field in editable:
                value = (ev.get(field) or "").strip()
                if value:
                    merged[field] = value
            if ev.get("source_doi"):
                merged["source_doi"] = ev["source_doi"].strip()
            if ev.get("source_url"):
                merged["source_url"] = ev["source_url"].strip()
        out.append(merged)
    return out


def write_rows(rows: List[Dict[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: List[str] = []
    seen = set()
    for row in rows:
        for field in row:
            if field not in seen:
                fields.append(field)
                seen.add(field)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--census", default="data/angiosperm_census.csv")
    parser.add_argument("--evidence", default="data/evidence_records.csv")
    parser.add_argument("--out", default="data/angiosperm_census_enriched.csv")
    args = parser.parse_args()
    merged = merge(read_rows(Path(args.census)), read_rows(Path(args.evidence)))
    write_rows(merged, Path(args.out))
    print(f"wrote {len(merged)} enriched census rows")


if __name__ == "__main__":
    main()

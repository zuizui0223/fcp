#!/usr/bin/env python3
"""Build a priority-independent blinded review universe from the OpenAlex v2.2 corpus.

This is the canonical upstream candidate-construction path for the JBI re-audit.
It deliberately does NOT read the legacy P1/P2/P3/P4/P5 screen priority and does NOT
assign within/among/mixed biological states.

Pipeline:
  v2.2 title/abstract retrieval -> focal binomial attribution -> strict GBIF Plantae
  species validation -> species x source blinded review rows.

OpenAlex query membership and taxonomic-resolution metadata are coordinator-only.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

from build_systematic_spatial_evidence_axes import resolve_gbif_name

BINOMIAL_RE = re.compile(r"\b([A-Z][a-z-]{2,})\s+([a-z][a-z-]{2,})\b")


def clean(value: Any) -> str:
    return " ".join(str(value or "").split()).strip()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def source_id(row: dict[str, str]) -> str:
    return clean(row.get("doi") or row.get("record_id") or row.get("dedup_key") or row.get("title"))


def candidates(text: str) -> list[str]:
    names: list[str] = []
    seen: set[str] = set()
    for genus, epithet in BINOMIAL_RE.findall(text or ""):
        name = f"{genus} {epithet}"
        if name not in seen:
            seen.add(name)
            names.append(name)
    return names


def focal_candidates(row: dict[str, str], abstract_chars: int = 2000) -> list[str]:
    """Return names plausibly attributable to the focal source, independent of priority.

    Prefer binomials in the title. If none occur there, use the first 2,000 abstract
    characters. Limit to five to avoid review papers turning into broad taxon censuses.
    """
    title_names = candidates(clean(row.get("title")))
    if title_names:
        return title_names[:5]
    return candidates(clean(row.get("abstract"))[:abstract_chars])[:5]


def resolve_one(name: str, timeout: int, retries: int) -> tuple[str, dict[str, Any]]:
    try:
        row = resolve_gbif_name(name, timeout, retries)
    except Exception as exc:  # noqa: BLE001
        row = {
            "input_name": name,
            "accepted": False,
            "match_type": "",
            "rank": "",
            "kingdom": "",
            "confidence": 0,
            "accepted_name": "",
            "family": "",
            "usage_key": "",
            "reason": f"resolution_error:{type(exc).__name__}",
        }
    return name, row


def reviewer_fields() -> dict[str, str]:
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
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise SystemExit(f"No rows for {path}")
    fields = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--works", required=True)
    p.add_argument("--historical-manifest", default="docs/supporting/frozen_classification_manifest.csv")
    p.add_argument("--blind-out", required=True)
    p.add_argument("--key-out", required=True)
    p.add_argument("--taxon-audit-out", required=True)
    p.add_argument("--summary-out", required=True)
    p.add_argument("--workers", type=int, default=24)
    p.add_argument("--timeout", type=int, default=25)
    p.add_argument("--retries", type=int, default=3)
    args = p.parse_args()

    works = read_csv(Path(args.works))
    if len(works) < 10000:
        raise SystemExit(f"Unexpectedly small v2.2 corpus: {len(works)}")

    links: list[dict[str, Any]] = []
    names: set[str] = set()
    records_with_focal = 0
    for row in works:
        focal = focal_candidates(row)
        if not focal:
            continue
        records_with_focal += 1
        sid = source_id(row)
        for name in focal:
            names.add(name)
            links.append({
                "input_name": name,
                "source_id": sid,
                "doi": clean(row.get("doi")),
                "record_id": clean(row.get("record_id")),
                "title": clean(row.get("title")),
                "abstract": clean(row.get("abstract")),
                "year": clean(row.get("year")),
                "journal": clean(row.get("journal")),
                "query_ids": clean(row.get("query_ids")),
            })

    if not names:
        raise SystemExit("No focal candidate binomials extracted from v2.2 corpus")

    resolution: dict[str, dict[str, Any]] = {}
    workers = max(1, min(args.workers, 48))
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(resolve_one, name, args.timeout, args.retries): name for name in sorted(names)}
        for index, future in enumerate(as_completed(futures), start=1):
            name, resolved = future.result()
            resolution[name] = resolved
            if index % 250 == 0:
                print({"resolved": index, "total": len(names)}, flush=True)

    accepted = {name: row for name, row in resolution.items() if bool(row.get("accepted")) and clean(row.get("accepted_name"))}
    errors = [row for row in resolution.values() if clean(row.get("reason")).startswith("resolution_error:")]

    taxon_rows = [resolution[name] for name in sorted(resolution)]
    write_csv(Path(args.taxon_audit_out), taxon_rows)

    grouped: dict[tuple[str, str], list[dict[str, Any]]] = {}
    family: dict[str, str] = {}
    for link in links:
        resolved = accepted.get(str(link["input_name"]))
        if not resolved:
            continue
        canonical = clean(resolved.get("accepted_name"))
        sid = clean(link.get("source_id"))
        if not canonical or not sid:
            continue
        family[canonical] = clean(resolved.get("family"))
        grouped.setdefault((canonical, sid), []).append(link)

    blind_rows: list[dict[str, Any]] = []
    key_rows: list[dict[str, Any]] = []
    for idx, ((canonical, sid), rows) in enumerate(sorted(grouped.items()), start=1):
        best = max(rows, key=lambda x: (len(clean(x.get("abstract"))), len(clean(x.get("title")))))
        rid = f"JV22-{idx:05d}"
        blind_rows.append({
            "source_review_id": rid,
            "canonical_name": canonical,
            "family": family.get(canonical, ""),
            "source_id": sid,
            "doi": clean(best.get("doi")),
            "record_id": clean(best.get("record_id")),
            "title": clean(best.get("title")),
            "year": clean(best.get("year")),
            "journal": clean(best.get("journal")),
            "navigation_excerpt": clean(best.get("abstract"))[:2400],
            **reviewer_fields(),
        })
        key_rows.append({
            "source_review_id": rid,
            "canonical_name": canonical,
            "source_id": sid,
            "input_names": ";".join(sorted({clean(x.get("input_name")) for x in rows})),
            "query_ids": ";".join(sorted({q for x in rows for q in clean(x.get("query_ids")).split(";") if q})),
            "gbif_family": family.get(canonical, ""),
            "coordinator_guard": "Query membership and taxon-resolution metadata are hidden from independent reviewers.",
        })

    write_csv(Path(args.blind_out), blind_rows)
    write_csv(Path(args.key_out), key_rows)

    historical = read_csv(Path(args.historical_manifest))
    historical_names = {clean(row.get("canonical_name")) for row in historical}
    candidate_species = {row["canonical_name"] for row in blind_rows}
    historical_recovered = sorted(historical_names.intersection(candidate_species))
    historical_missing = sorted(historical_names - candidate_species)

    summary = {
        "status": "complete",
        "input_v22_deduplicated_works": len(works),
        "records_with_primary_binomial_candidate": records_with_focal,
        "unique_primary_binomial_strings": len(names),
        "gbif_valid_input_names": len(accepted),
        "gbif_resolution_errors": len(errors),
        "gbif_valid_candidate_species_after_synonym_collapse": len(candidate_species),
        "species_source_blind_review_rows": len(blind_rows),
        "historical_34_species_represented_in_priority_independent_universe": len(historical_recovered),
        "historical_34_species_missing_from_priority_independent_universe": historical_missing,
        "reviewer_facing_query_membership_columns": 0,
        "reviewer_facing_legacy_priority_columns": 0,
        "semantic_guard": "This is a priority-independent taxon/source review universe. It does not infer natural eligibility or within/among/mixed state.",
    }
    Path(args.summary_out).write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)

    if len(errors) > max(20, int(0.03 * len(names))):
        raise SystemExit(f"Too many GBIF resolution errors: {len(errors)}/{len(names)}")
    if len(candidate_species) < 100:
        raise SystemExit(f"Unexpectedly few GBIF-valid candidate species: {len(candidate_species)}")
    if len(historical_missing) > 5:
        raise SystemExit(f"Unexpected benchmark species loss in v2.2 taxon/source universe: {historical_missing}")


if __name__ == "__main__":
    main()

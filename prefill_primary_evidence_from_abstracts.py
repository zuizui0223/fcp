#!/usr/bin/env python3
"""Conservatively prefill primary-evidence fields from title/abstract evidence.

This script never upgrades a record to fully verified. It only fills fields when
explicit phrases are present, and records every inferred field for later audit.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path
from typing import Dict, List, Tuple

COLOUR_WORDS = (
    "white", "pink", "purple", "blue", "yellow", "red", "orange", "green",
    "magenta", "scarlet", "cream", "violet", "black",
)

BOOL_RULES = {
    "within_population_polymorphism": [
        r"polymorphic population", r"within (?:a |the )?population", r"coexist in some populations",
        r"natural populations? .* (?:two|three|four|five|six) .* morph",
    ],
    "geographic_mosaic": [
        r"geographic(?:al)? (?:pattern|variation|distribution)", r"spatial variation",
        r"distinct geographic regions", r"disjunct distribution",
    ],
    "morph_frequency_reported": [
        r"morph frequenc", r"ratio of (?:flower )?morphs", r"frequency of .* morph",
    ],
    "frequency_dependent_selection": [
        r"frequency[- ]dependent selection", r"negative frequency[- ]dependent",
    ],
    "directional_selection": [r"directional selection"],
    "spatially_varying_selection": [
        r"spatial(?:ly)? var(?:ying|iation) selection", r"divergent selection .* (?:space|geograph)",
        r"regional selection",
    ],
    "temporally_varying_selection": [
        r"temporal(?:ly)? var(?:ying|iation) selection", r"spatiotemporal variation",
        r"fluctuating temperatures",
    ],
    "abiotic_selection": [
        r"abiotic", r"temperature", r"aridity", r"soil characteristics", r"climate",
    ],
    "herbivory_selection": [r"herbivor", r"seed predator"],
    "genetic_basis": [
        r"genetically controlled", r"genetic basis", r"loci controlling", r"segregation at .* loci",
    ],
    "environmental_basis": [
        r"environmental correlat", r"environmental basis", r"plasticity", r"temperature",
    ],
}


def read_rows(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def infer_colours(text: str) -> List[str]:
    found = []
    lower = text.lower()
    for colour in COLOUR_WORDS:
        if re.search(rf"\b{re.escape(colour)}(?:-flowered| flowers?| floral)?\b", lower):
            found.append(colour)
    return found


def infer_pollination(text: str) -> Tuple[str, str, str]:
    lower = text.lower()
    if "self-pollinating" in lower or "selfing" in lower:
        return "self", "none_or_incidental", "unknown"
    if "bumblebee" in lower or "bombus" in lower:
        return "insect", "bumblebee", "present"
    if "pollinator" in lower or "pollination" in lower:
        return "insect_or_mixed", "unknown", "unknown"
    return "unknown", "unknown", "unknown"


def infer_row(row: Dict[str, str]) -> Tuple[Dict[str, str], List[str]]:
    out = dict(row)
    text = " ".join((row.get("best_title", ""), row.get("review_notes", ""))).strip()
    lower = text.lower()
    inferred: List[str] = []

    for field, patterns in BOOL_RULES.items():
        if out.get(field, "unknown") not in {"", "unknown"}:
            continue
        if any(re.search(pattern, lower) for pattern in patterns):
            out[field] = "yes"
            inferred.append(field)

    colours = infer_colours(text)
    if colours and not out.get("colour_morphs"):
        out["colour_morphs"] = ";".join(colours)
        inferred.append("colour_morphs")
    if len(colours) >= 2 and not out.get("n_colour_morphs"):
        out["n_colour_morphs"] = str(len(colours))
        inferred.append("n_colour_morphs")

    pollination_mode, dominant_group, bumblebee_status = infer_pollination(text)
    for field, value in (
        ("pollination_mode", pollination_mode),
        ("dominant_pollinator_group", dominant_group),
        ("bumblebee_status", bumblebee_status),
    ):
        if value != "unknown" and out.get(field, "unknown") in {"", "unknown"}:
            out[field] = value
            inferred.append(field)

    if re.search(r"island|archipelago|islands", lower) and out.get("island_status", "unknown") in {"", "unknown"}:
        out["island_status"] = "island_or_archipelago"
        inferred.append("island_status")

    if re.search(r"wild population|natural population", lower) and out.get("native_or_introduced", "unknown") in {"", "unknown"}:
        out["native_or_introduced"] = "wild_population_reported"
        inferred.append("native_or_introduced")

    if inferred:
        if out.get("evidence_level") == "title_abstract_candidate":
            out["evidence_level"] = "abstract_prefilled_candidate"
        prior = out.get("review_notes", "")
        audit = "abstract_prefill=" + ";".join(sorted(set(inferred)))
        out["review_notes"] = (prior + " | " + audit).strip(" |")

    return out, sorted(set(inferred))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ledger", default="data/primary_evidence_ledger.csv")
    parser.add_argument("--out", default="data/primary_evidence_ledger_prefilled.csv")
    parser.add_argument("--qc-out", default="data/primary_evidence_prefill_qc.json")
    args = parser.parse_args()

    rows = read_rows(Path(args.ledger))
    if not rows:
        raise RuntimeError("Primary evidence ledger is empty")

    enriched: List[Dict[str, str]] = []
    field_counts: Counter[str] = Counter()
    species_with_prefill = 0
    for row in rows:
        new_row, fields = infer_row(row)
        enriched.append(new_row)
        if fields:
            species_with_prefill += 1
            field_counts.update(fields)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(enriched)

    qc = {
        "ledger_species": len(rows),
        "species_with_abstract_prefill": species_with_prefill,
        "species_still_requiring_primary_review": sum(
            row.get("verification_status") == "pending_primary_review" for row in enriched
        ),
        "prefilled_field_counts": dict(sorted(field_counts.items())),
        "verification_upgrades": 0,
        "mode": "conservative_title_abstract_prefill_v1",
    }
    Path(args.qc_out).write_text(json.dumps(qc, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(qc, ensure_ascii=False), flush=True)

    if len(enriched) != len(rows):
        raise RuntimeError("Prefill changed ledger row count")
    if species_with_prefill < 10:
        raise RuntimeError(f"Abstract prefill unexpectedly weak: {species_with_prefill} species")


if __name__ == "__main__":
    main()

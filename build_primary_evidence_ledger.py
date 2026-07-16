#!/usr/bin/env python3
"""Build a structured primary-evidence ledger from the natural-polymorphism review queue."""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from typing import Dict, List

FIELDS = [
    "review_priority", "canonical_name", "family", "best_title", "best_doi",
    "best_openalex_id", "evidence_class", "verification_status",
    "outcome_state", "within_population_polymorphism", "geographic_mosaic",
    "n_colour_morphs", "colour_morphs", "morph_frequency_reported",
    "morph_frequency_data", "population_count", "sample_size",
    "pollination_mode", "dominant_pollinator_group", "bumblebee_status",
    "bumblebee_visitation_share", "pollinator_richness",
    "frequency_dependent_selection", "directional_selection",
    "spatially_varying_selection", "temporally_varying_selection",
    "abiotic_selection", "herbivory_selection", "mating_system",
    "self_compatibility", "genetic_basis", "environmental_basis",
    "study_region", "latitude", "longitude", "island_status",
    "native_or_introduced", "evidence_level", "full_text_status",
    "review_notes", "reviewer", "review_date",
]

KNOWN_POSITIVES = {
    "Dactylorhiza sambucina",
    "Gymnadenia rhellicani",
    "Linanthus parryae",
    "Silene littorea",
}


def read_rows(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def inferred_row(row: Dict[str, str]) -> Dict[str, str]:
    name = row.get("canonical_name", "").strip()
    evidence = row.get("evidence_class", "natural_polymorphism").strip() or "natural_polymorphism"
    known = name in KNOWN_POSITIVES
    return {
        "review_priority": row.get("review_priority", ""),
        "canonical_name": name,
        "family": row.get("family", ""),
        "best_title": row.get("best_title", ""),
        "best_doi": row.get("best_doi", ""),
        "best_openalex_id": row.get("best_openalex_id", ""),
        "evidence_class": evidence,
        "verification_status": "seed_verified" if known else "pending_primary_review",
        "outcome_state": "unknown",
        "within_population_polymorphism": "yes" if known else "unknown",
        "geographic_mosaic": "unknown",
        "n_colour_morphs": "",
        "colour_morphs": "",
        "morph_frequency_reported": "unknown",
        "morph_frequency_data": "",
        "population_count": "",
        "sample_size": "",
        "pollination_mode": "unknown",
        "dominant_pollinator_group": "unknown",
        "bumblebee_status": "unknown",
        "bumblebee_visitation_share": "",
        "pollinator_richness": "",
        "frequency_dependent_selection": "unknown",
        "directional_selection": "unknown",
        "spatially_varying_selection": "unknown",
        "temporally_varying_selection": "unknown",
        "abiotic_selection": "unknown",
        "herbivory_selection": "unknown",
        "mating_system": "unknown",
        "self_compatibility": "unknown",
        "genetic_basis": "unknown",
        "environmental_basis": "unknown",
        "study_region": "",
        "latitude": "",
        "longitude": "",
        "island_status": "unknown",
        "native_or_introduced": "unknown",
        "evidence_level": "verified_seed" if known else "title_abstract_candidate",
        "full_text_status": "needed",
        "review_notes": row.get("best_match_evidence", ""),
        "reviewer": "",
        "review_date": "",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--queue", default="data/global_flower_colour_review_queue.csv")
    parser.add_argument("--out", default="data/primary_evidence_ledger.csv")
    parser.add_argument("--qc-out", default="data/primary_evidence_ledger_qc.json")
    args = parser.parse_args()

    rows = read_rows(Path(args.queue))
    if not rows:
        raise RuntimeError("Natural-polymorphism review queue is empty")

    ledger = [inferred_row(row) for row in rows]
    names = [row["canonical_name"] for row in ledger]
    if len(names) != len(set(names)):
        raise RuntimeError("Duplicate species in primary evidence ledger")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(ledger)

    status_counts = Counter(row["verification_status"] for row in ledger)
    full_text_counts = Counter(row["full_text_status"] for row in ledger)
    required_review = sum(row["verification_status"] == "pending_primary_review" for row in ledger)
    qc = {
        "ledger_species": len(ledger),
        "verified_seed_species": status_counts.get("seed_verified", 0),
        "pending_primary_review": required_review,
        "verification_status_counts": dict(sorted(status_counts.items())),
        "full_text_status_counts": dict(sorted(full_text_counts.items())),
        "core_fields_per_species": 30,
        "estimated_review_minutes_per_species": {"minimum": 8, "likely": 15, "complex": 25},
        "ledger_mode": "species_primary_evidence_v1",
    }
    Path(args.qc_out).write_text(json.dumps(qc, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(qc, ensure_ascii=False), flush=True)

    if len(ledger) < 10:
        raise RuntimeError(f"Ledger unexpectedly small: {len(ledger)}")
    if status_counts.get("seed_verified", 0) < 4:
        raise RuntimeError("Known positive controls missing from ledger")


if __name__ == "__main__":
    main()

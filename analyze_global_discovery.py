#!/usr/bin/env python3
"""Validate and prioritize context-filtered flower-colour discovery outputs offline."""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from typing import Dict, List

KNOWN_POSITIVES = {
    "Dactylorhiza sambucina",
    "Gymnadenia rhellicani",
    "Silene littorea",
    "Linanthus parryae",
}
SHORTLIST_FIELDS = [
    "priority_tier", "rank", "canonical_name", "family", "n_works",
    "max_score", "total_score", "n_title_matches", "n_context_matches",
    "best_title", "best_doi", "best_openalex_id", "best_match_evidence",
    "review_status",
]


def read_csv(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        raise RuntimeError(f"Missing required file: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def to_int(value: str) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return 0


def priority_tier(row: Dict[str, str]) -> str:
    n_works = to_int(row.get("n_works", "0"))
    max_score = to_int(row.get("max_score", "0"))
    title_matches = to_int(row.get("n_title_matches", "0"))
    context_matches = to_int(row.get("n_context_matches", "0"))
    if title_matches >= 1 and max_score >= 18:
        return "A"
    if title_matches >= 1 or context_matches >= 2 or (n_works >= 2 and max_score >= 14):
        return "B"
    return "C"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--works", default="data/global_flower_colour_works.csv")
    parser.add_argument("--species", default="data/global_flower_colour_species_ranked.csv")
    parser.add_argument("--shortlist-out", default="data/global_flower_colour_shortlist.csv")
    parser.add_argument("--qc-out", default="data/global_flower_colour_qc.json")
    parser.add_argument("--top", type=int, default=250)
    args = parser.parse_args()

    works = read_csv(Path(args.works))
    species = read_csv(Path(args.species))
    if not works or not species:
        raise RuntimeError("Discovery outputs are empty; refusing to report success")

    work_ids = [row.get("openalex_id", "").strip() for row in works]
    duplicate_work_ids = sum(count - 1 for count in Counter(work_ids).values() if count > 1)
    species_names = [row.get("canonical_name", "").strip() for row in species]
    duplicate_species = sum(count - 1 for count in Counter(species_names).values() if count > 1)
    if duplicate_work_ids or duplicate_species:
        raise RuntimeError(
            f"Duplicate rows detected: work_ids={duplicate_work_ids}, species={duplicate_species}"
        )

    enriched: List[Dict[str, str]] = []
    for row in species:
        enriched.append({"priority_tier": priority_tier(row), **row})
    enriched.sort(key=lambda row: (
        {"A": 0, "B": 1, "C": 2}[row["priority_tier"]],
        -to_int(row.get("n_title_matches", "0")),
        -to_int(row.get("max_score", "0")),
        -to_int(row.get("total_score", "0")),
        row.get("canonical_name", ""),
    ))
    shortlist = enriched[: max(1, args.top)]

    shortlist_path = Path(args.shortlist_out)
    shortlist_path.parent.mkdir(parents=True, exist_ok=True)
    with shortlist_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SHORTLIST_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(shortlist)

    families = Counter(row.get("family", "").strip() or "unknown" for row in species)
    tiers = Counter(row["priority_tier"] for row in enriched)
    positive_controls = sorted(KNOWN_POSITIVES.intersection(species_names))
    raw_mentions = [to_int(row.get("raw_species_mentions", "0")) for row in works]
    title_supported_species = sum(to_int(row.get("n_title_matches", "0")) > 0 for row in species)
    context_supported_species = sum(to_int(row.get("n_context_matches", "0")) > 0 for row in species)
    qc = {
        "retained_works": len(works),
        "candidate_species": len(species),
        "unique_families": len(families),
        "priority_tiers": dict(sorted(tiers.items())),
        "title_supported_species": title_supported_species,
        "context_supported_species": context_supported_species,
        "works_with_more_than_5_raw_species_mentions": sum(value > 5 for value in raw_mentions),
        "max_raw_species_mentions_in_retained_work": max(raw_mentions, default=0),
        "top_families": families.most_common(20),
        "duplicate_work_ids": duplicate_work_ids,
        "duplicate_species": duplicate_species,
        "known_positive_controls_recovered": positive_controls,
        "known_positive_controls_recovered_count": len(positive_controls),
        "shortlist_rows": len(shortlist),
        "ranking_mode": "species_specific_context_v2",
    }
    qc_path = Path(args.qc_out)
    qc_path.write_text(json.dumps(qc, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(qc, ensure_ascii=False), flush=True)
    if len(positive_controls) < 3:
        raise RuntimeError(
            f"Context filtering recovered only {len(positive_controls)}/4 positive controls; thresholds need review"
        )


if __name__ == "__main__":
    main()

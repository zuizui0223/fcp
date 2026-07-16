#!/usr/bin/env python3
"""Validate and prioritize global flower-colour discovery outputs offline."""
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
    "max_score", "total_score", "best_title", "best_doi",
    "best_openalex_id", "review_status",
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
    total_score = to_int(row.get("total_score", "0"))
    if n_works >= 2 or max_score >= 14 or total_score >= 20:
        return "A"
    if max_score >= 10 or total_score >= 12:
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
        to_int(row.get("rank", "0")) or 10**9,
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
    qc = {
        "mapped_works": len(works),
        "candidate_species": len(species),
        "unique_families": len(families),
        "priority_tiers": dict(sorted(tiers.items())),
        "top_families": families.most_common(20),
        "duplicate_work_ids": duplicate_work_ids,
        "duplicate_species": duplicate_species,
        "known_positive_controls_recovered": positive_controls,
        "known_positive_controls_recovered_count": len(positive_controls),
        "shortlist_rows": len(shortlist),
    }
    qc_path = Path(args.qc_out)
    qc_path.write_text(json.dumps(qc, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(qc, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()

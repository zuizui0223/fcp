#!/usr/bin/env python3
"""Aggregate raw mass-screen hits into a species-level ranked review table."""
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path
from typing import Dict, List

FIELDS = [
    "rank", "canonical_name", "family", "max_score", "total_score",
    "n_unique_hits", "n_unique_dois", "best_title", "best_doi",
    "best_source", "review_status",
]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hits", default="data/mass_screen_hits.csv")
    parser.add_argument("--out", default="data/mass_screen_species_ranked.csv")
    parser.add_argument("--top", type=int, default=10000)
    args = parser.parse_args()

    grouped: Dict[str, List[dict]] = defaultdict(list)
    with Path(args.hits).open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            name = (row.get("canonical_name") or "").strip()
            if name:
                grouped[name].append(row)

    ranked = []
    for name, rows in grouped.items():
        unique = {}
        for row in rows:
            key = (row.get("doi") or row.get("title") or "").strip().lower()
            if key and key not in unique:
                unique[key] = row
        vals = list(unique.values())
        vals.sort(key=lambda r: int(r.get("candidate_score") or 0), reverse=True)
        best = vals[0]
        scores = [int(r.get("candidate_score") or 0) for r in vals]
        dois = {r.get("doi", "").strip().lower() for r in vals if r.get("doi", "").strip()}
        ranked.append({
            "canonical_name": name,
            "family": best.get("family", ""),
            "max_score": max(scores),
            "total_score": sum(scores),
            "n_unique_hits": len(vals),
            "n_unique_dois": len(dois),
            "best_title": best.get("title", ""),
            "best_doi": best.get("doi", ""),
            "best_source": best.get("source", ""),
            "review_status": "unreviewed_machine_candidate",
        })

    ranked.sort(
        key=lambda r: (
            int(r["max_score"]), int(r["total_score"]), int(r["n_unique_hits"])
        ),
        reverse=True,
    )
    ranked = ranked[: max(1, args.top)]
    for i, row in enumerate(ranked, 1):
        row["rank"] = i

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(ranked)

    print(f"wrote {len(ranked)} ranked species to {out}")


if __name__ == "__main__":
    main()

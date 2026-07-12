#!/usr/bin/env python3
"""Prioritize an all-angiosperm census for evidence review without changing inclusion.

The census remains complete. This script only assigns review priority so scarce
manual effort can be focused on informative taxa while preserving a known
sampling universe for unbiased analyses.
"""
from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path
from typing import Dict, List

WIND_FAMILIES = {
    "Poaceae", "Cyperaceae", "Juncaceae", "Betulaceae", "Fagaceae",
    "Juglandaceae", "Urticaceae", "Cannabaceae", "Casuarinaceae",
    "Myricaceae", "Nothofagaceae",
}
KNOWN_RICH_FAMILIES = {
    "Orchidaceae", "Asteraceae", "Fabaceae", "Lamiaceae", "Ericaceae",
    "Caryophyllaceae", "Polemoniaceae", "Plantaginaceae", "Primulaceae",
    "Brassicaceae", "Solanaceae", "Ranunculaceae", "Iridaceae",
    "Onagraceae", "Boraginaceae", "Rubiaceae", "Phrymaceae",
}
SPECIALIZED_POLLINATION_FAMILIES = {
    "Bromeliaceae", "Heliconiaceae", "Strelitziaceae", "Cactaceae",
    "Gesneriaceae", "Bignoniaceae", "Proteaceae", "Melastomataceae",
    "Aristolochiaceae", "Magnoliaceae", "Nymphaeaceae",
}


def priority(row: Dict[str, str]) -> str:
    fam = row.get("family", "")
    if fam in KNOWN_RICH_FAMILIES:
        return "A_colour_polymorphism_rich"
    if fam in WIND_FAMILIES:
        return "B_wind_control"
    if fam in SPECIALIZED_POLLINATION_FAMILIES:
        return "C_specialized_pollination"
    return "D_background"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="data/angiosperm_census.csv")
    parser.add_argument("--out", default="data/angiosperm_census_prioritized.csv")
    args = parser.parse_args()

    src, dst = Path(args.input), Path(args.out)
    with src.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows: List[Dict[str, str]] = list(reader)
        fields = list(reader.fieldnames or [])
    if "review_priority" not in fields:
        fields.append("review_priority")

    counts = Counter()
    for row in rows:
        row["review_priority"] = priority(row)
        counts[row["review_priority"]] += 1

    dst.parent.mkdir(parents=True, exist_ok=True)
    with dst.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    for key, value in sorted(counts.items()):
        print(f"{key}: {value}")


if __name__ == "__main__":
    main()

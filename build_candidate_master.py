"""Build a deduplicated candidate master table from screening waves.

This script does not upgrade evidence status. It only combines candidate tables,
keeps the strongest available evidence label per species, and reports conflicts
for manual review.
"""
from __future__ import annotations

import csv
from pathlib import Path

INPUTS = [
    Path("data/candidate_screening_expanded.csv"),
    Path("data/candidate_screening_wave2.csv"),
]
OUTPUT = Path("data/candidate_screening_master.csv")
CONFLICTS = Path("data/candidate_screening_conflicts.csv")

EVIDENCE_RANK = {"verified": 3, "probable": 2, "screening": 1, "excluded": 0}
PRIORITY_RANK = {"high": 3, "medium": 2, "low": 1}


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def choose_preferred(a: dict[str, str], b: dict[str, str]) -> dict[str, str]:
    ea = EVIDENCE_RANK.get(a.get("evidence_status", ""), -1)
    eb = EVIDENCE_RANK.get(b.get("evidence_status", ""), -1)
    if eb > ea:
        return b
    if ea > eb:
        return a
    pa = PRIORITY_RANK.get(a.get("priority", ""), -1)
    pb = PRIORITY_RANK.get(b.get("priority", ""), -1)
    return b if pb > pa else a


def main() -> None:
    rows: list[dict[str, str]] = []
    for path in INPUTS:
        if not path.exists():
            raise FileNotFoundError(path)
        rows.extend(read_rows(path))

    by_species: dict[str, dict[str, str]] = {}
    conflicts: list[dict[str, str]] = []

    for row in rows:
        species = row["species"].strip()
        if species not in by_species:
            by_species[species] = row
            continue
        old = by_species[species]
        compared_fields = [
            "family",
            "broad_clade",
            "evidence_status",
            "pollination_mode_to_verify",
            "outcome_to_verify",
        ]
        differing = [f for f in compared_fields if old.get(f, "") != row.get(f, "")]
        if differing:
            conflicts.append(
                {
                    "species": species,
                    "differing_fields": ";".join(differing),
                    "record_a": repr(old),
                    "record_b": repr(row),
                }
            )
        by_species[species] = choose_preferred(old, row)

    master = sorted(by_species.values(), key=lambda r: (r.get("family", ""), r["species"]))
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=master[0].keys())
        writer.writeheader()
        writer.writerows(master)

    with CONFLICTS.open("w", newline="", encoding="utf-8") as handle:
        fields = ["species", "differing_fields", "record_a", "record_b"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(conflicts)

    counts: dict[str, int] = {}
    for row in master:
        status = row.get("evidence_status", "unknown")
        counts[status] = counts.get(status, 0) + 1

    print(f"input rows: {len(rows)}")
    print(f"unique species: {len(master)}")
    print(f"conflicts: {len(conflicts)}")
    print("evidence counts:", counts)
    print(f"wrote {OUTPUT}")
    print(f"wrote {CONFLICTS}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Build deterministic review batches from the full angiosperm census.

The queue keeps known/candidate colour-polymorphism taxa high priority while also
sampling background taxa so that review does not collapse into case-only
ascertainment. Unknown is preserved as unknown until evidence is reviewed.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path
from typing import Dict, Iterable, List, Set


def read_rows(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def canonical(text: str) -> str:
    return " ".join((text or "").strip().split()).lower()


def stable_score(text: str, seed: str) -> int:
    digest = hashlib.sha256(f"{seed}|{text}".encode("utf-8")).hexdigest()
    return int(digest[:16], 16)


def candidate_names(rows: Iterable[Dict[str, str]]) -> Set[str]:
    return {canonical(row.get("species", "")) for row in rows if row.get("species")}


def priority(row: Dict[str, str], candidates: Set[str]) -> int:
    name = canonical(row.get("canonical_name") or row.get("scientific_name", ""))
    if name in candidates:
        return 0
    status = (row.get("evidence_status") or "").strip().lower()
    if status not in {"", "unreviewed", "unknown"}:
        return 1
    family = (row.get("family") or "").strip()
    # Families with classic colour-polymorphism systems receive an early pass,
    # but all other families remain eligible through the background lane.
    focal = {
        "Orchidaceae", "Caryophyllaceae", "Polemoniaceae", "Ericaceae",
        "Iridaceae", "Primulaceae", "Plantaginaceae", "Brassicaceae",
        "Fabaceae", "Asteraceae", "Ranunculaceae", "Solanaceae",
    }
    if family in focal:
        return 2
    return 3


def build_queue(census: List[Dict[str, str]], candidates: Set[str], seed: str,
                batch_size: int, max_batches: int | None) -> List[Dict[str, str]]:
    ranked = []
    for row in census:
        name = row.get("canonical_name") or row.get("scientific_name", "")
        ranked.append((priority(row, candidates), stable_score(name, seed), row))
    ranked.sort(key=lambda x: (x[0], x[1]))

    limit = len(ranked) if max_batches is None else min(len(ranked), batch_size * max_batches)
    output: List[Dict[str, str]] = []
    for i, (tier, _, row) in enumerate(ranked[:limit]):
        out = dict(row)
        out["review_priority"] = str(tier)
        out["review_batch"] = str(i // batch_size + 1)
        out["review_status"] = "pending"
        output.append(out)
    return output


def write_rows(rows: List[Dict[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError("review queue is empty")
    fields = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--census", default="data/angiosperm_census.csv")
    parser.add_argument("--candidates", default="data/candidate_screening_expanded.csv")
    parser.add_argument("--out", default="data/review_queue.csv")
    parser.add_argument("--seed", default="fcp-review-v1")
    parser.add_argument("--batch-size", type=int, default=250)
    parser.add_argument("--max-batches", type=int, default=None)
    args = parser.parse_args()

    census = read_rows(Path(args.census))
    candidates = candidate_names(read_rows(Path(args.candidates)))
    queue = build_queue(census, candidates, args.seed, args.batch_size, args.max_batches)
    write_rows(queue, Path(args.out))
    print(f"wrote {len(queue)} taxa in {(len(queue) + args.batch_size - 1) // args.batch_size} batches")


if __name__ == "__main__":
    main()

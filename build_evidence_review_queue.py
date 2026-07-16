#!/usr/bin/env python3
"""Build a manually reviewable evidence queue from context-filtered discovery outputs."""
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
    "Linanthus parryae",
    "Silene littorea",
}

FIELDS = [
    "review_priority", "canonical_name", "family", "n_works", "n_title_matches",
    "n_context_matches", "max_score", "total_score", "best_title", "best_doi",
    "best_openalex_id", "best_match_evidence", "review_status", "review_reason",
]


def read_rows(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def to_int(value: object) -> int:
    try:
        return int(str(value or "0"))
    except ValueError:
        return 0


def classify(row: Dict[str, str]) -> tuple[str, str]:
    title = to_int(row.get("n_title_matches"))
    context = to_int(row.get("n_context_matches"))
    works = to_int(row.get("n_works"))
    max_score = to_int(row.get("max_score"))
    name = row.get("canonical_name", "")

    if name in KNOWN_POSITIVES:
        return "P0", "known_positive_control"
    if title >= 1 and (works >= 2 or context >= 1 or max_score >= 18):
        return "P1", "direct_title_support_plus_replication_or_context"
    if title >= 1 and max_score >= 14:
        return "P2", "direct_title_support"
    if context >= 2 and works >= 2 and max_score >= 14:
        return "P3", "replicated_local_context_support"
    return "DEFER", "insufficient_direct_support"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", default="data/global_flower_colour_species_ranked.csv")
    parser.add_argument("--out", default="data/global_flower_colour_review_queue.csv")
    parser.add_argument("--qc-out", default="data/global_flower_colour_review_queue_qc.json")
    parser.add_argument("--max-review", type=int, default=150)
    args = parser.parse_args()

    rows = read_rows(Path(args.species))
    if not rows:
        raise RuntimeError("No candidate species available for review queue")

    selected: List[Dict[str, str]] = []
    class_counts: Counter[str] = Counter()
    for row in rows:
        priority, reason = classify(row)
        class_counts[priority] += 1
        if priority == "DEFER":
            continue
        selected.append({
            "review_priority": priority,
            **row,
            "review_reason": reason,
        })

    selected.sort(key=lambda row: (
        {"P0": 0, "P1": 1, "P2": 2, "P3": 3}[row["review_priority"]],
        -to_int(row.get("n_title_matches")),
        -to_int(row.get("n_context_matches")),
        -to_int(row.get("n_works")),
        -to_int(row.get("max_score")),
        row.get("canonical_name", ""),
    ))
    selected = selected[: max(1, args.max_review)]

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(selected)

    recovered = sorted(KNOWN_POSITIVES.intersection({r.get("canonical_name", "") for r in selected}))
    qc = {
        "input_candidate_species": len(rows),
        "review_queue_rows": len(selected),
        "classification_counts": dict(sorted(class_counts.items())),
        "known_positive_controls_in_review_queue": recovered,
        "known_positive_controls_in_review_queue_count": len(recovered),
        "review_mode": "direct_title_and_replicated_context_v1",
    }
    Path(args.qc_out).write_text(json.dumps(qc, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(qc, ensure_ascii=False), flush=True)

    if len(selected) < 10:
        raise RuntimeError(f"Review queue unexpectedly small: {len(selected)}")
    if len(recovered) < 4:
        raise RuntimeError(f"Known-positive review recall failed: {len(recovered)}/4")


if __name__ == "__main__":
    main()

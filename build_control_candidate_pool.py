#!/usr/bin/env python3
"""Generate taxonomically matched candidate controls without assuming monomorphism.

These are screening candidates only. Unknown evidence is never recoded as monomorphic.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Set

FIELDS = [
    "focal_species", "focal_family", "focal_genus", "control_candidate",
    "control_family", "control_genus", "match_level", "match_rank",
    "control_verification_status", "monomorphism_evidence", "exclusion_reason",
]


def read_rows(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def genus_of(name: str) -> str:
    parts = name.split()
    return parts[0] if parts else ""


def candidate_names(paths: Iterable[Path]) -> Set[str]:
    names: Set[str] = set()
    for path in paths:
        if not path.exists():
            continue
        for row in read_rows(path):
            name = row.get("canonical_name", "").strip()
            if name:
                names.add(name)
    return names


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--census", required=True)
    parser.add_argument("--focals", default="data/global_flower_colour_review_queue.csv")
    parser.add_argument("--deferred", default="data/global_flower_colour_deferred.csv")
    parser.add_argument("--out", default="data/control_candidate_pool.csv")
    parser.add_argument("--qc-out", default="data/control_candidate_pool_qc.json")
    parser.add_argument("--per-focal", type=int, default=5)
    args = parser.parse_args()

    focals = read_rows(Path(args.focals))
    census = read_rows(Path(args.census))
    excluded = candidate_names([Path(args.focals), Path(args.deferred)])

    by_genus: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    by_family: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    for row in census:
        name = row.get("canonical_name", "").strip()
        family = row.get("family", "").strip()
        genus = row.get("genus", "").strip() or genus_of(name)
        if not name or " " not in name or name in excluded:
            continue
        clean = {"name": name, "family": family, "genus": genus}
        if genus:
            by_genus[genus].append(clean)
        if family:
            by_family[family].append(clean)

    for rows in by_genus.values():
        rows.sort(key=lambda r: r["name"])
    for rows in by_family.values():
        rows.sort(key=lambda r: (r["genus"], r["name"]))

    output: List[Dict[str, str]] = []
    seen_pairs: Set[tuple[str, str]] = set()
    level_counts: Counter[str] = Counter()

    for focal in focals:
        focal_name = focal.get("canonical_name", "").strip()
        family = focal.get("family", "").strip()
        genus = genus_of(focal_name)
        selected: List[tuple[Dict[str, str], str]] = []

        for row in by_genus.get(genus, []):
            selected.append((row, "same_genus"))
            if len(selected) >= args.per_focal:
                break

        if len(selected) < args.per_focal:
            used = {row["name"] for row, _ in selected}
            for row in by_family.get(family, []):
                if row["name"] in used:
                    continue
                selected.append((row, "same_family"))
                if len(selected) >= args.per_focal:
                    break

        for rank, (row, level) in enumerate(selected, start=1):
            pair = (focal_name, row["name"])
            if pair in seen_pairs:
                continue
            seen_pairs.add(pair)
            level_counts[level] += 1
            output.append({
                "focal_species": focal_name,
                "focal_family": family,
                "focal_genus": genus,
                "control_candidate": row["name"],
                "control_family": row["family"],
                "control_genus": row["genus"],
                "match_level": level,
                "match_rank": str(rank),
                "control_verification_status": "pending_monomorphism_review",
                "monomorphism_evidence": "unknown",
                "exclusion_reason": "",
            })

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(output)

    focal_with_controls = len({row["focal_species"] for row in output})
    qc = {
        "focal_species": len(focals),
        "focal_species_with_controls": focal_with_controls,
        "control_candidate_rows": len(output),
        "unique_control_candidates": len({row["control_candidate"] for row in output}),
        "match_level_counts": dict(sorted(level_counts.items())),
        "per_focal_target": args.per_focal,
        "semantic_guard": "unknown_is_not_monomorphic",
        "pool_mode": "taxonomically_matched_screening_candidates_v1",
    }
    Path(args.qc_out).write_text(json.dumps(qc, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(qc, ensure_ascii=False), flush=True)

    if len(focals) < 10:
        raise RuntimeError("Focal review queue unexpectedly small")
    if focal_with_controls < int(0.9 * len(focals)):
        raise RuntimeError(
            f"Too few focal species received controls: {focal_with_controls}/{len(focals)}"
        )


if __name__ == "__main__":
    main()

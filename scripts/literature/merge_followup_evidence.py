#!/usr/bin/env python3
"""Merge and aggregate targeted OpenAlex evidence before reclassification."""
from __future__ import annotations

import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path

DIRECT_RE = re.compile(
    r"(?:flower|floral|petal|corolla|bract|labellum).{0,80}(?:colou?r|pigment).{0,80}(?:polymorph|morph|variation|dimorph|cline)"
    r"|(?:polymorph|dimorph|morph|variation|cline).{0,80}(?:flower|floral|petal|corolla|bract|labellum).{0,80}(?:colou?r|pigment)",
    re.I,
)
NATURAL_RE = re.compile(
    r"natural population|wild population|in the wild|within[- ]population|morph frequenc|"
    r"geographic variation|spatial variation|geographic cline|population differentiation|"
    r"balancing selection|frequency[- ]dependent|pollinator[- ]mediated|local adaptation",
    re.I,
)
ARTIFICIAL_RE = re.compile(
    r"cultivar|breeding|ornamental|transgenic|mutagenesis|horticultur|floricultur|crop improvement",
    re.I,
)


def read(path: str) -> list[dict[str, str]]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def as_int(value: object) -> int:
    try:
        return int(str(value or 0))
    except ValueError:
        return 0


def evidence_strength(row: dict[str, str]) -> tuple[int, bool, bool, bool]:
    text = f"{row.get('title', '')} {row.get('evidence_snippet', '')}"
    direct = bool(DIRECT_RE.search(text))
    natural = bool(NATURAL_RE.search(text))
    artificial = bool(ARTIFICIAL_RE.search(text))
    strength = as_int(row.get("score")) + 18 * direct + 14 * natural - 16 * artificial
    return strength, direct, natural, artificial


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", default="data/global_flower_colour_species_ranked.csv")
    parser.add_argument("--followup", action="append", default=[])
    parser.add_argument("--out", default="data/global_flower_colour_species_ranked_resolved.csv")
    parser.add_argument("--qc-out", default="data/followup_merge_qc.json")
    args = parser.parse_args()

    species = read(args.species)
    by_name = {row["canonical_name"]: row for row in species}
    before = {
        name: (as_int(row.get("n_works")), as_int(row.get("n_context_matches")), as_int(row.get("max_score")))
        for name, row in by_name.items()
    }

    evidence_by_name: dict[str, list[dict[str, str]]] = defaultdict(list)
    seen: set[tuple[str, str]] = set()
    for path in args.followup:
        if not Path(path).exists():
            continue
        for evidence in read(path):
            name = evidence.get("canonical_name", "")
            work_id = evidence.get("openalex_id", "")
            if name not in by_name or not work_id or (name, work_id) in seen:
                continue
            seen.add((name, work_id))
            evidence_by_name[name].append(evidence)

    for name, evidence_rows in evidence_by_name.items():
        row = by_name[name]
        scored = [(evidence_strength(e), e) for e in evidence_rows]
        direct_count = sum(int(meta[1]) for meta, _ in scored)
        natural_count = sum(int(meta[2]) for meta, _ in scored)
        artificial_count = sum(int(meta[3]) for meta, _ in scored)
        title_count = sum(int(name.lower() in e.get("title", "").lower()) for _, e in scored)

        row["n_works"] = str(as_int(row.get("n_works")) + len(evidence_rows))
        row["n_context_matches"] = str(as_int(row.get("n_context_matches")) + len(evidence_rows))
        row["n_title_matches"] = str(as_int(row.get("n_title_matches")) + title_count)
        row["total_score"] = str(as_int(row.get("total_score")) + sum(as_int(e.get("score")) for _, e in scored))
        row["followup_evidence_count"] = str(len(evidence_rows))
        row["followup_direct_count"] = str(direct_count)
        row["followup_natural_count"] = str(natural_count)
        row["followup_artificial_count"] = str(artificial_count)

        best_meta, best = max(scored, key=lambda item: item[0][0])
        best_score = as_int(best.get("score"))
        if best_meta[0] > as_int(row.get("max_score")) or direct_count or natural_count:
            row["max_score"] = str(max(best_score, as_int(row.get("max_score"))))
            row["best_title"] = best.get("title", "")
            row["best_doi"] = best.get("doi", "")
            row["best_openalex_id"] = best.get("openalex_id", "")
            top = sorted(scored, key=lambda item: item[0][0], reverse=True)[:3]
            row["best_match_evidence"] = " || ".join(
                f"{e.get('title', '')}. {e.get('evidence_snippet', '')}" for _, e in top
            )[:1800]

    extra_fields = [
        "followup_evidence_count", "followup_direct_count",
        "followup_natural_count", "followup_artificial_count",
    ]
    fieldnames = list(species[0].keys())
    for field in extra_fields:
        if field not in fieldnames:
            fieldnames.append(field)
    for row in species:
        for field in extra_fields:
            row.setdefault(field, "0")

    species.sort(key=lambda row: (
        -as_int(row.get("followup_direct_count")),
        -as_int(row.get("followup_natural_count")),
        -as_int(row.get("n_title_matches")),
        -as_int(row.get("max_score")),
        -as_int(row.get("total_score")),
        row.get("canonical_name", ""),
    ))
    for rank, row in enumerate(species, 1):
        row["rank"] = str(rank)

    with Path(args.out).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(species)

    touched = set(evidence_by_name)
    strengthened = sum(
        1 for name in touched
        if (
            as_int(by_name[name].get("n_works")),
            as_int(by_name[name].get("n_context_matches")),
            as_int(by_name[name].get("max_score")),
        ) != before[name]
    )
    qc = {
        "followup_files": args.followup,
        "unique_followup_pairs": len(seen),
        "merged_rows": len(seen),
        "species_touched": len(touched),
        "species_strengthened": strengthened,
        "species_with_direct_followup": sum(as_int(by_name[n].get("followup_direct_count")) > 0 for n in touched),
        "species_with_natural_followup": sum(as_int(by_name[n].get("followup_natural_count")) > 0 for n in touched),
        "mode": "merge_followup_evidence_aggregated_v2",
    }
    Path(args.qc_out).write_text(json.dumps(qc, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(qc, ensure_ascii=False))


if __name__ == "__main__":
    main()

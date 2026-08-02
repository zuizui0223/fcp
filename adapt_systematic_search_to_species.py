#!/usr/bin/env python3
"""Convert the systematic publication queue into the legacy species-evidence schema.

This is a prioritisation adapter, not biological adjudication. Candidate species are
restricted to exact canonical names in the versioned angiosperm census.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

PRIORITY_SCORE = {
    "P1_high_natural_itv": 30,
    "P1_high_population_itv": 28,
    "P2_possible_itv": 18,
    "P3_likely_exclusion_review": 8,
    "P4_flower_colour_context_only": 5,
    "P5_low_relevance": 1,
}
HIGH = {"P1_high_natural_itv", "P1_high_population_itv"}
REVIEWABLE = HIGH | {"P2_possible_itv"}


def read_rows(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        yield from csv.DictReader(handle)


def write_rows(path: Path, fields: list[str], rows: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def split_species(value: str) -> list[str]:
    return [x.strip() for x in (value or "").split(";") if x.strip()]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--queue", required=True)
    ap.add_argument("--census", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--max-species-per-record", type=int, default=12)
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    census_rows = list(read_rows(Path(args.census)))
    census = {r["canonical_name"].strip(): r for r in census_rows if r.get("canonical_name")}

    per_species: dict[str, list[dict]] = defaultdict(list)
    priority_counts = Counter()
    database_counts = Counter()
    query_counts = Counter()
    input_records = 0
    records_with_valid_species = 0
    rejected_overbroad = 0

    works_fields = [
        "canonical_name", "family", "database", "query_ids", "record_id", "doi",
        "title", "abstract", "year", "journal", "url", "cited_by_count",
        "screen_priority", "provisional_scale_signal", "provisional_variation_form",
        "display_colour_signal", "variation_signal", "within_signal", "among_signal",
        "natural_signal", "cultivated_signal", "induced_signal", "ontogenetic_signal",
        "matched_terms", "review_status",
    ]
    works: list[dict] = []

    for row in read_rows(Path(args.queue)):
        input_records += 1
        priority = row.get("screen_priority", "")
        priority_counts[priority] += 1
        database_counts[row.get("database", "unknown")] += 1
        for query_id in (row.get("query_ids") or "").split(";"):
            if query_id.strip():
                query_counts[query_id.strip()] += 1

        names = []
        seen = set()
        for name in split_species(row.get("candidate_species", "")):
            if name in census and name not in seen:
                seen.add(name)
                names.append(name)
        if not names:
            continue
        if len(names) > args.max_species_per_record:
            rejected_overbroad += 1
            continue
        records_with_valid_species += 1
        for name in names:
            item = dict(row)
            item["canonical_name"] = name
            item["family"] = census[name].get("family", "")
            per_species[name].append(item)
            works.append({field: item.get(field, "") for field in works_fields})

    ranked = []
    for name, rows in per_species.items():
        def row_score(r: dict) -> int:
            score = PRIORITY_SCORE.get(r.get("screen_priority", ""), 0)
            score += 4 if str(r.get("natural_signal", "")).lower() == "true" else 0
            score += 3 if str(r.get("within_signal", "")).lower() == "true" else 0
            score += 2 if str(r.get("among_signal", "")).lower() == "true" else 0
            score += min(5, int(float(r.get("cited_by_count") or 0)) // 100)
            return score

        ordered = sorted(rows, key=lambda r: (row_score(r), int(float(r.get("cited_by_count") or 0))), reverse=True)
        best = ordered[0]
        title_matches = sum(name.lower() in (r.get("title") or "").lower() for r in rows)
        context_matches = sum(r.get("screen_priority") in REVIEWABLE for r in rows)
        scores = [row_score(r) for r in rows]
        ranked.append({
            "rank": 0,
            "canonical_name": name,
            "family": census[name].get("family", ""),
            "n_works": len(rows),
            "max_score": max(scores),
            "total_score": sum(scores),
            "n_title_matches": title_matches,
            "n_context_matches": context_matches,
            "best_title": best.get("title", ""),
            "best_doi": best.get("doi", ""),
            "best_openalex_id": best.get("record_id", "") if "openalex" in best.get("database", "") else "",
            "best_match_evidence": (best.get("abstract") or "")[:1000],
            "review_status": "unreviewed",
            "p1_works": sum(r.get("screen_priority") in HIGH for r in rows),
            "p2_works": sum(r.get("screen_priority") == "P2_possible_itv" for r in rows),
            "within_signal_works": sum(str(r.get("within_signal", "")).lower() == "true" for r in rows),
            "among_signal_works": sum(str(r.get("among_signal", "")).lower() == "true" for r in rows),
            "natural_signal_works": sum(str(r.get("natural_signal", "")).lower() == "true" for r in rows),
        })

    ranked.sort(key=lambda r: (-r["p1_works"], -r["max_score"], -r["n_context_matches"], -r["n_works"], r["canonical_name"]))
    for index, row in enumerate(ranked, 1):
        row["rank"] = index

    ranked_fields = [
        "rank", "canonical_name", "family", "n_works", "max_score", "total_score",
        "n_title_matches", "n_context_matches", "best_title", "best_doi",
        "best_openalex_id", "best_match_evidence", "review_status", "p1_works",
        "p2_works", "within_signal_works", "among_signal_works", "natural_signal_works",
    ]
    write_rows(out / "systematic_flower_colour_works_by_species.csv", works_fields, works)
    write_rows(out / "global_flower_colour_species_ranked_systematic.csv", ranked_fields, ranked)
    write_rows(out / "global_flower_colour_species_shortlist_systematic.csv", ranked_fields, ranked[:1000])

    family_counts = Counter(r["family"] or "Unknown" for r in ranked)
    family_rows = [
        {"family": family, "candidate_species": count,
         "p1_supported_species": sum(r["family"] == family and r["p1_works"] > 0 for r in ranked)}
        for family, count in family_counts.most_common()
    ]
    write_rows(out / "systematic_family_candidate_summary.csv",
               ["family", "candidate_species", "p1_supported_species"], family_rows)

    qc = {
        "status": "complete",
        "input_publication_records": input_records,
        "records_with_valid_census_species": records_with_valid_species,
        "records_rejected_as_overbroad": rejected_overbroad,
        "species_candidates": len(ranked),
        "species_with_p1_evidence": sum(r["p1_works"] > 0 for r in ranked),
        "species_with_multiple_p1_works": sum(r["p1_works"] >= 2 for r in ranked),
        "publication_priority_counts": dict(sorted(priority_counts.items())),
        "database_record_counts": dict(sorted(database_counts.items())),
        "top_query_counts": dict(query_counts.most_common(30)),
        "semantic_guard": "candidate_species_and_screening_priorities_are_not_adjudicated_natural_itv_or_strict_fcp",
        "next_required_stage": "duplicate_source_level_human_review_and_adjudication_before_macroecological_case_control_models",
    }
    (out / "systematic_species_adapter_qc.json").write_text(json.dumps(qc, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(qc, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

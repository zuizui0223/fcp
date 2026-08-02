#!/usr/bin/env python3
"""Build an exploratory within/among classification from systematic-search P1 records.

Candidate binomials are accepted when their genus occurs in the angiosperm census.
Exact census species matches are preferred; species absent from the census can still
proceed when the genus maps unambiguously to a family. The downstream GBIF backbone
step provides strict taxonomic validation.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

BINOMIAL = re.compile(r"^[A-Z][a-z]{2,} [a-z][a-z-]{2,}$")
TRUTHY = {"1", "true", "yes"}


def truthy(value: object) -> bool:
    return str(value or "").strip().lower() in TRUTHY


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--queue", required=True)
    parser.add_argument("--census", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--qc-out", required=True)
    args = parser.parse_args()

    out_path = Path(args.out)
    qc_path = Path(args.qc_out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    exact_family: dict[str, str] = {}
    genus_families: dict[str, set[str]] = defaultdict(set)
    with Path(args.census).open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            name = str(row.get("canonical_name") or "").strip()
            genus = str(row.get("genus") or (name.split()[0] if name else "")).strip()
            family = str(row.get("family") or "").strip()
            if name and family:
                exact_family[name] = family
            if genus and family:
                genus_families[genus].add(family)

    unique_genus_family = {
        genus: next(iter(families))
        for genus, families in genus_families.items()
        if len(families) == 1
    }

    diagnostics: Counter[str] = Counter()
    evidence: dict[str, list[dict[str, str]]] = defaultdict(list)
    family_for_name: dict[str, str] = {}
    source_for_name: dict[str, str] = {}

    with Path(args.queue).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        diagnostics["queue_columns"] = len(reader.fieldnames or [])
        for row in reader:
            diagnostics["records_total"] += 1
            if row.get("screen_priority") not in {
                "P1_high_natural_itv",
                "P1_high_population_itv",
            }:
                continue
            diagnostics["records_p1"] += 1

            # Keep P1 records unless they are explicitly artificial/induced or
            # ontogenetic. P1 already encodes flower-colour variation relevance.
            if any(
                truthy(row.get(column))
                for column in ("cultivated_signal", "induced_signal", "ontogenetic_signal")
            ):
                diagnostics["records_hard_exclusion"] += 1
                continue

            within = truthy(row.get("within_signal"))
            among = truthy(row.get("among_signal"))
            if not within and not among:
                diagnostics["records_without_direction"] += 1
                continue

            accepted_in_record = 0
            raw_names = str(row.get("candidate_species") or "")
            for token in raw_names.split(";"):
                name = " ".join(token.strip().split())
                if not BINOMIAL.fullmatch(name):
                    diagnostics["candidate_tokens_non_binomial"] += 1
                    continue
                genus = name.split()[0]
                if name in exact_family:
                    family = exact_family[name]
                    source = "exact_census_species"
                    diagnostics["candidate_exact_census"] += 1
                elif genus in unique_genus_family:
                    family = unique_genus_family[genus]
                    source = "census_genus_family_then_gbif_validation"
                    diagnostics["candidate_genus_rescued"] += 1
                else:
                    diagnostics["candidate_unknown_or_ambiguous_genus"] += 1
                    continue

                evidence[name].append(row)
                family_for_name[name] = family
                source_for_name[name] = source
                accepted_in_record += 1

            if accepted_in_record:
                diagnostics["records_with_accepted_plant_candidate"] += 1
            else:
                diagnostics["records_without_accepted_plant_candidate"] += 1

    rows: list[dict[str, object]] = []
    excluded: Counter[str] = Counter()
    for name, records in evidence.items():
        within_records = [row for row in records if truthy(row.get("within_signal"))]
        among_records = [row for row in records if truthy(row.get("among_signal"))]
        if within_records and among_records:
            excluded["mixed_species"] += 1
            continue
        if not within_records and not among_records:
            excluded["no_direction_species"] += 1
            continue

        directional = within_records or among_records
        scale = "within_population" if within_records else "among_population"
        best = max(
            directional,
            key=lambda row: (
                int(float(row.get("cited_by_count") or 0)),
                len(row.get("abstract") or ""),
            ),
        )
        rows.append(
            {
                "canonical_name": name,
                "family": family_for_name[name],
                "spatial_scale": scale,
                "classification_source": "baseline_unambiguous",
                "automated_classification_source": "systematic_P1_one_direction_plant_genus_validated",
                "taxon_candidate_source": source_for_name[name],
                "n_supporting_p1_records": len(records),
                "within_signal_records": len(within_records),
                "among_signal_records": len(among_records),
                "natural_signal_records": sum(truthy(r.get("natural_signal")) for r in records),
                "discrete_signal_records": sum(
                    r.get("provisional_variation_form") == "discrete_signal" for r in records
                ),
                "best_title": best.get("title", ""),
                "best_doi": best.get("doi", ""),
                "best_openalex_id": best.get("record_id", ""),
                "best_match_evidence": (best.get("abstract") or "")[:1200],
                "review_status": "unreviewed",
                "role": "focal",
            }
        )

    rows.sort(
        key=lambda row: (
            str(row["spatial_scale"]),
            -int(row["n_supporting_p1_records"]),
            str(row["canonical_name"]),
        )
    )
    classes = {str(row["spatial_scale"]) for row in rows}
    status = "complete" if len(rows) >= 20 and len(classes) == 2 else "insufficient"
    qc = {
        "status": status,
        "eligible_species_before_gbif_validation": len(rows),
        "families": len({str(row["family"]) for row in rows}),
        "within_species": sum(row["spatial_scale"] == "within_population" for row in rows),
        "among_species": sum(row["spatial_scale"] == "among_population" for row in rows),
        "record_diagnostics": dict(diagnostics),
        "excluded_species": dict(excluded),
        "semantic_guard": (
            "Exploratory automated P1 one-direction labels. Candidate binomials are "
            "plant-genus screened and require downstream strict GBIF validation; they "
            "are not source-level adjudicated biological classifications."
        ),
    }
    qc_path.write_text(json.dumps(qc, indent=2) + "\n", encoding="utf-8")

    if rows:
        with out_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

    print(json.dumps(qc, indent=2))
    if status != "complete":
        raise SystemExit(json.dumps(qc))


if __name__ == "__main__":
    main()

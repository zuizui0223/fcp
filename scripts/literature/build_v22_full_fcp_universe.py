#!/usr/bin/env python3
"""Build an all-species natural flower-colour polymorphism universe before C/S outcomes.

This corrects the main selection problem in the earlier replacement analysis: species
must enter the ecological universe because they have source-supported natural discrete
intraspecific floral-display colour variation, not because a C or S positive phrase was
found. C (local coexistence) and S (spatial segregation) are then measured as separate
positive documented-evidence axes inside that independently constructed universe.

Historical information is used only as an exact source->taxon rescue map when automated
source-text taxon extraction fails. Historical labels, climate values and model results
are never read.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import build_v22_coexistence_segregation_refined_v5 as rules
from build_systematic_spatial_evidence_axes import resolve_gbif_name

base = rules.base

# A conservative positive context that makes the source relevant to naturally occurring
# population-level variation. Generic greenhouse/genetic/cultivar material is handled by
# HARD_CONFLICT_RE and therefore cannot enter the universe from this signal alone.
NATURAL_POPULATION_CONTEXT_RE = re.compile(
    r"(?is)\b(?:wild(?:\s+populations?)?|natural(?:ly\s+occurring)?|natural\s+populations?|"
    r"field(?:\s+study|\s+studies|\s+populations?)?|in\s+situ|populations?|sites?|localit(?:y|ies)|"
    r"geographic(?:al)?\s+(?:range|regions?)|regions?|islands?|stands?)\b"
)


def clean(value: Any) -> str:
    return base.clean(value)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def resolve_one(name: str, timeout: int, retries: int) -> tuple[str, dict[str, Any]]:
    try:
        resolved = resolve_gbif_name(name, timeout, retries)
    except Exception as exc:  # noqa: BLE001
        resolved = {
            "input_name": name,
            "accepted": False,
            "accepted_name": "",
            "family": "",
            "usage_key": "",
            "reason": f"resolution_error:{type(exc).__name__}",
        }
    return name, resolved


def source_eligibility(title: str, abstract: str, work_type: str) -> dict[str, Any]:
    text = clean(f"{title} {abstract}")
    primary = work_type in base.PRIMARY_WORK_TYPES and not bool(base.NONPRIMARY_TEXT_RE.search(text))
    display = bool(base.DISPLAY_RE.search(text))
    discrete = bool(base.DISCRETE_RE.search(text))
    intraspecific = bool(base.INTRASPECIFIC_RE.search(text))
    population_context = bool(NATURAL_POPULATION_CONTEXT_RE.search(text))
    hard_conflict = bool(base.HARD_CONFLICT_RE.search(text))
    community_only = bool(base.COMMUNITY_RE.search(text)) and not intraspecific
    eligible = (
        primary
        and display
        and discrete
        and intraspecific
        and population_context
        and not hard_conflict
        and not community_only
    )
    return {
        "text": text,
        "primary": primary,
        "display": display,
        "discrete": discrete,
        "intraspecific": intraspecific,
        "population_context": population_context,
        "hard_conflict": hard_conflict,
        "community_only": community_only,
        "eligible": eligible,
    }


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--blind", required=True)
    p.add_argument("--key", required=True)
    p.add_argument("--historical-manifest", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--workers", type=int, default=32)
    p.add_argument("--timeout", type=int, default=20)
    p.add_argument("--retries", type=int, default=2)
    args = p.parse_args()

    blind = read_csv(Path(args.blind))
    key = read_csv(Path(args.key))
    if len(blind) != 12064 or len(key) != 12064:
        raise SystemExit(f"Expected 12064 blind/key rows, found {len(blind)}/{len(key)}")
    bmap = {r["record_review_id"]: r for r in blind}
    kmap = {r["record_review_id"]: r for r in key}
    if set(bmap) != set(kmap):
        raise SystemExit("Blind/key record IDs differ")

    historical_map = base.load_historical_taxon_map(Path(args.historical_manifest))
    historical_sources_seen = 0
    source_rows: list[dict[str, Any]] = []
    names_to_resolve: set[str] = set()

    for rid, row in bmap.items():
        hidden = kmap[rid]
        title = clean(row.get("title"))
        abstract = clean(row.get("abstract"))
        work_type = clean(row.get("work_type")).lower()
        source_key = base.norm_source(row.get("source_id"))
        historical_taxon = historical_map.get(source_key, "")
        if historical_taxon:
            historical_sources_seen += 1

        sig = source_eligibility(title, abstract, work_type)
        text = sig["text"]
        c_excerpt = base.positive_context(base.C_PATTERNS, text, axis="C") if sig["eligible"] else ""
        s_excerpt = base.positive_context(base.S_PATTERNS, text, axis="S") if sig["eligible"] else ""
        c_positive = bool(c_excerpt)
        s_positive = bool(s_excerpt)

        # Critical correction: taxon resolution is triggered by independent FCP
        # eligibility, NOT by C/S positivity.
        candidates = base.ranked_candidates(hidden, title, abstract, historical_taxon) if sig["eligible"] else []
        names_to_resolve.update(candidates)

        if sig["eligible"]:
            if c_positive and s_positive:
                status = "eligible_C_and_S_positive"
            elif c_positive:
                status = "eligible_C_positive"
            elif s_positive:
                status = "eligible_S_positive"
            else:
                status = "eligible_spatial_unresolved"
        elif sig["hard_conflict"]:
            status = "artificial_or_conflict_excluded"
        elif not sig["primary"]:
            status = "nonprimary_excluded"
        elif sig["community_only"]:
            status = "community_level_excluded"
        else:
            status = "not_strict_FCP_eligible"

        source_rows.append({
            "record_review_id": rid,
            "source_id": clean(row.get("source_id")),
            "title": title,
            "year": clean(row.get("year")),
            "work_type": work_type,
            "FCP_source_status": status,
            "FCP_eligible_source": int(sig["eligible"]),
            "primary_source_signal": int(sig["primary"]),
            "display_colour_signal": int(sig["display"]),
            "discrete_polymorphism_signal": int(sig["discrete"]),
            "intraspecific_signal": int(sig["intraspecific"]),
            "natural_population_context_signal": int(sig["population_context"]),
            "hard_conflict_signal": int(sig["hard_conflict"]),
            "community_level_signal": int(sig["community_only"]),
            "C_local_coexistence_documented_strict": int(c_positive),
            "S_spatial_segregation_documented_strict": int(s_positive),
            "C_evidence_excerpt": c_excerpt,
            "S_evidence_excerpt": s_excerpt,
            "historical_source_taxon_rescue": historical_taxon,
            "candidate_taxon_strings": ";".join(candidates),
        })

    workers = max(1, min(args.workers, 48))
    cache: dict[str, dict[str, Any]] = {}
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(resolve_one, name, args.timeout, args.retries): name for name in sorted(names_to_resolve)}
        for i, fut in enumerate(as_completed(futs), start=1):
            name, resolved = fut.result()
            cache[name] = resolved
            if i % 100 == 0:
                print({"gbif_resolved": i, "gbif_total": len(futs)}, flush=True)

    resolved_sources: list[dict[str, Any]] = []
    for row in source_rows:
        accepted: dict[str, dict[str, Any]] = {}
        for name in [x for x in row["candidate_taxon_strings"].split(";") if x]:
            resolved = cache.get(name, {})
            if resolved.get("accepted") and resolved.get("accepted_name"):
                accepted[str(resolved["accepted_name"])] = resolved
        if len(accepted) == 1:
            accepted_name, resolved = next(iter(accepted.items()))
            taxon_status = "resolved_unique"
            family = clean(resolved.get("family"))
            usage_key = clean(resolved.get("usage_key"))
        elif len(accepted) > 1:
            accepted_name = family = usage_key = ""
            taxon_status = "multiple_accepted_taxa_unresolved"
        else:
            accepted_name = family = usage_key = ""
            taxon_status = "no_accepted_taxon_resolved"
        completed = dict(row)
        completed.update({
            "accepted_name": accepted_name,
            "family": family,
            "gbif_usage_key": usage_key,
            "taxon_resolution_status": taxon_status,
        })
        resolved_sources.append(completed)

    # Species universe is defined by FCP eligibility first.
    by_species: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in resolved_sources:
        if int(row["FCP_eligible_source"]) == 1 and row["accepted_name"] and row["taxon_resolution_status"] == "resolved_unique":
            by_species[row["accepted_name"]].append(row)

    species_rows: list[dict[str, Any]] = []
    for species, rows in sorted(by_species.items()):
        c_sources = [r for r in rows if int(r["C_local_coexistence_documented_strict"]) == 1]
        s_sources = [r for r in rows if int(r["S_spatial_segregation_documented_strict"]) == 1]
        c = bool(c_sources)
        s = bool(s_sources)
        if c and s:
            state = "coexistence_and_segregation"
        elif c:
            state = "local_coexistence_only"
        elif s:
            state = "spatial_segregation_only"
        else:
            state = "organization_unresolved"
        years = [int(r["year"]) for r in rows if str(r.get("year", "")).isdigit()]
        family = next((r["family"] for r in rows if r["family"]), "")
        species_rows.append({
            "canonical_name": species,
            "family": family,
            "organization_state": state,
            "C_local_coexistence_documented": int(c),
            "S_spatial_segregation_documented": int(s),
            "n_FCP_eligible_sources": len(rows),
            "n_C_positive_sources": len(c_sources),
            "n_S_positive_sources": len(s_sources),
            "earliest_source_year": min(years) if years else "",
            "latest_source_year": max(years) if years else "",
            "source_year_span": (max(years) - min(years)) if years else "",
            "FCP_source_ids": ";".join(r["source_id"] for r in rows),
            "C_source_ids": ";".join(r["source_id"] for r in c_sources),
            "S_source_ids": ";".join(r["source_id"] for r in s_sources),
        })

    informative = [r for r in species_rows if int(r["C_local_coexistence_documented"]) or int(r["S_spatial_segregation_documented"])]
    unresolved = [r for r in species_rows if r["organization_state"] == "organization_unresolved"]
    unresolved_eligible_sources = [
        r for r in resolved_sources
        if int(r["FCP_eligible_source"]) == 1 and r["taxon_resolution_status"] != "resolved_unique"
    ]

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    write_csv(outdir / "v22_full_fcp_source_audit.csv", resolved_sources, list(resolved_sources[0].keys()))
    write_csv(outdir / "v22_full_fcp_species_universe.csv", species_rows, list(species_rows[0].keys()) if species_rows else ["canonical_name"])
    write_csv(outdir / "v22_full_fcp_informative_states.csv", informative, list(species_rows[0].keys()) if species_rows else ["canonical_name"])
    write_csv(outdir / "v22_full_fcp_organization_unresolved.csv", unresolved, list(species_rows[0].keys()) if species_rows else ["canonical_name"])
    write_csv(outdir / "v22_full_fcp_taxon_unresolved_sources.csv", unresolved_eligible_sources, list(resolved_sources[0].keys()))
    (outdir / "v22_full_fcp_gbif_cache.json").write_text(json.dumps(cache, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    states: dict[str, int] = defaultdict(int)
    for row in species_rows:
        states[row["organization_state"]] += 1
    source_status: dict[str, int] = defaultdict(int)
    for row in resolved_sources:
        source_status[row["FCP_source_status"]] += 1

    summary = {
        "status": "complete",
        "input_records": len(resolved_sources),
        "historical_exact_sources_seen_for_taxon_rescue_only": historical_sources_seen,
        "eligible_FCP_source_records": sum(int(r["FCP_eligible_source"]) for r in resolved_sources),
        "eligible_FCP_sources_taxonomically_unresolved": len(unresolved_eligible_sources),
        "gbif_names_queried": len(cache),
        "FCP_species_universe": len(species_rows),
        "C_positive_species": sum(int(r["C_local_coexistence_documented"]) for r in species_rows),
        "S_positive_species": sum(int(r["S_spatial_segregation_documented"]) for r in species_rows),
        "informative_C_or_S_species": len(informative),
        "organization_unresolved_species": len(unresolved),
        "species_state_counts": dict(states),
        "source_status_counts": dict(source_status),
        "universe_definition": (
            "Species enter because at least one primary source supports discrete intraspecific floral-display colour variation "
            "in a natural/population context without an artificial/conflict exclusion; C/S positivity is not an inclusion criterion."
        ),
        "axis_definition": {
            "C": "explicit local coexistence of discrete natural floral-colour variants in the same population/site",
            "S": "explicit spatial segregation/structuring of colour variants or morph frequencies among geographic units",
            "zero_semantics": "not documented by this strict pass; not biological absence",
        },
        "historical_34_role": "taxon rescue and later historical sensitivity only; never universe membership",
    }
    (outdir / "v22_full_fcp_universe_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)

    if len(resolved_sources) != 12064:
        raise SystemExit(f"Expected 12064 source rows, found {len(resolved_sources)}")
    if historical_sources_seen != 34:
        raise SystemExit(f"Expected 34 exact historical source records for rescue audit, found {historical_sources_seen}")
    if len(species_rows) <= 34:
        raise SystemExit(f"All-species universe did not expand beyond historical scale: {len(species_rows)} species")
    if not unresolved:
        raise SystemExit("Expected some organization-unresolved species in an independently constructed universe")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Refined single-pass audit of local coexistence (C) and spatial segregation (S).

Both axes are positive-evidence variables. C=0 never means biological non-coexistence,
and S=0 never means biological homogeneity; zero means not documented by this strict
pass. Ambiguous, non-primary, continuous-only, artificial, cultivated, ontogenetic,
community-level, or taxonomically ambiguous records remain unresolved.

Historical labels are never used. The historical manifest is used only as an exact
source->taxon rescue map for its 34 known source records.
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

from build_systematic_spatial_evidence_axes import resolve_gbif_name

PRIMARY_WORK_TYPES = {"article", "dissertation", "preprint", "report"}
COLORS = r"(?:white|yellow|pink|purple|blue|red|orange|green|cream|violet)"


def clean(value: Any) -> str:
    return " ".join(str(value or "").split()).strip()


def norm_source(value: Any) -> str:
    x = clean(value).lower().rstrip("./")
    x = re.sub(r"^https?://(?:dx\.)?doi\.org/", "", x)
    x = re.sub(r"^doi:\s*", "", x)
    if "openalex.org/" in x:
        return "openalex:" + x.rsplit("/", 1)[-1].upper()
    if re.fullmatch(r"w\d+", x, re.I):
        return "openalex:" + x.upper()
    if x.startswith("10."):
        return "doi:" + x
    return x


DISPLAY_RE = re.compile(
    rf"(?is)\b(?:flower|floral|corolla|petal|tepal|perianth|labellum|showy\s+bract)[-\s]*(?:colou?r|pigment(?:ation)?|reflectance)\b"
    rf"|\b(?:colou?r|pigment(?:ation)?|reflectance)\b.{{0,30}}\b(?:flower|floral|corolla|petal|tepal|perianth|labellum)\b"
    rf"|\b(?:two|multiple|several)\s+colou?rs?\s+of\s+flowers?\b"
    rf"|\b{COLORS}\b.{{0,35}}\b(?:and|or)\b.{{0,35}}\b{COLORS}\b.{{0,25}}\bflowers?\b"
)

DISCRETE_RE = re.compile(
    rf"(?is)\b(?:"
    rf"(?:flower|floral|petal|corolla|perianth)[-\s]*colou?r\s+(?:morphs?|forms?|phenotypes?|polymorph\w*)"
    rf"|colou?r\s+(?:morphs?|forms?|phenotypes?|polymorph\w*)"
    rf"|(?:flower|floral)\s+polymorph\w*"
    rf"|heterocyanic\s+(?:population|populations)"
    rf"|(?:two|multiple|several|three|four|five)\s+(?:distinct\s+)?(?:flower|floral|petal)?\s*colou?rs?(?:ed)?\s+(?:morphs?|forms?|phenotypes?|flowers?|individuals?)"
    rf"|(?:two|multiple|several)\s+colou?rs?\s+of\s+flowers?"
    rf"|{COLORS}\b.{{0,35}}\b(?:and|or)\b.{{0,35}}\b{COLORS}\b.{{0,35}}\b(?:flowers?|flowered\s+individuals?|morphs?|forms?|phenotypes?)"
    rf"|(?:flower|floral|colou?r).{{0,50}}polytypism"
    rf"|polytypism.{{0,50}}(?:flower|floral|colou?r)"
    rf")\b"
)

HARD_CONFLICT_RE = re.compile(
    r"(?is)\b(?:cultivars?|horticultural\s+lines?|commercial\s+variet(?:y|ies)|breeding\s+lines?|breeding\s+(?:program|process)|"
    r"mapping\s+population|transgenic|gene[- ]edited|CRISPR|induced\s+mutat(?:ion|ions|ed)?|irradiat(?:ed|ion)?|mutagenesis|"
    r"tissue\s+culture|somaclonal|ontogenetic\s+colou?r\s+change|flower\s+colou?r\s+change\s+with\s+age|"
    r"colou?r\s+change\s+after\s+pollination|changes?\s+colou?r\s+after\s+pollination|"
    r"experimental\s+populations?|artificial\s+populations?|artificial\s+flowers?|arrays?\s+of\s+artificial\s+flowers?)\b"
)

NONPRIMARY_TEXT_RE = re.compile(
    r"(?is)\b(?:systematic\s+review|meta[- ]analysis|review\s+article|we\s+review|we\s+also\s+review|"
    r"here\s+we\s+summari[sz]e|this\s+review|review\s+the\s+incidence)\b"
)
COMMUNITY_RE = re.compile(
    r"(?is)\b(?:species\s+proportion\s+of\s+different\s+flower\s+colou?rs?|community[- ]level|"
    r"across\s+species|flowering\s+plants\s+of\s+\w+\s+mountain)\b"
)
INTRASPECIFIC_RE = re.compile(
    r"(?is)\b(?:intraspecific|within[- ]population|within\s+(?:a|the|some|each)?\s*population|among\s+populations|"
    r"colou?r\s+morphs?|flower[- ]?colou?r\s+polymorph|floral[- ]?colou?r\s+polymorph|polymorphic\s+populations?)\b"
)
DEFINITION_RE = re.compile(r"(?is)\b(?:is\s+defined\s+as|defined\s+as|is\s+the\s+occurrence\s+of|refers?\s+to)\b")
NEG_C_RE = re.compile(
    r"(?is)\b(?:do(?:es)?\s+not|did\s+not|not|never|no)\b.{0,50}\b(?:co[- ]?occur|coexist|within\s+(?:the\s+same\s+)?population|same\s+population)"
)
NEG_S_RE = re.compile(
    r"(?is)(?:\b(?:no|not|without)\b.{0,70}\b(?:geographic|spatial|among[- ]population|between[- ]population|morph\s+frequenc|flower[- ]?colou?r)\b.{0,40}\b(?:variation|difference|structure|differentiation|pattern|association)"
    r"|\bconsistent\s+(?:ratio|ratios|frequency|frequencies)\s+of\s+(?:colou?r\s+)?morphs?\s+among\s+populations\b)"
)
OTHER_SPECIES_RE = re.compile(
    r"(?is)\b(?:other\s+species|across\s+species|comparison\s+with\s+other\s+species|relative\s+to\s+those\s+that\s+fix\s+between\s+populations)\b"
)

DET = r"(?:a|the|same|some|several|many|certain|both|each|one|single|study|these|those)?"
POP = r"(?:population|populations|site|sites|stand|stands)"

C_PATTERNS = [
    re.compile(rf"(?is)\b(?:co[- ]?occur(?:s|red|ring)?|coexist(?:s|ed|ing)?)\b.{{0,100}}\b(?:within|in)\s+{DET}\s*{POP}\b"),
    re.compile(rf"(?is)\b(?:within|in)\s+{DET}\s*{POP}\b.{{0,100}}\b(?:co[- ]?occur(?:s|red|ring)?|coexist(?:s|ed|ing)?)\b"),
    re.compile(rf"(?is)\b{COLORS}\b.{{0,60}}\b(?:and|or)\b.{{0,60}}\b{COLORS}\b.{{0,80}}\b(?:within\s+{DET}\s*{POP}|same\s+{POP})\b"),
    re.compile(rf"(?is)\b(?:within\s+{DET}\s*{POP}|same\s+{POP})\b.{{0,100}}\b{COLORS}\b.{{0,60}}\b(?:and|or)\b.{{0,60}}\b{COLORS}\b"),
    re.compile(r"(?is)\b(?:morphs?|forms?|phenotypes?)\b.{0,90}\bin\s+each\s+of\b.{0,60}\b(?:wild\s+)?populations?\b"),
    re.compile(r"(?is)\b(?:study\s+population|population)\b.{0,100}\b(?:contains?|contained|has|had|with|bears?|three|two|multiple|several)\b.{0,100}\b(?:floral|flower)?\s*colou?r\s+(?:morphs?|forms?|phenotypes?|polymorphism)\b"),
    re.compile(r"(?is)\b(?:both|two|multiple|several)\b.{0,90}\b(?:flowered\s+individuals?|colou?r\s+morphs?|colou?rs?\s+of\s+flowers?)\b.{0,100}\b(?:both\s+)?within\s+(?:a|the|same|some)?\s*populations?\b"),
    re.compile(r"(?is)\b(?:bears?|has|have)\s+(?:two|multiple|several)\s+colou?rs?\s+of\s+flowers?\b.{0,100}\bwithin\s+populations?\b"),
    re.compile(r"(?is)\b(?:flower|floral)?\s*colou?r\s+(?:morphs?|forms?|phenotypes?)\b.{0,100}\b(?:var(?:y|ies|ied)|occur|present)\b.{0,50}\bwithin\s+(?:and\s+among\s+)?(?:the\s+)?(?:\d+\s+)?populations?\b"),
    re.compile(r"(?is)\b(?:morphs?|forms?|phenotypes?)\b.{0,70}\b(?:var(?:y|ies|ied)|occur|present)\b.{0,50}\bwithin\s+and\s+among\b.{0,35}\b(?:the\s+)?(?:\d+\s+)?populations?\b"),
    re.compile(r"(?is)\bthere\s+(?:is|are|were)\b.{0,40}\bpolymorph\w*\b.{0,60}\bwithin\s+(?:the\s+)?(?:\d+\s+)?populations?\b"),
    re.compile(rf"(?is)\b(?:perianth|flowers?|floral\s+display)\b.{{0,80}}\b(?:either\s+)?{COLORS}\b.{{0,45}}\b(?:or|and)\b.{{0,45}}\b{COLORS}(?:/\w+)?\b.{{0,80}}\bwithin\s+(?:some\s+|the\s+|a\s+)?populations?\b"),
    re.compile(r"(?is)\b(?:segregating\s+variation|variation)\s+in\s+(?:flower|floral)\s+colou?r\b.{0,80}\bwithin\s+(?:a|the|one)\s+population\b"),
    re.compile(r"(?is)\bcolou?r[- ]polymorphic\s+populations?\b.{0,100}\b(?:comprised|composed|consisting)\s+of\b.{0,120}\b(?:flowers?|morphs?|forms?)\b"),
    re.compile(r"(?is)\bnatural\s+population\b.{0,180}\b(?:two|three|multiple|several)\s+(?:types?\s+of\s+)?(?:floral|flower)\s+colou?r\s+morphs?\b"),
    re.compile(r"(?is)\b(?:\d+|one|two|three|four|five|six|seven|eight|nine|ten|some|several)\s+populations?\s+(?:were|are)\s+mixed\b"),
    re.compile(r"(?is)\bpopulations?\s+(?:were|are)\s+uniformly\s+colou?red\b.{0,100}\b(?:\d+|one|two|three|four|five|six|seven|eight|nine|ten|some|several)\s+(?:were|are)\s+mixed\b"),
]

S_PATTERNS = [
    re.compile(r"(?is)\b(?:flower[- ]?colou?r|floral[- ]?colou?r|colou?r\s+morphs?|morph\s+frequenc(?:y|ies))\s+(?:variation|differentiation|differences?|distribution|frequenc(?:y|ies))\b.{0,80}\b(?:among|between|across)\s+(?:populations|sites|localities|regions|islands|stands)\b"),
    re.compile(r"(?is)\b(?:among|between|across)\s+(?:populations|sites|localities|regions|islands|stands)\b.{0,90}\b(?:flower[- ]?colou?r|floral[- ]?colou?r|colou?r\s+morphs?|morph\s+frequenc(?:y|ies))\b"),
    re.compile(r"(?is)\b(?:geographic|geographical|spatial|latitudinal|longitudinal|altitudinal|elevational)\s+(?:variation|differentiation|structure|structuring|cline|pattern|gradient|distribution)\s+(?:in|of|for)\s+(?:flower[- ]?colou?r|floral[- ]?colou?r|colou?r\s+morphs?|morph\s+frequenc(?:y|ies)|(?:flower|floral)\s+colou?r\s+polymorph\w*)\b"),
    re.compile(r"(?is)\b(?:flower[- ]?colou?r|floral[- ]?colou?r|colou?r\s+morphs?|morph\s+frequenc(?:y|ies)|(?:flower|floral)\s+colou?r\s+polymorph\w*)\b.{0,50}\b(?:geographic|geographical|spatial|latitudinal|longitudinal|altitudinal|elevational)\s+(?:variation|differentiation|structure|structuring|cline|pattern|gradient|distribution)\b"),
    re.compile(r"(?is)\b(?:flower|floral)?\s*colou?r\s+(?:morphs?|forms?|phenotypes?)\b.{0,100}\b(?:var(?:y|ies|ied)|differ(?:s|ed)?|occur)\b.{0,50}\bwithin\s+and\s+among\b.{0,35}\b(?:the\s+)?(?:\d+\s+)?populations?\b"),
    re.compile(r"(?is)\b(?:morphs?|forms?|phenotypes?)\b.{0,70}\b(?:var(?:y|ies|ied)|differ(?:s|ed)?)\b.{0,50}\bwithin\s+and\s+among\b.{0,35}\b(?:the\s+)?(?:\d+\s+)?populations?\b"),
    re.compile(r"(?is)\bpolytypism\b.{0,90}\bamong\s+populations?\b"),
    re.compile(r"(?is)\b(?:white|yellow|pink|purple|blue|red|orange|green)\s*[- ]?flowered\s+(?:plants?|individuals?)\b.{0,120}\b(?:restricted|confined|limited)\s+to\b"),
    re.compile(r"(?is)\b(?:white|yellow|pink|purple|blue|red|orange|green)\s+(?:morph|form|phenotype)s?\b.{0,120}\b(?:restricted|confined|limited)\s+to\b"),
    re.compile(r"(?is)\b(?:spatial|geographic|geographical)\s+autocorrelation\b.{0,100}\b(?:morph\s+frequency|colou?r\s+morph|flower[- ]?colou?r)\b"),
    re.compile(r"(?is)\b(?:morph\s+frequency|colou?r\s+morph|flower[- ]?colou?r)\b.{0,100}\b(?:spatial|geographic|geographical)\s+autocorrelation\b"),
    re.compile(r"(?is)\b(?:flower|floral|colou?r).{0,50}polymorph\w*\b.{0,60}\bboth\s+within\s+and\s+among\s+(?:the\s+)?(?:\d+\s+)?populations?\b"),
    re.compile(r"(?is)\b(?:\d+|one|two|three|four|five|six|seven|eight|nine|ten|some|several)\s+populations?\s+(?:were|are)\s+uniformly\s+colou?red\b.{0,100}\b(?:\d+|one|two|three|four|five|six|seven|eight|nine|ten|some|several)\s+(?:were|are)\s+mixed\b"),
    re.compile(r"(?is)\b(?:morph\s+frequenc(?:y|ies)|frequency\s+of\s+(?:the\s+)?(?:colou?r\s+)?morphs?)\b.{0,120}\b(?:var(?:y|ies|ied)|differ(?:s|ed)?|range[sd]?)\b.{0,80}\b(?:among|between|across)\s+(?:populations|sites|regions|localities|islands|stands)\b"),
]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def load_historical_taxon_map(path: Path) -> dict[str, str]:
    return {norm_source(row.get("source_id")): clean(row.get("canonical_name")) for row in read_csv(path)}


def positive_context(patterns: list[re.Pattern[str]], text: str, *, axis: str) -> str:
    for pattern in patterns:
        for match in pattern.finditer(text):
            start = max(0, match.start() - 130)
            end = min(len(text), match.end() + 130)
            context = text[start:end]
            if axis == "C":
                if NEG_C_RE.search(context) or DEFINITION_RE.search(context):
                    continue
            else:
                if NEG_S_RE.search(context) or OTHER_SPECIES_RE.search(context):
                    continue
            return clean(context)
    return ""


def ranked_candidates(hidden: dict[str, str], title: str, abstract: str, historical_taxon: str) -> list[str]:
    if historical_taxon:
        return [historical_taxon]
    detected = [clean(x) for x in clean(hidden.get("detected_binomial_strings")).split(";") if clean(x)]
    if not detected:
        return []
    title_lower = title.lower()
    title_hits = [name for name in detected if name.lower() in title_lower]
    if title_hits:
        return title_hits[:3]
    abstract_head = abstract[:2500].lower()
    abstract_hits = [name for name in detected if name.lower() in abstract_head]
    return abstract_hits[:3]


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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--blind", required=True)
    parser.add_argument("--key", required=True)
    parser.add_argument("--historical-manifest", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--workers", type=int, default=32)
    parser.add_argument("--timeout", type=int, default=25)
    parser.add_argument("--retries", type=int, default=3)
    args = parser.parse_args()

    blind = read_csv(Path(args.blind))
    key = read_csv(Path(args.key))
    if len(blind) != 12064 or len(key) != 12064:
        raise SystemExit(f"Expected 12064 blind/key rows, found {len(blind)}/{len(key)}")
    bmap = {row["record_review_id"]: row for row in blind}
    kmap = {row["record_review_id"]: row for row in key}
    if set(bmap) != set(kmap):
        raise SystemExit("Blind/key record IDs differ")

    historical_map = load_historical_taxon_map(Path(args.historical_manifest))
    names_to_resolve: set[str] = set()
    source_rows: list[dict[str, Any]] = []
    historical_sources_seen = 0

    for rid, row in bmap.items():
        hidden = kmap[rid]
        title = clean(row.get("title"))
        abstract = clean(row.get("abstract"))
        text = clean(f"{title} {abstract}")
        work_type = clean(row.get("work_type")).lower()
        source_key = norm_source(row.get("source_id"))
        historical_taxon = historical_map.get(source_key, "")
        if historical_taxon:
            historical_sources_seen += 1

        primary = work_type in PRIMARY_WORK_TYPES and not bool(NONPRIMARY_TEXT_RE.search(text))
        display = bool(DISPLAY_RE.search(text))
        discrete = bool(DISCRETE_RE.search(text))
        hard_conflict = bool(HARD_CONFLICT_RE.search(text))
        community_only = bool(COMMUNITY_RE.search(text)) and not bool(INTRASPECIFIC_RE.search(text))
        eligible_for_positive = primary and display and discrete and not hard_conflict and not community_only

        c_excerpt = positive_context(C_PATTERNS, text, axis="C") if eligible_for_positive else ""
        s_excerpt = positive_context(S_PATTERNS, text, axis="S") if eligible_for_positive else ""
        c_positive = bool(c_excerpt)
        s_positive = bool(s_excerpt)

        if c_positive and s_positive:
            evidence_status = "C_and_S_positive"
        elif c_positive:
            evidence_status = "C_positive"
        elif s_positive:
            evidence_status = "S_positive"
        elif hard_conflict:
            evidence_status = "conflict_unresolved"
        elif not primary:
            evidence_status = "nonprimary_unresolved"
        elif community_only:
            evidence_status = "community_level_unresolved"
        elif display and discrete:
            evidence_status = "polymorphism_relevant_unresolved"
        else:
            evidence_status = "no_strict_positive_evidence"

        candidates = ranked_candidates(hidden, title, abstract, historical_taxon) if (c_positive or s_positive or historical_taxon) else []
        names_to_resolve.update(candidates)

        source_rows.append({
            "record_review_id": rid,
            "source_id": clean(row.get("source_id")),
            "title": title,
            "year": clean(row.get("year")),
            "work_type": work_type,
            "strict_evidence_status": evidence_status,
            "primary_source_signal": int(primary),
            "display_colour_signal": int(display),
            "discrete_polymorphism_signal": int(discrete),
            "hard_conflict_signal": int(hard_conflict),
            "community_level_signal": int(community_only),
            "C_local_coexistence_documented_strict": int(c_positive),
            "S_spatial_segregation_documented_strict": int(s_positive),
            "C_evidence_excerpt": c_excerpt,
            "S_evidence_excerpt": s_excerpt,
            "historical_source_taxon_rescue": historical_taxon,
            "candidate_taxon_strings": ";".join(candidates),
        })

    workers = max(1, min(args.workers, 48))
    cache: dict[str, dict[str, Any]] = {}
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(resolve_one, name, args.timeout, args.retries): name for name in sorted(names_to_resolve)}
        for index, future in enumerate(as_completed(futures), start=1):
            name, resolved = future.result()
            cache[name] = resolved
            if index % 100 == 0:
                print({"gbif_resolved": index, "gbif_total": len(futures)}, flush=True)

    source_evidence: list[dict[str, Any]] = []
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
        source_evidence.append(completed)

    by_species: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in source_evidence:
        if row["accepted_name"]:
            by_species[row["accepted_name"]].append(row)

    species_rows: list[dict[str, Any]] = []
    for species, rows in sorted(by_species.items()):
        c_sources = [row for row in rows if int(row["C_local_coexistence_documented_strict"]) == 1]
        s_sources = [row for row in rows if int(row["S_spatial_segregation_documented_strict"]) == 1]
        conflicts = [row for row in rows if row["strict_evidence_status"] == "conflict_unresolved"]
        c_positive = bool(c_sources)
        s_positive = bool(s_sources)
        if c_positive and s_positive:
            state = "coexistence_and_segregation"
        elif c_positive:
            state = "local_coexistence_only"
        elif s_positive:
            state = "spatial_segregation_only"
        else:
            state = "unresolved"
        family = next((row["family"] for row in rows if row["family"]), "")
        species_rows.append({
            "canonical_name": species,
            "family": family,
            "organization_state": state,
            "C_local_coexistence_documented": int(c_positive),
            "S_spatial_segregation_documented": int(s_positive),
            "n_resolved_sources": len(rows),
            "n_C_positive_sources": len(c_sources),
            "n_S_positive_sources": len(s_sources),
            "n_conflict_sources": len(conflicts),
            "C_source_ids": ";".join(row["source_id"] for row in c_sources),
            "S_source_ids": ";".join(row["source_id"] for row in s_sources),
        })

    freeze_rows = [row for row in species_rows if row["organization_state"] != "unresolved"]
    unresolved_sources = [
        row for row in source_evidence
        if row["strict_evidence_status"] not in {"C_positive", "S_positive", "C_and_S_positive"}
        or (row["C_local_coexistence_documented_strict"] == 1 or row["S_spatial_segregation_documented_strict"] == 1)
        and row["taxon_resolution_status"] != "resolved_unique"
    ]

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    write_csv(outdir / "v22_cs_source_evidence.csv", source_evidence, list(source_evidence[0].keys()))
    write_csv(outdir / "v22_cs_species_states.csv", species_rows, list(species_rows[0].keys()) if species_rows else ["canonical_name"])
    write_csv(outdir / "v22_cs_strict_freeze.csv", freeze_rows, list(species_rows[0].keys()) if species_rows else ["canonical_name"])
    write_csv(outdir / "v22_cs_unresolved_sources.csv", unresolved_sources, list(source_evidence[0].keys()))
    (outdir / "v22_cs_gbif_cache.json").write_text(json.dumps(cache, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    state_counts: dict[str, int] = defaultdict(int)
    for row in species_rows:
        state_counts[row["organization_state"]] += 1
    evidence_counts: dict[str, int] = defaultdict(int)
    for row in source_evidence:
        evidence_counts[row["strict_evidence_status"]] += 1
    positive_source_rows = [row for row in source_evidence if row["C_local_coexistence_documented_strict"] == 1 or row["S_spatial_segregation_documented_strict"] == 1]
    unresolved_positive_taxa = [row for row in positive_source_rows if row["taxon_resolution_status"] != "resolved_unique"]

    summary = {
        "status": "complete",
        "input_records": len(source_evidence),
        "historical_exact_sources_seen": historical_sources_seen,
        "gbif_names_queried": len(cache),
        "strict_positive_source_records": len(positive_source_rows),
        "C_positive_source_records": sum(int(row["C_local_coexistence_documented_strict"]) for row in source_evidence),
        "S_positive_source_records": sum(int(row["S_spatial_segregation_documented_strict"]) for row in source_evidence),
        "positive_sources_taxonomically_unresolved": len(unresolved_positive_taxa),
        "resolved_species_with_any_source": len(species_rows),
        "strict_freeze_species": len(freeze_rows),
        "species_state_counts": dict(state_counts),
        "source_evidence_status_counts": dict(evidence_counts),
        "definition": {
            "C": "explicit local coexistence of discrete natural floral-colour variants in the same population/site",
            "S": "explicit spatial segregation/structuring of discrete floral-colour variants or morph frequencies among geographic units",
            "zero_semantics": "not documented by the strict pass; not biological absence",
        },
        "semantic_guard": (
            "C and S are independent positive-evidence axes. S is not defined as the negation of C. "
            "A species may be positive on both axes at different spatial scales."
        ),
    }
    (outdir / "v22_cs_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False), flush=True)

    if len(source_evidence) != 12064:
        raise SystemExit(f"Expected 12064 source rows, found {len(source_evidence)}")
    if historical_sources_seen != 34:
        raise SystemExit(f"Expected 34 historical exact source rows, found {historical_sources_seen}")
    if not freeze_rows:
        raise SystemExit("No strict C/S freeze species were produced")


if __name__ == "__main__":
    main()

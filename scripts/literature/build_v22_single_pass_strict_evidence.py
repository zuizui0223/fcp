#!/usr/bin/env python3
"""Build a conservative single-pass evidence freeze from the canonical v2.2 record universe.

This path intentionally replaces the duplicate-review gate. It does NOT treat generic
morph terminology or absence of wording as a biological state. Only explicit positive
source wording can set a spatial evidence axis. Ambiguous records remain unresolved.

Historical labels are never used. The historical manifest is used only as a source->taxon
rescue map so exact benchmark papers are not lost when title/abstract taxon extraction
fails.
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
    r"(?is)(?:flower|floral|corolla|petal|tepal|labellum|showy\s+bract).{0,90}"
    r"(?:colou?r|pigment|reflectance|anthocyanin)|"
    r"(?:colou?r|pigment|reflectance|anthocyanin).{0,90}"
    r"(?:flower|floral|corolla|petal|tepal|labellum|showy\s+bract)"
)
DISCRETE_RE = re.compile(
    r"(?is)(?:colou?r\s+morphs?|flower[- ]?colou?r\s+polymorph|floral[- ]?colou?r\s+polymorph|"
    r"colou?r\s+polymorph|discrete\s+(?:flower|floral).{0,40}colou?r|"
    r"(?:two|multiple|several|both)\s+(?:flower|floral)?\s*colou?r\s+(?:forms?|morphs?))"
)
CONTINUOUS_RE = re.compile(
    r"(?is)(?:continuous\s+(?:flower|floral).{0,50}colou?r|"
    r"(?:hue|chroma|lightness|reflectance).{0,60}(?:gradient|continuous|variation))"
)
NATURAL_CONTEXT_RE = re.compile(
    r"(?is)\b(?:wild|natural\s+population|field\s+(?:population|site|study)|"
    r"population|populations|localit(?:y|ies)|sites?|geograph(?:ic|ical))\b"
)
CONFLICT_RE = re.compile(
    r"(?is)\b(?:cultivar|cultivars|horticultural\s+line|horticultural\s+lines|"
    r"commercial\s+variet|breeding\s+line|transgenic|gene[- ]edited|CRISPR|"
    r"induced\s+mutat|irradiat|mutagenesis|tissue\s+culture|somaclonal|"
    r"ontogenetic\s+colou?r\s+change|flower\s+colou?r\s+change\s+with\s+age)\b"
)
LOCAL_PATTERNS = [
    re.compile(r"(?is)\b(?:co[- ]?occur|coexist)(?:s|ed|ing)?\b.{0,90}\b(?:within|in)\b.{0,35}\b(?:population|site)\b"),
    re.compile(r"(?is)\b(?:within|in)\b.{0,35}\b(?:population|site)\b.{0,90}\b(?:co[- ]?occur|coexist)(?:s|ed|ing)?\b"),
    re.compile(r"(?is)\bsame\s+(?:natural\s+)?(?:population|site)\b.{0,100}\b(?:colou?r\s+morph|colou?r\s+form|polymorph)"),
    re.compile(r"(?is)\b(?:colou?r\s+morph|colou?r\s+form|polymorph).{0,100}\bsame\s+(?:natural\s+)?(?:population|site)\b"),
]
GEO_PATTERNS = [
    re.compile(r"(?is)\b(?:among|between)\s+(?:populations|sites|localities|regions|islands)\b"),
    re.compile(r"(?is)\b(?:geographic|geographical|spatial|latitudinal|longitudinal|altitudinal|elevational)\s+(?:variation|differentiation|structure|structuring|cline|pattern|gradient)\b"),
    re.compile(r"(?is)\b(?:variation|frequency|frequencies|differs?|varies?)\b.{0,70}\b(?:among|between)\b.{0,35}\b(?:populations|sites|localities|regions|islands)\b"),
    re.compile(r"(?is)\b(?:regional|geographic|geographical)\s+(?:restriction|replacement|segregation)\b"),
]


def match_snippet(patterns: list[re.Pattern[str]], text: str) -> str:
    for pat in patterns:
        m = pat.search(text)
        if m:
            a = max(0, m.start() - 90)
            b = min(len(text), m.end() + 90)
            return clean(text[a:b])
    return ""


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as h:
        return list(csv.DictReader(h))


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as h:
        w = csv.DictWriter(h, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def load_historical_taxon_map(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for row in read_csv(path):
        out[norm_source(row.get("source_id"))] = clean(row.get("canonical_name"))
    return out


def resolve_one(name: str, timeout: int, retries: int) -> tuple[str, dict[str, Any]]:
    try:
        result = resolve_gbif_name(name, timeout, retries)
    except Exception as exc:  # noqa: BLE001
        result = {
            "input_name": name,
            "accepted": False,
            "accepted_name": "",
            "family": "",
            "usage_key": "",
            "reason": f"resolution_error:{type(exc).__name__}",
        }
    return name, result


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--blind", required=True)
    p.add_argument("--key", required=True)
    p.add_argument("--historical-manifest", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--workers", type=int, default=24)
    p.add_argument("--timeout", type=int, default=25)
    p.add_argument("--retries", type=int, default=3)
    args = p.parse_args()

    blind = read_csv(Path(args.blind))
    key = read_csv(Path(args.key))
    if len(blind) != 12064 or len(key) != 12064:
        raise SystemExit(f"Expected 12064 blind/key rows, found {len(blind)}/{len(key)}")
    bmap = {r["record_review_id"]: r for r in blind}
    kmap = {r["record_review_id"]: r for r in key}
    if set(bmap) != set(kmap):
        raise SystemExit("Blind/key record IDs differ")

    historical_map = load_historical_taxon_map(Path(args.historical_manifest))
    source_rows: list[dict[str, Any]] = []
    names_to_resolve: set[str] = set()
    historical_sources_seen = 0

    for rid, row in bmap.items():
        hidden = kmap[rid]
        text = clean(f"{row.get('title','')} {row.get('abstract','')}")
        source = norm_source(row.get("source_id"))
        historical_taxon = historical_map.get(source, "")
        if historical_taxon:
            historical_sources_seen += 1

        display = bool(DISPLAY_RE.search(text))
        discrete = bool(DISCRETE_RE.search(text))
        continuous = bool(CONTINUOUS_RE.search(text))
        natural_context = bool(NATURAL_CONTEXT_RE.search(text))
        conflict = bool(CONFLICT_RE.search(text))
        local_snippet = match_snippet(LOCAL_PATTERNS, text) if display and discrete else ""
        geo_snippet = match_snippet(GEO_PATTERNS, text) if display else ""
        local = bool(local_snippet)
        geo = bool(geo_snippet)

        if conflict:
            eligibility = "conflict_unresolved"
        elif display and natural_context and (discrete or continuous or local or geo):
            eligibility = "eligible_high_confidence"
        elif display:
            eligibility = "display_relevant_unresolved"
        else:
            eligibility = "unresolved_or_outside_display_domain"

        detected = [clean(x) for x in clean(hidden.get("detected_binomial_strings")).split(";") if clean(x)]
        candidates: list[str] = []
        if historical_taxon:
            candidates.append(historical_taxon)
        for name in detected:
            if name not in candidates:
                candidates.append(name)
        # Only resolve names for potentially informative records. All other records stay visible in the audit.
        if historical_taxon or local or geo or eligibility in {"eligible_high_confidence", "conflict_unresolved"}:
            names_to_resolve.update(candidates[:8])

        if discrete and continuous:
            variation_form = "both"
        elif discrete:
            variation_form = "discrete"
        elif continuous:
            variation_form = "continuous"
        else:
            variation_form = "unclear"

        source_rows.append({
            "record_review_id": rid,
            "source_id": clean(row.get("source_id")),
            "title": clean(row.get("title")),
            "year": clean(row.get("year")),
            "eligibility_status": eligibility,
            "display_colour_signal": int(display),
            "natural_context_signal": int(natural_context),
            "eligibility_conflict_signal": int(conflict),
            "variation_form": variation_form,
            "local_coexistence_documented_strict": int(local),
            "geographic_structure_documented_strict": int(geo),
            "local_evidence_excerpt": local_snippet,
            "geographic_evidence_excerpt": geo_snippet,
            "historical_source_taxon_rescue": historical_taxon,
            "candidate_taxon_strings": ";".join(candidates[:8]),
        })

    workers = max(1, min(args.workers, 48))
    cache: dict[str, dict[str, Any]] = {}
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futures = {ex.submit(resolve_one, n, args.timeout, args.retries): n for n in sorted(names_to_resolve)}
        for i, fut in enumerate(as_completed(futures), 1):
            name, result = fut.result()
            cache[name] = result
            if i % 250 == 0:
                print({"gbif_resolved": i, "gbif_total": len(futures)}, flush=True)

    # Attach a taxon only when the source maps unambiguously to one accepted species.
    source_evidence: list[dict[str, Any]] = []
    for row in source_rows:
        candidates = [x for x in row["candidate_taxon_strings"].split(";") if x]
        accepted: dict[str, dict[str, Any]] = {}
        for name in candidates:
            r = cache.get(name, {})
            if r.get("accepted") and r.get("accepted_name"):
                accepted[str(r["accepted_name"])] = r
        if len(accepted) == 1:
            accepted_name, resolved = next(iter(accepted.items()))
            taxon_status = "resolved_unique"
            family = clean(resolved.get("family"))
            usage_key = clean(resolved.get("usage_key"))
        elif len(accepted) > 1:
            accepted_name = ""
            family = ""
            usage_key = ""
            taxon_status = "multiple_accepted_taxa_unresolved"
        else:
            accepted_name = ""
            family = ""
            usage_key = ""
            taxon_status = "no_accepted_taxon_resolved"
        x = dict(row)
        x.update({
            "accepted_name": accepted_name,
            "family": family,
            "gbif_usage_key": usage_key,
            "taxon_resolution_status": taxon_status,
        })
        source_evidence.append(x)

    species_sources: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in source_evidence:
        if row["accepted_name"]:
            species_sources[row["accepted_name"]].append(row)

    species_rows: list[dict[str, Any]] = []
    for species, rows in sorted(species_sources.items()):
        eligible = [r for r in rows if r["eligibility_status"] == "eligible_high_confidence"]
        local_sources = [r for r in eligible if int(r["local_coexistence_documented_strict"]) == 1]
        geo_sources = [r for r in eligible if int(r["geographic_structure_documented_strict"]) == 1]
        conflicts = [r for r in rows if r["eligibility_status"] == "conflict_unresolved"]
        local = bool(local_sources)
        geo = bool(geo_sources)
        if local and geo:
            state = "mixed_evidence"
        elif local:
            state = "within_evidence_only"
        elif geo:
            state = "among_evidence_only"
        else:
            state = "unresolved"
        family = next((r["family"] for r in rows if r["family"]), "")
        species_rows.append({
            "canonical_name": species,
            "family": family,
            "documented_state": state,
            "local_coexistence_documented": int(local),
            "geographic_structure_documented": int(geo),
            "n_resolved_sources": len(rows),
            "n_high_confidence_eligible_sources": len(eligible),
            "n_local_positive_sources": len(local_sources),
            "n_geographic_positive_sources": len(geo_sources),
            "n_eligibility_conflict_sources": len(conflicts),
            "local_source_ids": ";".join(r["source_id"] for r in local_sources),
            "geographic_source_ids": ";".join(r["source_id"] for r in geo_sources),
        })

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    source_fields = list(source_evidence[0].keys())
    species_fields = list(species_rows[0].keys()) if species_rows else ["canonical_name"]
    write_csv(outdir / "v22_single_pass_source_evidence.csv", source_evidence, source_fields)
    write_csv(outdir / "v22_single_pass_species_states.csv", species_rows, species_fields)
    (outdir / "v22_single_pass_gbif_cache.json").write_text(json.dumps(cache, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    state_counts: dict[str, int] = defaultdict(int)
    for row in species_rows:
        state_counts[row["documented_state"]] += 1
    source_status_counts: dict[str, int] = defaultdict(int)
    for row in source_evidence:
        source_status_counts[row["eligibility_status"]] += 1
    summary = {
        "status": "complete",
        "input_records": len(blind),
        "historical_exact_sources_seen": historical_sources_seen,
        "candidate_names_queried_to_gbif": len(names_to_resolve),
        "resolved_species_with_any_source": len(species_rows),
        "species_state_counts": dict(sorted(state_counts.items())),
        "source_eligibility_counts": dict(sorted(source_status_counts.items())),
        "strict_local_positive_sources": sum(int(r["local_coexistence_documented_strict"]) for r in source_evidence),
        "strict_geographic_positive_sources": sum(int(r["geographic_structure_documented_strict"]) for r in source_evidence),
        "semantic_guard": (
            "Single-pass strict documented-evidence audit. Generic morph terminology cannot set local coexistence. "
            "Missing wording is never biological absence. Ambiguous eligibility/taxon attribution remains unresolved. "
            "Historical spatial labels and climatic outcomes are not used."
        ),
    }
    (outdir / "v22_single_pass_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)

    if historical_sources_seen != 34:
        raise SystemExit(f"Expected all 34 historical source records, found {historical_sources_seen}")
    errors = [r for r in cache.values() if str(r.get("reason", "")).startswith("resolution_error:")]
    if len(errors) > max(10, int(0.05 * max(1, len(cache)))):
        raise SystemExit(f"Too many GBIF resolution errors: {len(errors)}/{len(cache)}")


if __name__ == "__main__":
    main()

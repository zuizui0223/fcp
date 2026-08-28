#!/usr/bin/env python3
"""Derive the primary floral-display FCP core from the expanded all-species universe.

Core membership is independent of C/S positivity. It requires at least one already-FCP-
eligible source that directly documents discrete intraspecific *display* colour variation.
Sources restricted to stigma/gynoecium/anther/pollen colour, and review/synthesis-only
records, do not define core membership. Expanded-universe analyses retain them as a
prespecified sensitivity.

For a core species, C/S evidence is aggregated only across core-eligible display-colour
sources, so a non-display source cannot change a display-colour organization state.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Any

COLORS = r"(?:white|yellow|pink|purple|blue|red|orange|green|cream|violet|magenta)"
DIRECT_PATTERNS: list[tuple[str, re.Pattern[str]]] = [
    ("explicit_flower_colour_polymorphism", re.compile(
        r"(?is)\b(?:flower|floral|petal|corolla|perianth|labellum|bract|inflorescence)[- ]?colou?r\s+polymorph\w*\b"
    )),
    ("explicit_flower_colour_morphs", re.compile(
        r"(?is)\b(?:flower|floral|petal|corolla|perianth|labellum|bract|inflorescence)[- ]?colou?r\s+(?:morphs?|forms?|variants?|dimorph\w*)\b"
    )),
    ("colour_morphs_with_display_context", re.compile(
        r"(?is)\bcolou?r\s+morphs?\b.{0,120}\b(?:flower|floral|petal|corolla|perianth|labellum|bract|inflorescence|population)\b"
    )),
    ("display_context_with_colour_morphs", re.compile(
        r"(?is)\b(?:flower|floral|petal|corolla|perianth|labellum|bract|inflorescence|population)\b.{0,120}\bcolou?r\s+morphs?\b"
    )),
    ("geographic_flower_colour_variation", re.compile(
        r"(?is)\b(?:geograph\w*|spatial)\s+(?:variation|differentiation|distribution)\b.{0,80}\b(?:flower|floral|petal|corolla|perianth|labellum)\s+colou?r\b"
    )),
    ("intraspecific_flower_colour_differentiation", re.compile(
        r"(?is)\bintraspecific\s+(?:variation|differentiation\w*)\b.{0,100}\b(?:flower|floral|flowering|petal)\s+colou?r"
    )),
    ("petal_colour_morph_cline", re.compile(
        r"(?is)\b(?:morph[- ]ratio\s+cline|morph\s+frequenc\w*|variants?)\b.{0,120}\b(?:keel\s+)?petal\s+colou?r\b|\b(?:keel\s+)?petal\s+colou?r\b.{0,120}\b(?:morph[- ]ratio\s+cline|morph\s+frequenc\w*|variants?)\b"
    )),
    ("named_colour_flower_categories", re.compile(
        rf"(?is)\b{COLORS}\b.{{0,55}}\b(?:and|or|versus|vs\.?)\b.{{0,55}}\b{COLORS}\b.{{0,70}}\b(?:flowers?|flowered\s+individuals?|morphs?|forms?)\b"
    )),
    ("named_colour_reverse_categories", re.compile(
        rf"(?is)\b(?:flowers?|flowered\s+individuals?|morphs?|forms?)\b.{{0,70}}\b{COLORS}\b.{{0,55}}\b(?:and|or|versus|vs\.?)\b.{{0,55}}\b{COLORS}\b"
    )),
    ("heterocyanic_population", re.compile(r"(?is)\bheterocyanic\s+populations?\b")),
    ("within_population_flower_colour", re.compile(
        r"(?is)\bwithin[- ]population\b.{0,100}\b(?:flower|floral|petal|corolla|perianth)\s+colou?r\b"
    )),
    ("flower_colour_within_population", re.compile(
        r"(?is)\b(?:flower|floral|petal|corolla|perianth)\s+colou?r\b.{0,100}\bwithin[- ]population\b"
    )),
    ("wild_population_display_colour", re.compile(
        r"(?is)\bwild\s+populations?\b.{0,150}\b(?:petal|flower|floral|bract|inflorescence)\s+colou?r\s+(?:variants?|morphs?|polymorph\w*)\b"
    )),
]

REVIEW_SYNTHESIS_RE = re.compile(
    r"(?is)^\s*(?:review\b|research\s+progress\b|paradigm\s+shifts\b.*\bintroduction\b|"
    r"systematic\s+review\b|meta[- ]analysis\b)|\bcompiled\s+(?:information|data)\s+from\s+(?:hundreds|many)\s+of\s+species\b"
)
NONDISPLAY_ONLY_RE = re.compile(
    r"(?is)\b(?:stigma|gynoecium|anther|pollen|androecium)\s+colou?r\s+polymorph\w*\b"
)
DISPLAY_STRUCTURE_RE = re.compile(
    r"(?is)\b(?:petal|corolla|perianth|tepal|labellum|bract|inflorescence)\s+colou?r\b|"
    r"\b(?:flower|floral)[- ]?colou?r\s+(?:morphs?|forms?|variants?|polymorph\w*)\b"
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as h:
        return list(csv.DictReader(h))


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as h:
        w=csv.DictWriter(h, fieldnames=fields); w.writeheader(); w.writerows(rows)


def direct_reason(text: str, title: str) -> str:
    if REVIEW_SYNTHESIS_RE.search(title or "") or REVIEW_SYNTHESIS_RE.search(text or ""):
        return ""
    if NONDISPLAY_ONLY_RE.search(text or "") and not DISPLAY_STRUCTURE_RE.search(text or ""):
        return ""
    for name, pattern in DIRECT_PATTERNS:
        if pattern.search(text or ""):
            return name
    return ""


def main() -> None:
    p=argparse.ArgumentParser()
    p.add_argument("--source-audit", required=True)
    p.add_argument("--record-screen", required=True)
    p.add_argument("--outdir", required=True)
    args=p.parse_args()

    sources=read_csv(Path(args.source_audit))
    records=read_csv(Path(args.record_screen))
    rmap={r['record_review_id']:r for r in records}

    core_source_rows=[]
    expanded_eligible_by_species: dict[str,list[dict[str,str]]]=defaultdict(list)
    core_eligible_by_species: dict[str,list[dict[str,str]]]=defaultdict(list)
    exclusion_rows=[]

    for row in sources:
        if str(row.get('FCP_eligible_source'))!='1' or row.get('taxon_resolution_status')!='resolved_unique' or not row.get('accepted_name'):
            continue
        species=row['accepted_name']
        expanded_eligible_by_species[species].append(row)
        rr=rmap.get(row['record_review_id'],{})
        title=str(rr.get('title','') or '')
        text=' '.join(str(rr.get(k,'') or '') for k in ('title','abstract'))
        reason=direct_reason(text,title)
        if reason:
            x=dict(row); x['core_membership_reason']=reason
            core_source_rows.append(x)
            core_eligible_by_species[species].append(row)
        else:
            exclusion_rows.append({
                'accepted_name':species,
                'source_id':row.get('source_id',''),
                'title':title,
                'review_or_synthesis_flag':int(bool(REVIEW_SYNTHESIS_RE.search(title) or REVIEW_SYNTHESIS_RE.search(text))),
                'nondisplay_only_flag':int(bool(NONDISPLAY_ONLY_RE.search(text) and not DISPLAY_STRUCTURE_RE.search(text))),
                'reason':'no_high_specificity_display_FCP_membership_evidence',
            })

    core_species_names=set(core_eligible_by_species)
    species_rows=[]
    for species in sorted(core_species_names):
        rows=core_eligible_by_species[species]
        c=[r for r in rows if str(r.get('C_local_coexistence_documented_strict'))=='1']
        s=[r for r in rows if str(r.get('S_spatial_segregation_documented_strict'))=='1']
        C=bool(c); S=bool(s)
        state='coexistence_and_segregation' if C and S else 'local_coexistence_only' if C else 'spatial_segregation_only' if S else 'organization_unresolved'
        family=next((r.get('family','') for r in rows if r.get('family')), '')
        years=[int(r['year']) for r in rows if str(r.get('year','')).isdigit()]
        species_rows.append({
            'canonical_name':species,
            'family':family,
            'organization_state':state,
            'C_local_coexistence_documented':int(C),
            'S_spatial_segregation_documented':int(S),
            'n_FCP_eligible_sources':len(rows),
            'n_expanded_FCP_eligible_sources':len(expanded_eligible_by_species[species]),
            'n_core_membership_sources':len(rows),
            'n_C_positive_sources':len(c),
            'n_S_positive_sources':len(s),
            'earliest_source_year':min(years) if years else '',
            'latest_source_year':max(years) if years else '',
            'source_year_span':max(years)-min(years) if years else '',
            'core_source_ids':';'.join(r['source_id'] for r in rows),
            'C_source_ids':';'.join(r['source_id'] for r in c),
            'S_source_ids':';'.join(r['source_id'] for r in s),
        })

    out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True)
    write_csv(out/'v22_core_fcp_membership_sources.csv',core_source_rows,list(core_source_rows[0].keys()) if core_source_rows else ['accepted_name'])
    write_csv(out/'v22_core_fcp_species_universe.csv',species_rows,list(species_rows[0].keys()) if species_rows else ['canonical_name'])
    write_csv(out/'v22_core_fcp_excluded_expanded_sources.csv',exclusion_rows,list(exclusion_rows[0].keys()) if exclusion_rows else ['accepted_name'])
    states=defaultdict(int)
    for r in species_rows: states[r['organization_state']]+=1
    expanded_species=set(expanded_eligible_by_species)
    summary={
        'status':'complete',
        'core_protocol_version':'display-core-v2',
        'expanded_FCP_species_with_resolved_eligible_source':len(expanded_species),
        'expanded_eligible_resolved_sources':sum(len(v) for v in expanded_eligible_by_species.values()),
        'core_membership_source_records':len(core_source_rows),
        'core_FCP_species_universe':len(species_rows),
        'expanded_species_excluded_from_core':len(expanded_species-core_species_names),
        'expanded_species_excluded_from_core_names':sorted(expanded_species-core_species_names),
        'species_state_counts':dict(states),
        'C_positive_species':sum(int(r['C_local_coexistence_documented']) for r in species_rows),
        'S_positive_species':sum(int(r['S_spatial_segregation_documented']) for r in species_rows),
        'organization_unresolved_species':sum(r['organization_state']=='organization_unresolved' for r in species_rows),
        'membership_rule':'At least one resolved eligible source directly supports discrete intraspecific floral-display colour variation; review/synthesis-only and non-display sexual-organ-only sources cannot define core membership. C/S positivity is not used for membership.',
        'axis_rule':'C/S evidence is aggregated only across core-eligible display-colour sources for a core-member species.',
        'expanded_role':'prespecified sensitivity retaining the broader eligibility gate',
    }
    (out/'v22_core_fcp_universe_summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(summary,indent=2))
    if len(species_rows)<20:
        raise SystemExit(f'Core universe unexpectedly small: {len(species_rows)}')
    if not any(r['organization_state']=='organization_unresolved' for r in species_rows):
        raise SystemExit('Core membership leaked C/S positivity requirement')
    if not expanded_species-core_species_names:
        raise SystemExit('Core/expanded sensitivity boundary collapsed; review exclusions may not be applied')

if __name__=='__main__':
    main()

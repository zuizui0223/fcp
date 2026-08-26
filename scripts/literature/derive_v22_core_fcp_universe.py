#!/usr/bin/env python3
"""Derive a high-specificity core FCP universe from the expanded all-species audit.

Core membership is based on direct flower-colour polymorphism wording in title/abstract,
not on C or S positivity. Once a species has at least one core membership source, all of
its already-eligible FCP sources may contribute positive C/S evidence.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Any

COLORS = r"(?:white|yellow|pink|purple|blue|red|orange|green|cream|violet)"
DIRECT_PATTERNS: list[tuple[str, re.Pattern[str]]] = [
    ("explicit_flower_colour_polymorphism", re.compile(
        r"(?is)\b(?:flower|floral|petal|corolla|perianth|labellum)[- ]?colou?r\s+polymorph\w*\b"
    )),
    ("explicit_flower_colour_morphs", re.compile(
        r"(?is)\b(?:flower|floral|petal|corolla|perianth|labellum)[- ]?colou?r\s+(?:morphs?|forms?|variants?|dimorph\w*)\b"
    )),
    ("colour_morphs_with_flower_context", re.compile(
        r"(?is)\bcolou?r\s+morphs?\b.{0,100}\b(?:flower|floral|petal|corolla|perianth|labellum|population)\b"
    )),
    ("flower_context_with_colour_morphs", re.compile(
        r"(?is)\b(?:flower|floral|petal|corolla|perianth|labellum|population)\b.{0,100}\bcolou?r\s+morphs?\b"
    )),
    ("named_colour_flower_categories", re.compile(
        rf"(?is)\b{COLORS}\b.{{0,55}}\b(?:and|or|versus|vs\.?)\b.{{0,55}}\b{COLORS}\b.{{0,60}}\b(?:flowers?|flowered\s+individuals?|morphs?|forms?)\b"
    )),
    ("named_colour_reverse_categories", re.compile(
        rf"(?is)\b(?:flowers?|flowered\s+individuals?|morphs?|forms?)\b.{{0,60}}\b{COLORS}\b.{{0,55}}\b(?:and|or|versus|vs\.?)\b.{{0,55}}\b{COLORS}\b"
    )),
    ("heterocyanic_population", re.compile(r"(?is)\bheterocyanic\s+populations?\b")),
    ("within_population_flower_colour", re.compile(
        r"(?is)\bwithin[- ]population\b.{0,80}\b(?:flower|floral|petal)\s+colou?r\b"
    )),
    ("flower_colour_within_population", re.compile(
        r"(?is)\b(?:flower|floral|petal)\s+colou?r\b.{0,80}\bwithin[- ]population\b"
    )),
    ("wild_population_petals", re.compile(
        r"(?is)\bwild\s+populations?\b.{0,120}\b(?:petal|flower|floral)\s+colou?r\s+(?:variants?|morphs?|polymorph\w*)\b"
    )),
]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as h:
        return list(csv.DictReader(h))


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as h:
        w=csv.DictWriter(h, fieldnames=fields); w.writeheader(); w.writerows(rows)


def direct_reason(text: str) -> str:
    for name, pattern in DIRECT_PATTERNS:
        if pattern.search(text):
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
    eligible_by_species: dict[str,list[dict[str,str]]]=defaultdict(list)
    for row in sources:
        if str(row.get('FCP_eligible_source'))!='1' or row.get('taxon_resolution_status')!='resolved_unique' or not row.get('accepted_name'):
            continue
        species=row['accepted_name']
        eligible_by_species[species].append(row)
        rr=rmap.get(row['record_review_id'],{})
        text=' '.join(str(rr.get(k,'') or '') for k in ('title','abstract'))
        reason=direct_reason(text)
        if reason:
            x=dict(row); x['core_membership_reason']=reason
            core_source_rows.append(x)

    core_species_names={r['accepted_name'] for r in core_source_rows}
    species_rows=[]
    for species in sorted(core_species_names):
        rows=eligible_by_species[species]
        membership=[r for r in core_source_rows if r['accepted_name']==species]
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
            'n_core_membership_sources':len(membership),
            'n_C_positive_sources':len(c),
            'n_S_positive_sources':len(s),
            'earliest_source_year':min(years) if years else '',
            'latest_source_year':max(years) if years else '',
            'source_year_span':max(years)-min(years) if years else '',
            'core_source_ids':';'.join(r['source_id'] for r in membership),
            'FCP_source_ids':';'.join(r['source_id'] for r in rows),
            'C_source_ids':';'.join(r['source_id'] for r in c),
            'S_source_ids':';'.join(r['source_id'] for r in s),
        })

    out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True)
    write_csv(out/'v22_core_fcp_membership_sources.csv',core_source_rows,list(core_source_rows[0].keys()) if core_source_rows else ['accepted_name'])
    write_csv(out/'v22_core_fcp_species_universe.csv',species_rows,list(species_rows[0].keys()) if species_rows else ['canonical_name'])
    states=defaultdict(int)
    for r in species_rows: states[r['organization_state']]+=1
    summary={
        'status':'complete',
        'expanded_eligible_resolved_sources':sum(1 for r in sources if str(r.get('FCP_eligible_source'))=='1' and r.get('taxon_resolution_status')=='resolved_unique'),
        'core_membership_source_records':len(core_source_rows),
        'core_FCP_species_universe':len(species_rows),
        'species_state_counts':dict(states),
        'C_positive_species':sum(int(r['C_local_coexistence_documented']) for r in species_rows),
        'S_positive_species':sum(int(r['S_spatial_segregation_documented']) for r in species_rows),
        'organization_unresolved_species':sum(r['organization_state']=='organization_unresolved' for r in species_rows),
        'membership_rule':'At least one already-eligible source contains a high-specificity direct flower-colour polymorphism/morph statement. C/S positivity is not used for membership.',
        'axis_rule':'C/S evidence is aggregated across all eligible sources for a core-member species.',
    }
    (out/'v22_core_fcp_universe_summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(summary,indent=2))
    if len(species_rows)<20:
        raise SystemExit(f'Core universe unexpectedly small: {len(species_rows)}')
    if not any(r['organization_state']=='organization_unresolved' for r in species_rows):
        raise SystemExit('Core membership leaked C/S positivity requirement')

if __name__=='__main__':
    main()

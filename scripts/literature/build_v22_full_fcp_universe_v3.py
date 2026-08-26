#!/usr/bin/env python3
"""All-species FCP universe with post-GBIF hierarchical focal-taxon fallback.

Species membership is independent of C/S evidence. Focal taxon attribution proceeds by
source structure *after* GBIF resolution:
  exact historical-source rescue -> title -> abstract -> navigation hint.
A tier is accepted only when it resolves to exactly one accepted plant species. If a
higher tier resolves to zero species, the next tier is tried. If it resolves to multiple
species, the source is left ambiguous rather than guessing.
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

import build_v22_full_fcp_universe as v1

base=v1.base
BINOMIAL_RE=re.compile(r"([A-Z][a-z][A-Za-z-]{1,})\s+([a-z][a-z-]{2,})")


def clean(x: Any)->str:
    return base.clean(x)


def extract(text: str, limit: int)->list[str]:
    out=[]; seen=set()
    for genus, epithet in BINOMIAL_RE.findall(text or ""):
        for name in [f"{genus} {epithet}"] + ([f"{genus} {epithet[:-2]}"] if epithet.endswith('isis') and len(epithet)>6 else []):
            if name not in seen:
                seen.add(name); out.append(name)
                if len(out)>=limit: return out
    return out


def nav_candidates(hidden: dict[str,str])->list[str]:
    return [clean(x) for x in clean(hidden.get('detected_binomial_strings')).split(';') if clean(x)][:12]


def read_csv(path: Path)->list[dict[str,str]]:
    with path.open(newline='',encoding='utf-8') as h: return list(csv.DictReader(h))


def write_csv(path: Path, rows: list[dict[str,Any]], fields: list[str])->None:
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='',encoding='utf-8') as h:
        w=csv.DictWriter(h,fieldnames=fields); w.writeheader(); w.writerows(rows)


def accepted_in_tier(names: list[str], cache: dict[str,dict[str,Any]])->dict[str,dict[str,Any]]:
    out={}
    for name in names:
        r=cache.get(name,{})
        if r.get('accepted') and r.get('accepted_name'):
            out[str(r['accepted_name'])]=r
    return out


def main()->None:
    p=argparse.ArgumentParser()
    p.add_argument('--blind',required=True); p.add_argument('--key',required=True)
    p.add_argument('--historical-manifest',required=True); p.add_argument('--outdir',required=True)
    p.add_argument('--workers',type=int,default=32); p.add_argument('--timeout',type=int,default=20); p.add_argument('--retries',type=int,default=2)
    args=p.parse_args()

    blind=read_csv(Path(args.blind)); key=read_csv(Path(args.key))
    if len(blind)!=12064 or len(key)!=12064: raise SystemExit('Expected 12064 rows')
    bmap={r['record_review_id']:r for r in blind}; kmap={r['record_review_id']:r for r in key}
    if set(bmap)!=set(kmap): raise SystemExit('Blind/key IDs differ')
    historical=base.load_historical_taxon_map(Path(args.historical_manifest))

    rows=[]; names=set(); hist_seen=0
    for rid,row in bmap.items():
        hidden=kmap[rid]; title=clean(row.get('title')); abstract=clean(row.get('abstract')); wt=clean(row.get('work_type')).lower()
        source_key=base.norm_source(row.get('source_id')); hist=historical.get(source_key,'')
        if hist: hist_seen+=1
        sig=v1.source_eligibility(title,abstract,wt); text=sig['text']
        cex=base.positive_context(base.C_PATTERNS,text,axis='C') if sig['eligible'] else ''
        sex=base.positive_context(base.S_PATTERNS,text,axis='S') if sig['eligible'] else ''
        cp=bool(cex); sp=bool(sex)
        if sig['eligible']:
            status='eligible_C_and_S_positive' if cp and sp else 'eligible_C_positive' if cp else 'eligible_S_positive' if sp else 'eligible_spatial_unresolved'
        elif sig['hard_conflict']: status='artificial_or_conflict_excluded'
        elif not sig['primary']: status='nonprimary_excluded'
        elif sig['community_only']: status='community_level_excluded'
        else: status='not_strict_FCP_eligible'

        hist_t=[hist] if hist and sig['eligible'] else []
        title_t=extract(title,12) if sig['eligible'] and not hist_t else []
        abstract_t=extract(abstract[:5000],24) if sig['eligible'] and not hist_t else []
        nav_t=nav_candidates(hidden) if sig['eligible'] and not hist_t else []
        for tier in (hist_t,title_t,abstract_t,nav_t): names.update(tier)
        rows.append({
            'record_review_id':rid,'source_id':clean(row.get('source_id')),'title':title,'year':clean(row.get('year')),'work_type':wt,
            'FCP_source_status':status,'FCP_eligible_source':int(sig['eligible']),
            'primary_source_signal':int(sig['primary']),'display_colour_signal':int(sig['display']),
            'discrete_polymorphism_signal':int(sig['discrete']),'intraspecific_signal':int(sig['intraspecific']),
            'natural_population_context_signal':int(sig['population_context']),'hard_conflict_signal':int(sig['hard_conflict']),
            'community_level_signal':int(sig['community_only']),
            'C_local_coexistence_documented_strict':int(cp),'S_spatial_segregation_documented_strict':int(sp),
            'C_evidence_excerpt':cex,'S_evidence_excerpt':sex,'historical_source_taxon_rescue':hist,
            'historical_candidate_taxa':';'.join(hist_t),'title_candidate_taxa':';'.join(title_t),
            'abstract_candidate_taxa':';'.join(abstract_t),'navigation_candidate_taxa':';'.join(nav_t),
        })

    cache={}; workers=max(1,min(args.workers,48))
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs={ex.submit(v1.resolve_one,n,args.timeout,args.retries):n for n in sorted(names)}
        for i,fut in enumerate(as_completed(futs),1):
            n,r=fut.result(); cache[n]=r
            if i%100==0: print({'gbif_resolved':i,'gbif_total':len(futs)},flush=True)

    resolved=[]; tier_counts=defaultdict(int)
    for row in rows:
        chosen_name=family=usage=''; status='not_applicable'; chosen_tier=''
        if int(row['FCP_eligible_source'])==1:
            tiers=[
                ('historical_exact',[x for x in row['historical_candidate_taxa'].split(';') if x]),
                ('title',[x for x in row['title_candidate_taxa'].split(';') if x]),
                ('abstract',[x for x in row['abstract_candidate_taxa'].split(';') if x]),
                ('navigation_hint',[x for x in row['navigation_candidate_taxa'].split(';') if x]),
            ]
            for tier,names_t in tiers:
                if not names_t: continue
                acc=accepted_in_tier(names_t,cache)
                if len(acc)==1:
                    chosen_name,r=next(iter(acc.items())); family=clean(r.get('family')); usage=clean(r.get('usage_key'))
                    status='resolved_unique'; chosen_tier=tier; break
                if len(acc)>1:
                    status=f'multiple_accepted_taxa_unresolved_at_{tier}'; chosen_tier=tier; break
            if status=='not_applicable': status='no_accepted_taxon_resolved'
        tier_counts[chosen_tier or status]+=1
        x=dict(row); x.update({'accepted_name':chosen_name,'family':family,'gbif_usage_key':usage,
                               'taxon_resolution_status':status,'taxon_resolution_tier':chosen_tier})
        resolved.append(x)

    by=defaultdict(list)
    for r in resolved:
        if int(r['FCP_eligible_source'])==1 and r['taxon_resolution_status']=='resolved_unique' and r['accepted_name']:
            by[r['accepted_name']].append(r)
    species=[]
    for name,rr in sorted(by.items()):
        cs=[r for r in rr if int(r['C_local_coexistence_documented_strict'])==1]
        ss=[r for r in rr if int(r['S_spatial_segregation_documented_strict'])==1]
        C=bool(cs); S=bool(ss); state='coexistence_and_segregation' if C and S else 'local_coexistence_only' if C else 'spatial_segregation_only' if S else 'organization_unresolved'
        years=[int(r['year']) for r in rr if str(r.get('year','')).isdigit()]
        fam=next((r['family'] for r in rr if r['family']),'')
        species.append({'canonical_name':name,'family':fam,'organization_state':state,
                        'C_local_coexistence_documented':int(C),'S_spatial_segregation_documented':int(S),
                        'n_FCP_eligible_sources':len(rr),'n_C_positive_sources':len(cs),'n_S_positive_sources':len(ss),
                        'earliest_source_year':min(years) if years else '','latest_source_year':max(years) if years else '',
                        'source_year_span':max(years)-min(years) if years else '',
                        'FCP_source_ids':';'.join(r['source_id'] for r in rr),'C_source_ids':';'.join(r['source_id'] for r in cs),
                        'S_source_ids':';'.join(r['source_id'] for r in ss)})
    informative=[r for r in species if int(r['C_local_coexistence_documented']) or int(r['S_spatial_segregation_documented'])]
    unresolved=[r for r in species if r['organization_state']=='organization_unresolved']
    taxon_unresolved=[r for r in resolved if int(r['FCP_eligible_source'])==1 and r['taxon_resolution_status']!='resolved_unique']

    out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True)
    write_csv(out/'v22_full_fcp_source_audit.csv',resolved,list(resolved[0].keys()))
    fields=list(species[0].keys()) if species else ['canonical_name']
    write_csv(out/'v22_full_fcp_species_universe.csv',species,fields)
    write_csv(out/'v22_full_fcp_informative_states.csv',informative,fields)
    write_csv(out/'v22_full_fcp_organization_unresolved.csv',unresolved,fields)
    write_csv(out/'v22_full_fcp_taxon_unresolved_sources.csv',taxon_unresolved,list(resolved[0].keys()))
    (out/'v22_full_fcp_gbif_cache.json').write_text(json.dumps(cache,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    states=defaultdict(int); source_status=defaultdict(int)
    for r in species: states[r['organization_state']]+=1
    for r in resolved: source_status[r['FCP_source_status']]+=1
    summary={'status':'complete','protocol_version':'full-fcp-universe-v3-post-gbif-hierarchical-taxon',
             'input_records':len(resolved),'historical_exact_sources_seen_for_taxon_rescue_only':hist_seen,
             'eligible_FCP_source_records':sum(int(r['FCP_eligible_source']) for r in resolved),
             'eligible_FCP_sources_taxonomically_unresolved':len(taxon_unresolved),'gbif_names_queried':len(cache),
             'FCP_species_universe':len(species),'C_positive_species':sum(int(r['C_local_coexistence_documented']) for r in species),
             'S_positive_species':sum(int(r['S_spatial_segregation_documented']) for r in species),'informative_C_or_S_species':len(informative),
             'organization_unresolved_species':len(unresolved),'species_state_counts':dict(states),'source_status_counts':dict(source_status),
             'taxon_resolution_tier_counts':dict(tier_counts),
             'universe_definition':'Species enter because at least one primary source supports discrete intraspecific floral-display colour variation in a natural/population context without an artificial/conflict exclusion; C/S positivity is not an inclusion criterion.',
             'taxon_attribution':'post-GBIF hierarchical fallback: historical exact source -> title -> abstract -> navigation hint; zero-accepted tiers fall through, multi-accepted tiers remain unresolved',
             'axis_definition':{'C':'explicit local coexistence of discrete natural floral-colour variants in the same population/site','S':'explicit spatial segregation/structuring of colour variants or morph frequencies among geographic units','zero_semantics':'not documented by this strict pass; not biological absence'},
             'historical_34_role':'taxon rescue and later historical sensitivity only; never universe membership'}
    (out/'v22_full_fcp_universe_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(summary,ensure_ascii=False,indent=2),flush=True)
    if len(resolved)!=12064 or hist_seen!=34: raise SystemExit('retrieval/rescue boundary failed')
    if len(species)<=34 or not unresolved: raise SystemExit('independent universe boundary failed')

if __name__=='__main__': main()

#!/usr/bin/env python3
"""Parallel species-level wrapper for paginated GBIF occurrence sampling.

Biological/QC rules are identical to build_gbif_occurrence_sample_paginated.py. Only the
outer loop over species is parallelized; each species still pages its own GBIF query
sequentially. Results are merged in canonical-name order for deterministic output.
"""
from __future__ import annotations

import argparse, csv, json
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import build_gbif_occurrence_sample_paginated as base


def process_one(item: dict[str,str], args: argparse.Namespace) -> dict[str,Any]:
    name=item['canonical_name']
    local_exclusions=Counter()
    try:
        match=base.resolve_taxon(name,timeout=args.timeout,retries=args.retries)
        raw,requests,reported=base.occurrence_pages(
            match['taxon_key'],max_records=args.max_records,page_size=args.page_size,
            timeout=args.timeout,retries=args.retries,delay=args.delay)
    except Exception as exc:
        return {'name':name,'error':{'canonical_name':name,'error':f'{type(exc).__name__}: {exc}'[:500]},
                'rows':[],'audit':None,'requests':0,'reported':0,'exclusions':local_exclusions}
    seen=set(); rows=[]; basis_counts=Counter()
    for record in raw:
        basis=str(record.get('basisOfRecord') or '')
        if basis not in base.ACCEPTED_BASIS:
            local_exclusions['basis_of_record']+=1; continue
        coordinate=base.finite_coordinate(record)
        if coordinate is None:
            local_exclusions['invalid_coordinate']+=1; continue
        lat,lon=coordinate
        if not base.uncertainty_ok(record,args.max_coordinate_uncertainty_m):
            local_exclusions['coordinate_uncertainty']+=1; continue
        key=base.record_key(record,lat,lon,args.dedup_decimals)
        if key in seen:
            local_exclusions['spatial_duplicate']+=1; continue
        seen.add(key); basis_counts[basis]+=1
        rows.append({**item,'role':'focal','focal_species':name,'match_level':'self',
                     'gbif_key':record.get('key',''),'gbif_taxon_key':match['taxon_key'],
                     'gbif_accepted_name':match['accepted_name'],'decimalLatitude':lat,
                     'decimalLongitude':lon,'coordinateUncertaintyInMeters':record.get('coordinateUncertaintyInMeters',''),
                     'year':record.get('year',''),'basisOfRecord':basis,'datasetKey':record.get('datasetKey','')})
    audit={'canonical_name':name,'family':item['family'],**match,'gbif_reported_count':reported,
           'records_retrieved_before_filters':len(raw),'records_retained':len(rows),
           'basis_counts':json.dumps(dict(sorted(basis_counts.items())),sort_keys=True)}
    return {'name':name,'error':None,'rows':rows,'audit':audit,'requests':requests+1,
            'reported':reported,'exclusions':local_exclusions}


def main() -> None:
    p=argparse.ArgumentParser()
    p.add_argument('--species',required=True); p.add_argument('--out',required=True); p.add_argument('--qc-out',required=True); p.add_argument('--taxon-audit-out',required=True)
    p.add_argument('--max-records',type=int,default=3000); p.add_argument('--page-size',type=int,default=300); p.add_argument('--dedup-decimals',type=int,default=3)
    p.add_argument('--max-coordinate-uncertainty-m',type=float,default=20000); p.add_argument('--timeout',type=int,default=45); p.add_argument('--retries',type=int,default=4); p.add_argument('--delay',type=float,default=.10); p.add_argument('--species-workers',type=int,default=6)
    args=p.parse_args()
    raw_rows=base.read_rows(args.species); species={}
    for row in raw_rows:
        name=str(row.get('canonical_name') or '').strip()
        if name: species[name]={'canonical_name':name,'family':str(row.get('family') or row.get('family_class') or '').strip()}
    if not species: raise SystemExit('No species found in input')
    workers=max(1,min(args.species_workers,8)); results=[]
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs={ex.submit(process_one,item,args):name for name,item in sorted(species.items())}
        for i,fut in enumerate(as_completed(futs),1):
            results.append(fut.result())
            if i%10==0 or i==len(futs): print({'species_completed':i,'species_total':len(futs)},flush=True)
    results.sort(key=lambda r:r['name'])
    output=[]; audits=[]; errors=[]; exclusions=Counter(); requests=0; reported=0
    for r in results:
        output.extend(r['rows']); requests+=r['requests']; reported+=r['reported']; exclusions.update(r['exclusions'])
        if r['audit'] is not None: audits.append(r['audit'])
        if r['error'] is not None: errors.append(r['error'])
    out=Path(args.out); out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('w',newline='',encoding='utf-8') as h:
        w=csv.DictWriter(h,fieldnames=base.OUTPUT_FIELDS); w.writeheader(); w.writerows(output)
    fields=['canonical_name','family','taxon_key','match_type','rank','scientific_name','accepted_name','status','gbif_reported_count','records_retrieved_before_filters','records_retained','basis_counts']
    with Path(args.taxon_audit_out).open('w',newline='',encoding='utf-8') as h:
        w=csv.DictWriter(h,fieldnames=fields); w.writeheader(); w.writerows(audits)
    per_species=Counter(str(r['canonical_name']) for r in output)
    qc={'species_requested':len(species),'species_resolved':len(audits),'species_with_retained_records':len(per_species),
        'species_ge20_retained_records':sum(v>=20 for v in per_species.values()),'retained_coordinate_rows':len(output),
        'gbif_requests_including_name_resolution':requests,'sum_gbif_reported_occurrence_counts':reported,
        'max_records_retrieved_per_species':args.max_records,'page_size':min(args.page_size,300),'hasGeospatialIssue_query_value':False,
        'accepted_basis_of_record':sorted(base.ACCEPTED_BASIS),'maximum_coordinate_uncertainty_m':args.max_coordinate_uncertainty_m,
        'coordinate_deduplication':f'one record per rounded latitude/longitude cell at {args.dedup_decimals} decimals',
        'species_workers':workers,'exclusions':dict(sorted(exclusions.items())),'failed_species':errors,
        'method_guard':'Same paginated GBIF occurrence-search QC as sequential sampler; only species-level outer loop is parallelized.'}
    Path(args.qc_out).write_text(json.dumps(qc,indent=2)+'\n',encoding='utf-8'); print(json.dumps(qc,indent=2))

if __name__=='__main__': main()

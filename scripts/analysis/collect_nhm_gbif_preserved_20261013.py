#!/usr/bin/env python3
"""Bounded provenance-preserving GBIF PRESERVED_SPECIMEN site collector.

The manifest is frozen independently of GBIF occurrence geography. This script
records specimen coordinates but NEVER reads FCP photo occurrence coordinates or
its flower colour. It does not estimate climate/phenotype effects.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,time,urllib.parse,urllib.request,urllib.error
from collections import Counter
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
from typing import Any

MANIFEST_SHA='f859bd83deb8eda50cd6527fb399cecf821a38718775fd1bdd5883f7d2be6758'
SEED='fcp_20261013_gbif_specimen_fixed_pages_v1'
GBIF='https://api.gbif.org/v1'
INAT_DATASET='50c9509d-22c7-4a22-a47d-8c48425ef4a7'
FIELDS=['canonical_name','genus','botanical_colour','gbif_id','gbif_taxon_key','datasetKey',
        'institutionCode','collectionCode','catalogNumber','decimalLatitude','decimalLongitude',
        'coordinateUncertaintyInMeters','year','basisOfRecord','gbif_issues']

def hashint(x:str)->int:return int(hashlib.sha256((SEED+'|'+x).encode()).hexdigest()[:15],16)

def urls(taxonKey:int,n:int,window:int=120)->list[tuple[int,int]]:
    """Two precommitted independently offset bounded occurrence windows."""
    if n<=0:return []
    if n<=window*2:return [(0,min(n,window*2))]
    kmax=n-window
    # Span both halves of all accessible search results rather than newest only.
    split=n//2
    max_left=max(0,split-window)
    max_right=max(0,n-split-window)
    left=hashint('first:'+str(taxonKey))%(max_left+1)
    right=split+hashint('second:'+str(taxonKey))%(max_right+1)
    assert 0<=left<right and left+window<=right and right+window<=n
    return [(left,window),(right,window)]

def request(endpoint:str,params:dict[str,object],*,timeout:int=25,retries:int=4)->dict:
    url=GBIF+endpoint+'?'+urllib.parse.urlencode(params)
    for k in range(retries):
        try:
            req=urllib.request.Request(url,headers={'User-Agent':'fcp-preserved-botanical-20261013/1.0 (zuizui0223/fcp)',
                    'Accept':'application/json'})
            with urllib.request.urlopen(req,timeout=timeout) as resp:dat=json.load(resp)
            if not isinstance(dat,dict):raise ValueError('GBIF unexpected response')
            time.sleep(0.10)  # respect public GBIF API, bound request concurrency elsewhere
            return dat
        except (TimeoutError,urllib.error.HTTPError,urllib.error.URLError,json.JSONDecodeError) as e:
            if k==retries-1:raise RuntimeError(f'{type(e).__name__}: {str(e)[:160]}') from e
            time.sleep(min(12,2**k))
    raise AssertionError

def resolve(name:str)->tuple[int,str]:
    r=request('/species/match',{'name':name,'kingdom':'Plantae','strict':'true'})
    rank=str(r.get('rank',''))
    canonical=str(r.get('canonicalName',''))
    if r.get('matchType')!='EXACT' or rank!='SPECIES' or canonical.lower()!=name.lower():
        raise ValueError('HOLD_TAXON_NOT_EXACT_SPECIES:'+str(r.get('matchType'))+':'+rank+':'+canonical)
    accepted=str(r.get('acceptedScientificName','') or r.get('scientificName',''))
    # No free synonym-acceptance from another binomial: target MUST be the GBIF accepted species
    if r.get('status')!='ACCEPTED':raise ValueError('HOLD_NOT_ACCEPTED_SPECIES:'+str(r.get('status')))
    key=r.get('speciesKey') or r.get('usageKey')
    if key is None:raise ValueError('HOLD_MISSING_SPECIES_KEY')
    return int(key),accepted

def valid(row:dict[str,object],name:str,key:int):
    if str(row.get('basisOfRecord'))!='PRESERVED_SPECIMEN':return None,'not_preserved_specimen'
    if str(row.get('datasetKey'))==INAT_DATASET:return None,'inat_dataset'
    if 'CitizenScience' in (row.get('datasetCategory') or []):return None,'citizen_science_category'
    if str(row.get('occurrenceStatus','')).upper() not in ('PRESENT',''):return None,'not_present'
    # Species exact by interpreted GBIF canonical species, no genus-only fallback or synonym.
    if str(row.get('species','')).strip().lower()!=name.lower():return None,'interpreted_species_mismatch'
    if row.get('taxonRank') and str(row.get('taxonRank')) not in ('SPECIES','SUBSPECIES','VARIETY','FORM'):
        return None,'wrong_record_rank'
    issues=[str(x) for x in (row.get('issues') or [])]
    if 'COORDINATE_INVALID' in issues or 'ZERO_COORDINATE' in issues:return None,'geospatial_flag'
    try:lat=float(row.get('decimalLatitude'));lon=float(row.get('decimalLongitude'))
    except (ValueError,TypeError):return None,'no_coordinate'
    if not (-90<=lat<=90 and -180<=lon<=180) or (abs(lat)<1e-8 and abs(lon)<1e-8):return None,'invalid_coordinate'
    unc=row.get('coordinateUncertaintyInMeters')
    if unc not in (None,''):
        try:uu=float(unc)
        except (ValueError,TypeError):return None,'uncertainty_invalid'
        if not math.isfinite(uu) or uu<0 or uu>20000:return None,'uncertainty_gt_20km'
    # Explicit same-site dedup is performed after this quality test.
    try:yr=int(row.get('year'))
    except (ValueError,TypeError):yr=None
    if yr is not None and (yr<1850 or yr>2026):return None,'impossible_collection_year'
    record={'canonical_name':name,'gbif_id':str(row.get('key')),'gbif_taxon_key':key,
            'datasetKey':str(row.get('datasetKey','')),'institutionCode':str(row.get('institutionCode') or ''),
            'collectionCode':str(row.get('collectionCode') or ''),'catalogNumber':str(row.get('catalogNumber') or ''),
            'decimalLatitude':lat,'decimalLongitude':lon,
            'coordinateUncertaintyInMeters':unc if unc is not None else '',
            'year':yr if yr is not None else '', 'basisOfRecord':'PRESERVED_SPECIMEN',
            'gbif_issues':'|'.join(issues)}
    return record,None

def collect_one(spec:dict[str,str])->dict:
    name=spec['canonical_name'];k=None; n=None; request_n=0;raw=[];reason=Counter();accepted=''
    try:
        k,accepted=resolve(name);request_n+=1
        p={'taxonKey':k,'basisOfRecord':'PRESERVED_SPECIMEN','hasCoordinate':'true',
           'hasGeospatialIssue':'false','occurrenceStatus':'present'}
        count=request('/occurrence/search',{**p,'limit':0});request_n+=1
        n=int(count.get('count',0))
        if n>100000:return {'name':name,'status':'HOLD_GBIF_TOO_MANY_RESULTS','count':n,'requests':request_n,'rows':[]}
        for offset,limit in urls(k,n):
            d=request('/occurrence/search',{**p,'limit':limit,'offset':offset});request_n+=1
            if not isinstance(d.get('results'),list):raise RuntimeError('Non-list GBIF occurrence results')
            raw.extend(d['results'])
    except Exception as e:
        return {'name':name,'status':'HOLD_SOURCE_OR_RESOLUTION','error':str(e)[:300],
                'count':n,'requests':request_n,'rows':[]}
    seen=set();cells=set();rows=[];source_datasets=set()
    for r in raw:
        v,bad=valid(r,name,k)
        if bad:reason[bad]+=1;continue
        if not v['gbif_id'] or v['gbif_id'] in seen:reason['duplicate_gbif_id']+=1;continue
        seen.add(v['gbif_id'])
        # coarsen repeated samples of same 5km-ish latitude/longitude grid cell
        xy=(round(float(v['decimalLatitude'])*20),round(float(v['decimalLongitude'])*20))
        if xy in cells:reason['geographic_cell_duplicate']+=1;continue
        cells.add(xy);v['genus']=spec['genus'];v['botanical_colour']=spec['botanical_colour']
        rows.append(v);source_datasets.add(v['datasetKey'])
    return {'name':name,'status':'SOURCE_COLLECTED' if rows else 'HOLD_NO_RETAINED_SPECIMEN',
      'taxonKey':k,'accepted_name':accepted,'count':n,'requests':request_n,
      'fetched_records':len(raw),'retained_unique_sites':len(rows),'datasets':len(source_datasets),
      'exclusions':dict(reason),'rows':rows}

def run(manifest:Path,outdir:Path,workers:int=4)->dict:
    if hashlib.sha256(manifest.read_bytes()).hexdigest()!=MANIFEST_SHA:
        raise RuntimeError('frozen NHM taxon identity manifest SHA drift')
    with manifest.open(newline='',encoding='utf-8') as f:spec=list(csv.DictReader(f))
    if len(spec)!=705 or len({v['canonical_name'] for v in spec})!=705:
        raise RuntimeError('frozen 705 strict botanical binomials mutated')
    outdir.mkdir(parents=True,exist_ok=True)
    out=[]
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs={ex.submit(collect_one,s):s['canonical_name'] for s in spec}
        for i,fut in enumerate(as_completed(futs),1):
            out.append(fut.result())
            if i%50==0 or i==len(spec):print(f'GBIF frozen taxa completed {i}/{len(spec)}',flush=True)
    out.sort(key=lambda z:z['name'])
    rows=[r for sp in out for r in sp['rows']]; counters=Counter(x['status'] for x in out)
    # Persist full source dataset+institution+record provenance in collector only.
    with (outdir/'specimen_sites.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=FIELDS);w.writeheader();w.writerows(rows)
    audits=[]
    for sp in out:
        x={k:v for k,v in sp.items() if k!='rows'};audits.append(x)
    (outdir/'taxon_source_audit.json').write_text(json.dumps(audits,ensure_ascii=False,indent=2)+'\n')
    qc={'schema':'fcp_external_NHM_GBIF_preserved_specimen_source_v1','manifest_sha256':MANIFEST_SHA,
      'requested_n_species':len(spec),'source_status_counts':dict(counters),
      'gbif_requests_total':sum(t.get('requests',0) for t in out),
      'species_with_at_least_five_unique_specimen_cells':sum(t.get('retained_unique_sites',0)>=5 for t in out),
      'species_with_any_specimen_sites':sum(t.get('retained_unique_sites',0)>=1 for t in out),
      'specimen_site_records':len(rows),'unique_specimen_dataset_keys':len({x['datasetKey'] for x in rows}),
      'east_90_180_specimen_rows':sum(90<=float(x['decimalLongitude'])<=180 for x in rows),
      'species_east':len({x['canonical_name'] for x in rows if 90<=float(x['decimalLongitude'])<=180}),
      'species_elsewhere':len({x['canonical_name'] for x in rows if not 90<=float(x['decimalLongitude'])<=180}),
      'basis_of_record':['PRESERVED_SPECIMEN'],'inat_photo_dataset_excluded':INAT_DATASET,
      'GBIF_basis_and_taxonomic_identity_checked':True,'max_records_per_species':240,
      'sampling':'fixed SHA offset two length120 index windows; 100000 API hard cap; indexed/temporal bias remains',
      'geographic_dedup':'1 record per species 0.05-degree cell, not uniform-km area',
      'data_not_independent_of_herbarium_collection_effort':True,
      'citable_GBIF_download_DOI_created':False}
    qc['specimen_csv_sha256']=hashlib.sha256((outdir/'specimen_sites.csv').read_bytes()).hexdigest()
    (outdir/'qc.json').write_text(json.dumps(qc,indent=2)+'\n')
    print(json.dumps(qc,indent=2),flush=True)
    return qc
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--manifest',type=Path,required=True);p.add_argument('--outdir',type=Path,required=True);p.add_argument('--workers',type=int,default=4);a=p.parse_args();run(a.manifest,a.outdir,a.workers)

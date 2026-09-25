#!/usr/bin/env python3
from __future__ import annotations
import argparse, io, json, math, time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from urllib.parse import urlencode
import numpy as np
import pandas as pd
import requests

API='https://api.globalbioticinteractions.org/interaction.csv'
QUERY_SPECS=[
    ('pollinator_to_plant','targetTaxon','pollinates'),
    ('pollinator_to_plant','targetTaxon','visitsFlowersOf'),
    ('plant_to_pollinator','sourceTaxon','pollinatedBy'),
    ('plant_to_pollinator','sourceTaxon','flowersVisitedBy'),
]
PRIMARY_GUILDS={'bee','other_hymenoptera','lepidoptera','diptera','coleoptera','aves','chiroptera'}
BEE_MARKERS=('apoidea','anthophila','apidae','megachilidae','halictidae','andrenidae','colletidae','melittidae','stenotritidae')

def textcol(df,*names):
    for n in names:
        if n in df.columns:
            return df[n].fillna('').astype(str)
    return pd.Series(['']*len(df),index=df.index,dtype=str)

def classify(path,name):
    s=(str(path)+' '+str(name)).lower()
    if 'sphingidae' in s:
        return 'lepidoptera','hawkmoth'
    if 'lepidoptera' in s:
        return 'lepidoptera','other_lepidoptera'
    if any(x in s for x in BEE_MARKERS):
        return 'bee','bee'
    if 'hymenoptera' in s:
        return 'other_hymenoptera','other_hymenoptera'
    if 'diptera' in s:
        return 'diptera','diptera'
    if 'coleoptera' in s:
        return 'coleoptera','coleoptera'
    if 'aves' in s or 'trochilidae' in s or 'nectariniidae' in s or 'meliphagidae' in s:
        return 'aves','aves'
    if 'chiroptera' in s or 'pteropodidae' in s or 'phyllostomidae' in s:
        return 'chiroptera','chiroptera'
    return 'other_unknown','other_unknown'

def request_one(species,orientation,param,itype,limit,timeout,retries):
    params={param:species,'interactionType':itype,'includeObservations':'true','limit':str(limit)}
    last=None
    for attempt in range(retries):
        try:
            r=requests.get(API,params=params,timeout=timeout,headers={'User-Agent':'fcp-white-pollinator-coverage/1.0'})
            if r.status_code!=200:
                raise RuntimeError(f'HTTP {r.status_code}')
            raw=r.text
            try:
                d=pd.read_csv(io.StringIO(raw))
            except pd.errors.EmptyDataError:
                d=pd.DataFrame()
            return {
                'species':species,'orientation':orientation,'param':param,'interaction_type_requested':itype,
                'http_status':r.status_code,'technical_valid':True,'rows':int(len(d)),'truncated':bool(len(d)>=limit),
                'url':r.url,'data':d,'error':''
            }
        except Exception as e:
            last=e
            if attempt+1<retries: time.sleep(min(2**attempt,8))
    return {
        'species':species,'orientation':orientation,'param':param,'interaction_type_requested':itype,
        'http_status':None,'technical_valid':False,'rows':0,'truncated':False,
        'url':API+'?'+urlencode(params),'data':pd.DataFrame(),'error':str(last)[:500]
    }

def standardize(res):
    d=res['data']
    if d.empty:
        return pd.DataFrame(columns=['plant_species','pollinator_name','pollinator_path','guild','subguild','interaction_type','decimalLatitude','decimalLongitude','referenceCitation','query_orientation'])
    if res['orientation']=='pollinator_to_plant':
        pname=textcol(d,'sourceTaxonName','source_taxon_name')
        ppath=textcol(d,'sourceTaxonPathNames','sourceTaxonPath','source_taxon_path')
    else:
        pname=textcol(d,'targetTaxonName','target_taxon_name')
        ppath=textcol(d,'targetTaxonPathNames','targetTaxonPath','target_taxon_path')
    interaction=textcol(d,'interactionTypeName','interaction_type')
    lat=pd.to_numeric(textcol(d,'decimalLatitude','latitude'),errors='coerce')
    lon=pd.to_numeric(textcol(d,'decimalLongitude','longitude'),errors='coerce')
    ref=textcol(d,'referenceCitation','citation','studySourceCitation')
    rows=[]
    for i in d.index:
        guild,sub=classify(ppath.loc[i],pname.loc[i])
        rows.append({
            'plant_species':res['species'],
            'pollinator_name':pname.loc[i].strip(),
            'pollinator_path':ppath.loc[i],
            'guild':guild,'subguild':sub,
            'interaction_type':interaction.loc[i] if i in interaction.index else res['interaction_type_requested'],
            'decimalLatitude':float(lat.loc[i]) if np.isfinite(lat.loc[i]) else np.nan,
            'decimalLongitude':float(lon.loc[i]) if np.isfinite(lon.loc[i]) else np.nan,
            'referenceCitation':ref.loc[i] if i in ref.index else '',
            'query_orientation':res['orientation'],
        })
    return pd.DataFrame(rows)

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--measured',required=True)
    p.add_argument('--outdir',required=True)
    p.add_argument('--limit',type=int,default=512)
    p.add_argument('--workers',type=int,default=8)
    p.add_argument('--timeout',type=int,default=60)
    p.add_argument('--retries',type=int,default=4)
    a=p.parse_args()
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)

    species=sorted(pd.read_csv(a.measured,usecols=['species']).species.dropna().astype(str).unique())
    if len(species)!=499:
        raise SystemExit(f'expected 499 third-cohort species, found {len(species)}')

    jobs=[(sp,*q) for sp in species for q in QUERY_SPECS]
    results=[]
    with ThreadPoolExecutor(max_workers=max(1,a.workers)) as ex:
        futs={ex.submit(request_one,sp,ori,param,itype,a.limit,a.timeout,a.retries):(sp,ori,param,itype) for sp,ori,param,itype in jobs}
        for n,f in enumerate(as_completed(futs),1):
            results.append(f.result())
            if n%100==0 or n==len(futs):
                print(json.dumps({'queries_done':n,'queries_total':len(futs),'technical_failures':sum(not r['technical_valid'] for r in results)}),flush=True)

    receipts=pd.DataFrame([{k:v for k,v in r.items() if k!='data'} for r in results]).sort_values(['species','orientation','interaction_type_requested'])
    receipts.to_csv(out/'query_receipts.csv',index=False)

    pieces=[standardize(r) for r in results if r['technical_valid'] and r['rows']>0]
    rec=pd.concat(pieces,ignore_index=True) if pieces else pd.DataFrame(columns=['plant_species','pollinator_name','pollinator_path','guild','subguild','interaction_type','decimalLatitude','decimalLongitude','referenceCitation','query_orientation'])
    rec=rec[rec.pollinator_name.astype(str).str.len()>0].copy()
    rec=rec.drop_duplicates(['plant_species','pollinator_name','interaction_type','decimalLatitude','decimalLongitude','referenceCitation'])
    rec.to_csv(out/'standardized_interaction_records.csv.gz',index=False,compression='gzip')

    coverage=[]
    for sp in species:
        g=rec[rec.plant_species.eq(sp)].copy()
        taxa=set(g.pollinator_name.astype(str))
        guilds=set(x for x in g.guild.astype(str) if x in PRIMARY_GUILDS)
        has_hawk=bool((g.subguild=='hawkmoth').any()) if len(g) else False
        geo=g.dropna(subset=['decimalLatitude','decimalLongitude']).copy()
        if len(geo):
            geo['cell1deg']=np.floor(geo.decimalLatitude).astype(int).astype(str)+':'+np.floor(geo.decimalLongitude).astype(int).astype(str)
            nguild=int(len(set(x for x in geo.guild.astype(str) if x in PRIMARY_GUILDS)))
            ncells=int(geo.cell1deg.nunique())
        else:
            nguild=0; ncells=0
        sr=receipts[receipts.species.eq(sp)]
        coverage.append({
            'species':sp,
            'all_four_queries_valid':bool(len(sr)==4 and sr.technical_valid.all()),
            'any_query_truncated':bool(sr.truncated.any()) if len(sr) else False,
            'n_standardized_records':int(len(g)),
            'n_distinct_pollinator_taxa':int(len(taxa)),
            'n_primary_guilds':int(len(guilds)),
            'guilds':'|'.join(sorted(guilds)),
            'has_sphingidae':has_hawk,
            'n_georeferenced_records':int(len(geo)),
            'n_geo_cells_1deg':ncells,
            'n_geo_primary_guilds':nguild,
        })
    cov=pd.DataFrame(coverage)
    cov.to_csv(out/'coverage_by_plant.csv',index=False)

    nvalid=int(cov.all_four_queries_valid.sum())
    p0=(nvalid/len(cov)>=0.95)
    covered=cov.n_distinct_pollinator_taxa.ge(1)
    p1_n1=int(covered.sum())
    p1_n3=int(cov.n_distinct_pollinator_taxa.ge(3).sum())
    p1=bool(p1_n1>=100 and p1_n3>=60)
    p2_n=int(cov.n_primary_guilds.ge(2).sum())
    p2=bool(p2_n>=40)
    hawk=int(cov.has_sphingidae.sum())
    comparison=int((covered & ~cov.has_sphingidae).sum())
    p3=bool(hawk>=20 and comparison>=20)
    spatial=(cov.n_georeferenced_records.ge(10)&cov.n_geo_cells_1deg.ge(2))
    p4_n=int(spatial.sum())
    p4_multiguild=int((spatial&cov.n_geo_primary_guilds.ge(2)).sum())
    p4=bool(p4_n>=50 and p4_multiguild>=30)

    status='coverage_complete' if p0 else 'coverage_not_evaluable_api_integrity'
    summary={
        'schema':'fcp_white_pollinator_coverage_v1',
        'status':status,
        'role':'outcome_blind_coverage_only',
        'plant_frame_species':len(cov),
        'live_api_warning':'GloBI live/API results are unstable and are not authorized for final biological inference; passing dimensions must be rebuilt from a stable versioned release.',
        'query_cap_per_species_orientation':a.limit,
        'queries_total':len(receipts),
        'queries_technically_valid':int(receipts.technical_valid.sum()),
        'queries_truncated':int(receipts.truncated.sum()),
        'P0_query_integrity':{'pass':p0,'species_all_four_valid':nvalid,'required_fraction':0.95},
        'P1_broad_pollinator_identifiability':{'pass':p1,'species_with_any_pollinator':p1_n1,'species_with_ge3_pollinator_taxa':p1_n3,'requirements':'>=100 any; >=60 with >=3 taxa'},
        'P2_multiguild_identifiability':{'pass':p2,'species_with_ge2_primary_guilds':p2_n,'requirement':'>=40'},
        'P3_hawkmoth_identifiability':{'pass':p3,'species_with_sphingidae':hawk,'covered_species_without_sphingidae':comparison,'requirements':'>=20 each; non-detection is not biological absence'},
        'P4_spatial_identifiability':{'pass':p4,'species_ge10_georeferenced_and_ge2_cells':p4_n,'of_these_species_with_ge2_guilds':p4_multiguild,'requirements':'>=50 spatial; >=30 spatial multi-guild'},
        'guild_species_counts':{g:int(cov.guilds.fillna('').str.split('|').apply(lambda xs:g in xs).sum()) for g in sorted(PRIMARY_GUILDS)},
        'hard_nonclaims':['coverage does not establish pollinator preference','coverage does not establish pollinator-driven evolution','coverage does not identify pigmented-to-white transition direction','missing GloBI records are not pollinator absence']
    }
    (out/'coverage_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
    if not p0:
        raise SystemExit('P0 API integrity gate failed')
if __name__=='__main__':
    main()

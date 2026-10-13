#!/usr/bin/env python3
"""Select NHM botanical colour taxa BEFORE inspecting GBIF specimen locations.

Only independently described botanical labels and a deterministic taxonomy hash are
used. Absolutely no FCP iNaturalist name/position/outcome filters are applied.
This is retrospective relative to the east-white-blue photo hypothesis.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,sys
from pathlib import Path
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parent))
# Parent PR #191 is tracked without change and provides independently guarded
# strict botanical flower-colour definitions (no FCP photograph outcomes).
from audit_independent_nhm_colour_rain_20261013 import (
    blob_sha,extract_botanical,NHM_BLOB_SHA,NHM_BYTES
)
SALT='FCP_EXTERNAL_GBIF_BOTANICAL_PRESERVED_20261013_v1'
MAX_GENERA=96
PER_COLOUR=4
MIN_PER_HUE_GENUS=3

def token(s:str)->str: return hashlib.sha256((SALT+'|'+str(s)).encode()).hexdigest()

def choose(labels: pd.DataFrame,max_genera:int=MAX_GENERA)->tuple[pd.DataFrame,dict]:
    d=labels[labels.label.isin(['white','blue_purple'])].copy()
    d['genus']=d.taxon.str.split().str[0]
    n=d.groupby(['genus','label']).taxon.nunique().unstack(fill_value=0)
    allowed=n.index[(n.get('white',0)>=MIN_PER_HUE_GENUS)&(n.get('blue_purple',0)>=MIN_PER_HUE_GENUS)].tolist()
    selected=sorted(allowed,key=lambda g:token('genus:'+g))[:max_genera]
    rows=[]
    for g in selected:
        for color in ['white','blue_purple']:
            names=d.loc[(d.genus==g)&(d.label==color),'taxon'].drop_duplicates().tolist()
            for name in sorted(names,key=lambda v:token('species:'+v))[:PER_COLOUR]:
                rows.append({'canonical_name':name,'genus':g,'botanical_colour':color,
                  'selection_order_genus':1+selected.index(g),'selection_hash':token('species:'+name)})
    out=pd.DataFrame(rows,columns=['canonical_name','genus','botanical_colour','selection_order_genus','selection_hash'])
    assert out.canonical_name.is_unique
    assert out.groupby('genus').botanical_colour.nunique().eq(2).all() if len(out) else True
    return out,{'selection_salt':SALT,'candidate_genera_with_at_least_three_each_hue':len(allowed),
        'selected_genera':len(selected),'selected_species':len(out),
        'selected_botanical_colour_counts':out.botanical_colour.value_counts().to_dict(),
        'max_genera':max_genera,'per_colour_cap':PER_COLOUR,
        'minimum_each_hue_per_genus':MIN_PER_HUE_GENUS,
        'uses_fcp_occurrence_or_photo_colour':False,'selection_uses_gbif_occurrences':False,
        'prior_gbif_api_feasibility_test':'Campanula punctata only; no candidate identity, environment or hue used to select manifest'}

def run(nhm:Path,outdir:Path,max_genera:int=MAX_GENERA)->dict:
    b=nhm.stat().st_size
    if b!=NHM_BYTES or blob_sha(nhm)!=NHM_BLOB_SHA:raise RuntimeError('NHM exact pinned source file mismatch')
    raw=pd.read_csv(nhm,low_memory=False)
    labels,counts=extract_botanical(raw)
    out,summary=choose(labels,max_genera)
    outdir.mkdir(parents=True,exist_ok=True)
    out.to_csv(outdir/'frozen_selected_taxa.csv',index=False)
    h=hashlib.sha256((outdir/'frozen_selected_taxa.csv').read_bytes()).hexdigest()
    result={'schema':'fcp_independent_preserved_specimen_preoutcome_taxa_v1',
            'nhm_source_size_bytes':b,'nhm_source_git_blob_sha1':NHM_BLOB_SHA,
            'nhm_trait_coverage':counts,'manifest_sha256':h,**summary,
            'source_independence':'NHM botanical description; no FCP photo identities, observed sites or palette labels',
            'inference_ceiling':'Not prospective relative to the earlier outcome-exposed east white-blue hypothesis'}
    (outdir/'preflight.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    return result
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--nhm-source',type=Path,required=True);p.add_argument('--outdir',type=Path,required=True);p.add_argument('--max-genera',type=int,default=MAX_GENERA);x=p.parse_args();print(json.dumps(run(x.nhm_source,x.outdir,x.max_genera),indent=2))

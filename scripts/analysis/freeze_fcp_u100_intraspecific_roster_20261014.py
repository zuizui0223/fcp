#!/usr/bin/env python3
"""Reconstruct frozen U100 from history WITHOUT opening phenotype/climate values.

Do NOT infer capacity from the 42,111 metadata source: that historical CSV has
only taxon ID, species, and discovery provenance. The U100 parents are two
independent identity surfaces: P100 3,730 (observer-cap available) and the
original 1,000 discovery/validation species, whose capacity was attested in
the 2026-09-13 original P100 artifact but not stored per legacy taxon.
"""
from __future__ import annotations

import argparse, hashlib, json
from pathlib import Path
import pandas as pd

VERSION="fcp_u100_within_species_outcome_blind_roster_v2"
SOURCE={
 "frame":"5d871bd1f1190c68fd513bb527e6c93de1606dc35ff8cfff5f6187850dc382cc",
 "P100":"1473aad680fe2fa84903c5e11ee104fd0eceef828957f16c2ef1a35f9dd6993c",
 "P500":"f54d07fb2338a20a3f92808b58f0ebc948ab0d5af3cb01fb26702f861184b7a4",
 "third_selected":"16db6a2fc265f0ab20e7dae4e6f9058b907f5c19af3ddab5b6077797dd1def59",
 "legacy_discovery":"ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4",
 "legacy_validation":"0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6",
}
EXPECTED={
 "legacy_discovery_or_validation":1000,
 "earlier_P500_selected_not_new_holdout":500,
 "third_cohort_selected_500_of_499_measured":500,
 "still_never_selected_for_high_depth":2730,
}
ID=["inat_taxon_id","species"]
POOL=ID+["after_observer_cap"]


def sha(path:Path)->str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def idframe(df:pd.DataFrame,cap=False)->pd.DataFrame:
    cols=POOL if cap else ID
    missing=set(cols)-set(df.columns)
    if missing:raise ValueError('missing source identity/capacity columns '+str(sorted(missing)))
    d=df[cols].copy()
    d.inat_taxon_id=pd.to_numeric(d.inat_taxon_id,errors='raise').astype('int64')
    d.species=d.species.astype('string').str.strip()
    if d.inat_taxon_id.duplicated().any() or d.species.duplicated().any():raise ValueError('source taxonomy is not one-to-one')
    if d.species.isna().any() or d.species.str.len().eq(0).any():raise ValueError('blank species')
    if cap:d.after_observer_cap=pd.to_numeric(d.after_observer_cap,errors='raise').astype('int64')
    return d.sort_values('inat_taxon_id',kind='mergesort').reset_index(drop=True)


def check_subset(d:pd.DataFrame, parent:pd.DataFrame, name:str):
    mapping=dict(zip(parent.inat_taxon_id,parent.species))
    for i,s in zip(d.inat_taxon_id,d.species):
        if mapping.get(int(i))!=s:raise ValueError(f'{name} missing/mismatched identity {i} {s}')


def roster_from_identities(frame,p100,p500,third,discovery,validation,expected:dict|None=None):
    f=idframe(frame);p=idframe(p100,cap=True);old=idframe(p500,cap=True)
    third=idframe(third);a=idframe(discovery);b=idframe(validation)
    for name,obj,parent in [('P100',p,f),('P500',old,p),('third',third,p),('discovery',a,f),('validation',b,f)]:check_subset(obj,parent,name)
    legacy=pd.concat([a,b],ignore_index=True)
    if legacy.inat_taxon_id.duplicated().any() or legacy.species.duplicated().any():raise ValueError('legacy discovery/validation overlap')
    idsP=set(p.inat_taxon_id);idsLegacy=set(legacy.inat_taxon_id)
    idsP500=set(old.inat_taxon_id);idsThird=set(third.inat_taxon_id)
    if idsLegacy & idsP:raise ValueError('legacy/P100 overlap')
    if idsThird & idsP500:raise ValueError('old P500 / third overlap')
    if p.after_observer_cap.lt(100).any():raise ValueError('P100 below frozen observer-cap gate')
    if not old.set_index('inat_taxon_id').after_observer_cap.eq(p.set_index('inat_taxon_id').after_observer_cap.loc[old.inat_taxon_id].to_numpy()).all():
        raise ValueError('P500 capacity diverged from P100')
    p=p.copy();p['capacity_provenance']='frozen_P100_observer_cap'
    l=legacy.copy();l['after_observer_cap']=pd.NA
    l['capacity_provenance']='legacy_U100_attested_parent_selection_20260913_per_species_value_not_archived'
    d=pd.concat([p,l],ignore_index=True)
    d['historical_membership']=[('legacy_discovery_or_validation' if tid in idsLegacy else
        'earlier_P500_selected_not_new_holdout' if tid in idsP500 else
        'third_cohort_selected_500_of_499_measured' if tid in idsThird else
        'still_never_selected_for_high_depth') for tid in d.inat_taxon_id]
    d['fresh_outcome_unexposed_candidate']=d.historical_membership.eq('still_never_selected_for_high_depth')
    d['new_colour_measurement_status']='NOT_PERFORMED_BY_THIS_FREEZE'
    d=d.sort_values('inat_taxon_id',kind='mergesort').reset_index(drop=True)
    found={k:int(d.historical_membership.eq(k).sum()) for k in EXPECTED}
    if expected is not None and found!=expected:raise ValueError('unexpected cohort attrition '+str(found))
    if not d.inat_taxon_id.is_unique or len(d)!=len(p)+len(legacy):raise ValueError('U100 identity double counted')
    return d,found


def run(inputs:dict[str,Path],out:Path):
    for key,h in SOURCE.items():
        if sha(inputs[key])!=h:raise RuntimeError(f'{key} historical SHA256 differs from frozen source: {sha(inputs[key])}')
    df=pd.read_csv(inputs['frame'],usecols=ID)
    p=pd.read_csv(inputs['P100'],usecols=POOL)
    q=pd.read_csv(inputs['P500'],usecols=POOL)
    t=pd.read_csv(inputs['third_selected'],sep='\t',usecols=ID)
    a=pd.read_csv(inputs['legacy_discovery'],usecols=ID).drop_duplicates()
    b=pd.read_csv(inputs['legacy_validation'],usecols=ID).drop_duplicates()
    if tuple(map(len,[df,p,q,t,a,b]))!=(42111,3730,500,500,500,500):raise RuntimeError('legacy high-depth/photo opportunity source denominators drift')
    roster,counts=roster_from_identities(df,p,q,t,a,b,EXPECTED)
    if len(roster)!=4730:raise RuntimeError('not historical U100 4730')
    out.mkdir(parents=True,exist_ok=True)
    full=out/'all_u100_4730_species_status.csv';un=out/'never_selected_2730_preoutcome_species.csv'
    roster.to_csv(full,index=False,lineterminator='\n')
    roster.loc[roster.fresh_outcome_unexposed_candidate,POOL].to_csv(un,index=False,lineterminator='\n')
    receipt={
      'schema':VERSION,'status':'U100_IDENTITIES_FROZEN_NO_NEW_COLOUR_OR_ENVIRONMENT_OPENED',
      'historical_opportunity_species':len(df),'high_depth_opportunity_species':len(roster),
      'historical_memberships':counts,'never_selected_high_depth':2730,
      'legacy_capacity_value_not_in_42111_csv':True,
      'historic_legacy_U100_membership_attested_in_parent':True,
      'photo_colours_retrieved_by_this_step':0,'new_W_tests_performed':0,'new_within_species_environment_tests_performed':0,
      'original_q_white_changed':False,'source_columns_used':['inat_taxon_id','species','after_observer_cap (only P100/P500)'],
      'source_SHA256':SOURCE,
      'output_SHA256':{full.name:sha(full),un.name:sha(un)},
      'claim_ceiling':'original U100 opportunity, not actual 4730 W eligible measurements or genetic polymorphism'
    }
    (out/'outcome_blind_roster_receipt.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
    return receipt


def main():
    parser=argparse.ArgumentParser()
    for key in SOURCE:parser.add_argument('--'+key.lower().replace('_','-'),dest=key,type=Path,required=True)
    parser.add_argument('--outdir',type=Path,required=True)
    a=parser.parse_args()
    inputs={k:getattr(a,k) for k in SOURCE}
    res=run(inputs,a.outdir)
    print(json.dumps({k:res[k] for k in ('status','historical_memberships','output_SHA256')},ensure_ascii=False))

if __name__=='__main__':main()

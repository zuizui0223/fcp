#!/usr/bin/env python3
"""Post-outcome *descriptive* leave-one-site-year-out audit. NOT prespecified evidence.

No source observations modified. Per-stratum hypergeometric draws recover the same
original 9,999-draw source-locked null stream, so zero omitted sites are imputed.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
from field_partition import (SEEDS, R, load_sisyrinchium,load_boechera,
                             select_strata,hold_reason,within_site,permute_conditional)

def audit_one(d, species, R=R):
    rng=np.random.default_rng(SEEDS[species])
    site_records=[]
    for z in select_strata(d,species):
        n,v=np.asarray(z['n'],dtype=int),np.asarray(z['v'],dtype=int)
        n_occupied=int(np.count_nonzero(n>0));N=int(n.sum());V=int(v.sum())
        if hold_reason(N,n_occupied,V,N-V) is not None:
            continue
        observed=within_site(n,v)
        null=permute_conditional(n,v,R,rng=rng)['B']-observed['expected_B']
        site_records.append({'siteyear':z['siteyear'],'site':z['site'],'observed_delta':observed['delta'],'null_delta':null,'N':N,'V':V})
    if not site_records:
        return [],None
    matrix=np.stack([z['null_delta'] for z in site_records])
    observed=np.array([z['observed_delta'] for z in site_records])
    full=matrix.mean(axis=0)
    original_mean=float(observed.mean())
    out=[]
    for j,rec in enumerate(site_records):
        if len(site_records)==1: continue
        minus_obs=float((observed.sum()-observed[j])/(len(site_records)-1))
        minus_null=(matrix.sum(axis=0)-matrix[j])/(len(site_records)-1)
        p=(1+np.count_nonzero(minus_null>=minus_obs-1e-12))/(R+1)
        out.append({'species':species,'excluded_siteyear':rec['siteyear'],'excluded_source_N':rec['N'],
                    'excluded_delta_percentage_points':100*rec['observed_delta'],
                    'full_mean_delta_percentage_points':100*original_mean,
                    'after_leave_one_out_delta_percentage_points':100*minus_obs,
                    'after_leave_one_out_p_one_sided':p,
                    'null_mean_after_exclusion':float(minus_null.mean()),
                    'fragility_scope':'POSTHOC_DESCRIPTION_NOT_REGISTERED_PRIMARY'})
    full_p=(1+np.count_nonzero(full>=original_mean-1e-12))/(R+1)
    return out,{'species':species,'site_equal_original':original_mean,'original_p_recomputed':full_p,
                'per_site_null_mean':float(full.mean()),'sites_used':len(site_records)}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--source-dir',type=Path,default=Path(__file__).resolve().parents[1]/'source');parser.add_argument('--outdir',type=Path,default=Path(__file__).resolve().parents[1]/'results');parser.add_argument('--reps',type=int,default=R)
    args=parser.parse_args()
    sis=load_sisyrinchium(args.source_dir/'Sisyrinchium_data.zip')
    _,boe,_=load_boechera(args.source_dir/'Boechera_natural_populations_individual.txt',args.source_dir/'Boechera_natural_populations_plot.txt')
    allrows=[]; checks={}
    for name,d in [('sisyrinchium',sis),('boechera',boe)]:
        rows,check=audit_one(d,name,R=args.reps)
        allrows.extend(rows);checks[name]=check
    args.outdir.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(allrows).to_csv(args.outdir/'posthoc_leave_one_siteyear_out.csv',index=False)
    (args.outdir/'posthoc_fragility_receipt.json').write_text(json.dumps({'status':'COMPUTED_POSTHOC_NOT_CONFIRMATORY','checks':checks,'n_left_out_rows':len(allrows)},indent=2)+'\n')
    print(json.dumps({'checks':checks,'boechera_loo': [r for r in allrows if r['species']=='boechera']},indent=2))
if __name__=='__main__':main()

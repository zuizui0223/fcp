#!/usr/bin/env python3
"""Exploratory species-conditional equal-area colour-turnover permutation, source locked.

Directly tests observed global PHOTO colour turnover vs fixed cell/photo counts.
No genetic IBD/IBE, local adaptation, environmental sorting or past distribution
is inferable from this statistical operation.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd
from fcp_global_photo_atlas_20261020 import (
    read_original, geographic_units, high_purity, group_stats, check, STATES
)

REPS=199
N_STATES=4


def eligible_units(d):
    opportunity=group_stats(d,['inat_taxon_id','cell'],10,3)
    keep=opportunity[['inat_taxon_id','cell']]
    out=d.merge(keep,on=['inat_taxon_id','cell'],how='inner',validate='many_to_one')
    return out,opportunity


def score(frequencies):
    a,b=np.triu_indices(len(frequencies),k=1)
    if len(a)==0:return np.nan,np.nan
    x=frequencies[a]-frequencies[b]
    return float(np.abs(x).sum(axis=1).mean()/2),float(np.abs(x[:,0]).mean())


def taxon_test(g,tag):
    g=g.sort_values(['cell','photo_id'],kind='mergesort')
    sites,n=np.unique(g.cell.to_numpy(),return_counts=True)
    if len(sites)<2:return None
    labels=pd.Categorical(g.morph, categories=STATES).codes
    check(np.min(labels)>=0,'bad morph code')
    starts=np.insert(np.cumsum(n)[:-1],0,0)
    denominators=n[:,None]
    obs_matrix=np.add.reduceat(np.eye(N_STATES,dtype=int)[labels],starts,axis=0)
    actual=score(obs_matrix/denominators)
    tid=int(g.inat_taxon_id.iloc[0]); seed=int.from_bytes(
       hashlib.sha256(f'20261020_{tag}_{tid}'.encode()).digest()[:8], 'little')
    rng=np.random.default_rng(seed)
    null=np.empty((REPS,2))
    for k in range(REPS):
        x=np.add.reduceat(np.eye(N_STATES,dtype=int)[rng.permutation(labels)],starts,axis=0)
        null[k]=score(x/denominators)
    return {'taxon_id':tid, 'n_photos_in_eligible_cells':int(len(g)),
           'n_eligible_cells':int(len(sites)), 'n_cell_pairs':int(len(sites)*(len(sites)-1)//2),
           'observed_tv':actual[0], 'observed_abs_white_delta':actual[1],
           'null_expected_tv':float(null[:,0].mean()),
           'null_expected_white_delta':float(null[:,1].mean()),
           'excess_tv':actual[0]-null[:,0].mean(),
           'excess_white_delta':actual[1]-null[:,1].mean()}, null


def calculate(d,tag):
    rows=[];nulls=[]
    for _,g in d.groupby('inat_taxon_id',sort=True):
        ret=taxon_test(g,tag)
        if ret is None:continue
        row,vals=ret;rows.append(row);nulls.append(vals)
    frame=pd.DataFrame(rows)
    check(len(frame)>0,'no eligible species in geographic turnover')
    sims=np.mean(np.stack(nulls,axis=0),axis=0)
    observed=frame[['observed_tv','observed_abs_white_delta']].mean().to_numpy()
    nullmean=sims.mean(axis=0)
    p=[float((1+np.count_nonzero(sims[:,col]>=observed[col]))/(REPS+1))
       for col in range(2)]
    return frame,{
        'n_species':int(len(frame)),
        'n_eligible_cell_pairs_not_independent':int(frame.n_cell_pairs.sum()),
        'mean_observed_between_cell_tv':float(observed[0]),
        'mean_null_tv':float(nullmean[0]),
        'mean_excess_tv':float(observed[0]-nullmean[0]),
        'mean_observed_white_delta':float(observed[1]),
        'mean_null_white_delta':float(nullmean[1]),
        'mean_excess_white_delta':float(observed[1]-nullmean[1]),
        'vertex_permutation_upper_p_tv_exploratory':p[0],
        'vertex_permutation_upper_p_white_exploratory':p[1],
        'repetitions':REPS,
        'species_equal_not_photo_pair_equal':True,
        'fixed_original_taxon_cell_photo_opportunities':True,
        'posthoc_not_prospective':True,
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--photos',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    a=parser.parse_args()
    data=read_original(a.photos)
    geo=geographic_units(data)
    allunits,opp=eligible_units(geo)
    original,summary_full=calculate(allunits,'all')
    pure=high_purity(data)
    qgeo=geographic_units(pure)
    pureunits,qopportunity=eligible_units(qgeo)
    qspecies=qopportunity.groupby('inat_taxon_id').cell.nunique()
    qids=set(qspecies[qspecies>=2].index)
    # Keep same exact cell membership in the two PHOTO quality conditions.
    qkeep=qopportunity[qopportunity.inat_taxon_id.isin(qids)][['inat_taxon_id','cell']]
    pure_paired=pureunits.merge(qkeep,on=['inat_taxon_id','cell'],validate='many_to_one')
    full_paired=allunits.merge(qkeep,on=['inat_taxon_id','cell'],validate='many_to_one')
    matched_orig,summary_matched_all=calculate(full_paired,'matched_all')
    matched_pure,summary_matched_pure=calculate(pure_paired,'matched_pure')
    check(set(matched_orig.taxon_id)==set(matched_pure.taxon_id),'purity taxa do not match')
    result={
      'schema':'fcp_actual_447500_species_grid_vertex_null_v1',
      'n_total_original_classifiable':int(len(data)),
      'n_total_pure_colour_roi':int(len(pure)),
      'all_classifiable_species_grid':summary_full,
      'matched_same_species_same_cells_all_classifiable':summary_matched_all,
      'matched_same_species_same_cells_high_purity':summary_matched_pure,
      'quality_rule':'assigned ROI >=65%, all competing coarse groups <=20%; not morph validation',
      'inference_ceiling':'geographic sorting of photographed colour labels only; not genetic population structure or dispersal history',
    }
    a.out.mkdir(parents=True,exist_ok=True)
    (a.out/'equal_area_spatial_photo_label_nulls.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('FCP_REAL_GRID_SPATIAL_NULL_RESULT',json.dumps(result,sort_keys=True,allow_nan=False))

if __name__=='__main__':main()

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


def bootstrap_component_share(frame:pd.DataFrame, reps:int=999):
    """Resample species, not photo-pairs. A ratio CI may be unbounded near zero."""
    a=frame[['excess_tv','excess_white_delta']].to_numpy(float)
    n=len(a)
    if n<20:return [None,None]
    rng=np.random.default_rng(2026102039)
    draws=rng.integers(0,n,size=(reps,n))
    means=a[draws].mean(axis=1)
    ratios=np.divide(means[:,1],means[:,0],
        out=np.full(reps,np.nan),where=np.abs(means[:,0])>1e-12)
    finite=ratios[np.isfinite(ratios)]
    return [float(x) for x in np.quantile(finite,[.025,.975])] if len(finite)>reps*.95 else [None,None]


def selftest():
    arr=np.array([[.5,.3,.2,0],[.2,.3,.3,.2]],float)
    tv,white,res=score(arr)
    check(abs(tv-.3)<1e-12 and abs(white-.3)<1e-12 and abs(res)<1e-12,
          "pure white/nonwhite compositional switch")
    arr=np.array([[.5,.5,0,0],[.5,0,.5,0]],float)
    tv,white,res=score(arr)
    check(abs(tv-.5)<1e-12 and abs(white)<1e-12 and abs(res-.5)<1e-12,
          "pure nonwhite hue switch")
    print("FCP_SPATIAL_COLOUR_DECOMPOSITION_TESTS_PASS")


def eligible_units(d):
    opportunity=group_stats(d,['inat_taxon_id','cell'],10,3)
    keep=opportunity[['inat_taxon_id','cell']]
    out=d.merge(keep,on=['inat_taxon_id','cell'],how='inner',validate='many_to_one')
    return out,opportunity


def score(frequencies):
    a,b=np.triu_indices(len(frequencies),k=1)
    if len(a)==0:return np.nan,np.nan
    x=frequencies[a]-frequencies[b]
    tv=float(np.abs(x).sum(axis=1).mean()/2)
    white=float(np.abs(x[:,0]).mean())
    residual=tv-white
    check(residual >= -1e-12,"invalid TV decomposition")
    return tv, white, max(0.0,residual)


def taxon_test(g,tag):
    g=g.sort_values(['cell','photo_id'],kind='mergesort')
    sites,n=np.unique(g.cell.to_numpy(),return_counts=True)
    if len(sites)<2:return None
    check(g.fcp_source_cohort.nunique()==1,'source-cohort taxon identity inconsistent')
    cohort=str(g.fcp_source_cohort.iloc[0])
    labels=pd.Categorical(g.morph, categories=STATES).codes
    n_white=int(np.count_nonzero(labels==0))
    n_nonwhite=int(len(labels)-n_white)
    n_nonwhite_hues=int(np.unique(labels[labels!=0]).size)
    check(np.min(labels)>=0,'bad morph code')
    starts=np.insert(np.cumsum(n)[:-1],0,0)
    denominators=n[:,None]
    obs_matrix=np.add.reduceat(np.eye(N_STATES,dtype=int)[labels],starts,axis=0)
    actual=score(obs_matrix/denominators)
    tid=int(g.inat_taxon_id.iloc[0]); seed=int.from_bytes(
       hashlib.sha256(f'20261020_{tag}_{tid}'.encode()).digest()[:8], 'little')
    rng=np.random.default_rng(seed)
    null=np.empty((REPS,3))
    for k in range(REPS):
        x=np.add.reduceat(np.eye(N_STATES,dtype=int)[rng.permutation(labels)],starts,axis=0)
        null[k]=score(x/denominators)
    return {'taxon_id':tid, 'cohort':cohort,
           'n_white_photo':n_white, 'n_nonwhite_photo':n_nonwhite,
           'n_nonwhite_hues':n_nonwhite_hues,
           'n_photos_in_eligible_cells':int(len(g)),
           'n_eligible_cells':int(len(sites)), 'n_cell_pairs':int(len(sites)*(len(sites)-1)//2),
           'observed_tv':actual[0], 'observed_abs_white_delta':actual[1],
           'observed_nonwhite_residual_tv':actual[2],
           'null_expected_tv':float(null[:,0].mean()),
           'null_expected_white_delta':float(null[:,1].mean()),
           'null_expected_nonwhite_residual':float(null[:,2].mean()),
           'excess_tv':actual[0]-null[:,0].mean(),
           'excess_white_delta':actual[1]-null[:,1].mean(),
           'excess_nonwhite_residual':actual[2]-null[:,2].mean()}, null


def mean_trimmed(a, fraction=.10, axis=0):
    """Trim exactly floor(f*n) from both tails, applied identically to observed/null."""
    z=np.sort(np.asarray(a,float),axis=axis)
    n=z.shape[axis]; k=int(n*fraction)
    t=np.take(z,np.arange(k,n-k),axis=axis)
    return np.mean(t,axis=axis)


def robust_white_summary(frame:pd.DataFrame, stacked_null:np.ndarray):
    """Recalibrate every robust statistic using the same original 199 vertex-null worlds.

    Null worlds are centered at their within-taxon permutation mean; both
    observed and null statistics receive the same transformations. Positive
    selection gates depend only on each taxon's source composition.
    All tests are explicitly exploratory and not independent confirmation.
    """
    cases={
        'all_spatially_evaluable_taxa':np.ones(len(frame),bool),
        'photo_white_and_nonwhite_at_least_5_each':(
            (frame.n_white_photo>=5)&(frame.n_nonwhite_photo>=5)).to_numpy(),
        'photo_two_nonwhite_hues_present':(frame.n_nonwhite_hues>=2).to_numpy(),
    }
    for c in ('discovery','validation','prospective','wave4','new_U100'):
        cases['leave_out_'+c]=frame.cohort.ne(c).to_numpy()
        cases['white5_'+c]=(
          ((frame.n_white_photo>=5)&(frame.n_nonwhite_photo>=5))
          & frame.cohort.eq(c)).to_numpy()
    cases['white5_leave_out_new_U100']=(
       (frame.n_white_photo>=5)&(frame.n_nonwhite_photo>=5)
       & frame.cohort.ne('new_U100')).to_numpy()
    output={}
    for name,mask in cases.items():
        part=frame.loc[mask]
        n=len(part)
        if n<20:
            output[name]={'n_species':int(n),'status':'insufficient_species'}
            continue
        obswhite=part.excess_white_delta.to_numpy(float)
        obsres=part.excess_nonwhite_residual.to_numpy(float)
        sims=stacked_null[mask]
        simswhite=sims[:,:,1]-sims[:,:,1].mean(axis=1,keepdims=True)
        simsres=sims[:,:,2]-sims[:,:,2].mean(axis=1,keepdims=True)
        white_minus_other=obswhite-obsres
        null_white_minus_other=simswhite-simsres
        def p_upper(a,worlds):
            return float((1+np.count_nonzero(worlds>=a))/(REPS+1))
        observed={
          'mean':float(obswhite.mean()),
          'median':float(np.median(obswhite)),
          'trimmed10':float(mean_trimmed(obswhite)),
          'positive_fraction':float(np.mean(obswhite>1e-12)),
          'white_minus_nonwhite_mean':float(white_minus_other.mean()),
        }
        nulls={
          'mean':simswhite.mean(axis=0),
          'median':np.median(simswhite,axis=0),
          'trimmed10':mean_trimmed(simswhite,axis=0),
          'positive_fraction':np.mean(simswhite>1e-12,axis=0),
          'white_minus_nonwhite_mean':null_white_minus_other.mean(axis=0),
        }
        rng=np.random.default_rng(2026102037+n)
        boot_idx=rng.integers(0,n,size=(999,n))
        boot=white_minus_other[boot_idx].mean(axis=1)
        output[name]={
          'n_species':int(n),
          'n_taxa_with_at_least_one_white_label':int((part.n_white_photo>0).sum()),
          'n_taxa_with_white_and_nonwhite_labels':int(((part.n_white_photo>0)&(part.n_nonwhite_photo>0)).sum()),
          'observed':observed,
          'null_matched_upper_p_exploratory':{key:p_upper(observed[key],nulls[key]) for key in observed},
          'white_minus_nonwhite_species_bootstrap_CI95':[float(z) for z in np.quantile(boot,[.025,.975])],
          'nonwhite_residual_mean_excess':float(obsres.mean()),
          'n_null_worlds':REPS,
          'posthoc_multiple_tests_not_confirmatory':True,
        }
    return output


def summarize_by_source_cohort(frame:pd.DataFrame, stacked_null:np.ndarray):
    """Same original species-equal nulls, now shown within disjoint source cohorts."""
    out={}
    for cohort in ('discovery','validation','prospective','wave4','new_U100'):
        chosen=frame.cohort.eq(cohort).to_numpy()
        if not chosen.any():
            out[cohort]={'n_species':0,'status':'underidentified'}
            continue
        part=frame.loc[chosen]
        means=stacked_null[chosen].mean(axis=0)
        obs=part[['observed_tv','observed_abs_white_delta','observed_nonwhite_residual_tv']].mean().to_numpy()
        nullmeans=means.mean(axis=0)
        excess=obs-nullmeans
        out[cohort]={
            'n_species':int(len(part)),
            'mean_excess_total_tv':float(excess[0]),
            'mean_excess_white_delta':float(excess[1]),
            'mean_excess_nonwhite_residual':float(excess[2]),
            'white_share_of_mean_excess_turnover':float(excess[1]/excess[0]) if abs(excess[0])>1e-12 else None,
            'species_bootstrap_95CI_white_share':bootstrap_component_share(part),
            'p_upper_white_axis_199':float((1+np.count_nonzero(means[:,1]>=obs[1]))/(REPS+1)),
            'p_upper_nonwhite_residual_199':float((1+np.count_nonzero(means[:,2]>=obs[2]))/(REPS+1)),
            'fraction_of_species_positive_white_excess':float((part.excess_white_delta>0).mean()),
            'exploratory_posthoc':True,
        }
    return out


def relative_component_excess_test(frame:pd.DataFrame, null_stack:np.ndarray, min_baseline:float):
    """Compare proportional geographic excess after equalizing each colour component's own opportunity.

    For species i:
       Rwhite_i = observed(abs white-frequency difference)_i / E_null(white)_i - 1
       Rother_i = observed(nonwhite residual TV)_i / E_null(other)_i - 1
    Only species with BOTH null expected components >= predeclared floor are eligible.
    Permutation distribution uses the identical expected denominators in each null
    world. This is a POST HOC relative-spatial-sorting diagnostic, not selection.
    """
    wn=frame.null_expected_white_delta.to_numpy(float)
    hn=frame.null_expected_nonwhite_residual.to_numpy(float)
    eligible=(wn>=min_baseline)&(hn>=min_baseline)&np.isfinite(wn)&np.isfinite(hn)
    n=int(eligible.sum())
    result={'min_component_null_baseline':float(min_baseline),'n_eligible_species':n,
            'posthoc_component_scale_sensitivity':True}
    if n<20:
        result['status']='not_estimable_under_20_species'
        return result
    x=frame.loc[eligible].copy()
    w=wn[eligible]
    h=hn[eligible]
    obs_w=x.observed_abs_white_delta.to_numpy(float)/w-1.0
    obs_h=x.observed_nonwhite_residual_tv.to_numpy(float)/h-1.0
    diff=obs_w-obs_h
    null_w=null_stack[eligible,:,1]/w[:,None]-1.0
    null_h=null_stack[eligible,:,2]/h[:,None]-1.0
    null_diff=null_w-null_h
    def p_upper(val, null):
        return float((1+np.count_nonzero(null>=val))/(REPS+1))
    nullmean=null_diff.mean(axis=0)
    nullmedian=np.median(null_diff,axis=0)
    nulltrim=mean_trimmed(null_diff,axis=0)
    rng=np.random.default_rng(2026102051+int(min_baseline*10000))
    boot=rng.integers(0,n,size=(999,n))
    bootmeans=diff[boot].mean(axis=1)
    result.update({
       'status':'estimable',
       'n_taxa_with_white_and_nonwhite_labels':int(((x.n_white_photo>0)&(x.n_nonwhite_photo>0)).sum()),
       'mean_relative_white_excess':float(obs_w.mean()),
       'mean_relative_nonwhite_excess':float(obs_h.mean()),
       'mean_white_minus_nonwhite_relative_excess':float(diff.mean()),
       'median_white_minus_nonwhite_relative_excess':float(np.median(diff)),
       'trimmed10_white_minus_nonwhite_relative_excess':float(mean_trimmed(diff)),
       'species_bootstrap_CI95_mean_white_minus_nonwhite':[float(z) for z in np.quantile(bootmeans,[.025,.975])],
       'matched_null_p_upper_mean_white_advantage':p_upper(float(diff.mean()),nullmean),
       'matched_null_p_upper_median_white_advantage':p_upper(float(np.median(diff)),nullmedian),
       'matched_null_p_upper_trim10_white_advantage':p_upper(float(mean_trimmed(diff)),nulltrim),
       'fraction_species_white_relative_excess_gt_nonwhite':float(np.mean(diff>0)),
       'n_per_cohort':{str(key):int(value) for key,value in x.cohort.value_counts().items()},
    })
    subset={}
    for cohort in ('discovery','validation','prospective','wave4','new_U100','leave_out_new_U100'):
        sel=(x.cohort.ne('new_U100') if cohort=='leave_out_new_U100' else x.cohort.eq(cohort)).to_numpy()
        if int(sel.sum())<20:
            subset[cohort]={'n':int(sel.sum()),'status':'under_20'}
            continue
        d=diff[sel]; nmat=null_diff[sel].mean(axis=0)
        subset[cohort]={'n':int(len(d)), 'mean_relative_white_minus_nonwhite':float(d.mean()),
                        'matched_null_p_upper_exploratory':p_upper(float(d.mean()),nmat)}
    result['cohort_sensitivity']=subset
    return result


def calculate(d,tag):
    rows=[];nulls=[]
    for _,g in d.groupby('inat_taxon_id',sort=True):
        ret=taxon_test(g,tag)
        if ret is None:continue
        row,vals=ret;rows.append(row);nulls.append(vals)
    frame=pd.DataFrame(rows)
    check(len(frame)>0,'no eligible species in geographic turnover')
    stack=np.stack(nulls,axis=0)
    sims=np.mean(stack,axis=0)
    observed=frame[['observed_tv','observed_abs_white_delta','observed_nonwhite_residual_tv']].mean().to_numpy()
    nullmean=sims.mean(axis=0)
    p=[float((1+np.count_nonzero(sims[:,col]>=observed[col]))/(REPS+1))
       for col in range(3)]
    return frame,{
        'n_species':int(len(frame)),
        'source_cohort_replication':summarize_by_source_cohort(frame,stack),
        'robust_white_effect_distribution':robust_white_summary(frame,stack),
        'relative_component_baseline_normalized':{
          'primary_min_null_baseline_0_01':relative_component_excess_test(frame,stack,.01),
          'sensitivity_min_null_baseline_0_005':relative_component_excess_test(frame,stack,.005),
          'sensitivity_min_null_baseline_0_02':relative_component_excess_test(frame,stack,.02),
        },
        'n_eligible_cell_pairs_not_independent':int(frame.n_cell_pairs.sum()),
        'mean_observed_between_cell_tv':float(observed[0]),
        'mean_null_tv':float(nullmean[0]),
        'mean_excess_tv':float(observed[0]-nullmean[0]),
        'mean_observed_white_delta':float(observed[1]),
        'mean_null_white_delta':float(nullmean[1]),
        'mean_excess_white_delta':float(observed[1]-nullmean[1]),
        'mean_observed_nonwhite_residual':float(observed[2]),
        'mean_null_nonwhite_residual':float(nullmean[2]),
        'mean_excess_nonwhite_residual':float(observed[2]-nullmean[2]),
        'white_share_of_mean_excess_turnover':float((observed[1]-nullmean[1])/(observed[0]-nullmean[0])) if abs(observed[0]-nullmean[0])>1e-12 else None,
        'species_bootstrap_95CI_white_share':bootstrap_component_share(frame),
        'vertex_permutation_upper_p_tv_exploratory':p[0],
        'vertex_permutation_upper_p_white_exploratory':p[1],
        'vertex_permutation_upper_p_nonwhite_residual_exploratory':p[2],
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
    selftest()
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

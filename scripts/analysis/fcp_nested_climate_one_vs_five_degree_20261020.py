#!/usr/bin/env python3
"""Decompose photo-colour/BIO5 association WITHIN 5-degree sampled geography.

Same source-verified photos and long-term WorldClim raster as the adjacent
within-cell/month analysis. Every eligible taxon x 5deg cell x calendar-month
has an exact algebraic decomposition of the white-nonwhite thermal contrast
into: (i) differences *within* its nested 1deg subcells and
(ii) variation *among* its nested 1deg subcells.

This is not genetic IBE, selection, annual weather, nor distribution history.
No coordinates, observer IDs, photos or individual row outputs are exported.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd

from fcp_withincell_month_climate_20261020 import (
    PHOTO_SHA, BIO_ZIP_SHA, PHOTO_ROWS, PHOTO_TAXA,
    photo_read, climate_sample, high_purity_mask, make_cell,
    COHORTS, SEED, REPS, stable_seed, require,
    MIN_LABEL_PER_STRATUM, MIN_PHOTOS_PER_STRATUM,
    MIN_DISTINCT_OBSERVERS_PER_LABEL, MIN_BIO5_RANGE,
)
from fcp_withincell_month_climate_20261020 import eligible_strata

AXES=("BIO1","BIO5")
COMPONENTS=("total_5deg_month","within_1deg_month","between_1deg_month")
MASTER="FCP_NESTED_1DEG_WITHIN_5DEG_SPLIT_20261020"

def local_cells(g):
    lon=g.longitude.to_numpy(float)
    lat=g.latitude.to_numpy(float)
    x=np.clip(np.floor(lon+180).astype(int),0,359)
    y=np.clip(np.floor(lat+90).astype(int),0,179)
    return y*360+x

def prepare_fixed(values,subcells):
    n=len(subcells)
    ids, inv=np.unique(subcells,return_inverse=True)
    nk=np.bincount(inv).astype(float)
    sumx=np.zeros((len(ids),2),float)
    np.add.at(sumx,inv,values)
    return inv,nk,sumx,values.sum(axis=0)

def exact_decomp(x,y,params):
    """E[T|white]-E[T|nonwhite] = within1 + between1 exactly."""
    inverse,counts,xsum,xsumall=params
    y=np.asarray(y,float)
    n=len(y);k=float(y.sum())
    require(0<k<n,"single-colour split in contrast")
    weighted=np.sum(x*y[:,None],axis=0)
    centered_cov=(weighted-k*xsumall/n)/n
    subwhite=np.bincount(inverse,weights=y,minlength=len(counts))
    covwithin=(weighted-np.sum(xsum*(subwhite/counts)[:,None],axis=0))/n
    scale=(k/n)*(1-k/n)
    total=centered_cov/scale
    within=covwithin/scale
    between=total-within
    if not np.allclose(total,within+between,rtol=0,atol=1e-8):
        raise RuntimeError("FCP_LOCAL_SCALE_COVARIANCE_IDENTITY_FAILED")
    return np.stack((total,within,between),axis=0)

def local_permutation(x,y,params,*,taxon,cell,month,purity):
    rng=np.random.default_rng(stable_seed(MASTER,taxon,cell,month,purity))
    obs=exact_decomp(x,y,params)
    inv=params[0]
    subgroups=[np.flatnonzero(inv==i) for i in range(len(params[1]))]
    null5=np.empty((REPS,3,2),float)
    null1=np.empty((REPS,3,2),float)
    for r in range(REPS):
        y5=y[rng.permutation(len(y))]
        null5[r]=exact_decomp(x,y5,params)
        y1=y.copy()
        for ix in subgroups:
            if len(ix)>1:
                y1[ix]=rng.permutation(y[ix])
        null1[r]=exact_decomp(x,y1,params)
    require(np.allclose(null1[:,2,:],obs[2,:],atol=1e-8,rtol=0),
            "within-1degree label shuffle changed between-subcell colour share")
    return obs,null5,null1

def analyse(d,purity,only_strata=None):
    g=make_cell(d,5)
    if purity:
        g=g.loc[high_purity_mask(g)].copy()
    diagnostics={
      "strata_5deg_month_any":0,
      "strata_two_colour_photo_opportunity":0,
      "strata_distinct_observers":0,
      "strata_bio5_range_ge_0p5C":0,
      "strata_5deg_qualified":0,
      "strata_5deg_with_2plus_nested_1deg_cells":0,
      "strata_with_1deg_subcell_containing_both_colour_labels":0,
      "strata_with_any_within1deg_temperature_variation":0,
    }
    groups={}
    valid_strata=set()
    for (tid,cell,month),part in g.groupby(["inat_taxon_id","spatialcell","month"],sort=True):
        key=(int(tid),int(cell),int(month))
        if only_strata is not None and key not in only_strata:continue
        diagnostics["strata_5deg_month_any"]+=1
        y=part.morph.eq("white").to_numpy()
        n=len(y)
        if n<MIN_PHOTOS_PER_STRATUM or min(int(y.sum()),n-int(y.sum()))<MIN_LABEL_PER_STRATUM:continue
        diagnostics["strata_two_colour_photo_opportunity"]+=1
        observers=part.observer_id.to_numpy()
        if (len(set(observers[y]))<MIN_DISTINCT_OBSERVERS_PER_LABEL or
            len(set(observers[~y]))<MIN_DISTINCT_OBSERVERS_PER_LABEL):continue
        diagnostics["strata_distinct_observers"]+=1
        x=part[list(AXES)].to_numpy(float)
        if np.ptp(x[:,1])<MIN_BIO5_RANGE:continue
        diagnostics["strata_bio5_range_ge_0p5C"]+=1
        diagnostics["strata_5deg_qualified"]+=1
        valid_strata.add(key)
        sub=local_cells(part)
        uk=np.unique(sub)
        if len(uk)>=2:diagnostics["strata_5deg_with_2plus_nested_1deg_cells"]+=1
        mixed=any(np.unique(y[sub==u]).size==2 for u in uk)
        if mixed:diagnostics["strata_with_1deg_subcell_containing_both_colour_labels"]+=1
        if any(np.ptp(x[sub==u,1])>1e-10 for u in uk):
            diagnostics["strata_with_any_within1deg_temperature_variation"]+=1
        fix=prepare_fixed(x,sub)
        obs,null5,null1=local_permutation(x,y,fix,taxon=int(tid),cell=int(cell),
                                         month=int(month),purity=purity)
        grp=groups.setdefault(int(tid),{
            "cohort":str(part.fcp_source_cohort.iloc[0]),
            "all_obs":[],"all_null5":[],"all_null1":[],
            "multi_obs":[],"multi_null5":[],"multi_null1":[],
        })
        require(grp["cohort"]==str(part.fcp_source_cohort.iloc[0]),"cohort non-disjoint")
        for key,v in (("all_obs",obs),("all_null5",null5),("all_null1",null1)):
            grp[key].append(v)
        if len(uk)>=2:
            for key,v in (("multi_obs",obs),("multi_null5",null5),("multi_null1",null1)):
                grp[key].append(v)
    result={}
    for scope in ("all","multi"):
        valid=[v for v in groups.values() if v[scope+"_obs"]]
        summary={}
        for cohort in ("pooled",*COHORTS):
            entries=valid if cohort=="pooled" else [v for v in valid if v["cohort"]==cohort]
            summary[cohort]=summarise_taxa(entries,scope,cohort,purity)
        result[scope+"_eligible_5deg_month_strata"]=summary
    taxon_means={tid:np.stack(v["all_obs"]).mean(axis=0)
                 for tid,v in groups.items() if v["all_obs"]}
    return ({"original_geocoded_observation_rows":int(len(g)),
             "quality":"high_purity" if purity else "all_classifiable",
             "conditional_exact_taxon_cell_month_support":only_strata is not None,
             "diagnostics":diagnostics,
             "results":result},
            valid_strata,taxon_means)

def summarise_taxa(entries,scope,cohort,purity):
    n=len(entries)
    if not n:return {"n_species":0,"status":"not_estimable"}
    obs=np.stack([np.stack(r[scope+"_obs"]).mean(axis=0) for r in entries])
    null5=np.stack([np.stack(r[scope+"_null5"]).mean(axis=0) for r in entries])
    null1=np.stack([np.stack(r[scope+"_null1"]).mean(axis=0) for r in entries])
    require(obs.shape==(n,3,2) and null5.shape==(n,REPS,3,2)
            and null1.shape==(n,REPS,3,2),"invalid source-permutation panel")
    mean=obs.mean(axis=0)
    worlds5=null5.mean(axis=0)
    worlds1=null1.mean(axis=0)
    require(np.allclose(mean[0],mean[1]+mean[2],atol=1e-9,rtol=0),
            "pooled white-temperature split inconsistent")
    one_sided_p=lambda arr, val:float((1+np.count_nonzero(arr>=val))/(REPS+1))
    if n>=20:
        p5={key:{axis:one_sided_p(worlds5[:,i,k],mean[i,k]) for k,axis in enumerate(AXES)}
              for i,key in enumerate(COMPONENTS)}
        p1={axis:one_sided_p(worlds1[:,1,k],mean[1,k]) for k,axis in enumerate(AXES)}
        rng=np.random.default_rng(stable_seed(MASTER,"taxon_boot",scope,cohort,purity))
        draws=rng.integers(0,n,size=(999,n))
        boot=obs[draws].mean(axis=1)
        ci={key:{axis:[float(x) for x in np.quantile(boot[:,i,k],[.025,.975])]
                 for k,axis in enumerate(AXES)} for i,key in enumerate(COMPONENTS)}
    else:
        p5=None;p1=None;ci=None
    return{
      "n_species":int(n),
      "n_5deg_month_strata_not_independent":int(sum(len(v[scope+"_obs"]) for v in entries)),
      "source_cohort":cohort,
      "mean_signed_white_minus_nonwhite_C":{component:{
             axis:float(mean[i,k]) for k,axis in enumerate(AXES)}
             for i,component in enumerate(COMPONENTS)},
      "positive_species_fraction":{component:{
             axis:float(np.mean(obs[:,i,k]>0)) for k,axis in enumerate(AXES)}
             for i,component in enumerate(COMPONENTS)},
      "species_bootstrap_95CI":ci,
      "five_degree_vertex_null_upper_p":p5,
      "one_degree_nested_conditioned_within_component_upper_p":p1,
      "pointwise_exploratory_no_multiple_testing_control":True,
      "if_no_within_cell_morph_climate_opportunity_within_component_is_zero":True,
      "n_original_null_worlds":REPS,
      "status":"exploratory" if n>=20 else "underpowered_no_p_value",
    }

def matched_quality_delta(all_photo,pure_photo):
    require(set(all_photo)==set(pure_photo),"strict versus full taxa not identical")
    ids=sorted(pure_photo)
    n=len(ids)
    require(n>=20,"matched quality subset has insufficient independent taxa")
    broad=np.stack([all_photo[k] for k in ids])
    strict=np.stack([pure_photo[k] for k in ids])
    diff=strict-broad
    rng=np.random.default_rng(stable_seed(MASTER,"paired_quality",n))
    draws=rng.integers(0,n,size=(999,n))
    boot=np.mean(diff[draws],axis=1)
    return {
      "n_same_taxa":n,
      "quality_subset_is_nested_not_independent":True,
      "species_equal_strict_minus_all_degrees_C":{
        component:{axis:float(diff[:,i,k].mean())
              for k,axis in enumerate(AXES)}
          for i,component in enumerate(COMPONENTS)},
      "species_bootstrap95_strict_minus_all_degrees_C":{
        component:{axis:[float(v) for v in np.quantile(boot[:,i,k],[.025,.975])]
              for k,axis in enumerate(AXES)}
          for i,component in enumerate(COMPONENTS)},
      "all_analysis_uses_identical_taxon_5deg_cell_calendar_month_strata":True,
      "strict_quality_may_change_the_number_of_1deg_subcells_within_those_strata":True,
      "posthoc_unequal_photo_count_after_filtering":True,
    }


def synthetic_test():
    # Four 10C and four 20C photos; white biased to warm subcell,
    # but temperature is constant INSIDE each 1degree subcell.
    x=np.array([[10.,10.]]*4+[[20.,20.]]*4)
    y=np.array([1,0,0,0,1,1,1,0],bool)
    params=prepare_fixed(x,np.array([0]*4+[1]*4))
    v=exact_decomp(x,y,params)
    require(np.allclose(v[0],[5.,5.]),"between-subcell total incorrect")
    require(np.allclose(v[1],0.),"within subcell must be zero")
    require(np.allclose(v[2],[5.,5.]),"between subcell must equal total")
    # One subcell with a white high temperature contrast.
    x2=np.array([[10.,10.],[10.,10.],[15.,15.],[15.,15.]])
    y2=np.array([False,False,True,True])
    k=exact_decomp(x2,y2,prepare_fixed(x2,np.array([0,0,0,0])))
    require(np.allclose(k[0],[5.,5.]) and np.allclose(k[1],[5.,5.])
            and np.allclose(k[2],0.),"within only synthetic case incorrect")
    obs,n5,n1=local_permutation(x,y,params,taxon=1,cell=2,month=6,purity=False)
    require(np.allclose(n1[:,2,:],obs[2,:]),"1deg-conditioned permutation broken")
    print("FCP_NESTED_1DEG_5DEG_CLIMATE_IDENTITY_AND_NULL_PASS")

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--photos",type=Path)
    ap.add_argument("--worldclim-zip",type=Path)
    ap.add_argument("--out",type=Path)
    ap.add_argument("--self-test",action="store_true")
    a=ap.parse_args()
    if a.self_test:synthetic_test();return
    require(a.photos is not None and a.worldclim_zip is not None and a.out is not None,
            "photos, worldclim and output mandatory")
    synthetic_test()
    d=climate_sample(photo_read(a.photos),a.worldclim_zip)
    broad,broad_strata,_=analyse(d,False)
    pure,strict_strata,pure_taxon_means=analyse(d,True)
    matched_broad,matched_strata,matched_taxon_means=analyse(d,False,
                  only_strata=strict_strata)
    require(strict_strata==matched_strata,
            "original photos do not retain all strict photo qualifying strata")
    matched=matched_quality_delta(matched_taxon_means,pure_taxon_means)
    require(matched["n_same_taxa"]==
            pure["results"]["all_eligible_5deg_month_strata"]["pooled"]["n_species"],
            "strict/original matched taxon denominator drift")
    result={
       "schema":"fcp_real_447500_within5deg_vs_between1deg_thermal_difference_v1",
       "source_photo_sha256":PHOTO_SHA,
       "worldclim_bio_zip_sha256":BIO_ZIP_SHA,
       "source_original_photo_rows":PHOTO_ROWS,
       "source_original_taxa":PHOTO_TAXA,
       "date_and_classifiable_geo_climate_rows":int(len(d)),
       "five_degree_all":broad,
       "five_degree_strict_purity":pure,
       "same_taxon_5degcell_month_original_opportunity":matched_broad,
       "exact_quality_opportunity_matched_delta":matched,
       "n_same_taxon_cell_month_strata_strict_and_original":len(strict_strata),
       "nested_spatial_decomposition_is_exact":True,
       "cal_month_match_across_years":True,
       "no_genetic_IBD_IBE_demographic_history_or_causal_adaptation":True,
       "data_quality":"classifiable original ROI photo labels; strict ROI purity is compositional sensitivity only",
       "climate_type":"WorldClim 1970-2000 long-term normals, not year-specific temperatures",
    }
    a.out.mkdir(parents=True,exist_ok=True)
    (a.out/"nested_five_degree_one_degree_climate_components.json").write_text(
         json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+"\n")
    compact=lambda z:{key:z["results"][key]["pooled"] for key in z["results"]}
    print("FCP_REAL_NESTED_CLIMATE_RESULTS",json.dumps({
         "n_source":len(d),"broad":compact(broad),"strict":compact(pure),
         "matched_original":compact(matched_broad),
         "matched_delta":matched,
         "opportunity_broad":broad["diagnostics"],
         "opportunity_strict":pure["diagnostics"]},sort_keys=True,allow_nan=False))

if __name__=="__main__":
    main()

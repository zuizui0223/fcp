#!/usr/bin/env python3
"""FCP all 4475 whole-source observer D repeatability + univariate BIO5 Firth.

All 4475 taxon IDs remain in the public species-level analysis table; 447500
authentic source photos and 227107 classifiable photo labels are mandatory.
A species with only one observed photo colour has NO estimable within-species
effect on white-vs-coloured photographs; no pseudocount turns it into a genotype
or a climatic slope. Univariate slopes here are not the earlier restricted
three-predictor 1574-species multivariate fit or a test of local adaptation.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))

import numpy as np
import pandas as pd

from fcp_full4475_hypothesis_census_20261021 import (
    SOURCE_SHA,BIO_ARCHIVE_SHA,TOTAL_TAXA,CLASSIFIABLE_ROWS,
    load_original,sample_climate,sampled_diversity,check,
)
from scripts.analysis.run_fcp_447500_environment_model_20261018 import fit_firth

SEED=20261021
ORIGINAL_PRIOR_SAMPLES=100


def rank(values):
    a=pd.Series(values,dtype=float)
    return a.rank(method="average").to_numpy(float)


def pearson(a,b):
    a=np.asarray(a,float)
    b=np.asarray(b,float)
    a=a-a.mean();b=b-b.mean()
    den=float(np.sqrt(np.sum(a*a)*np.sum(b*b)))
    return float(np.sum(a*b)/den) if den>0 else None


def observer_split(g):
    n=len(g)
    obs=g.observer_id.dropna().unique()
    if len(obs)<4:
        return {"n_photo":n,"n_distinct_observer":int(len(obs)),
                "n_split_a":0,"n_split_b":0,"D_split_a":None,"D_split_b":None,
                "D_split_status":"lt_four_independent_observers"}
    order=sorted(obs,key=lambda v:hashlib.sha256((str(SEED)+":"+str(v)).encode()).hexdigest())
    a=set(order[::2])
    if not a or len(a)==len(order):
        raise RuntimeError("observer-disjoint split error")
    first=g.loc[g.observer_id.isin(a)]
    second=g.loc[g.observer_id.notna() & ~g.observer_id.isin(a)]
    if len(first)<2 or len(second)<2:
        status="lt_two_classified_in_half"
        da=db=None
    else:
        status="estimable"
        da=sampled_diversity(first)
        db=sampled_diversity(second)
    return {"n_photo":n,"n_distinct_observer":int(len(obs)),
            "n_split_a":int(len(first)),"n_split_b":int(len(second)),
            "D_split_a":da,"D_split_b":db,
            "D_split_status":status}


def fit_bio5(g):
    n=len(g)
    y=g.morph.eq("white").to_numpy(float)
    w=int(y.sum());c=int(n-w)
    base={
     "n_classified_photo":n,
     "n_white":w,"n_nonwhite":c,
     "n_BIO5_usable_photo":0,
     "n_WHITE_BIO5_usable":0,"n_NONWHITE_BIO5_usable":0,
     "bio5_observed_SD_C":None,"bio5_observed_span_C":None,
     "BIO5_Firth_beta_logodds_per_C":None,
     "BIO5_Firth_Wald_se_approx_per_C":None,
     "BIO5_Firth_observer_association_status":"not_estimable",
     "no_genotype_or_fitness":True,
    }
    temp=pd.to_numeric(g.BIO5,errors="coerce").to_numpy(float)
    ok=np.isfinite(temp)
    xx=temp[ok];yy=y[ok]
    base["n_BIO5_usable_photo"]=len(xx)
    base["n_WHITE_BIO5_usable"]=int(yy.sum())
    base["n_NONWHITE_BIO5_usable"]=int(len(yy)-yy.sum())
    if len(xx)<3:
        base["BIO5_Firth_observer_association_status"]="n_climate_lt3";return base
    if min(base["n_WHITE_BIO5_usable"],base["n_NONWHITE_BIO5_usable"])==0:
        base["BIO5_Firth_observer_association_status"]="one_photo_colour_state_only";return base
    sd=float(np.std(xx))
    base["bio5_observed_SD_C"]=sd
    base["bio5_observed_span_C"]=float(np.ptp(xx))
    if sd<.05 or np.ptp(xx)<.1:
        base["BIO5_Firth_observer_association_status"]="insufficient_climate_contrast";return base
    x=np.column_stack([np.ones(len(xx)),(xx-np.mean(xx))/sd])
    fit=fit_firth(x,yy)
    base["BIO5_Firth_observer_association_status"]=str(fit["status"])
    if fit["status"]=="ok":
        base["BIO5_Firth_beta_logodds_per_C"]=float(fit["beta"][1]/sd)
        base["BIO5_Firth_Wald_se_approx_per_C"]=float(np.sqrt(fit["cov"][1,1])/sd)
    return base


def bootstrap_genera(df,col,seed,n=999):
    vals={}
    for row in df.itertuples(index=False):
        v=getattr(row,col)
        if v is None or not np.isfinite(v):
            continue
        key=str(row.species).split()[0]
        vals.setdefault(key,[]).append(float(v))
    if len(vals)<5:
        return {"n_genera":len(vals),"confidence_interval_95":None}
    genera=sorted(vals)
    totals=np.array([sum(vals[x]) for x in genera],float)
    counts=np.array([len(vals[x]) for x in genera],int)
    rng=np.random.default_rng(seed)
    draw=rng.integers(0,len(genera),size=(n,len(genera)))
    means=totals[draw].sum(axis=1)/counts[draw].sum(axis=1)
    return {"n_genera":len(genera),
            "confidence_interval_95":[float(x) for x in np.quantile(means,[.025,.975])],
            "n_resampling_worlds":n}


def stats(rows):
    df=pd.DataFrame(rows)
    check(len(df)==TOTAL_TAXA and df.inat_taxon_id.is_unique,
          "species-level observer-climate ledger has omitted taxa")
    s=df.loc[df.D_split_status.eq("estimable")].copy()
    r=pearson(rank(s.D_split_a),rank(s.D_split_b))
    if r is None:
        raise RuntimeError("species-level observer-disjoint D ranks degenerate")
    fits=df.loc[df.BIO5_Firth_observer_association_status.eq("ok")].copy()
    if len(fits):
        vals=fits.BIO5_Firth_beta_logodds_per_C.to_numpy(float)
        trim=np.sort(vals)[int(.1*len(vals)):len(vals)-int(.1*len(vals))]
        ci=bootstrap_genera(fits,"BIO5_Firth_beta_logodds_per_C",SEED)
    else:
        vals=np.array([]);trim=np.array([]);ci=None
    bycohort={}
    for cohort,g in df.groupby("source_cohort",sort=True):
        subset=g.loc[g.BIO5_Firth_observer_association_status.eq("ok")]
        bycohort[str(cohort)]={
         "n_taxa_total":len(g),
         "n_observer_D_split_estimable":int(g.D_split_status.eq("estimable").sum()),
         "n_species_BIO5_Firth_estimable":int(len(subset)),
         "mean_Firth_logodds_per_C":float(subset.BIO5_Firth_beta_logodds_per_C.mean())
            if len(subset) else None,
         "fraction_positive_given_estimation":float((subset.BIO5_Firth_beta_logodds_per_C>0).mean())
            if len(subset) else None,
        }
    out={
      "schema":"fcp_all4475_observer_split_D_and_univariate_Firth_BIO5_posthoc_v1",
      "complete_source_447500_photos":True,
      "complete_source_4475_taxa":True,
      "source_classifiable_photo_rows":CLASSIFIABLE_ROWS,
      "observer_D_split": {
         "n_species_estimable":len(s),
         "spearman_observer_disjoint_split_D":r,
         "n_species_missing_but_kept_in_roster":TOTAL_TAXA-len(s),
         "split_rule":"hash-sort the original observer IDs, alternating observers into independent two groups without using colour labels",
         "one_split_only_not_200_frozen_splits":True,
      },
      "BIO5_univariate_photo_classification_Firth": {
        "n_species_estimable":int(len(fits)),
        "n_species_missing_but_kept_in_roster":TOTAL_TAXA-int(len(fits)),
        "mean_beta_logodds_per_C":float(np.mean(vals)) if len(vals) else None,
        "median_beta_logodds_per_C":float(np.median(vals)) if len(vals) else None,
        "trimmed_10pct_mean_beta_logodds_per_C":float(np.mean(trim)) if len(trim) else None,
        "positive_fraction":float(np.mean(vals>0)) if len(vals) else None,
        "genus_bootstrap_mean_95CI":ci,
        "status":"exploratory univariate Firth; not earlier conditional multivariate Firth",
        "conditioned_BIO5_photo_time":"1970-2000 long-term WorldClim; no flowering weather or fitness",
        "strict_bio5_SD_minimum_C":0.05,
        "strict_bio5_span_minimum_C":0.1,
      },
      "by_source_cohort":bycohort,
      "no_colour_morph_genetic_alleles_measured":True,
      "all_taxon_source_row_statuses_in_saved_csv":True,
    }
    return out


def synthetic():
    g=pd.DataFrame({"morph":["white"]*3+["blue_purple"]*3,
                    "observer_id":[1,2,3,4,5,6],
                    "BIO5":[20.,22.,24.,26.,28.,30.]})
    row=fit_bio5(g)
    check(row["BIO5_Firth_observer_association_status"]=="ok",
          "original Firth package estimator not fitted in basic fixture")
    p=observer_split(g)
    check(p["D_split_status"]=="estimable" and p["n_split_a"]==3 and p["n_split_b"]==3,
          "observer split not valid")
    print("FCP_ALL4475_OBSERVER_AND_BIO5_FIRTH_SYNTHETIC_PASS")


def main():
    a=argparse.ArgumentParser()
    a.add_argument("--photos",type=Path)
    a.add_argument("--worldclim-zip",type=Path)
    a.add_argument("--out",type=Path)
    a.add_argument("--self-test",action="store_true")
    o=a.parse_args()
    synthetic()
    if o.self_test:return
    check(o.photos and o.worldclim_zip and o.out,"sources not set")
    old,d=load_original(o.photos)
    d=sample_climate(d,o.worldclim_zip)
    rows=[]
    for index,(tid,g) in enumerate(d.groupby("inat_taxon_id",sort=True)):
        rows.append({
          "inat_taxon_id":int(tid),
          "species":str(g.species.iloc[0]),
          "source_cohort":str(g.fcp_source_cohort.iloc[0]),
          **observer_split(g),**fit_bio5(g)})
        if (index+1)%500==0:
            print("FCP_ALL4475_FIRTH_AND_OBSERVER_PARTITIONS",index+1,flush=True)
    result=stats(rows)
    o.out.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(rows).to_csv(o.out/"all4475_species_observer_BIO5_Firth_NO_COORDINATES.csv",
                               index=False,float_format="%.12g")
    (o.out/"all4475_observer_D_and_Firth_BIO5_result.json").write_text(
       json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("FCP_REAL_ALL4475_OBSERVER_FIRTH_COMPLETE",json.dumps(result,sort_keys=True,allow_nan=False))

if __name__=="__main__":
    main()

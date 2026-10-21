#!/usr/bin/env python3
"""Exploratory 4,475-species joint phenotype-space vs geographical-climatic fate.

Read exactly one original source-backed, all-4,475-taxon signed species matrix.
Primary question: does W axis strength predict geographical photo-colour sorting
or climatic photo-colour associations after classifiability/observer/colour
composition and source-cohort controls?

No new photo-level permutation: all three derived response variables come from
the SAME original iNaturalist images and are statistically non-independent.
Genus-cluster intervals are descriptive, conditional on selected-photo
opportunity and on fixed residual nuisance fits, not causal/genetic inference.
The original species without estimable local or climate contrasts stay NA.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata

EXPECTED_SOURCE="23d3daadc9799c2e7b5867d54ccc1ed98369ba79a81b72dc7f831e5197e63c49"
EXPECTED_N=4475
SEED=20261021
BOOT=1999
COHORTS=("discovery","validation","prospective","wave4","new_U100")

COMPARE={
 "W_geometry_vs_geographic_sorting":("W_axis_sq_alignment","global_vs_local_50km_colour_mismatch_delta"),
 "W_geometry_vs_photo_BIO5_slope":("W_axis_sq_alignment","ob_BIO5_Firth_beta_logodds_per_C"),
 "W_geometry_vs_photo_temperature_seasonality":("W_axis_sq_alignment","clim_BIO4_signed_white_minus_nonwhite_SD"),
 "geographic_sorting_vs_photo_BIO5_slope":("global_vs_local_50km_colour_mismatch_delta","ob_BIO5_Firth_beta_logodds_per_C"),
}


def check(ok,message):
    if not ok:raise RuntimeError("FCP_ALL4475_JOINT_AXIS_FATE_EVIDENCE_INVALID: "+message)


def digest(p):
    h=hashlib.sha256()
    with Path(p).open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):h.update(chunk)
    return h.hexdigest()


def ranks(z):
    a=np.asarray(z,float)
    check(np.isfinite(a).all(),"nonfinite values passed to ranks")
    return rankdata(a,method="average")


def corr(a,b):
    x=np.asarray(a,float);y=np.asarray(b,float)
    cx=x-x.mean();cy=y-y.mean()
    den=np.sqrt((cx*cx).sum()*(cy*cy).sum())
    return float((cx*cy).sum()/den) if den>1e-14 else None


def rank_controls(frame, pair):
    # These photography-source controls affect both true colour biology
    # and ascertainment and should be viewed as diagnostic, not adjustment
    # for independent genetic ancestry/history.
    cols=["D_photo_pairwise_unbiased",
          "p_white_photo_label",
          "n_classified_photo",
          "n_distinct_observers_source"]
    if "geographic_sorting" in pair:
        cols+=["n_local50km_colour_photo_pairs"]
    if "BIO5" in pair:
        cols+=["ob_bio5_observed_SD_C"]
    use=frame.copy()
    X=[np.ones(len(use))]
    for col in cols:
        check(col in use.columns,"missing source covariate "+col)
        v=pd.to_numeric(use[col],errors="coerce").to_numpy(float)
        if not np.isfinite(v).all():
            # Never fill a missing covariate with a numeric zero.
            return None,cols,"missing_source_covariates"
        z=ranks(v)
        sd=z.std()
        if sd>1e-10:
            X.append((z-z.mean())/sd)
    for cohort in COHORTS[1:]:
        X.append(use.source_cohort.eq(cohort).to_numpy(float))
    M=np.column_stack(X)
    require_rank=np.linalg.matrix_rank(M)
    if require_rank<M.shape[1]:
        return None,cols,"rank_deficient_controls"
    return M,cols,"ok"


def fixed_nuisance_partial_r(frame,name,xcol,ycol):
    frame=frame.copy()
    a=ranks(frame[xcol].to_numpy(float))
    b=ranks(frame[ycol].to_numpy(float))
    xctrl,cols,status=rank_controls(frame,name)
    if status!="ok":
        return None,None,None,{"status":status,"controls":cols}
    resid_a=a-xctrl@np.linalg.lstsq(xctrl,a,rcond=None)[0]
    resid_b=b-xctrl@np.linalg.lstsq(xctrl,b,rcond=None)[0]
    out=corr(resid_a,resid_b)
    return out,resid_a,resid_b,{
       "status":"estimable" if out is not None else "degenerate_residual",
       "controls":cols,
       "n_effective_control_columns":int(xctrl.shape[1]),
       "rank_residual_sd_W_or_x":float(np.std(resid_a)),
       "rank_residual_sd_y":float(np.std(resid_b)),
    }


def genus_bootstrap_fixed_fit_residual_corr(genus,x,y,seed):
    g=pd.DataFrame({"genus":genus,"x":x,"y":y})
    agg=g.assign(xx=g.x*g.x,yy=g.y*g.y,xy=g.x*g.y).groupby("genus",sort=True).agg(
        n=("x","size"),sx=("x","sum"),sy=("y","sum"),
        sxx=("xx","sum"),syy=("yy","sum"),sxy=("xy","sum"))
    n_gen=len(agg)
    if n_gen<20:return None,n_gen
    vals=agg[["n","sx","sy","sxx","syy","sxy"]].to_numpy(float)
    rng=np.random.default_rng(seed)
    draws=rng.integers(0,n_gen,size=(BOOT,n_gen))
    # Each genus is one bootstrap cluster, even if the same genus contains
    # multiple included species.
    totals=vals[draws].sum(axis=1)
    N,sx,sy,sxx,syy,sxy=totals.T
    vx=sxx-sx*sx/N
    vy=syy-sy*sy/N
    cv=sxy-sx*sy/N
    rs=np.divide(cv,np.sqrt(np.maximum(vx*vy,0)),
                 out=np.full(BOOT,np.nan),where=(vx>1e-12)&(vy>1e-12))
    stable=rs[np.isfinite(rs)]
    if len(stable)<BOOT*.95:return None,n_gen
    return [float(x) for x in np.quantile(stable,[.025,.975])],n_gen


def compare(df,name,xcol,ycol):
    allowed=pd.to_numeric(df[xcol],errors="coerce").notna()&pd.to_numeric(df[ycol],errors="coerce").notna()
    # Additional BIO5-SD nuisance covariate may itself be missing.
    if "BIO5" in name:
        allowed &= pd.to_numeric(df.ob_bio5_observed_SD_C,errors="coerce").notna()
    f=df.loc[allowed].copy()
    n=len(f)
    if n<30:return {"n_species":n,"status":"insufficient_joint_opportunity"}
    check(f.inat_taxon_id.is_unique,"joined analytical data duplicates species IDs")
    rraw=corr(ranks(f[xcol]),ranks(f[ycol]))
    point,rx,ry,details=fixed_nuisance_partial_r(f,name,xcol,ycol)
    if details["status"]!="estimable":
        return {"n_species":n,"status":details["status"],"diagnostic":details}
    bci,n_gen=genus_bootstrap_fixed_fit_residual_corr(
       f.species.str.split().str[0].to_numpy(),rx,ry,SEED+sum(ord(c) for c in name))
    co=[]
    for cohort,g in f.groupby("source_cohort"):
        xa=g[xcol].to_numpy(float);yb=g[ycol].to_numpy(float)
        co.append({"cohort":str(cohort),"n_species":len(g),
            "spearman_raw":corr(ranks(xa),ranks(yb)) if len(g)>1 else None})
    loo=[]
    for leave in COHORTS:
        part=f.loc[~f.source_cohort.eq(leave)]
        if len(part)<100:continue
        fit,xx,yy,sens=fixed_nuisance_partial_r(part,name,xcol,ycol)
        loo.append({"left_out":leave,"n_remaining_species":len(part),
                    "partial_rank_r":fit,"status":sens["status"]})
    return{
        "n_species":n,
        "n_missing_taxa_retained_elsewhere":EXPECTED_N-n,
        "n_genera":n_gen,
        "raw_spearman":rraw,
        "partial_rank_r_controlling_photo_opportunity_and_source":point,
        "genus_bootstrap_95CI_conditional_fixed_residuals":bci,
        "control_model":details,
        "by_disjoint_species_cohort_raw":co,
        "leave_one_source_cohort_out_diagnostic":loo,
        "exploratory_posthoc":True,
        "potential_shared_source_photo_measurement_artifact":True,
        "not_test_of_genetic_IBD_IBE_history_or_selection":True,
    }


def synthetic():
    rng=np.random.default_rng(20261021)
    x=np.arange(40,dtype=float)
    y=x*2+3
    check(abs(corr(ranks(x),ranks(y))-1)<1e-12,
          "rank-correlation fixture")
    z=np.array([5,4,3,2,1],float)
    check(corr(ranks(z),ranks(np.arange(5)))<-.99,
          "opposite monotone-rank fixture")
    print("FCP_NP_4475_GEOMETRY_FATE_JOINT_DIAGNOSTIC_SYNTHETIC_PASS")


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--species-ledger",type=Path)
    p.add_argument("--receipt",type=Path)
    p.add_argument("--out",type=Path)
    p.add_argument("--self-test",action="store_true")
    a=p.parse_args()
    synthetic()
    if a.self_test:return
    check(a.species_ledger and a.receipt and a.out,"real input paths are required")
    check(digest(a.species_ledger)==EXPECTED_SOURCE,
          "unified 4475 source matrix SHA does not match frozen complete archive")
    receipt=json.loads(a.receipt.read_text())
    check(receipt["merged_species_ledger_sha256"]==EXPECTED_SOURCE
          and receipt["all_species_preserved"] is True
          and receipt["n_total_taxa"]==EXPECTED_N,
          "unified 4475 source receipt not authentic")
    df=pd.read_csv(a.species_ledger,low_memory=False)
    check(len(df)==EXPECTED_N and df.inat_taxon_id.is_unique,
          "analysis did not retain all 4475 original taxa")
    check(set(df.source_cohort)==set(COHORTS),"not all five original source cohorts")
    r={name:compare(df,name,*pair) for name,pair in COMPARE.items()}
    output={
      "schema":"fcp_entire4475_photo_W_geometry_geographic_climatic_fate_joint_exploratory_v1",
      "all_species_in_common_parent_roster":EXPECTED_N,
      "actual_source_species_matrix_sha256":EXPECTED_SOURCE,
      "source_of_447500_photos_sha256":receipt["source_actual_photo_sha256"],
      "primary_hypothesis":"similar photographed W axes do not automatically imply identical spatial or climatic associations among species",
      "results":r,
      "no_independent_photo_provider_or_genetic_data":True,
      "no_genotype_phenotype_fitness_history_causal_estimation":True,
      "strictly_posthoc_descriptive_association_from_same_classified_photos":True,
      "H4_local_genetic_morphs_testable":False,
      "H5_same_photographic_individual_pigment_mechanisms_testable":False,
      "note":"CIs resample genera of existing nuisance-model residuals without refitting; not a prospective permutation, not a mechanistic no-effect proof. Negative or null co-variation may reflect measurement error.",
    }
    a.out.mkdir(parents=True,exist_ok=True)
    (a.out/"all4475_W_vs_geography_climate_fate_result.json").write_text(
        json.dumps(output,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("FCP_REAL_4475_GEOMETRY_VS_FATE_JOINT_RESULT",json.dumps(output,sort_keys=True,allow_nan=False))

if __name__=="__main__":main()

#!/usr/bin/env python3
"""Independent original-photo macroclimate inference audit and cohort-deletion sensitivity.

The original source has five species-disjoint iNaturalist PHOTO cohorts; one
later high-volume cohort supplies many taxa. Check that the positive thermal
white-palette cross-species mean is not entirely created by one cohort and
compare the original/observer-thinned source on matched named taxon IDs.

A positive photo gradient is NOT evidence of inherited white petal morphs,
pigmentation costs, site-specific fitness, or independent botanical background.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

DOMAINS=("thermal_level","precip_amount","BIO4_seasonality",
         "BIO15_seasonality","SRADmean")
HUES=("white","yellow","orange","red","pink","magenta","purple","blue","bronze")
COMP=("total","between_original_photo_labels","within_original_photo_labels")
COHORTS=("discovery","validation","prospective","wave4","new_U100")
B=999
SEED=20261024
ARTIFACT_ID=12549840336
DIGEST="3f290b8028bbfcdb6fa2d1a4ae85f3155b21786d7c8d48175af4805c78c099a0"


def check(ok,desc):
    if not ok:raise RuntimeError("FCP_SHARED_CLIMATE_INDEPENDENT_SOURCE_AUDIT: "+desc)


def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):h.update(b)
    return h.hexdigest()


def load(source):
    p=source/"shared_thermal_moisture_nine_hue_source_result.json"
    j=json.loads(p.read_text())
    check(j["original_photos"]==447500 and j["original_species"]==4475
          and j["source_classifiable_photos"]==227107,
          "original 447500/4475/227107 photo denominator mismatch")
    check(tuple(j["fixed_five_climate_domains"])==DOMAINS,
          "original source climate domain variables changed")
    panels={}
    files={"source_result_json":sha(p)}
    for name in ("original","observer_thinned"):
        csv=source/f"fcp4475_{name}_thermal_moisture_support.csv"
        npz=source/f"fcp4475_{name}_thermal_moisture_nine_hue_coefficients.npz"
        files[name+"_support_csv"]=sha(csv)
        files[name+"_coefficient_npz"]=sha(npz)
        df=pd.read_csv(csv,low_memory=False)
        check(len(df)==4475 and df.inat_taxon_id.is_unique
              and df.n_original_classifiable.sum()==227107,
              "not full original 4475 source taxa in numeric support")
        check(set(df.source_cohort)==set(COHORTS),"not original five species source cohorts")
        with np.load(npz,allow_pickle=False) as z:
            ids=z["inat_taxon_id"]
            val=z["conditional_slope"].copy()
            raw=z["signed_photo_climate_covariance"].copy()
            check(tuple(z["domains"])==DOMAINS
                  and tuple(z["photo_hues"])==HUES
                  and tuple(z["photographed_components"])==COMP,
                  "original source climatic/photo colour coefficient order changed")
        check(np.array_equal(ids,df.inat_taxon_id.to_numpy()),
              "source taxon ID order differs")
        check(val.shape==(4475,5,3,9) and raw.shape==(4475,5,3,9),
              "full source 5 climate × 9 hue coefficient tensor shape drift")
        check(np.max(np.abs(np.nansum(val,axis=-1)[
             np.isfinite(val).all(axis=-1)]))<1e-7,
              "original source photo-climate signed colour fractions do not sum to zero")
        check(np.max(np.abs(np.nansum(raw,axis=-1)[
             np.isfinite(raw).all(axis=-1)]))<1e-7,
              "original source photographic covariance not compositional")
        for macro in ("thermal_level","precip_amount"):
            k=DOMAINS.index(macro)
            spec=j["results"][name]["components"]["total"]["per_domain"][macro]
            src=val[:,k,0,0]
            eligible=np.isfinite(src)
            check(spec["n_species"]==int(eligible.sum()),
                  "source thermal/precipitation eligibility changes")
            check(abs(float(src[eligible].mean())-
                    spec["white_photo_fraction_mean_partial_slope"])<1e-10,
                  "original macroclimate photo-white species-equal mean changed")
        panels[name]={"support":df,"coef":val,"cov":raw}
    return j,panels,files


def genus_boot(frame,col,seed,rep=B):
    values=frame[col].to_numpy(float)
    use=np.isfinite(values)
    x=values[use]
    if len(x)<35:return {"status":"not_estimable","n_species":int(len(x))}
    group,genes=pd.factorize(frame.loc[use,"genus"],sort=True)
    G=len(genes)
    sums=np.zeros(G,float)
    count=np.bincount(group,minlength=G)
    np.add.at(sums,group,x)
    choices=np.random.default_rng(seed).integers(G,size=(rep,G))
    bs=sums[choices].sum(axis=1)/count[choices].sum(axis=1)
    ci=np.quantile(bs,[.025,.975])
    return {
       "status":"whole_genus_cluster_original_source",
       "n_species":int(len(x)),
       "n_genera":int(G),
       "mean_original_photo_white_thermal_partial_slope":float(x.mean()),
       "fraction_taxa_positive":float(np.mean(x>0)),
       "genus_bootstrap_marginal95":[float(v) for v in ci],
       "marginal95_excludes_zero":bool(ci[0]>0 or ci[1]<0),
       "n_genus_draws":rep
    }


def audit(source):
    j,panels,files=load(source)
    output={}
    for idx,name in enumerate(("original","observer_thinned")):
        df=panels[name]["support"].copy()
        val=panels[name]["coef"]
        df["thermal_original_photo_white_slope"]=val[:,0,0,0]
        df["precip_original_photo_white_slope"]=val[:,1,0,0]
        cases={}
        for i,omit in enumerate(COHORTS):
            subset=df.loc[df.source_cohort!=omit]
            cases[omit]=genus_boot(subset,"thermal_original_photo_white_slope",
                                  seed=SEED+idx*703+i*73)
        bycohort={}
        for i,c in enumerate(COHORTS):
            subset=df.loc[df.source_cohort==c]
            bycohort[c]=genus_boot(subset,"thermal_original_photo_white_slope",
                                  seed=SEED+idx*701+i*79)
        vals=[v["mean_original_photo_white_thermal_partial_slope"]
              for v in bycohort.values() if v["status"]!="not_estimable"]
        output[name]={
          "source_cohort_exclusion_thermal_photo_hue":cases,
          "five_source_disjoint_cohorts":bycohort,
          "cohort_equal_weight_mean_of_five_photo_hue_slopes":float(np.mean(vals)),
          "five_cohort_mean_signs":[int(np.sign(v)) for v in vals],
          "main_species_equal_photo_hue_slope":float(
              df.thermal_original_photo_white_slope.mean()),
        }
    full=panels["original"]["support"].copy()
    thin=panels["observer_thinned"]["support"].copy()
    v0=panels["original"]["coef"][:,0,0,0]
    v1=panels["observer_thinned"]["coef"][:,0,0,0]
    check(full.inat_taxon_id.equals(thin.inat_taxon_id),"original/observer taxon panel mismatch")
    overlap=np.isfinite(v0)&np.isfinite(v1)
    overlap_df=full.loc[overlap,["genus"]].copy()
    overlap_df["difference"]=v0[overlap]-v1[overlap]
    result=genus_boot(overlap_df,"difference",seed=SEED+994)
    result["n_source_eligible_matched_taxa"]=int(np.sum(overlap))
    result["original_source_matched_taxa_mean"]=float(v0[overlap].mean())
    result["observer_thinned_matched_taxa_mean"]=float(v1[overlap].mean())
    result["pair_not_independent_source_photographs"]=True
    return {
       "schema":"fcp4475_macroclimate_original_photo_cohort_holdout_source_QA_v1",
       "source_github_actions_run_id":40066224301,
       "source_artifact_id":ARTIFACT_ID,
       "source_zip_sha256":DIGEST,
       "source_original_photos":447500,"source_original_taxa":4475,
       "source_original_nine_hue_classifiable":227107,
       "original_derived_input_file_sha256":files,
       "genus_cluster_resamples":B,
       "conclusions_describe_image_hue_not_petal_genetic_morphs":True,
       "source_taxa_not_imputed":True,
       "source_cohorts_disjoint_but_one_photo_provider":True,
       "five_cohort_original_photo_thermal_consistency":output,
       "one_photo_per_observer_matched_taxa_difference":result,
    }


def fixture():
    g=pd.DataFrame({"genus":["A","A","B","B","C","C"]*12,
                    "white":[1.,2.,2.,4.,4.,5.]*12})
    m=genus_boot(g,"white",seed=SEED)
    check(m["n_species"]==72 and m["mean_original_photo_white_thermal_partial_slope"]>2,
          "genus bootstrap synthetic mean wrong")
    print("FCP_MACRO_THERMAL_FIVE_COHORT_SOURCE_AUDIT_SYNTHETIC_PASS")


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path)
    p.add_argument("--out",type=Path)
    p.add_argument("--self-test",action="store_true")
    args=p.parse_args()
    fixture()
    if args.self_test:return
    check(args.source and args.out,"immutable original source artifact required")
    data=audit(args.source)
    args.out.mkdir(parents=True,exist_ok=True)
    (args.out/"five_cohort_thermal_leave_one_out_original_source.json").write_text(
        json.dumps(data,sort_keys=True,indent=2,allow_nan=False)+"\n")
    print("FCP_MACRO_THERMAL_SOURCE_COHORT_AUDIT_PASS",json.dumps({
        p:{
          "omit_each":[{"omitted":k,"n":v["n_species"],"mean":v["mean_original_photo_white_thermal_partial_slope"],"ci":v["genus_bootstrap_marginal95"]}
                       for k,v in val["source_cohort_exclusion_thermal_photo_hue"].items()],
          "mean_equal_five":val["cohort_equal_weight_mean_of_five_photo_hue_slopes"],
        } for p,val in data["five_cohort_original_photo_thermal_consistency"].items()
    },sort_keys=True,allow_nan=False))
    print("FCP_MACRO_THERMAL_SOURCE_AUDIT_FULL_JSON",
          json.dumps(data,sort_keys=True,separators=(',',':'),allow_nan=False))


if __name__=="__main__":
    main()

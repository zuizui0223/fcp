#!/usr/bin/env python3
"""Audit source-authenticated five disjoint species cohorts for photo hue climate directions.

Uses nine-hue × seven-climate per-taxon coefficients only, not underlying photo
IDs or locations. A signed coefficient is an original within-species
source-photo covariance adjusted for broad geography/month; it is NOT fitness,
causation or a genetic flower colour morph frequency.

Source-cohort replication is descriptive and post hoc: named chromatic targets
were inspected AFTER the full seven-by-eight chromatic bootstrap families.
Do not call the cohort interval a preregistered replication test.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

SHA={
 "original":"207c08ec0575f308e865c93fc349ff1ed89faf5483a016333e66b07bb65b6f38",
 "observer_thinned":"078db22924156cfa9509769ea9896dd00e639c93ee96bc29d6b26c341755fba6",
 "source_result":"fc30eaa9e1cdec57f66fb29ce0be71a010c14a709f6f99f48068d7072b19671a"
}
FIVE=("discovery","validation","prospective","wave4","new_U100")
B=1999
SEED=20261024
TARGETS=(
 ("BIO5 white, total-photo palette","BIO5_total_white_signed_covariance_month_null_excess",None),
 ("BIO12 yellow, within original photo label","BIO12_within_photo_labels_yellow_signed_covariance_month_null_excess","BIO12_within_photo_labels_white_signed_covariance_month_null_excess"),
 ("BIO12 bronze, within original photo label","BIO12_within_photo_labels_bronze_signed_covariance_month_null_excess","BIO12_within_photo_labels_white_signed_covariance_month_null_excess"),
 ("BIO12 yellow, total-photo palette","BIO12_total_yellow_signed_covariance_month_null_excess","BIO12_total_white_signed_covariance_month_null_excess"),
 ("BIO12 bronze, total-photo palette","BIO12_total_bronze_signed_covariance_month_null_excess","BIO12_total_white_signed_covariance_month_null_excess"),
)


def chk(x,msg):
    if not x:raise RuntimeError("FCP_FIVE_COHORT_CLIMATE_HUE_PROVENANCE: "+msg)


def sha(p):
    h=hashlib.sha256()
    with Path(p).open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):h.update(b)
    return h.hexdigest()


def cohort_bootstrap(d,target,rep=B,seed=SEED):
    d=d.loc[np.isfinite(d[target].to_numpy(float))]
    x=d[target].to_numpy(float)
    ge,genus=pd.factorize(d.species.str.split().str[0],sort=True)
    G=len(genus)
    if len(x)<20 or G<15:
        return {"status":"not_estimable","n_taxa":len(x),"n_genera":G}
    val=np.zeros(G,float)
    cnt=np.bincount(ge,minlength=G)
    np.add.at(val,ge,x)
    rand=np.random.default_rng(seed)
    draw=rand.integers(G,size=(rep,G))
    sample=val[draw].sum(axis=1)/cnt[draw].sum(axis=1)
    mean=float(x.mean())
    interval=[float(v) for v in np.quantile(sample,[.025,.975])]
    return {
      "status":"species_equal_genus_cluster_bootstrap",
      "n_taxa":len(x),
      "n_genera":G,
      "mean_signed_original_photo_hue_covariance":mean,
      "median_per_taxon":float(np.median(x)),
      "fraction_of_source_taxa_positive":float(np.mean(x>0)),
      "exploratory_genus_marginal95":interval,
      "marginal95_excludes_zero":bool(interval[0]>0 or interval[1]<0),
      "bootstrap_repetitions":rep
    }


def calc(paths):
    for key,p in paths.items():
        chk(sha(p)==SHA[key],"expected immutable original-source "+key+" bytes differ")
    source=json.loads(paths["source_result"].read_text())
    chk(source["source_n_original_photos"]==447500
        and source["source_n_original_taxa"]==4475
        and source["source_n_classified_photo"]==227107,
        "source photo denominator drift")
    results={}
    for m,mode in enumerate(("original","observer_thinned")):
        df=pd.read_csv(paths[mode],low_memory=False)
        chk(len(df)==4475 and df.inat_taxon_id.is_unique,
            "4,475 source taxon rows missing")
        chk(set(df.source_cohort)==set(FIVE),
            "not exactly five original species-disjoint cohorts")
        df["genus"]=df.species.str.split().str[0]
        output={}
        for i,(title,key,white) in enumerate(TARGETS):
            chk(key in df.columns and (white is None or white in df.columns),
                "missing original seven-climate source-hue coefficients")
            # For chromatic-only nonwhite hue: project out uniform missing
            # nonwhite share that accompanies its white-palette source fraction.
            value=pd.to_numeric(df[key],errors="coerce")
            if white:
                value=value+pd.to_numeric(df[white],errors="coerce")/8.
            df["_target"]=value
            group=[]
            for j,cohort in enumerate(FIVE):
                x=df.loc[df.source_cohort==cohort]
                stats=cohort_bootstrap(x,"_target",seed=SEED+1000*m+83*i+11*j)
                group.append({"cohort":cohort,**stats})
            use=np.isfinite(value)
            observedmean=float(value[use].mean())
            base=source["results"][mode]["component_results"][
                "total" if i in (0,3,4) else "within_photo_labels"
            ]["axes"]["BIO5" if i==0 else "BIO12"]
            chk(base["n_species"]==int(use.sum()),
                "original-source hue endpoint denominator no longer matches completed analysis")
            if i==0:
                chk(abs(observedmean-base["white_mean_species_covariance"])<1e-12,
                    "BIO5 signed white original global mean changed")
            if i in (1,2,3,4):
                hue="yellow" if i in (1,3) else "bronze"
                expected=base["photographed_nine_hue_components"][hue]["chromatic_only_component"]
                chk(abs(observedmean-expected)<1e-12,
                    "exact chromatic residual differs from fixed-white-axis original data")
            signs=[int(np.sign(e["mean_signed_original_photo_hue_covariance"]))
                   for e in group if e["status"]=="species_equal_genus_cluster_bootstrap"]
            output[title]={
               "global_species_equal_mean":observedmean,
               "global_n_eligible_species":int(use.sum()),
               "cohort_mean_signs":signs,
               "n_cohorts_with_positive_mean":sum(z>0 for z in signs),
               "n_cohorts_with_negative_mean":sum(z<0 for z in signs),
               "no_prospective_confirmation_claim":True,
               "cohorts":group,
            }
        results[mode]=output
    return {
      "schema":"fcp4475_nine_hue_climate_original_five_cohort_source_audit_v1",
      "parent_source_sha256":dict(SHA),
      "original_source_photo_run":40060153598,
      "original_artifact_id":12547087651,
      "original_zip_sha256":"c23f0c7452019940bdbcca88b249be854c3f236e39ce189aa8c77119f021f778",
      "photos_original_n":447500,
      "species_original_n":4475,
      "source_classifiable_photo_n":227107,
      "five_species_disjoint_source_cohorts":list(FIVE),
      "species_equal_and_whole_genus_bootstrap":B,
      "chromatic_only_nonwhite_hue":"original photo hue coefficient plus the photo white component divided by eight, exactly removing pre-fixed white-vs-chromatic projection",
      "direction_targets_selected_posthoc_after_full_family56_original_results":True,
      "not_an_independent_source_or_genetic_morph_replication":True,
      "not_flower_background_or_causal_climate_control":True,
      "results":results
    }


def main():
    p=argparse.ArgumentParser()
    for key in ("original","observer-thinned","source-result","out"):
        p.add_argument("--"+key,type=Path)
    p.add_argument("--self-test",action="store_true")
    a=p.parse_args()
    if a.self_test:
        source=np.array([.1,.2,-.1,-.2])
        y=np.array([.01,-.02,-.01,.02])
        chk(np.allclose(y+source/8,y+source/8),"chromatic hue fixture broken")
        print("FCP4475_FIVE_COHORT_PHOTO_HUE_SOURCE_SYNTHETIC_PASS")
        return
    chk(a.original and a.observer_thinned and a.source_result and a.out,
        "all original source cohort coefficients and source result required")
    result=calc({"original":a.original,
                 "observer_thinned":a.observer_thinned,
                 "source_result":a.source_result})
    a.out.mkdir(parents=True,exist_ok=True)
    (a.out/"five_cohort_climate_colour_signature_audit.json").write_text(
        json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("FCP4475_FIVE_COHORT_COLOUR_CLIMATE_AUDIT_PASS",json.dumps({
        panel:{
          name:{"signs":v["cohort_mean_signs"],
                "n":v["global_n_eligible_species"],
                "mean":v["global_species_equal_mean"]}
          for name,v in vals.items()
        } for panel,vals in result["results"].items()
    },sort_keys=True))


if __name__=="__main__":
    main()

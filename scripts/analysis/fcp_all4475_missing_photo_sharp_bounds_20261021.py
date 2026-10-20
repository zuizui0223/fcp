#!/usr/bin/env python3
"""Sharp, assumption-explicit missing-colour completion bounds for ALL 4,475 FCP taxa.

The immutable 447,500 original photos contain 227,107 with an original,
classified four-state photo label. For EACH of the original 100 observations
per species, hypothesize that every remaining unclassified observation could
be assigned ANY of the four original labels.

Compute sharp, attainable bounds on the TOTAL-100-photo unbiased Simpson pair
mismatch and 100-photo white proportion; all 4,475 taxa remain in the output.
These are purely algebraic sensitivity limits to missing photographed labels,
NOT claims that unclassified images depict comparable flowers, and NOT bounds
on genetically inherited natural polymorphism or selection/fitness.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd

SHA = "8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6"
N_PHOTO, N_TAXA, N_CLASSIFIABLE = 447500, 4475, 227107
STATES = ("white","yellow_orange","red_pink","blue_purple")
SEED = 20261021

def check(ok,message):
    if not ok:raise RuntimeError("FCP_ALL4475_MISSING_PHOTO_BOUNDS_INVALID: "+message)

def source_hash(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda:stream.read(1<<20),b""): h.update(chunk)
    return h.hexdigest()

def mismatch(counts):
    x=np.asarray(counts,dtype=np.int64)
    n=x.sum()
    if n<2:return None
    return float(1-(x*(x-1)).sum()/(n*(n-1)))

def sharp_completion(counts,missing):
    """Attain min/max unbiased pair-mismatch for 4-state assignments.

    Min: assign every missing colour to one of the largest existing states.
    Max: add each missing photo to a currently smallest colour state.
    This minimises/maximises sum(n_k^2) at fixed sum by convex exchange.
    """
    original=np.asarray(counts,dtype=np.int64)
    check(original.shape==(4,) and (original>=0).all() and missing>=0,
          "four nonnegative state counts and missing >=0 required")
    low=original.copy()
    low[int(np.argmax(low))]+=int(missing)
    high=original.copy()
    for _ in range(int(missing)):
        high[int(np.argmin(high))]+=1
    check(low.sum()==high.sum()==original.sum()+missing,"completion lost photo rows")
    a=mismatch(low);b=mismatch(high)
    check(a is not None and b is not None and a<=b+1e-14,
          "sharp lower bound exceeds upper bound")
    return a,b,low,high

def selftest():
    # All possible allocations of m=0..5 missing across 4 states.
    for obs in ([2,0,0,0],[1,1,0,0],[0,0,0,2],[1,1,1,1],[3,2,1,0]):
        for missing in range(6):
            lo,hi,_,_=sharp_completion(obs,missing)
            candidates=[]
            for i in range(missing+1):
                for j in range(missing-i+1):
                    for k in range(missing-i-j+1):
                        l=missing-i-j-k
                        candidates.append(mismatch(np.array(obs)+[i,j,k,l]))
            check(abs(lo-min(candidates))<1e-12 and abs(hi-max(candidates))<1e-12,
                  "sharp completion does not match exhaustive enumerated assignments")
    check(abs(mismatch([2,2,0,0])-2/3)<1e-12,"unbiased count mismatch arithmetic")
    print("FCP_ALL4475_MISSINGNESS_SHARP_EXTREMA_ENUMERATION_PASS")

def bootstrap_by_genus(table,columns,iterations=999):
    df=table.copy()
    df["genus"]=df.species.str.split().str[0]
    g=df.groupby("genus",sort=True)[list(columns)].agg(["sum","count"])
    groups=g.index.to_list()
    n=len(groups)
    check(n>100,"too few genera for hierarchical sensitivity")
    rng=np.random.default_rng(SEED)
    i=rng.integers(0,n,size=(iterations,n))
    output={}
    for col in columns:
        totals=g[col]["sum"].to_numpy(float)
        sizes=g[col]["count"].to_numpy(int)
        result=totals[i].sum(axis=1)/sizes[i].sum(axis=1)
        output[col]=[float(a) for a in np.quantile(result,[.025,.975])]
    return output

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--photos",type=Path)
    p.add_argument("--out",type=Path)
    p.add_argument("--self-test",action="store_true")
    z=p.parse_args()
    selftest()
    if z.self_test:return
    check(z.photos is not None and z.out is not None,"actual photo file/output required")
    check(source_hash(z.photos)==SHA,"authentic original photo CSV SHA mismatch")
    cols=["inat_taxon_id","species","photo_id","fcp_source_cohort","global_classifiable","morph"]
    source=pd.read_csv(z.photos,usecols=cols,low_memory=False)
    check(len(source)==N_PHOTO and source.inat_taxon_id.nunique()==N_TAXA,
          "complete 4475x100 original source denominator invalid")
    check(source.photo_id.is_unique,"duplicate photo ID in original source")
    photo_per_taxon=source.groupby("inat_taxon_id").size()
    check(photo_per_taxon.eq(100).all(),"source taxon lacks 100 original photos")
    check(source.groupby("inat_taxon_id").fcp_source_cohort.nunique().eq(1).all(),
          "taxon belongs to two sourcing cohorts")
    flag=(source.global_classifiable.fillna(False).astype(str).str.lower()
          .isin(["true","yes","1"]))
    check(int(flag.sum())==N_CLASSIFIABLE,"original classified total changed")
    classified=source.loc[flag]
    check(classified.morph.isin(STATES).all(),"unexpected classified state")
    coarse=(classified.groupby(["inat_taxon_id","morph"],observed=True)
            .size().unstack(fill_value=0).reindex(columns=STATES,fill_value=0))
    header=(source.groupby("inat_taxon_id",sort=True)
            .agg(species=("species","first"),cohort=("fcp_source_cohort","first")))
    check(len(coarse)==len(header),"at least one taxon has no assigned photograph")
    table=header.join(coarse,how="left").fillna({k:0 for k in STATES})
    out=[]
    for tid,rec in table.iterrows():
        counts=rec[list(STATES)].to_numpy(dtype=int)
        classified_n=int(counts.sum());miss=100-classified_n
        lo,hi,_,_=sharp_completion(counts,miss)
        D=mismatch(counts)
        w=int(counts[0])
        out.append({
            "inat_taxon_id":int(tid),
            "species":str(rec.species),"source_cohort":str(rec.cohort),
            "original_photos":100,
            "classified_four_state_photos":classified_n,
            "missing_original_photo_colour_labels":miss,
            "coarse_photo_white":w,
            "coarse_photo_yellow_orange":int(counts[1]),
            "coarse_photo_red_pink":int(counts[2]),
            "coarse_photo_blue_purple":int(counts[3]),
            "D_given_classifiable_unbiased":D,
            "D_if_all_100_four_state_lower_sharp":lo,
            "D_if_all_100_four_state_upper_sharp":hi,
            "D_if_all_100_four_state_sharp_width":hi-lo,
            "white_fraction_if_all_100_lower":w/100.,
            "white_fraction_if_all_100_upper":(w+miss)/100.,
            "both_photo_states_directly_recorded":int(np.count_nonzero(counts)>1),
            "n_photo_colour_labels_per_species_below40":bool(classified_n<40),
            "all_four_state_completions_are_hypothetical_not_observed":True,
        })
    df=pd.DataFrame(out).sort_values("inat_taxon_id")
    check(len(df)==N_TAXA and df.inat_taxon_id.nunique()==N_TAXA,"4,475 taxa not retained")
    check(df.classified_four_state_photos.sum()==N_CLASSIFIABLE,
          "classified source photo count lost")
    check(np.all((df.D_given_classifiable_unbiased>=0)&(df.D_given_classifiable_unbiased<=1)),
          "D original classified out of bounds")
    bins=pd.cut(df.classified_four_state_photos,[0,19,39,59,79,100],
                labels=["8–19","20–39","40–59","60–79","80–100"],include_lowest=True)
    df["classifiable_photo_bin"]=bins.astype(str)
    summary={
        "schema":"fcp_all4475_adversarial_photo_colour_missingness_sharp_bounds_v1",
        "source_sha256":SHA,
        "n_taxa":N_TAXA,
        "n_original_photos":N_PHOTO,
        "n_observed_classifiable_four_state":N_CLASSIFIABLE,
        "n_missing_four_state_photo_labels":N_PHOTO-N_CLASSIFIABLE,
        "mean_D_conditional_on_classified":float(df.D_given_classifiable_unbiased.mean()),
        "mean_D_hypothetical_all100_sharp_lower":float(df.D_if_all_100_four_state_lower_sharp.mean()),
        "mean_D_hypothetical_all100_sharp_upper":float(df.D_if_all_100_four_state_upper_sharp.mean()),
        "mean_D_hypothetical_all100_bound_width":float(df.D_if_all_100_four_state_sharp_width.mean()),
        "median_species_D_hypothetical_lower":float(df.D_if_all_100_four_state_lower_sharp.median()),
        "median_species_D_hypothetical_upper":float(df.D_if_all_100_four_state_upper_sharp.median()),
        "n_original_species_with_at_least_two_coarse_colours":int((df.both_photo_states_directly_recorded).sum()),
        "n_species_with_D_all100_lower_above_0p1":int((df.D_if_all_100_four_state_lower_sharp>0.1).sum()),
        "n_species_with_D_all100_upper_below_0p1":int((df.D_if_all_100_four_state_upper_sharp<0.1).sum()),
        "cohort_n_taxa":df.groupby("source_cohort").size().to_dict(),
        "strata_by_classifiable_n":[{"photo_bin":str(k),"n_taxa":len(g),
              "mean_classified_D":float(g.D_given_classifiable_unbiased.mean()),
              "mean_all100_lower":float(g.D_if_all_100_four_state_lower_sharp.mean()),
              "mean_all100_upper":float(g.D_if_all_100_four_state_upper_sharp.mean()),
              "mean_missing_photo_labels":float(g.missing_original_photo_colour_labels.mean())}
             for k,g in df.groupby("classifiable_photo_bin",sort=True)],
        "genus_cluster_bootstrap_95CI":bootstrap_by_genus(df,[
          "D_given_classifiable_unbiased",
          "D_if_all_100_four_state_lower_sharp",
          "D_if_all_100_four_state_upper_sharp"
        ]),
        "bounds_assumption":"Every unclassified original source photo hypothetically has one of the same four original biological colour labels. If this assumption fails (flower organ absent, image fails ROI QC), these are NOT biological true-colour bounds.",
        "note":"Sharp counts are algebraic for an entirely hypothetical 100-photo labelling completion, without missing-at-random assumptions. They are NOT bounds on wild genetic FCP or true pigment cost.",
        "taxa_eliminated_or_missing_records_imputed":0,
        "no_unclassified_photo_label_reconstruction":True,
        "no_native_wild_population_genotype_evidence":True,
        "research_role":"post hoc whole-4475 source classification sensitivity and limitation of global phenotypic diversity interpretations",
    }
    z.out.mkdir(parents=True,exist_ok=True)
    df.to_csv(z.out/"all4475_source_label_missingness_sharp_bounds_no_coordinates.csv",
              index=False,float_format="%.12g")
    (z.out/"all4475_source_label_missingness_summary.json").write_text(
        json.dumps(summary,indent=2,sort_keys=True,allow_nan=False,ensure_ascii=False)+"\n")
    print("FCP_ALL4475_AUTHENTIC_CLASSIFIABILITY_BOUNDS_COMPLETE",
          json.dumps(summary,sort_keys=True,allow_nan=False,ensure_ascii=False))

if __name__=="__main__":
    main()

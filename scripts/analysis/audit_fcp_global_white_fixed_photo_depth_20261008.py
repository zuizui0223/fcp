#!/usr/bin/env python3
"""Standardize visible white/nonwhite detection at equal 5-photo depth.

This is retrospective photo-label support, not genetic white-morph frequency.
Input is the exact SHA-verified species-by-absolute-latitude geographic
source summary from completed 42,111/85,337 photo measurements.

For each sampled species with >=K classifiable taxon-cell anchors:
P(both states in K draws without replacement) =
1 - choose(n_white,K)/choose(n,K)
  - choose(n_nonwhite,K)/choose(n,K).

Also report a finite-photo multinomial/Bernoulli within-species variance
sensitivity; it requires exchangeable photo labels and is NOT an
evolutionary variance-component estimate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

SOURCE_SUMMARY_SHA="52d732b9bbd7e9c2292232e1328b48dca175dd34ca71e901d04b9e3d9fcc7fee"
REGIONS=("ALL_GLOBAL","low_0_30","middle_30_60","high_60_90")
MIN_DEPTHS=(5,10)
STANDARDIZED_DRAW=5
BOOTSTRAPS=999
SEED=20261008155


def file_sha(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as handle:
        for c in iter(lambda:handle.read(1<<20),b""):
            h.update(c)
    return h.hexdigest()


def all_same_prob(n:np.ndarray,count:int|np.ndarray,k:int)->np.ndarray:
    n=np.asarray(n,dtype=np.int64)
    count=np.asarray(count,dtype=np.int64)
    if n.shape!=count.shape or np.any(n<k) or np.any(count<0) or np.any(count>n):
        raise ValueError("Invalid per-species photo colour counts for exact finite sampling")
    p=np.ones(len(n),float)
    valid=count>=k
    p[~valid]=0.0
    for i in range(k):
        if valid.any():
            p[valid]*=(count[valid]-i)/(n[valid]-i)
    return p


def mixed_prob_fixed_sample(n:np.ndarray,w:np.ndarray,k:int)->np.ndarray:
    n=np.asarray(n,dtype=np.int64)
    w=np.asarray(w,dtype=np.int64)
    return np.clip(1-all_same_prob(n,w,k)-all_same_prob(n,n-w,k),0.0,1.0)


def bootstrap_ci(values:np.ndarray,seed:int)->list[float]|None:
    x=np.asarray(values,float)
    if len(x)<2:return None
    rng=np.random.default_rng(seed)
    replicates=np.empty(BOOTSTRAPS,float)
    for b in range(BOOTSTRAPS):
        replicates[b]=float(x[rng.integers(0,len(x),len(x))].mean())
    return [float(v) for v in np.quantile(replicates,[.025,.975])]


def finite_photo_partition(n:np.ndarray,w:np.ndarray)->dict:
    n=np.asarray(n,float)
    w=np.asarray(w,float)
    if not len(n):
        return {"estimable":False}
    if np.any(n<2) or np.any(w<0) or np.any(w>n):
        raise ValueError("Finite photo variance estimator requires >=2 known photos per species")
    p=w/n
    global_white=float(w.sum()/n.sum())
    total=global_white*(1-global_white)
    weighted_within_observed=float(np.average(p*(1-p),weights=n))
    weighted_within_corrected=float(np.average((n/(n-1))*p*(1-p),weights=n))
    equal_white=float(p.mean())
    equal_total=equal_white*(1-equal_white)
    eq_observed=float(np.mean(p*(1-p)))
    eq_corrected=float(np.mean((n/(n-1))*p*(1-p)))
    return {
        "estimable":True,
        "n_species":int(len(n)),
        "mean_white_fraction_photo_weighted":global_white,
        "mean_white_fraction_species_equal":equal_white,
        "uncorrected_between_species_share_photo_weighted":float((total-weighted_within_observed)/total) if total else None,
        "finite_photo_corrected_between_species_share_photo_weighted":float((total-weighted_within_corrected)/total) if total else None,
        "uncorrected_between_species_share_species_equal":float((equal_total-eq_observed)/equal_total) if equal_total else None,
        "finite_photo_corrected_between_species_share_species_equal":float((equal_total-eq_corrected)/equal_total) if equal_total else None,
        "assumption":"finite-photo sensitivity assumes species-conditional photo observations approximately exchangeable independent labels; violates if locality/genotype/photo-error clustered",
        "not_fitness_population_variance":True
    }


def analyze(species_by_region:pd.DataFrame)->dict:
    required={"region","inat_taxon_id","n_classified","n_white","observed_cells"}
    if not required.issubset(species_by_region):
        raise ValueError("Original per-species photo opportunity columns missing")
    if species_by_region[["region","inat_taxon_id"]].duplicated().any():
        raise ValueError("Photo unit or species by region counted twice")
    if not set(species_by_region.region).issubset(REGIONS):
        raise ValueError("Unknown geographic scope")
    data=species_by_region.copy()
    for col in ("n_classified","n_white","observed_cells"):
        data[col]=pd.to_numeric(data[col],errors="raise").astype(int)
    if (data.n_white>data.n_classified).any() or (data.n_classified>data.observed_cells).any():
        raise ValueError("Colour classification counts exceed source photo opportunity")
    out=[]
    for ri,region in enumerate(REGIONS):
        part=data.loc[data.region==region]
        for threshold in MIN_DEPTHS:
            adequate=part.loc[part.n_classified>=threshold]
            n=adequate.n_classified.to_numpy(dtype=int)
            w=adequate.n_white.to_numpy(dtype=int)
            if not len(n):
                out.append({
                    "region":region,"minimum_classifiable_cell_photos":threshold,
                    "n_species":0,"status":"HOLD_NO_REPEAT_PHOTO_OPPORTUNITY"
                });continue
            present=(w>0)&(w<n)
            fixed=mixed_prob_fixed_sample(n,w,STANDARDIZED_DRAW)
            out.append({
                "region":region,
                "minimum_classifiable_cell_photos":threshold,
                "n_species":int(len(n)),
                "n_classifiable_original_taxon_cell_photos":int(n.sum()),
                "observed_fraction_species_with_both_colours_at_all_available_depth":float(present.mean()),
                "expected_fraction_species_showing_both_colours_in_exactly_5_of_their_existing_photos":float(fixed.mean()),
                "photo_standardized_expected_mixed_fraction_species_bootstrap_95CI":bootstrap_ci(fixed,SEED+ri*20+threshold),
                "observed_vs_standardized_mix_detection_difference":float(present.mean()-fixed.mean()),
                "exact_number_with_both_observed_states":int(present.sum()),
                "fraction_species_with_all_observed_labels_white":float(np.mean(w==n)),
                "fraction_species_with_all_observed_labels_nonwhite":float(np.mean(w==0)),
                "photo_variance_partition_finite_depth_sensitivity":finite_photo_partition(n,w),
                "status":"EXPLORATORY_FIXED_OBSERVED_PHOTO_LABELS_DEPTH_DIAGNOSTIC",
                "inference_boundary":"5-photo hypergeometric draws sample only already observed photo labels from these same taxa and cells; not hypothetical natural allele frequency or an independently evaluated population sample",
            })
    return {
        "schema":"fcp_all42111_white_50_50_fixed_photo_depth_sensitivity_v1",
        "date_jst":"2026-10-08",
        "status":"RETROSPECTIVE_EQUAL_FIVE_PHOTO_DETECTION_OPPORTUNITY_DIAGNOSTIC",
        "source_species_region_table_sha256":SOURCE_SUMMARY_SHA,
        "source_global_species_frame":42111,
        "photographic_measurement_frame":85337,
        "standardized_photo_depth":STANDARDIZED_DRAW,
        "results":out,
        "hard_nonclaims":[
            "five-photo standardized mixture probability measures subsampling of ALREADY OBSERVED visible photo colours, not probability of genetic FCP",
            "no extra original images were accessed, species photo IDs were not changed",
            "all-zero or all-white photo categories are not proven genetic monomorphisms",
            "finite photo correction assumes approximate independent conditional photographs; geographical clustering and image colour error can violate it",
            "sampling-depth thresholds select different subsets of plant species; differences are not direct longitudinal changes in the same samples",
            "the global 50:50 white/coloured photo share does not imply balanced lifetime reproduction or pollen/herbivore cost equilibrium"
        ],
        "confirmatory_decisions_changed":False
    }


def main()->None:
    p=argparse.ArgumentParser()
    p.add_argument("--original-species-regions",type=Path,required=True)
    p.add_argument("--outdir",type=Path,required=True)
    args=p.parse_args()
    if file_sha(args.original_species_regions)!=SOURCE_SUMMARY_SHA:
        raise RuntimeError("Historical 42k world species/regional observed colour counts have changed")
    tab=pd.read_csv(args.original_species_regions,low_memory=False)
    report=analyze(tab)
    args.outdir.mkdir(parents=True,exist_ok=True)
    (args.outdir/"result.json").write_text(
        json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False))


if __name__=="__main__":
    main()

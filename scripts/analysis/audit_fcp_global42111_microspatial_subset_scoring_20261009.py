#!/usr/bin/env python3
"""Outcome-blind locally exchangeable PHOTO subset scoring for FCP global null.

This post-hoc sensitivity evaluates original 250/500km congeneric photo models
ONLY on original photo IDs belonging to >=2-species geographic microgroups
at <=50/100km diameter, regardless of colour class. Models are still fitted
on all original 250/500km source photographs, as in the frozen main assay.
Only evaluation weights differ; fixed taxon-ID CV and photo-colour count-
preserving spatial microgroup null refitting remain unchanged.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from audit_fcp_global42111_local_congener_moisture_null_20261009 import (
    fixed_source_cohort, EXPECTED_ORIGINAL_N, EXPECTED_OBS_GAIN, SEED,
)
from audit_fcp_global42111_local_congeners_20261009 import (
    photo_populations, fold_plan, train_only_group_predict, GEO, TEMP, RAIN, CLASSES,
)
from audit_fcp_global42111_microspatial_photo_shuffle_20261009 import microgroups

GROUP_CAPS=(250,500)
MICRO_CAPS=(50,100)
PERMUTATIONS=199
SOURCE_SCHEMA="fcp_global42111_original_nearby_exchangeable_photo_subset_null_v1"
MIN_EVAL_COUNT=100
MIN_EVAL_FRACTION=.1


def original_fixed_fold_photo_gain(source:pd.DataFrame,fold:np.ndarray)->np.ndarray:
    """Exactly original Brier gain per held-out photo, before any scoring mask."""
    y=pd.Categorical(source.morph,categories=CLASSES).codes
    if (y<0).any():raise ValueError("Unclassified photo colour cannot enter frozen cohort")
    truth=np.eye(4)[y]
    contribution=np.full(len(source),np.nan,float)
    for k in range(5):
        tr=np.flatnonzero(fold!=k)
        te=np.flatnonzero(fold==k)
        if not len(te):continue
        train=source.iloc[tr];test=source.iloc[te]
        if not test.genus_cell_id.isin(train.genus_cell_id.value_counts().loc[lambda x:x>=2].index).all():
            raise RuntimeError("Local original photo genus-cell training supports fewer than two source species")
        basic=train_only_group_predict(train,test,GEO+TEMP)
        plus=train_only_group_predict(train,test,GEO+TEMP+RAIN)
        contribution[te]=(np.square(basic-truth[te]).sum(axis=1)-
                          np.square(plus-truth[te]).sum(axis=1))
    if not np.isfinite(contribution).all():
        raise RuntimeError("FCP source photo heldout prediction unassigned")
    return contribution


def evaluate_subset(source:pd.DataFrame, micro_distance:int,nperm:int,seed:int)->dict:
    groups=microgroups(source,float(micro_distance))
    # Primary evaluation eligibility uses ONLY original photographed-site
    # geometry and species counts, NEVER the source colour outcome.
    eligible=np.unique(np.concatenate([g for g in groups if len(g)>=2])) if any(len(g)>=2 for g in groups) else np.array([],int)
    n=len(eligible);ratio=n/len(source)
    base=source.morph.to_numpy(str)
    stat={
        "microgroup_photo_diameter_km":micro_distance,
        "source_original_photo_species":len(source),
        "n_original_evaluation_photo_ids_in_multispecies_microgroups":n,
        "source_species_fraction_scored":ratio,
        "n_original_geographical_microgroups":len(groups),
        "n_original_multispecies_microgroups":sum(len(g)>=2 for g in groups),
        "test_species_selection_uses_only_original_photo_coordinates_and_species_ID":True,
        "no_photo_colour_class_used_to_select_scoring_subset":True,
        "same_original_complete_cohort_fit_for_all_subset_evaluations":True,
        "microgroup_four_colour_composition_preserved_per_null":True,
        "n_fixed_null_model_refits":0,
        "source_label_shuffles":nperm,
    }
    if n<MIN_EVAL_COUNT or ratio<MIN_EVAL_FRACTION:
        stat.update({"status":"HOLD_INSUFFICIENT_COLOUR_BLIND_MULTISPECIES_PHOTO_NEIGHBORHOODS",
                     "permutation_p":None})
        return stat
    fold=fold_plan(source)
    scores=original_fixed_fold_photo_gain(source,fold)
    obs=float(scores[eligible].mean())
    rng=np.random.default_rng(seed)
    null=[]
    for _ in range(nperm):
        shuffled=base.copy()
        for micro in groups:
            if len(micro)>=2:
                shuffled[micro]=rng.permutation(base[micro])
        variant=source.copy()
        variant["morph"]=shuffled
        score=original_fixed_fold_photo_gain(variant,fold)
        null.append(float(score[eligible].mean()))
    arr=np.asarray(null,float)
    stat.update({
        "status":"SOURCE_COLOR_BLIND_MICRONEIGHBORHOOD_SCORING_AND_NULL",
        "n_fixed_null_model_refits":nperm,
        "original_heldout_brier_gain_all_source_species":float(scores.mean()),
        "original_heldout_brier_gain_geographic_scoring_subset":obs,
        "original_photo_subset_null_mean":float(arr.mean()),
        "original_photo_subset_null_sd":float(arr.std(ddof=1)) if nperm>1 else None,
        "original_photo_subset_null_quantiles_2p5_50_97p5":[float(v) for v in np.quantile(arr,[.025,.5,.975])],
        "null_exceedances":int(np.sum(arr>=obs-1e-12)),
        "permutation_p":float((1+np.sum(arr>=obs-1e-12))/(nperm+1)),
    })
    return stat


def run(source:pd.DataFrame,*,strict:bool=True,nperm:int=PERMUTATIONS)->dict:
    pops,coverage=photo_populations(source,strict=strict)
    d=pops["CLIMATE_ALL"]
    allout={}
    for cap in GROUP_CAPS:
        x=fixed_source_cohort(d,cap)
        if strict and len(x)!=EXPECTED_ORIGINAL_N[str(cap)]:
            raise RuntimeError("Immutable source-photo taxon sample changed")
        full=original_fixed_fold_photo_gain(x,fold_plan(x)).mean()
        if strict and abs(full-EXPECTED_OBS_GAIN[str(cap)])>1e-8:
            raise RuntimeError("Independent source all-photo previous published Brier statistic drift")
        inner={}
        for micro in MICRO_CAPS:
            inner[str(micro)]=evaluate_subset(x,micro,nperm,SEED+cap*100+micro)
            item=inner[str(micro)]
            if item["status"]=="SOURCE_COLOR_BLIND_MICRONEIGHBORHOOD_SCORING_AND_NULL":
                if abs(item["original_heldout_brier_gain_all_source_species"]-full)>1e-12:
                    raise RuntimeError("Source model fitted on different photos")
        allout[str(cap)]={"source_species":len(x),
                          "full_cohort_original_moisture_gain":float(full),
                          "outcome_blind_microgroup_scoring":inner}
    return {
        "schema":SOURCE_SCHEMA,
        "date_jst":"2026-10-09",
        "status":"POSTHOC_SOURCE_PHOTO_SPATIAL_EXCHANGEABLE_SUBSET_CALIBRATION",
        "original_source_species":len(source),
        "original_colour_classifiable_source_species":coverage["source_classified_photo_taxa"],
        "original_climate_complete_photo_species":coverage["climate_eligible_classified_source_taxa"],
        "original_distance_caps_km":list(GROUP_CAPS),
        "original_nearby_complete_link_caps_km":list(MICRO_CAPS),
        "source_four_state_labels_are_original_and_never_remeasured":True,
        "colour_blind_microgroup_eval_eligibility_fixed_before_permutation":True,
        "same_observed_model_as_previously_frozen_250_500_cohort":True,
        "source_original_cohorts":allout,
        "nonclaims":[
            "This is a new post-hoc photo subset estimand, not the original full-cohort p-value",
            "A 50-100km genus-cell photo neighbourhood is not a common population or genetic family tree",
            "The source photo group data are selected on geography, observed genus and original colour classifiability",
            "A rare source conditional label shuffle does not establish climatic selection or evolutionary adaptation",
            "The two cohort thresholds and microgroups overlap; none independently validates the other",
        ],
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-breadth-abiotic",required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    args=p.parse_args()
    source=pd.read_csv(args.original_breadth_abiotic,low_memory=False)
    result=run(source)
    args.outdir.mkdir(parents=True,exist_ok=True)
    (args.outdir/"result.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,indent=2,sort_keys=True),flush=True)


if __name__=="__main__":
    main()

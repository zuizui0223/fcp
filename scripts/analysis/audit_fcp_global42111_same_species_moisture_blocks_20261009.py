#!/usr/bin/env python3
"""Ablate temperature vs precipitation in ORIGINAL same-species FCP photo pairs.

Exact source 13,416 observer-disjoint pairs, 3,196 doubly photographed
four-colour outcomes, ~3,182 actual two-site climate-complete source pairs.
NO genetic, individual plant, pollen/pollinator or adaptive fitness inference.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from analyze_fcp_global42111_within_species_pair_abiotic_20261009 import (
    N_PAIRS,BOTH_CLASSIFIABLE,DISCORDANT,N_FOLDS,N_BOOT,SEED,
    BASE,read_original_pairs,join_origins
)

THERMAL=("delta_wc_bio1","delta_wc_bio5")
MOISTURE=("delta_wc_bio12","delta_wc_bio15")
MODELS={
    "GEOGRAPHY":BASE,
    "GEOGRAPHY_TEMPERATURE":BASE+THERMAL,
    "GEOGRAPHY_MOISTURE":BASE+MOISTURE,
    "GEOGRAPHY_ALL_CLIMATE":BASE+THERMAL+MOISTURE,
}
MIN_PAIRS=300
ASSERT_SOURCE_CLIMATE_PAIRS=3182
ASSERT_SOURCE_CLIMATE_DISCORDANT=793


def model_group(d:pd.DataFrame,group_name:str)->dict:
    from sklearn.model_selection import GroupKFold
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score,log_loss

    if group_name not in ("genus","pair_midpoint_cell_162"):
        raise ValueError("Use only historical source genus or geographic midpoint group")
    if d.inat_taxon_id.duplicated().any():
        raise ValueError("One source photographed pair per species required")
    y=d.pair_state.eq("discordant").to_numpy(int)
    grouping=d[group_name].to_numpy()
    if len(d)<MIN_PAIRS or len(set(grouping))<N_FOLDS or y.sum()<20:
        return {"status":"HOLD_PAIR_GEOGRAPHIC_OR_TAXON_SUPPORT"}
    folds=list(GroupKFold(n_splits=N_FOLDS).split(d,y,grouping))
    if any(len(np.unique(y[tr]))!=2 or int(y[tr].sum())<20 for tr,te in folds):
        return {"status":"HOLD_SOURCE_PAIR_TRAINING_CLASS_COVERAGE"}
    score={}
    loss={}
    for name,cols in MODELS.items():
        x=d[list(cols)].to_numpy(float)
        if not np.isfinite(x).all():raise ValueError("Source climate matching contains missing values")
        probs=np.full(len(d),np.nan)
        for tr,te in folds:
            model=make_pipeline(StandardScaler(),LogisticRegression(
                max_iter=750,C=1.0,random_state=SEED))
            model.fit(x[tr],y[tr])
            probs[te]=np.clip(model.predict_proba(x[te])[:,1],1e-6,1-1e-6)
        if np.isnan(probs).any():raise RuntimeError("Source climate model failed geographic holdout")
        loss[name]=-(y*np.log(probs)+(1-y)*np.log1p(-probs))
        score[name]={
            "source_pairs_evaluated":len(d),
            "heldout_log_loss":float(log_loss(y,probs,labels=[0,1])),
            "heldout_AUC":float(roc_auc_score(y,probs)),
        }
    comp={
        "total_climate_beyond_geography":("GEOGRAPHY","GEOGRAPHY_ALL_CLIMATE"),
        "unique_temperature_beyond_geography_moisture":("GEOGRAPHY_MOISTURE","GEOGRAPHY_ALL_CLIMATE"),
        "unique_moisture_beyond_geography_temperature":("GEOGRAPHY_TEMPERATURE","GEOGRAPHY_ALL_CLIMATE"),
    }
    unique,gi=np.unique(grouping,return_inverse=True)
    n=np.bincount(gi)
    increments={}
    for name,(base,full) in comp.items():
        delta=loss[base]-loss[full]
        group_sums=np.bincount(gi,weights=delta)
        rng=np.random.default_rng(SEED+len(name))
        sel=rng.integers(0,len(unique),size=(N_BOOT,len(unique)))
        boot=group_sums[sel].sum(axis=1)/n[sel].sum(axis=1)
        increments[name]={
            "heldout_log_loss_reduction":float(delta.mean()),
            "original_genus_or_cell_block_bootstrap_95CI":[float(v) for v in np.quantile(boot,[.025,.975])],
            "positive_increment_supported_by_conditional_bootstrap":bool(np.quantile(boot,.025)>0),
        }
    return {
        "status":"ORIGINAL_SAME_SPECIES_PHOTOGRAPHED_PAIR_CLIMATE_BLOCK_MODEL",
        "heldout_group":group_name,
        "n_original_source_species_pairs":len(d),
        "n_photo_colour_discordant":int(y.sum()),
        "n_heldout_original_groups":len(unique),
        "same_species_pairs_and_split_in_all_models":True,
        "models":score,
        "climate_increments":increments,
        "inference_boundary":"One observer-disjoint photographed pair per same species; no genotype/pigment/local selection measured",
    }


def evaluate(pairs:pd.DataFrame,env:pd.DataFrame,*,strict:bool=True)->dict:
    # Reuse the original source-identity/year-independent matching code, but
    # deliberately DO NOT impose any SoilGrids complete-case qualification.
    subset,coverage=join_origins(pairs,env,require_soil=False)
    if strict and (len(subset)!=ASSERT_SOURCE_CLIMATE_PAIRS or
                   int(subset.pair_state.eq("discordant").sum())!=ASSERT_SOURCE_CLIMATE_DISCORDANT):
        raise RuntimeError("Previously checked source 3182/793 original photo-colour pairs drift")
    report={
        "schema":"fcp_global42111_3182_same_species_photo_pairs_temperature_vs_moisture_v1",
        "status":"SOURCE_CLIMATE_ONLY_ORIGINAL_PHOTO_PAIR_ABLATION_NOT_CAUSAL",
        "source_original_fixed_species_pairs":N_PAIRS,
        "source_original_four_colour_classifiable_pairs":BOTH_CLASSIFIABLE,
        "source_original_discordant_photo_pairs":DISCORDANT,
        "source_two_site_climate_evaluable_pairs":len(subset),
        "source_two_site_climate_evaluable_discordant":int(subset.pair_state.eq("discordant").sum()),
        "source_pair_missing_photo_colour_preserved":coverage["n_original_pairs_missing_classified_response"],
        "source_pair_environment_coverage":coverage,
        "source_temperature_features":list(THERMAL),
        "source_moisture_features":list(MOISTURE),
        "fixed_source_genus_and_midpoint_cell_fivefold":True,
        "no_source_photo_pixels_redownloaded":True,
        "no_new_photo_classification_or_prospective_source_opened":True,
        "original_highdepth_and_FCP_manuscript_not_reinterpreted":True,
        "models":[model_group(subset,"genus"),model_group(subset,"pair_midpoint_cell_162")],
        "nonclaims":[
            "This tests macroclimatic contrasts between photographs of different localities, not genetically identified same-population morphs",
            "Climate-only original pair selection differs from genus-conditioned BETWEEN-SPECIES sampled taxa",
            "Correlated BIO1/BIO5 and BIO12/BIO15 blocks do not provide independent causal contributions",
            "The source 10220 undetermined photographic pairs are not counted as matching flower-colour outcomes",
            "Blocked intervals condition on fixed folds and are retrospective, not external genetic or experimental confirmation",
        ],
    }
    return report


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--fixed-pairs",required=True,type=Path)
    p.add_argument("--original-photo-environment",required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    a=p.parse_args()
    pairs=read_original_pairs(a.fixed_pairs)
    env=pd.read_csv(a.original_photo_environment,low_memory=False)
    result=evaluate(pairs,env)
    a.outdir.mkdir(parents=True,exist_ok=True)
    (a.outdir/"result.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2,sort_keys=True),flush=True)
if __name__=="__main__":
    main()

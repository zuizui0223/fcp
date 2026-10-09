#!/usr/bin/env python3
"""Photo-colour permutation null for original global FCP local congeneric rain.

Preserve each original genus×162-cell photographed species group and its four
coarse photographed colour counts, original photo coordinates and WorldClim.
Shuffle original photographed colour labels only BETWEEN the distinct local
congeneric species, with a fixed source-only test-fold assignment. Refit the
exact same training-only within-group models. This calibrates the small
MOISTURE increment in the 250 and 500km source photo groups; does not establish
climate selection or genetic/within-species effects.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd

from audit_fcp_global42111_local_congeners_20261009 import (
    photo_populations,fold_plan,train_only_group_predict,CLASSES,GEO,TEMP,RAIN
)
from audit_fcp_global42111_local_congener_diameters_20261009 import group_geometry

DISTANCE_LIMITS=(250.,500.)
PERMUTATIONS=199
SEED=20261009
EXPECTED_ORIGINAL_N={"250":872,"500":1761}
EXPECTED_OBS_GAIN={"250":0.008598034415982893,"500":0.005115260801465297}


def fixed_source_cohort(pop:pd.DataFrame,maxkm:float)->pd.DataFrame:
    geometry=group_geometry(pop)
    keep=geometry.loc[geometry.max_source_photo_distance_km.le(maxkm+1e-8),"genus_cell_id"]
    x=pop.loc[pop.genus_cell_id.isin(keep)].copy().reset_index(drop=True)
    if x.inat_taxon_id.duplicated().any() or x.genus_cell_id.nunique()!=len(keep):
        raise RuntimeError("The exact pre-outcome source photo group membership changed")
    if (x.genus_cell_id.value_counts()<3).any():
        raise RuntimeError("Insufficient original local source congeners")
    return x


def outofsource_photo_moisture_gain(d:pd.DataFrame,fold:np.ndarray)->float:
    y=pd.Categorical(d.morph,categories=CLASSES).codes
    if (y<0).any():raise RuntimeError("Illegal original flower photo classification in permuted sample")
    ys=np.eye(len(CLASSES))[y]
    model1=GEO+TEMP
    model2=GEO+TEMP+RAIN
    total=np.zeros(len(d),float)
    visited=np.zeros(len(d),bool)
    for k in range(5):
        tr_idx=np.flatnonzero(fold!=k);te_idx=np.flatnonzero(fold==k)
        if len(te_idx)==0:continue
        train=d.iloc[tr_idx];test=d.iloc[te_idx]
        if test.genus_cell_id.isin(train.genus_cell_id.value_counts().loc[
                lambda x:x>=2].index).eq(False).any():
            raise RuntimeError("Permuted original test photo has no train-only congeners")
        p0=train_only_group_predict(train,test,model1)
        p1=train_only_group_predict(train,test,model2)
        target=ys[te_idx]
        total[te_idx]=np.sum((p0-target)**2,axis=1)-np.sum((p1-target)**2,axis=1)
        visited[te_idx]=True
    if not visited.all():raise RuntimeError("Not every original photo species received test prediction")
    return float(total.mean())


def null_calibration(d:pd.DataFrame,nperm:int=PERMUTATIONS,seed:int=SEED)->dict:
    if len(d)<300:raise ValueError("Too few local congeners for permutation inference")
    groups=[np.asarray(v,dtype=int) for _,v in d.groupby("genus_cell_id",sort=True).indices.items()]
    if not all(len(g)>=3 for g in groups):raise RuntimeError("Original photo group had <3 species")
    fold=fold_plan(d)
    original=d.morph.to_numpy(str)
    obs=outofsource_photo_moisture_gain(d,fold)
    rng=np.random.default_rng(seed)
    vals=[]
    for _ in range(nperm):
        labels=original.copy()
        for ids in groups:
            labels[ids]=rng.permutation(original[ids])
        fake=d.copy()
        fake["morph"]=labels
        val=outofsource_photo_moisture_gain(fake,fold)
        vals.append(val)
    values=np.asarray(vals,float)
    p=float((1+np.sum(values>=obs-1e-12))/(nperm+1))
    return {
        "n_original_photo_species":len(d),
        "n_original_genus_cell_groups":len(groups),
        "n_original_genera":int(d.genus.nunique()),
        "n_source_photo_geographic_cells":int(d.photo_cell_162.nunique()),
        "n_colour_classes_in_source":int(len(np.unique(original))),
        "null_label_shuffle_scope":"EXACT_SAME_GENUS_AND_SAME_ORIGINAL_162_CELL_ONLY",
        "group_colour_composition_preserved_in_every_null":True,
        "test_fold_assignments_preserved_in_every_null":True,
        "same_fitted_method_refit_for_every_null":True,
        "observed_moisture_added_multiclass_brier_gain":obs,
        "n_null_permutations":nperm,
        "null_mean_gain":float(values.mean()),
        "null_sd_gain":float(values.std(ddof=1)),
        "null_2p5_50_97p5":[float(x) for x in np.quantile(values,[.025,.5,.975])],
        "null_maximum_gain":float(values.max()),
        "one_sided_permutation_p_source_conditional":p,
        "n_null_gain_greater_than_observed":int(np.sum(values>=obs-1e-12)),
    }


def run(source:pd.DataFrame,*,strict:bool=True,nperm:int=PERMUTATIONS)->dict:
    pops,coverage=photo_populations(source,strict=strict)
    climat=pops["CLIMATE_ALL"]
    out={}
    for km in DISTANCE_LIMITS:
        key=str(int(km))
        chosen=fixed_source_cohort(climat,km)
        if strict and len(chosen)!=EXPECTED_ORIGINAL_N[key]:
            raise RuntimeError(f"{key}: frozen original local photo support changed")
        receipt=null_calibration(chosen,nperm=nperm,seed=SEED+int(km))
        if strict and abs(receipt["observed_moisture_added_multiclass_brier_gain"]-
                          EXPECTED_OBS_GAIN[key])>1e-8:
            raise RuntimeError("Independent source-verified heldout photo Brier gain mismatch")
        out[key]=receipt
    return {
        "schema":"fcp_global42111_local_congener_moisture_permutation_null_v1",
        "date_jst":"2026-10-09",
        "status":"HISTORICAL_SOURCE_WITHINGROUP_LABEL_SHUFFLE_CALIBRATION",
        "original_global_source_taxa":42111 if strict else len(source),
        "original_photo_colour_classified_source_taxa":coverage["source_classified_photo_taxa"],
        "original_climate_complete_species":coverage["climate_eligible_classified_source_taxa"],
        "distance_caps_km_precommitted_all_tested":list(DISTANCE_LIMITS),
        "permutations_precommitted":nperm,
        "source_groups_chosen_before_photo_colour_permutation":True,
        "source_covid_photo_pixels_or_new_labels_opened":False,
        "frozen_climate_photo_moisture_increment_results":out,
        "hard_nonclaims":[
            "A positive permutation calibrated association among congeneric photos is not causal climate selection",
            "Permutation is conditional on the original local genus-cell group and its complete-case source-photo sampling",
            "Overlapping 250km and 500km original sources are NOT independent evidence replications",
            "One photo per species, geographic range and within-cell environment still confound phenotypic biology",
            "199 permutations yield only 0.005 resolution for one-sided source-conditional p-values",
        ],
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-breadth-abiotic",required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    a=p.parse_args()
    d=pd.read_csv(a.original_breadth_abiotic,low_memory=False)
    result=run(d)
    a.outdir.mkdir(parents=True,exist_ok=True)
    (a.outdir/"result.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,indent=2,sort_keys=True),flush=True)


if __name__=="__main__":
    main()

#!/usr/bin/env python3
"""Matched original-photo predictive contrast: within species vs between species.

Both levels predict the EXACT SAME original photographed four-colour outcome
in a held-out geographic cell. Same climate BIO1/BIO5/BIO12/BIO15 and site GEO.
Within training: geographically distinct photos of the same named species.
Between training: ONLY other named species in same genus. Species-fold filter
strictly prevents test species entering any between training mean or slope.
The different *training information sets* are explicitly retained; Brier
increment contrasts are not effects that can be causally pooled or equated.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold

from compare_fcp_original85337_species_fixed_multiabiotic_20261010 import (
    match_original,populations,train_predict
)
from test_fcp_global42111_within_genus_environment_20261009 import (
    GEO,TEMPERATURE,MOISTURE,CLASSES
)

SCHEMA="fcp_original85337_matched_within_between_precipitation_comparison_v1"
N_CELLS=5
N_SPECIES_FOLDS=5
SOURCE_SEED=202610101915
BOOTSTRAPS=999
RIDGE=10.
MIN_ELIGIBLE=500
MIN_EVAL_SPECIES=100
BASE_FEATURES=GEO+TEMPERATURE
FULL_FEATURES=GEO+TEMPERATURE+MOISTURE
STUDY_BLOCKS={"geography_plus_temperature":BASE_FEATURES,
              "geography_plus_temperature_plus_precipitation":FULL_FEATURES}


def phylogeny_free_species_fold(taxon_id:int)->int:
    b=f"fcp-same-original-photo-fold-20261010|{int(taxon_id)}".encode()
    return int.from_bytes(hashlib.sha256(b).digest()[:8],"little")%N_SPECIES_FOLDS


def source_genus(v:pd.Series)->pd.Series:
    g=v.astype(str).str.strip().str.split().str[0]
    if g.eq("").any():raise RuntimeError("No valid original genus name")
    return g


def species_equal_training_genus_summaries(train:pd.DataFrame,features:tuple[str,...]):
    """Other species ONLY, equal species weight, not photo or genus abundance."""
    if train.empty:raise ValueError("No between-species training photographs")
    if train.inat_taxon_id.duplicated().all():
        raise RuntimeError("Training must include distinct original taxa")
    coords=train[list(features)].to_numpy(float)
    y=pd.get_dummies(pd.Categorical(train.morph,categories=CLASSES)).to_numpy(float)
    if y.shape[1]!=4 or not np.isfinite(coords).all():raise RuntimeError("Source original feature or colour mismatch")
    means=coords.mean(axis=0);sds=coords.std(axis=0)
    sds[sds<1e-8]=1.
    x=(coords-means)/sds
    d=pd.DataFrame(x,columns=list(features)).assign(
       genus=train.genus.to_numpy(str),
       taxon=train.inat_taxon_id.to_numpy(int),
       **{f"colour_{i}":y[:,i] for i in range(4)})
    if d.groupby("taxon").genus.nunique().max()>1:raise RuntimeError("source taxon belongs to more than one genus")
    source=d.groupby(["taxon","genus"],sort=False).mean(numeric_only=True).reset_index()
    sizes=source.genus.value_counts()
    px=list(features);py=[f"colour_{i}" for i in range(4)]
    g=source.groupby("genus",sort=False)[px+py].mean()
    xx=source[px].to_numpy(float)
    yy=source[py].to_numpy(float)
    centered_x=xx-g.loc[source.genus,px].to_numpy(float)
    centered_y=yy-g.loc[source.genus,py].to_numpy(float)
    beta=np.linalg.solve(centered_x.T@centered_x+RIDGE*np.eye(len(features)),
                         centered_x.T@centered_y)
    return {"feature_mean":means,"feature_sd":sds,"genera":g,
            "counts":sizes,"beta":beta}


def predict_between_other_species(train:pd.DataFrame,target:pd.DataFrame,
                                   features:tuple[str,...])->np.ndarray:
    """Never use target taxon in BETWEEN training at any geographic cell."""
    if set(train.inat_taxon_id)&set(target.inat_taxon_id):
        raise RuntimeError("test species leaked into between-species training")
    summary=species_equal_training_genus_summaries(train,features)
    count=summary["counts"]
    if not target.genus.isin(count[count>=3].index).all():
        raise RuntimeError("heldout genus has fewer than 3 independent source species")
    g=summary["genera"].loc[target.genus]
    feat=list(features)
    z=(target[feat].to_numpy(float)-summary["feature_mean"])/summary["feature_sd"]
    out=g[[f"colour_{i}" for i in range(4)]].to_numpy(float)+(
        z-g[feat].to_numpy(float))@summary["beta"]
    out=np.maximum(0,out)
    total=out.sum(axis=1,keepdims=True)
    return np.divide(out,total,out=np.full_like(out,.25),where=total>0)


def cluster_ci(values:np.ndarray,ids:np.ndarray,*,label:str)->list[float]:
    _,inverse=np.unique(ids,return_inverse=True)
    sums=np.bincount(inverse,weights=values)
    counts=np.bincount(inverse)
    if len(sums)<5:raise RuntimeError("Insufficient independent source groups")
    seed=int.from_bytes(hashlib.sha256(f"{SOURCE_SEED}|{label}".encode()).digest()[:8],"little")
    rng=np.random.default_rng(seed)
    idx=rng.integers(len(sums),size=(BOOTSTRAPS,len(sums)))
    draws=sums[idx].sum(axis=1)/counts[idx].sum(axis=1)
    return [float(v) for v in np.quantile(draws,[.025,.975])]


def run(original_cells:pd.DataFrame,full_sites:pd.DataFrame,*,strict=True):
    joined,old=match_original(original_cells,full_sites,strict=strict)
    pools,ledger=populations(joined)
    data=pools["CLIMATE_SOURCE_ONLY"].copy()
    assert set(FULL_FEATURES).issubset(data.columns)
    data["genus"]=source_genus(data.species)
    data["species_fold"]=data.inat_taxon_id.map(phylogeny_free_species_fold)
    response=pd.Categorical(data.morph,categories=CLASSES).codes
    if (response<0).any():raise RuntimeError("source colour drift")
    pred={lvl:{m:np.full((len(data),4),np.nan,float) for m in STUDY_BLOCKS}
          for lvl in ["within_species","between_species"]}
    supported=np.zeros(len(data),bool)
    geometries=[]
    for spatial_fold,(itr,ite) in enumerate(GroupKFold(n_splits=N_CELLS).split(
        data,response,groups=data.cell_id)):
        tr=data.iloc[itr]
        te=data.iloc[ite]
        n_local=0
        for taxon_fold in range(N_SPECIES_FOLDS):
            ids=np.flatnonzero(te.species_fold.to_numpy(int)==taxon_fold)
            if len(ids)==0:continue
            target=te.iloc[ids]
            other=tr.loc[tr.species_fold.ne(taxon_fold)].copy()
            # Strictly held-out species, NOT some photo from it.
            if set(other.inat_taxon_id)&set(target.inat_taxon_id):
                raise RuntimeError("between-species target leakage")
            genus_counts=other[["genus","inat_taxon_id"]].drop_duplicates().genus.value_counts()
            species_train=tr.inat_taxon_id.value_counts()
            allowable=(target.genus.isin(genus_counts[genus_counts>=3].index)
                       &target.inat_taxon_id.isin(species_train.index))
            matched=ids[allowable.to_numpy(bool)]
            if len(matched)==0:continue
            eligible_target=te.iloc[matched]
            if set(eligible_target.cell_id)&set(tr.cell_id):
                raise RuntimeError("same geographical cell entered both sides")
            for model,features in STUDY_BLOCKS.items():
                pred["within_species"][model][ite[matched]]=train_predict(tr,te,features,matched)
                pred["between_species"][model][ite[matched]]=predict_between_other_species(
                    other,eligible_target,features)
            supported[ite[matched]]=True
            n_local+=len(matched)
        geometries.append({"geographic_fold":spatial_fold,
            "n_test_photos":len(ite),
            "n_common_supported_test_photos":n_local,
            "n_geographic_cells_heldout":te.cell_id.nunique()})
    n=int(supported.sum())
    species=data.loc[supported,"inat_taxon_id"].to_numpy(int)
    genus=data.loc[supported,"genus"].to_numpy(str)
    site=data.loc[supported,"cell_id"].to_numpy(int)
    if n<MIN_ELIGIBLE or len(np.unique(species))<MIN_EVAL_SPECIES:
        return {"status":"HOLD_INSUFFICIENT_EXACTLY_MATCHED_COMMON_SOURCE_PHOTOS",
                "n_common_test_photos":n,
                "n_common_test_species":int(len(np.unique(species)))}
    truth=np.eye(4)[response[supported]]
    results={}
    per_photo={}
    for level in pred:
        p=pred[level]
        for name,arr in p.items():
            if not np.isfinite(arr[supported]).all() or (arr[supported]<0).any():
                raise RuntimeError("A source model lacks a matched prediction")
        before=np.square(truth-p["geography_plus_temperature"][supported]).sum(axis=1)
        after=np.square(truth-p["geography_plus_temperature_plus_precipitation"][supported]).sum(axis=1)
        gain=before-after
        per_photo[level]=gain
        results[level]={
            "photo_equal_Brier_before":float(before.mean()),
            "photo_equal_Brier_after":float(after.mean()),
            "precipitation_unique_Brier_gain":float(gain.mean()),
            "species_cluster_95CI":cluster_ci(gain,species,label=level+"_species"),
            "genus_cluster_95CI":cluster_ci(gain,genus,label=level+"_genus"),
            "original_geographical_cell_cluster_95CI":cluster_ci(gain,site,label=level+"_cell"),
            "n_source_species":int(len(np.unique(species))),
            "n_source_genera":int(len(np.unique(genus))),
            "n_original_geographical_cells":int(len(np.unique(site))),
        }
    differential=per_photo["within_species"]-per_photo["between_species"]
    return {
       "schema":SCHEMA,"status":"POSTHOC_SOURCE_MATCHED_WITHIN_BETWEEN_PREDICTIVE_SCORE_COMPARISON",
       "n_source_cell_ledger":len(original_cells),
       "n_original_fourstate_classifiable_cells":old["source_four_state_colour_classifiable"],
       "n_original_climate_complete_colour_cells":ledger["n_original_classified_photo_cells_climate_complete"],
       "n_common_test_photos":n,
       "n_common_species":int(len(np.unique(species))),
       "n_common_genera":int(len(np.unique(genus))),
       "n_source_cells_used_in_heldout_evaluation":int(len(np.unique(site))),
       "n_source_photo_geography_folds":N_CELLS,"n_source_taxon_disjoint_between_training_folds":N_SPECIES_FOLDS,
       "source_train_test_key_guards":True,
       "source_colour_response_categories":list(CLASSES),
       "common_geographic_features":list(GEO),
       "common_temperature_features":list(TEMPERATURE),
       "added_precipitation_features":list(MOISTURE),
       "matched_test_photo_ids":True,
       "models":results,
       "within_minus_between_precipitation_gain":{
         "paired_gain_difference":float(differential.mean()),
         "species_cluster_95CI":cluster_ci(differential,species,label="diff_species"),
         "genus_cluster_95CI":cluster_ci(differential,genus,label="diff_genus"),
         "geographical_cell_cluster_95CI":cluster_ci(differential,site,label="diff_cell"),
         "claim":"difference in conditional predictive gains under unequal source-training-information sets; not a difference in causal ecological coefficients"},
       "fivefold_source_cell_support":geometries,
       "date_jst":"2026-10-10",
       "confirmatory_decisions_changed":False,
       "hard_nonclaims":[
          "Common photo scoring does NOT imply that separately estimated within- and between-species regressions share a causal estimator.",
          "Between species excludes the focal source test species from ALL training photos; within species explicitly includes that species in other geographic cells.",
          "Actual geography was entered as fixed lat/lon/elevation terms, NOT a full true-site stochastic spatial kernel or genetic model.",
          "Species/genus definitions are nominal source taxonomy; related congeners have unresolved within-genus tree distance.",
          "Colour outcomes are historical photograph classifications without human blinded expert verification.",
          "All Brier intervals condition on this fixed spatial/species crossvalidation and are post hoc.",
          "A positive rain predictor does not identify common ecological selection, inherited floral pigments or fitness.",
       ]
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-taxon-cell",type=Path,required=True)
    p.add_argument("--all-original-expanded-photo-sites",type=Path,required=True)
    p.add_argument("--outdir",type=Path,required=True)
    a=p.parse_args()
    src=pd.read_csv(a.original_taxon_cell,low_memory=False)
    sites=pd.read_csv(a.all_original_expanded_photo_sites,low_memory=False)
    z=run(src,sites)
    a.outdir.mkdir(parents=True,exist_ok=True)
    (a.outdir/"result.json").write_text(json.dumps(z,indent=2)+"\n")
    print(json.dumps(z,indent=2))
if __name__=="__main__":main()

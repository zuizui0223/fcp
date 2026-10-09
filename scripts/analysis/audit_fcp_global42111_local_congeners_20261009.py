#!/usr/bin/env python3
"""Local congeneric species contrast in the 42,111-original-photo FCP atlas.

Stricter spatial composition control than whole-genus training: only compare
DIFFERENT original source species WITHIN the SAME genus AND SAME original 162
equal-area geography cell. Group-specific colour means and environment means
come from the OTHER training species, never the held-out photo.
Not within-species, not genetic polymorphism, not causal climate selection.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

CLASSES=("white","yellow_orange","red_pink","blue_purple")
GEO=("abs_latitude","lon_sin","lon_cos","wc_elevation_m")
TEMP=("wc_bio1","wc_bio5")
RAIN=("wc_bio12","wc_bio15")
SOIL=("soil_pH","soil_SOC","soil_N","soil_clay","soil_available_water_proxy")
FAMILIES={
    "GENUS_CELL_BASELINE":(),
    "GENUS_CELL_GEO":GEO,
    "GENUS_CELL_GEO_TEMPERATURE":GEO+TEMP,
    "GENUS_CELL_GEO_MOISTURE":GEO+RAIN,
    "GENUS_CELL_GEO_ALL_CLIMATE":GEO+TEMP+RAIN,
    "GENUS_CELL_GEO_ALL_CLIMATE_SOIL":GEO+TEMP+RAIN+SOIL,
}
SOURCE_TOTAL=42111
SOURCE_CLASSIFIED=18457
CLIMATE_AVAILABLE=18413
SOIL_AVAILABLE=14136
FOLDS=5
MIN_GROUP_TAXA=3
MIN_EVAL=300
BOOT=999
RIDGE=10.0
SEED=20261009


def cell_id(lat:np.ndarray,lon:np.ndarray)->np.ndarray:
    lat=np.asarray(lat,dtype=float);lon=np.asarray(lon,dtype=float)
    out=np.full(len(lat),-1,dtype=int)
    good=np.isfinite(lat)&np.isfinite(lon)&(np.abs(lat)<=90)&(np.abs(lon)<=180)
    y=np.clip(np.floor((1+np.sin(np.deg2rad(lat[good])))*4.5).astype(int),0,8)
    x=np.clip(np.floor((lon[good]+180)/20).astype(int),0,17)
    out[good]=18*y+x
    return out


def photo_populations(source:pd.DataFrame,*,strict=True)->tuple[dict[str,pd.DataFrame],dict]:
    necessary={"inat_taxon_id","species","morph","measurement_status","latitude","longitude",
               "site_geo_status","environment_climate_complete","environment_soil_complete",
               "wc_elevation_m",*TEMP,*RAIN,*SOIL}
    if not necessary.issubset(source):
        raise ValueError("Frozen original photo environmental fields missing")
    if source.inat_taxon_id.duplicated().any():
        raise ValueError("Original one-photo species duplicated")
    if strict and len(source)!=SOURCE_TOTAL:
        raise ValueError("Frozen original source 42111 species denominator changed")
    z=source.copy()
    classified=z.morph.isin(CLASSES)&z.measurement_status.eq("classified_four_state_morph")
    if strict and classified.sum()!=SOURCE_CLASSIFIED:
        raise ValueError("Original source colour-classifiable 18457 count changed")
    lat=pd.to_numeric(z.latitude,errors="coerce")
    lon=pd.to_numeric(z.longitude,errors="coerce")
    z["abs_latitude"]=lat.abs()
    z["lon_sin"]=np.sin(np.deg2rad(lon))
    z["lon_cos"]=np.cos(np.deg2rad(lon))
    z["photo_cell_162"]=cell_id(lat.to_numpy(float),lon.to_numpy(float))
    z["genus"]=z.species.fillna("").astype(str).str.strip().str.split().str[0]
    z["genus_cell_id"]=z.genus+"|"+z.photo_cell_162.astype(str)
    pos=z.site_geo_status.eq("VALID_PUBLIC_ORIGINAL_PHOTO_POINT")&z.photo_cell_162.ge(0)
    clim=classified&pos&z.environment_climate_complete.fillna(False).astype(bool)&z[list(GEO+TEMP+RAIN)].notna().all(axis=1)
    soil=clim&z.environment_soil_complete.fillna(False).astype(bool)&z[list(SOIL)].notna().all(axis=1)
    pools={"CLIMATE_ALL":z.loc[clim].copy(),"CLIMATE_SOIL_COMPLETE":z.loc[soil].copy()}
    stats={
        "source_taxa":len(z),
        "source_classified_photo_taxa":int(classified.sum()),
        "source_unclassified_photo_taxa":int((~classified).sum()),
        "climate_eligible_classified_source_taxa":int(clim.sum()),
        "soil_eligible_classified_source_taxa":int(soil.sum()),
        "n_source_photo_cells_with_climate":int(z.loc[clim,"photo_cell_162"].nunique()),
        "n_source_genera_with_climate":int(z.loc[clim,"genus"].nunique()),
    }
    if strict and (stats["climate_eligible_classified_source_taxa"]!=CLIMATE_AVAILABLE or
                   stats["soil_eligible_classified_source_taxa"]!=SOIL_AVAILABLE):
        raise ValueError("Frozen 18413 and 14136 original photo environmental case counts drifted")
    for pop in pools.values():
        if pop.genus.eq("").any() or pop.genus_cell_id.str.contains(r"\|-1$").any():
            raise ValueError("Source genus or original cell absent")
    return pools,stats


def source_group_coverage(pop:pd.DataFrame)->dict:
    freq=pop.groupby("genus_cell_id",sort=False).size()
    enough=freq.ge(MIN_GROUP_TAXA)
    target=pop.genus_cell_id.isin(freq[enough].index)
    return {
        "n_original_complete_photo_species":len(pop),
        "n_original_genus_cell_groups":int(len(freq)),
        "n_groups_with_three_or_more_original_species":int(enough.sum()),
        "n_photo_species_in_groups_with_three_or_more":int(target.sum()),
        "n_original_genera_with_locally_supported_congeners":int(pop.loc[target,"genus"].nunique()),
        "n_original_photo_cells_with_locally_supported_congeners":int(pop.loc[target,"photo_cell_162"].nunique()),
        "n_one_or_two_species_genus_cells":int((~enough).sum()),
        "n_source_species_lacking_three_local_congeners":int((~target).sum()),
        "group_support_decided_without_opening_flower_colour_outcomes":True,
    }


def fold_plan(pop:pd.DataFrame)->np.ndarray:
    if pop.inat_taxon_id.duplicated().any():
        raise ValueError("Duplicate original photo-taxonomy identity")
    # Stable taxon-ID-only hash, independent of PHOTO COLOUR or abiotic values.
    fold=np.full(len(pop),-1,dtype=int)
    for group,indices in pop.groupby("genus_cell_id",sort=True).indices.items():
        order=sorted(indices,key=lambda i:hashlib.sha256(
            f"fcp_20261009_local_congeners|{int(pop.iloc[i].inat_taxon_id)}".encode()).hexdigest())
        for k,i in enumerate(order):
            fold[i]=k%FOLDS
    if (fold<0).any():raise RuntimeError("Missing local congener test fold")
    return fold


def train_only_group_predict(train:pd.DataFrame,test:pd.DataFrame,features:tuple[str,...])->np.ndarray:
    if train.inat_taxon_id.duplicated().any() or test.inat_taxon_id.duplicated().any():
        raise ValueError("One source photo per taxon is mandatory")
    if set(train.inat_taxon_id)&set(test.inat_taxon_id):
        raise ValueError("Test photographed species in source training group")
    group="genus_cell_id"
    groups=train[group].to_numpy(str)
    counts=pd.Series(groups).value_counts()
    if not test[group].isin(counts[counts>=2].index).all():
        raise ValueError("Unknown or singly sampled local genus-cell test group")
    y=pd.get_dummies(pd.Categorical(train.morph,categories=CLASSES)).to_numpy(float)
    if y.shape[1]!=4:raise ValueError("Unknown original source photo label")
    ym=pd.DataFrame(y,columns=CLASSES).assign(group=groups).groupby("group",sort=False)[list(CLASSES)].mean()
    baseline=ym.reindex(test[group].to_numpy(str)).to_numpy(float)
    if not features:return baseline
    xt=train[list(features)].to_numpy(float)
    xe=test[list(features)].to_numpy(float)
    if not np.isfinite(xt).all() or not np.isfinite(xe).all():
        raise ValueError("Invalid source environmental predictor in local comparison")
    mean=xt.mean(axis=0);sd=xt.std(axis=0);sd[sd<1e-10]=1.0
    xt=(xt-mean)/sd;xe=(xe-mean)/sd
    xm=pd.DataFrame(xt,columns=list(features)).assign(group=groups).groupby("group",sort=False)[list(features)].mean()
    centered_x=xt-xm.reindex(groups).to_numpy(float)
    centered_y=y-ym.reindex(groups).to_numpy(float)
    beta=np.linalg.solve(centered_x.T@centered_x+RIDGE*np.eye(len(features)),centered_x.T@centered_y)
    pred=baseline+(xe-xm.reindex(test[group].to_numpy(str)).to_numpy(float))@beta
    pred=np.clip(pred,0,None)
    totals=pred.sum(axis=1,keepdims=True)
    return np.divide(pred,totals,out=np.full_like(pred,0.25),where=totals>0)


def grouped_bootstrap_loss(gain:np.ndarray,groups:np.ndarray,seed:int)->list[float]:
    unique,ids=np.unique(groups,return_inverse=True)
    n=np.bincount(ids)
    sums=np.bincount(ids,weights=gain)
    rng=np.random.default_rng(seed)
    draw=rng.integers(0,len(unique),size=(BOOT,len(unique)))
    values=sums[draw].sum(axis=1)/n[draw].sum(axis=1)
    return [float(v) for v in np.quantile(values,[.025,.975])]


def test_local_congeners(pop:pd.DataFrame,*,include_soil:bool)->dict:
    info=source_group_coverage(pop)
    result={
        "schema":"fcp_original_42111_genus_cell_fixed_photo_climate_v1",
        "status":"SOURCE_GROUP_SUPPORT_PREFLIGHT",
        "original_photo_support":info,
        "min_source_species_per_genus_and_geographic_cell":MIN_GROUP_TAXA,
        "photo_colour_never_used_for_group_selection":True,
        "one_photo_per_source_species":True,
        "within_same_genus_and_same_source_162_cell_only":True,
        "geographic_cell_holdout_not_feasible_for_identical_local_group":True,
    }
    freq=pop.genus_cell_id.value_counts()
    take=pop.genus_cell_id.isin(freq.loc[freq>=MIN_GROUP_TAXA].index)
    d=pop.loc[take].copy().reset_index(drop=True)
    result["n_comparable_original_local_photo_species"]=len(d)
    if len(d)<MIN_EVAL or len(set(d.morph))<4:
        result["status"]="HOLD_INSUFFICIENT_SAME_GENUS_SAME_CELL_SPECIES"
        return result
    fold=fold_plan(d)
    family={k:v for k,v in FAMILIES.items() if include_soil or k!="GENUS_CELL_GEO_ALL_CLIMATE_SOIL"}
    y=pd.Categorical(d.morph,categories=CLASSES).codes
    onehot=np.eye(4)[y]
    predicted={name:np.full((len(d),4),np.nan,float) for name in family}
    tested=np.zeros(len(d),bool)
    fold_receipts=[]
    for k in range(FOLDS):
        train_idx=np.flatnonzero(fold!=k)
        test_idx=np.flatnonzero(fold==k)
        if len(test_idx)==0:continue
        train=d.iloc[train_idx];test=d.iloc[test_idx]
        counts=train.genus_cell_id.value_counts()
        good=test.genus_cell_id.isin(counts[counts>=2].index).to_numpy()
        test_idx=test_idx[good];test=d.iloc[test_idx]
        if len(test_idx)==0:continue
        if test.inat_taxon_id.isin(train.inat_taxon_id).any():
            raise RuntimeError("Same photographed species appeared on both sides")
        tested[test_idx]=True
        for name,features in family.items():
            predicted[name][test_idx]=train_only_group_predict(train,test,features)
        fold_receipts.append({"fold":k,"n_test_with_training_local_congeners":len(test_idx),
                              "n_source_genus_cells_in_test":int(test.genus_cell_id.nunique())})
    if int(tested.sum())<MIN_EVAL:
        result["status"]="HOLD_INSUFFICIENT_LOCAL_TRAIN_ONLY_GENUS_CELL_SUPPORT"
        return result
    if np.any(~np.isfinite(np.stack([v[tested] for v in predicted.values()]))):
        raise RuntimeError("Local genus cell test photo did not get all comparable model predictions")
    losses={}
    summary={}
    weights=[]
    for name,prob in predicted.items():
        err=np.sum((prob[tested]-onehot[tested])**2,axis=1)
        losses[name]=err
        summary[name]={"source_test_photo_species":int(tested.sum()),"heldout_multiclass_brier":float(err.mean())}
    extra={
        "climate_beyond_local_geography":("GENUS_CELL_GEO","GENUS_CELL_GEO_ALL_CLIMATE"),
        "temperature_unique_beyond_local_geography_rain":("GENUS_CELL_GEO_MOISTURE","GENUS_CELL_GEO_ALL_CLIMATE"),
        "moisture_unique_beyond_local_geography_temperature":("GENUS_CELL_GEO_TEMPERATURE","GENUS_CELL_GEO_ALL_CLIMATE"),
    }
    if include_soil:extra["soil_beyond_local_geography_climate"]=("GENUS_CELL_GEO_ALL_CLIMATE","GENUS_CELL_GEO_ALL_CLIMATE_SOIL")
    increment={}
    for key,(baseline,full) in extra.items():
        gain=losses[baseline]-losses[full]
        region=d.loc[tested,"photo_cell_162"].to_numpy(int)
        genera=d.loc[tested,"genus"].to_numpy(str)
        local=d.loc[tested,"genus_cell_id"].to_numpy(str)
        increment[key]={
            "mean_heldout_brier_reduction":float(gain.mean()),
            "original_cell_cluster_bootstrap_95CI":grouped_bootstrap_loss(gain,region,SEED+len(key)),
            "genus_cluster_bootstrap_95CI":grouped_bootstrap_loss(gain,genera,SEED+len(key)+100),
            "local_genus_cell_cluster_bootstrap_95CI":grouped_bootstrap_loss(gain,local,SEED+len(key)+200),
            "all_three_CIs_entirely_positive":False,
        }
        z=increment[key]
        z["all_three_CIs_entirely_positive"]=all(z[k][0]>0 for k in (
            "original_cell_cluster_bootstrap_95CI","genus_cluster_bootstrap_95CI",
            "local_genus_cell_cluster_bootstrap_95CI"))
    result.update({
        "status":"SOURCE_LOCAL_GENUS_CELL_HELDOUT_PHOTO_SPECIES_ANALYSIS",
        "n_original_photo_species_evaluated":int(tested.sum()),
        "n_original_taxonomic_genera_evaluated":int(d.loc[tested,"genus"].nunique()),
        "n_local_genus_cell_groups_evaluated":int(d.loc[tested,"genus_cell_id"].nunique()),
        "n_global_photo_cells_evaluated":int(d.loc[tested,"photo_cell_162"].nunique()),
        "five_photo_source_species_folds":fold_receipts,
        "heldout_4class_photo_brier":summary,
        "incremental_source_photo_prediction":increment,
        "predictors":{k:list(v) for k,v in family.items()},
        "inference_limits":[
            "The training/control compares DIFFERENT species WITHIN the same genus and geographic cell, not different sites for one species",
            "The original 162 grid has substantial geographic extent; within-cell spatial gradients and environmental covariance remain",
            "Local source groups were required to have >=3 taxa without using photo-colour outcome; excluded source genera/cells are not inferentially represented",
            "The same local group necessarily appears in training; this is local conditional out-of-species prediction, NOT out-of-cell extrapolation",
            "Resampling fixed predictions is not retrained or phylogenetically independent confirmation",
            "Predicted group means do not identify genetic polymorphism or selection fitness",
        ],
    })
    return result


def run(source:pd.DataFrame,*,strict=True)->dict:
    pools,coverage=photo_populations(source,strict=strict)
    return {
        "schema":"fcp_global42111_local_congeneric_photo_colour_environment_v1",
        "date_jst":"2026-10-09",
        "status":"HISTORIC_SOURCE_LOCAL_CONGENERIC_PHOTOGRAPHIC_GRADIENT_TEST",
        "original_42111_denominator":SOURCE_TOTAL if strict else len(source),
        "original_photo_source_coverage":coverage,
        "climate_only_local_congeners":test_local_congeners(pools["CLIMATE_ALL"],include_soil=False),
        "soil_complete_local_congeners":test_local_congeners(pools["CLIMATE_SOIL_COMPLETE"],include_soil=True),
        "old_highdepth_and_prospective_taxa_unmodified":True,
        "no_new_original_photo_colour_or_raster_read":True,
        "interpretation":"Original-photo congeneric species difference conditional on same 162-cell group; not causal local selection or within-species flower-colour evolution",
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-breadth-abiotic",required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    a=p.parse_args()
    source=pd.read_csv(a.original_breadth_abiotic,low_memory=False)
    result=run(source)
    a.outdir.mkdir(parents=True,exist_ok=True)
    (a.outdir/"result.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,indent=2,sort_keys=True),flush=True)


if __name__=="__main__":
    main()

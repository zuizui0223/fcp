#!/usr/bin/env python3
"""Original FCP 85337 repeated species-cell flower photos: true spatial kernel.

Compare environment after TRAIN-only species colour means and a low-rank
Nyström approximation of exp(-great-circle-distance/scale) covariance.
Original location/photo ID/morph fixed. Original equal-area cells held out,
and group/site-missing rows remain in the source denominator.

This is a predictive spatial covariance approximation, not proof residual
Moran's I is zero, heritable colour variation or environmental selection.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from compare_fcp_original85337_species_fixed_multiabiotic_20261010 import (
    match_original, populations, train_predict
)
from compare_fcp_global42111_full_climate_sun_soil_blocks_20261010 import (
    BLOCKS,GEO,CLASSES
)

SCHEMA="fcp_original_85337_species_fixed_geodesic_spatial_kernel_abiotic_v1"
SCALES_KM=(100.,500.)
N_LANDMARKS=128
SEED=20261010
BOOT=199
MIN_ELIGIBLE=200
N_FOLDS=5
EARTH_KM=6371.0088


def xyz(real_photo_lat:np.ndarray,real_photo_lon:np.ndarray)->np.ndarray:
    a=np.asarray(real_photo_lat,float)
    b=np.asarray(real_photo_lon,float)
    if a.ndim!=1 or b.ndim!=1 or len(a)!=len(b):
        raise ValueError("Original photo site arrays have inconsistent size")
    if not (np.isfinite(a).all() and np.isfinite(b).all() and
            np.all(np.abs(a)<=90) and np.all(np.abs(b)<=180)):
        raise ValueError("Source photo geographic coordinates are missing or impossible")
    a=np.deg2rad(a);b=np.deg2rad(b)
    return np.column_stack((np.cos(a)*np.cos(b),np.cos(a)*np.sin(b),np.sin(a)))


def geodesic_exponential_kernel(a:np.ndarray,b:np.ndarray,scale_km:float)->np.ndarray:
    if scale_km<=0:raise ValueError("Source spatial length scale positive")
    dot=np.clip(a@b.T,-1,1)
    return np.exp(-np.arccos(dot)*EARTH_KM/scale_km)


def nystrom_train_test(train:np.ndarray,test:np.ndarray,*,scale_km:float,
                       nlandmarks:int=N_LANDMARKS,seed:int=SEED)->tuple[np.ndarray,np.ndarray,dict]:
    if len(train)<nlandmarks or nlandmarks<4:
        raise ValueError("Insufficient TRAIN original photo site landmarks")
    # Landmarks chosen from TRAIN geographic photo positions, not colour outcome
    rng=np.random.default_rng(seed)
    indexes=rng.choice(len(train),nlandmarks,replace=False)
    land=train[indexes]
    landmark_matrix=geodesic_exponential_kernel(land,land,scale_km)
    landmark_matrix=(landmark_matrix+landmark_matrix.T)/2
    eig,vec=np.linalg.eigh(landmark_matrix)
    positive=np.maximum(eig,1e-6)
    inverse_root=(vec/np.sqrt(positive))@vec.T
    ftrain=geodesic_exponential_kernel(train,land,scale_km)@inverse_root
    ftest=geodesic_exponential_kernel(test,land,scale_km)@inverse_root
    if not(np.isfinite(ftrain).all() and np.isfinite(ftest).all()):
        raise RuntimeError("Original source location Nyström basis nonfinite")
    return ftrain,ftest,{
        "n_source_training_only_landmarks":nlandmarks,
        "geodesic_kernel_scale_km":scale_km,
        "min_original_source_landmark_gram_eigenvalue":float(eig.min()),
        "landmark_selection_colour_outcome_blind":True,
    }


def families(*,soil:bool,nlandmarks:int=N_LANDMARKS)->dict[str,tuple[str,...]]:
    included={k:v for k,v in BLOCKS.items() if soil or k!="soil"}
    env=tuple(v for block in included.values() for v in block)
    basis=tuple("spatial_kernel_feature_"+str(i) for i in range(nlandmarks))
    controls=GEO+basis
    result={
        "SPECIES_GEOGRAPHY_ONLY":GEO,
        "SPECIES_GEOGRAPHY_SPATIAL_KERNEL":controls,
        "SPECIES_GEOGRAPHY_SPATIAL_KERNEL_ALL_ENV":controls+env,
    }
    for name,variables in included.items():
        result["SPATIAL_ALL_MINUS_"+name.upper()]=controls+tuple(v for v in env if v not in variables)
    return result


def residual_neighbour_autocorrelation(
    eligible_source:pd.DataFrame,original_classes:np.ndarray,
    prediction:dict[str,np.ndarray],*,radius_km:tuple[int,...]=(100,500),
)->dict:
    """Descriptive directed k-nearest-photo residual Moran I by colour.

    No permutation p-value, no guarantee of source field stationarity, and
    no claim that a near-zero statistic eliminates all spatial confounding.
    """
    from sklearn.neighbors import BallTree
    n=len(eligible_source)
    if n<12:
        raise ValueError("Too few original photo sites for 8-nearest residual audit")
    # BallTree haversine requires latitude then longitude, in RADIANS.
    original_sites=eligible_source[["latitude","longitude"]].to_numpy(float)
    if not (np.isfinite(original_sites).all() and
            np.all(np.abs(original_sites[:,0])<=90) and
            np.all(np.abs(original_sites[:,1])<=180)):
        raise ValueError("Original source photographed coordinates missing or impossible")
    coordinates=np.deg2rad(original_sites)
    tree=BallTree(coordinates,metric="haversine")
    dist,other=tree.query(coordinates,k=min(n,9))
    dist=dist[:,1:]*EARTH_KM
    other=other[:,1:]
    classes=np.asarray(original_classes,int)
    if classes.shape!=(n,) or np.any((classes<0)|(classes>=4)):
        raise ValueError("Mismatched original photo-colour outcome in residual audit")
    yy=np.eye(4)[classes]
    out={}
    for scale in radius_km:
        mask=(dist<=scale)
        edge_count=int(mask.sum())
        if edge_count==0:
            out[str(scale)]={"status":"HOLD_NO_GENUINELY_NEARBY_SOURCE_PHOTO_PAIRS",
                             "n_nearby_directed_edges":0}
            continue
        summary={}
        for name,prob in prediction.items():
            if prob.shape!=(n,4) or not np.isfinite(prob).all():
                raise ValueError("One of the fixed heldout predicted four-photo states was unavailable")
            residual=yy-prob
            centered=residual-residual.mean(axis=0)
            denominator=np.sum(centered**2,axis=0)
            # Preserve four separate source photo-colour classes, and total
            # directed source-neighbor edge weight, no pseudoreplication p.
            separate=[]
            for c in range(4):
                weighted=(centered[:,c,None]*centered[other,c])
                val=float(n/edge_count * weighted[mask].sum()/denominator[c]) if denominator[c]>1e-12 else None
                separate.append(val)
            finite=[x for x in separate if x is not None]
            summary[name]={
                "source_four_colour_class_residual_Moran_I":dict(zip(CLASSES,separate)),
                "mean_absolute_class_Moran_I":float(np.mean(np.abs(finite))) if finite else None,
                "n_source_geo_photo_records":n,
            }
        out[str(scale)]={
            "status":"ORIGINAL_NEAREST_PHOTO_RESIDUAL_SPATIAL_AUTOCORRELATION_DIAGNOSTIC",
            "n_nearby_directed_edges":edge_count,
            "mean_neighbors_within_radius":edge_count/n,
            "max_nearest_neighbors_examined":min(n-1,8),
            "residual_methods":summary,
            "descriptive_not_a_significance_test":True,
        }
    return out


def fit_source(d:pd.DataFrame,*,soil:bool,scale_km:float,nboot:int=BOOT,
               nlandmarks:int=N_LANDMARKS)->dict:
    from sklearn.model_selection import GroupKFold
    if len(d)<MIN_ELIGIBLE or d.cell_id.nunique()<N_FOLDS:
        return {"status":"HOLD_SOURCE_SPECIES_REGION_SPATIAL_COVERAGE"}
    if not d.morph.isin(CLASSES).all() or d.inat_taxon_id.isna().any():
        raise ValueError("Original unclassifiable photo or source taxon disappeared")
    methods=families(soil=soil,nlandmarks=nlandmarks)
    source_labels=pd.Categorical(d.morph,categories=CLASSES).codes
    coords=xyz(d.latitude.to_numpy(float),d.longitude.to_numpy(float))
    predicted={k:np.full((len(d),4),np.nan,float) for k in methods}
    eligible=np.zeros(len(d),bool)
    fold_report=[]
    spatial_receipts=[]
    for fold,(it,ie) in enumerate(GroupKFold(n_splits=N_FOLDS).split(
        d,source_labels,groups=d.cell_id)):
        tr=d.iloc[it].copy();te=d.iloc[ie].copy()
        fit_tr,fit_te,ledger=nystrom_train_test(coords[it],coords[ie],
                                 scale_km=scale_km,nlandmarks=nlandmarks,seed=SEED+fold)
        spatial_receipts.append(ledger)
        for i in range(nlandmarks):
            col="spatial_kernel_feature_"+str(i)
            tr[col]=fit_tr[:,i];te[col]=fit_te[:,i]
        training_taxa=tr.inat_taxon_id.value_counts()
        eligible_in_test=np.flatnonzero(te.inat_taxon_id.isin(training_taxa.index).to_numpy())
        eligible[ie[eligible_in_test]]=True
        fold_report.append({
            "fold":fold,"n_test_original_source_photo_cell_records":len(ie),
            "n_eligible_test_original_photo_species_seen_in_other_regions":len(eligible_in_test),
            "n_original_photo_cells_heldout":int(te.cell_id.nunique())})
        for name,feature_cols in methods.items():
            if len(eligible_in_test):
                predicted[name][ie[eligible_in_test]]=train_predict(
                    tr,te,feature_cols,eligible_in_test)
    count=int(eligible.sum())
    if count<MIN_ELIGIBLE or len(set(source_labels[eligible]))<4:
        return {"status":"HOLD_INSUFFICIENT_OUTOFCELL_TRAIN_KNOWN_SOURCE_SPECIES",
                "n_test_original_photos":count,"fold_coverage":fold_report}
    if not all(np.isfinite(z[eligible]).all() for z in predicted.values()):
        raise RuntimeError("Original test species support varied across spatial/environment models")
    target=np.eye(4)[source_labels[eligible]]
    errors={k:np.square(p[eligible]-target).sum(axis=1) for k,p in predicted.items()}
    models={k:{"n_same_original_heldout_photo_rows":count,
               "mean_heldout_fourstate_brier":float(v.mean())} for k,v in errors.items()}
    comparisons={
        "spatial_kernel_beyond_linear_geography":("SPECIES_GEOGRAPHY_ONLY",
                                                  "SPECIES_GEOGRAPHY_SPATIAL_KERNEL"),
        "all_abiotic_beyond_species_real_geodesic_spatial_covariance":(
            "SPECIES_GEOGRAPHY_SPATIAL_KERNEL","SPECIES_GEOGRAPHY_SPATIAL_KERNEL_ALL_ENV"),
    }
    for block in BLOCKS:
        if soil or block!="soil":
            comparisons["unique_"+block+"_beyond_species_spatial_kernel_and_other_abiotic"]=(
                "SPATIAL_ALL_MINUS_"+block.upper(),"SPECIES_GEOGRAPHY_SPATIAL_KERNEL_ALL_ENV")
    groupings={"source_species":d.loc[eligible,"inat_taxon_id"].to_numpy(int),
               "source_geographical_cell":d.loc[eligible,"cell_id"].to_numpy(int)}
    gains={}
    for label,(a,b) in comparisons.items():
        loss=errors[a]-errors[b]
        intervals={}
        for i,(kind,group) in enumerate(groupings.items()):
            unique,inv=np.unique(group,return_inverse=True)
            sm=np.bincount(inv,weights=loss)
            nn=np.bincount(inv)
            rng=np.random.default_rng(SEED+len(label)+i)
            draw=rng.integers(0,len(unique),size=(nboot,len(unique)))
            samples=sm[draw].sum(axis=1)/nn[draw].sum(axis=1)
            intervals[kind]=[float(z) for z in np.quantile(samples,[.025,.975])]
        gains[label]={
            "outofregion_source_photo_brier_reduction":float(loss.mean()),
            "species_block_bootstrap_95CI":intervals["source_species"],
            "geographical_cell_block_bootstrap_95CI":intervals["source_geographical_cell"],
            "positive_in_both_fixed_prediction_group_resamplings":bool(
                intervals["source_species"][0]>0 and
                intervals["source_geographical_cell"][0]>0),
        }
    spatial_residual=residual_neighbour_autocorrelation(
        d.loc[eligible],source_labels[eligible],
        {k:predicted[k][eligible] for k in (
            "SPECIES_GEOGRAPHY_ONLY",
            "SPECIES_GEOGRAPHY_SPATIAL_KERNEL",
            "SPECIES_GEOGRAPHY_SPATIAL_KERNEL_ALL_ENV")},
    )
    return {
        "status":"SOURCE_SPECIES_INTERCEPT_GEODESIC_NYSTROM_SPATIAL_KERNEL_EXPLORATION",
        "n_original_source_fourstate_environment_complete":len(d),
        "n_original_heldout_source_photos_with_training_species":count,
        "n_distinct_original_source_species_in_heldout":int(d.loc[eligible,"inat_taxon_id"].nunique()),
        "n_heldout_source_geographic_cells":int(d.loc[eligible,"cell_id"].nunique()),
        "fixed_geodesic_covariance_scale_km":scale_km,
        "n_source_training_only_spatial_landmarks":nlandmarks,
        "approximation":"low-rank Nyström approximation of exp(-great_circle_distance/scale), outcome-blind TRAIN-only landmarks",
        "spatial_kernel_not_full_high_rank_GP":True,
        "source_colour_means_by_species_training_cells_only":True,
        "all_models_same_test_photo_ids_and_geocell_folds":True,
        "folds":fold_report,"spatial_landmark_receipts":spatial_receipts,
        "fourstate_model_scores":models,"environment_increment":gains,
        "fixed_prediction_nearby_residual_photo_spatial_Moran_diagnostic":spatial_residual,
        "confidence_intervals_conditional_on_fixed_crossvalidated_predictions":True,
    }


def run(taxon_cell:pd.DataFrame,original_sites:pd.DataFrame,*,strict:bool=True,
        nboot:int=BOOT,nlandmarks:int=N_LANDMARKS)->dict:
    merged,source=match_original(taxon_cell,original_sites,strict=strict)
    selections,coverage=populations(merged)
    output={}
    for cohort,soil in (("climate",False),("soil",True)):
        frame=selections["CLIMATE_AND_SOIL" if soil else "CLIMATE_SOURCE_ONLY"]
        output[cohort]={
            str(int(km)):fit_source(frame,soil=soil,scale_km=km,
                                    nboot=nboot,nlandmarks=nlandmarks)
            for km in SCALES_KM}
    return {
        "schema":SCHEMA,"date_jst":"2026-10-10",
        "status":"RETROSPECTIVE_SPECIES_INTERCEPT_LOW_RANK_REAL_PHOTO_SPATIAL_KERNEL_ABIOTIC",
        "original_source_repeated_photo_ledger":source,
        "original_source_env_coverage":coverage,
        "n_source_genus_or_tree_proxy_added_to_species_intercept":0,
        "source_time_invariant_phylogenetic_main_effect_absorbed_by_species_intercept":True,
        "all_original_photo_outcomes_and_source_sites_unchanged":True,
        "all_environment_blocks":{k:list(v) for k,v in BLOCKS.items()},
        "spatial_kernels_km":list(SCALES_KM),
        "results":output,
        "hard_nonclaims":[
            "Nyström is a finite low-rank approximation to geodesic spatial covariance; no exact full Gaussian spatial field is fitted",
            "Cross-location original taxon-cell photo colour is not measured inherited within-population pigment",
            "Species intercept absorbs species-invariant relatedness but species-dependent environment slopes could differ",
            "Group-bootstrap of fixed fit predictions is not a spatially constrained source label shuffle",
            "Residual Moran's I and geographic fine-scale matched null still need audit",
            "Original direct dated LCVP tip coverage is not sufficient for full 42111 species phylogenetic selection claim",
        ],
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-taxon-cell",type=Path,required=True)
    p.add_argument("--expanded-all-original-photo-sites",type=Path,required=True)
    p.add_argument("--outdir",type=Path,required=True)
    a=p.parse_args()
    d=run(pd.read_csv(a.original_taxon_cell,low_memory=False),
          pd.read_csv(a.expanded_all_original_photo_sites,low_memory=False))
    a.outdir.mkdir(parents=True,exist_ok=True)
    (a.outdir/"result.json").write_text(json.dumps(d,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"schema":SCHEMA,"source":d["original_source_env_coverage"],
          "results":{k:{scale:{"status":v["status"],"test_n":v.get("n_original_heldout_source_photos_with_training_species"),
            "gains":v.get("environment_increment")} for scale,v in record.items()}
            for k,record in d["results"].items()}},sort_keys=True),flush=True)


if __name__=="__main__":
    main()

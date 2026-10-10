#!/usr/bin/env python3
"""Original 85,337 FCP species×region flower photos: nonlinear + local spatial kernel.

Retrospective head-to-head source photo-colour prediction:
  M0: train-only species original colour baseline + spherical degree3 geographic basis
  M1: M0 + low-rank GREAT-CIRCLE spatial exponential covariance feature field,
      knots fit to TRAIN geographic sites only, no test colour/location fit.
  M2: M1 + seven WorldClim/SoilGrids environmental blocks on SAME photo cases.
  M3: M2 minus one of seven blocks (negative blocks all reported).

All 5 folds hold out entire original geographic source cells. A fixed species
baseline is estimated from OTHER regions, so a species-invariant LCVP effect
cannot separately enter this within-species comparison. Spatial Nyström-like
basis approximates an exponential kernel; a valid field basis, not an exact
infinite GP. Does not establish plant adaptation or genetically inherited colour.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from audit_fcp85337_species_nonlinear_spatial_abiotic_20261010 import (
    add_spatial_terms,SPATIAL,DEGREE
)
from compare_fcp_original85337_species_fixed_multiabiotic_20261010 import (
    match_original,populations,train_predict,
)
from compare_fcp_global42111_full_climate_sun_soil_blocks_20261010 import (
    BLOCKS,GEO,CLASSES,
)

SCHEMA="fcp_original_85337_species_fixed_geodesic_lowrank_space_all_environment_v1"
SPACE_SCALES_KM=(250.,1000.)
N_FIXED_TRAIN_KNOTS=96
FOLDS=5
MIN_EVAL=200
BOOT=199
SEED=20261010
EARTH_RADIUS=6371.0088
RBF=tuple(f"spatial_trainonly_rbf_{i:03d}" for i in range(N_FIXED_TRAIN_KNOTS))


def unit_sphere(lat:np.ndarray,lon:np.ndarray)->np.ndarray:
    lat=np.asarray(lat,float);lon=np.asarray(lon,float)
    if lat.ndim!=1 or lon.ndim!=1 or len(lat)!=len(lon) or not np.isfinite(lat).all() or not np.isfinite(lon).all():
        raise ValueError("Original photo latitude/longitude must be finite and source-identical")
    if (abs(lat)>90).any() or (abs(lon)>180).any():
        raise ValueError("Illegal original photo geographical latitude/longitude")
    phi=np.deg2rad(lat);theta=np.deg2rad(lon)
    return np.column_stack((np.cos(phi)*np.cos(theta),
                            np.cos(phi)*np.sin(theta),np.sin(phi)))


def fit_training_geo_knots(train:pd.DataFrame,n_knots:int=N_FIXED_TRAIN_KNOTS)->np.ndarray:
    from sklearn.cluster import MiniBatchKMeans
    if n_knots<2 or len(train)<n_knots:
        raise ValueError("Insufficient original TRAIN photographed source sites for spatial basis")
    xyz=unit_sphere(train.latitude.to_numpy(float),train.longitude.to_numpy(float))
    # Unsupervised training-coordinate density approximation; NEVER uses old
    # photo colours, training species labels, or heldout source coordinates.
    model=MiniBatchKMeans(n_clusters=n_knots,random_state=SEED,
                          batch_size=1024,n_init=3,max_iter=120)
    model.fit(xyz)
    out=model.cluster_centers_.astype(float)
    norm=np.linalg.norm(out,axis=1)
    if (norm<=1e-8).any() or not np.isfinite(out).all():
        raise RuntimeError("Unusable train-only photographed geographical centres")
    return out/norm[:,None]


def rbf_columns(d:pd.DataFrame,knots:np.ndarray,scale_km:float,
                names:tuple[str,...]=RBF)->pd.DataFrame:
    xyz=unit_sphere(d.latitude.to_numpy(float),d.longitude.to_numpy(float))
    if not np.isfinite(knots).all() or knots.shape!=(len(names),3):
        raise ValueError("Original training-only geographical basis dimensions drifted")
    if scale_km<=0:raise ValueError("Source spatial distance bandwidth must be positive")
    dots=np.clip(xyz@knots.T,-1.,1.)
    distance=EARTH_RADIUS*np.arccos(dots)
    phi=np.exp(-distance/scale_km)
    if (phi<0).any() or (phi>1.000001).any() or not np.isfinite(phi).all():
        raise RuntimeError("Invalid original-photo geodesic spatial covariance feature")
    out=d.copy()
    out.loc[:,list(names)]=phi
    return out


def candidate_models(*,soil:bool)->dict[str,tuple[str,...]]:
    blocks={k:v for k,v in BLOCKS.items() if soil or k!="soil"}
    physical=tuple(v for x in blocks.values() for v in x)
    smooth=GEO+SPATIAL
    spatial=smooth+RBF
    families={
        "SPECIES_NONGP_COORDINATES":smooth,
        "SPECIES_GEODESIC_SPATIAL_FIELD":spatial,
        "SPECIES_GEODESIC_SPATIAL_FIELD_ALL_ENV":spatial+physical,
    }
    for name,v in blocks.items():
        families["SPATIAL_FIELD_FULL_MINUS_"+name.upper()]=spatial+tuple(
            f for f in physical if f not in v)
    return families


def evaluate(d:pd.DataFrame,*,soil:bool,bandwidth:float,nboot:int=BOOT,
             n_knots:int=N_FIXED_TRAIN_KNOTS)->dict:
    from sklearn.model_selection import GroupKFold
    input_frame=add_spatial_terms(d)
    if len(input_frame)<MIN_EVAL or input_frame.cell_id.nunique()<FOLDS:
        return {"status":"HOLD_MISSING_REGION_FIXED_SOURCE_PHOTO_OPPORTUNITY"}
    y=pd.Categorical(input_frame.morph,categories=CLASSES).codes
    if (y<0).any():
        raise ValueError("Original source unclassifiable photo entered flower outcome")
    model_names=candidate_models(soil=soil)
    if n_knots!=len(RBF):
        raise ValueError("Fixed spatial RBF basis dimension must not be tuned after outcome")
    pred={k:np.full((len(input_frame),4),np.nan) for k in model_names}
    eligible=np.zeros(len(input_frame),bool)
    receipts=[]
    for fold,(train_index,test_index) in enumerate(GroupKFold(n_splits=FOLDS).split(
        input_frame,y,groups=input_frame.cell_id)):
        tr=input_frame.iloc[train_index]
        te=input_frame.iloc[test_index]
        source_train_species=tr.inat_taxon_id.value_counts()
        source_eval=np.flatnonzero(te.inat_taxon_id.isin(source_train_species.index).to_numpy())
        eligible[test_index[source_eval]]=True
        n=len(source_eval)
        receipts.append({
            "fold":fold,
            "n_region_heldout_original_photo_records":int(len(test_index)),
            "n_species_train_identifiable_original_photo_records":int(n),
            "n_original_independent_source_photo_cells":int(te.cell_id.nunique()),
        })
        if n==0:continue
        coords=fit_training_geo_knots(tr,n_knots=n_knots)
        # Same TRAIN-defined knot positions and same original test photo IDs
        # for every compared predictor, not a post-outcome spatial selection.
        train=rbf_columns(tr,coords,bandwidth)
        test=rbf_columns(te,coords,bandwidth)
        for model,fields in model_names.items():
            pred[model][test_index[source_eval]]=train_predict(
                train,test,fields,source_eval)
    count=int(eligible.sum())
    if count<MIN_EVAL or len(set(y[eligible]))<4:
        return {"status":"HOLD_MISSING_REPLICATED_SOURCE_SPECIES_PHOTO_COLOUR_CLASSES",
                "n_original_evaluable_heldout_photos":count}
    if any(not np.isfinite(v[eligible]).all() for v in pred.values()):
        raise RuntimeError("Original source species and geographic test photos changed between model families")
    yy=np.eye(4)[y[eligible]]
    loss={k:np.square(v[eligible]-yy).sum(axis=1) for k,v in pred.items()}
    compare={
        "distance_covariance_basis_beyond_nonlinear_geographic_trends":(
            "SPECIES_NONGP_COORDINATES","SPECIES_GEODESIC_SPATIAL_FIELD"),
        "all_abiotic_beyond_species_nonlinear_space_and_geodesic_covariance":(
            "SPECIES_GEODESIC_SPATIAL_FIELD","SPECIES_GEODESIC_SPATIAL_FIELD_ALL_ENV"),
    }
    blocks={k:v for k,v in BLOCKS.items() if soil or k!="soil"}
    for k in blocks:
        compare["unique_"+k+"_beyond_species_geodesic_spatial_and_other_abiotic"]=(
            "SPATIAL_FIELD_FULL_MINUS_"+k.upper(),
            "SPECIES_GEODESIC_SPATIAL_FIELD_ALL_ENV")
    holdout_species=input_frame.loc[eligible,"inat_taxon_id"].to_numpy(int)
    holdout_cell=input_frame.loc[eligible,"cell_id"].to_numpy(int)
    gains={}
    for name,(before,after) in compare.items():
        delta=loss[before]-loss[after]
        intervals={}
        for label,group in (("source_species",holdout_species),("original_geographic_cell",holdout_cell)):
            unique,ids=np.unique(group,return_inverse=True)
            sum_loss=np.bincount(ids,weights=delta)
            size=np.bincount(ids)
            rng=np.random.default_rng(SEED+len(name)+(0 if label=="source_species" else 1))
            draws=rng.integers(0,len(unique),size=(nboot,len(unique)))
            vals=sum_loss[draws].sum(axis=1)/size[draws].sum(axis=1)
            intervals[label]=[float(x) for x in np.quantile(vals,[.025,.975])]
        gains[name]={
            "heldout_multiclass_brier_reduction":float(delta.mean()),
            "source_species_cluster_95CI":intervals["source_species"],
            "original_region_cell_cluster_95CI":intervals["original_geographic_cell"],
            "both_intervals_positive":bool(intervals["source_species"][0]>0 and
                                      intervals["original_geographic_cell"][0]>0),
        }
    return {
        "status":"SOURCE_SPECIES_INTERCEPT_AND_GEODESIC_LOW_RANK_SPATIAL_FIELD_EXPLORATION",
        "original_source_observed_photo_records":len(d),
        "n_train_species_heldout_original_photos":count,
        "n_distinct_source_taxa_in_evaluation":int(input_frame.loc[eligible,"inat_taxon_id"].nunique()),
        "n_original_source_geographic_cells":int(input_frame.loc[eligible,"cell_id"].nunique()),
        "spatial_exponential_bandwidth_km":bandwidth,
        "spatial_inducing_knots_count":n_knots,
        "training_only_unsupervised_source_photo_geographic_knots":True,
        "spatial_covariance_interpretation":"finite-rank PSD approximation to geodesic exponential kernel; NOT full GP posterior/complete spatial deconfounding",
        "species_means_learned_only_from_other_source_geographic_regions":True,
        "every_source_photo_same_heldout_fold_and_model_comparison":True,
        "source_fourstate_photo_labels_never_reclassified":True,
        "fixed_source_cell_group_cv":receipts,
        "model_source_heldout_brier":{k:float(x.mean()) for k,x in loss.items()},
        "source_block_predictive_gains":gains,
        "group_cluster_bootstrap_conditional_on_fixed_spatial_knot_and_model_fits":True,
    }


def run(original_cells:pd.DataFrame,full_sites:pd.DataFrame,*,strict=True,nboot:int=BOOT)->dict:
    merged,source=match_original(original_cells,full_sites,strict=strict)
    populations_by_missingness,cover=populations(merged)
    modes={
        "climate_without_soil_selection":(populations_by_missingness["CLIMATE_SOURCE_ONLY"],False),
        "climate_and_soil_complete":(populations_by_missingness["CLIMATE_AND_SOIL"],True),
    }
    outputs={}
    for label,(d,soil) in modes.items():
        outputs[label]={
            str(int(scale))+"km":evaluate(d,soil=soil,bandwidth=scale,nboot=nboot)
            for scale in SPACE_SCALES_KM
        }
    return {
        "schema":SCHEMA,"date_jst":"2026-10-10",
        "status":"RETROSPECTIVE_SOURCE_SAME_SPECIES_WITH_GEO_DISTANCE_SPATIAL_FIELD_ENVIRONMENT",
        "original_85337_taxon_cell_source":source,
        "original_42111_species_environment_coverage":cover,
        "source_photo_site_spatial_kernel_scales_km":list(SPACE_SCALES_KM),
        "source_train_only_knot_count":N_FIXED_TRAIN_KNOTS,
        "source_spherical_spatial_polynomial_degree":DEGREE,
        "n_environment_blocks_separately_tested":len(BLOCKS),
        "source_original_photo_outcomes_unchanged":True,
        "source_original_independently_reserved_taxa_untouched":True,
        "source_species_invariant_phylogenetic_component_absorbed_in_species_intercepts":True,
        "source_spatial_field_is_low_rank_not_full_GP":True,
        "results":outputs,
        "hard_nonclaims":[
            "A distance-based finite-rank spatial covariance field is not complete causal spatial deconfounding",
            "Photo-classifiable original taxon-cell observations are source-observational, not genetic population morph proportions",
            "Source genus and species identity/photographer confounding and distant region sampling remain",
            "All seven correlated environment groups are post hoc and original spatial bandwidths are sensitivity assumptions",
            "A positive Brier increment after a spatial field does not imply environmental adaptation or selection",
            "Species intercept absorbs species-invariant phylogeny; it is not a separately identified true LCVP covariance effect",
            "An independent, repeated-genotype and fine-spatial residual Moran/photo-opportunity validation remains missing",
        ],
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-taxon-cell",required=True,type=Path)
    p.add_argument("--all-original-expanded-photo-sites",required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    a=p.parse_args()
    old=pd.read_csv(a.original_taxon_cell,low_memory=False)
    src=pd.read_csv(a.all_original_expanded_photo_sites,low_memory=False)
    out=run(old,src)
    a.outdir.mkdir(parents=True,exist_ok=True)
    (a.outdir/"result.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "schema":SCHEMA,"coverage":out["original_42111_species_environment_coverage"],
        "models":{k:{scale:r["source_block_predictive_gains"] if r["status"].startswith("SOURCE_") else {"hold":r["status"]}
                    for scale,r in v.items()} for k,v in out["results"].items()},
    },sort_keys=True),flush=True)


if __name__=="__main__":
    main()

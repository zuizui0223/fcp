#!/usr/bin/env python3
"""Within-species global photographed colour beyond nonlinear spatial trends.

Use original 85,337 taxon×cell photos; keep photo identity and source-colour
status immutable. True original photo lat/lon only for degree-1..3 spherical
polynomial spatial basis, with NO cell-centroid environmental substitution.
One species' four-colour baseline is learned ONLY from other training cells.
Compare source species+linear geography, +flexible spherical spatial basis,
then + all abiotic data and each drop-one ecological block on IDENTICAL source
photo IDs, training regions and complete cases. Exploratory, noncausal; a
spatial polynomial is not full spatial covariance or genetic evolution.
"""
from __future__ import annotations
import argparse
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd

from compare_fcp_original85337_species_fixed_multiabiotic_20261010 import (
    match_original,populations,train_predict,BOOT as ORIGINAL_BOOT
)
from compare_fcp_global42111_full_climate_sun_soil_blocks_20261010 import (
    BLOCKS,GEO,CLASSES
)

SCHEMA="fcp_original_85337_species_fixed_nonlinear_spatial_environment_v1"
DEGREE=3
FOLDS=5
SEED=20261010
BOOT=199
MIN_TEST=200
SPATIAL=tuple(f"sphere_monomial_{k}" for k in range(1,20))


def add_spatial_terms(source:pd.DataFrame)->pd.DataFrame:
    if not {"latitude","longitude","site_geo_status"}.issubset(source):
        raise ValueError("Real source photographed geographical positions unavailable")
    d=source.copy()
    latitude=pd.to_numeric(d.latitude,errors="coerce").to_numpy(float)
    longitude=pd.to_numeric(d.longitude,errors="coerce").to_numpy(float)
    lat=np.deg2rad(latitude);lon=np.deg2rad(longitude)
    X=np.column_stack([np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)])
    n=0
    for degree in range(1,DEGREE+1):
        for inds in itertools.combinations_with_replacement(range(3),degree):
            n+=1
            d[f"sphere_monomial_{n}"]=np.prod(X[:,list(inds)],axis=1)
    if n!=len(SPATIAL):
        raise RuntimeError("Source spatial polynomial basis has wrong frozen dimension")
    bad=d.site_geo_status.ne("VALID_PUBLIC_ORIGINAL_PHOTO_POINT")
    d.loc[bad,list(SPATIAL)]=np.nan
    if d.loc[bad,list(SPATIAL)].notna().any().any():
        raise RuntimeError("Source no-public-photo-location entered spatial basis")
    return d


def families(*,with_soil:bool)->dict[str,tuple[str,...]]:
    blocks={k:tuple(v) for k,v in BLOCKS.items() if with_soil or k!="soil"}
    environment=tuple(i for q in blocks.values() for i in q)
    controls=GEO+SPATIAL
    out={
        "SPECIES_GEO":GEO,
        "SPECIES_GEO_SPHERICAL_3":controls,
        "SPECIES_GEO_SPHERICAL_3_ALL_ENV":controls+environment,
    }
    for name,cols in blocks.items():
        out["SPATIAL_FULL_MINUS_"+name.upper()]=controls+tuple(
            q for q in environment if q not in cols)
    return out


def spatial_species_heldout(d:pd.DataFrame,*,soil:bool,nboot:int=BOOT)->dict:
    from sklearn.model_selection import GroupKFold
    frame=add_spatial_terms(d)
    required=set(GEO+SPATIAL+tuple(v for k,b in BLOCKS.items() if soil or k!="soil" for v in b))
    if not required.issubset(frame):
        raise ValueError("Source original environmental feature missing")
    if len(frame)<MIN_TEST or frame.cell_id.nunique()<FOLDS:
        return {"status":"HOLD_SOURCE_WITHIN_SPECIES_REGION_GROUP_SUPPORT"}
    y=pd.Categorical(frame.morph,categories=CLASSES).codes
    if (y<0).any():
        raise ValueError("Unclassifiable flower photograph leaked into source assay")
    methods=families(with_soil=soil)
    eligible=np.zeros(len(frame),bool)
    predictions={key:np.full((len(frame),4),np.nan) for key in methods}
    folds=[]
    for fold,(itr,ite) in enumerate(GroupKFold(n_splits=FOLDS).split(
        frame,y,groups=frame.cell_id)):
        train=frame.iloc[itr];test=frame.iloc[ite]
        source_species=train.inat_taxon_id.value_counts()
        original_photo_ids=np.flatnonzero(test.inat_taxon_id.isin(source_species.index).to_numpy())
        eligible[ite[original_photo_ids]]=True
        folds.append({"fold":fold,"n_source_heldout_photos":len(ite),
                      "n_original_photo_species_seen_in_other_cells":len(original_photo_ids),
                      "n_original_heldout_cells":int(test.cell_id.nunique())})
        for name,cols in methods.items():
            if len(original_photo_ids):
                predictions[name][ite[original_photo_ids]]=train_predict(
                    train,test,cols,original_photo_ids)
    n=int(eligible.sum())
    if n<MIN_TEST or len(set(y[eligible]))<4:
        return {"status":"HOLD_INSUFFICIENT_EVALUABLE_SOURCE_SPECIES_CROSS_REGION_PHOTOS",
                "n_original_photo_species_training_intercept_heldout":n}
    if not all(np.isfinite(proba[eligible]).all() for proba in predictions.values()):
        raise RuntimeError("Different source photo test support among spatial/environment models")
    target=np.eye(4)[y[eligible]]
    errors={name:np.square(value[eligible]-target).sum(axis=1)
            for name,value in predictions.items()}
    comparison={
        "nonlinear_space_over_linear_geography":("SPECIES_GEO","SPECIES_GEO_SPHERICAL_3"),
        "all_environment_beyond_species_plus_nonlinear_space":(
            "SPECIES_GEO_SPHERICAL_3","SPECIES_GEO_SPHERICAL_3_ALL_ENV"),
    }
    for name in BLOCKS:
        if name!="soil" or soil:
            comparison["unique_"+name+"_beyond_species_nonlinear_space_other_environment"]=(
                "SPATIAL_FULL_MINUS_"+name.upper(),"SPECIES_GEO_SPHERICAL_3_ALL_ENV")
    indices={"original_species":frame.loc[eligible,"inat_taxon_id"].to_numpy(int),
             "original_photo_cell":frame.loc[eligible,"cell_id"].to_numpy(int)}
    gains={}
    for name,(base,full) in comparison.items():
        dloss=errors[base]-errors[full]
        intervals={}
        for level,groups in indices.items():
            unique,group_idx=np.unique(groups,return_inverse=True)
            group_sums=np.bincount(group_idx,weights=dloss)
            group_counts=np.bincount(group_idx)
            rng=np.random.default_rng(SEED+len(name)+(0 if level=="original_species" else 1))
            draws=rng.integers(0,len(unique),size=(nboot,len(unique)))
            vals=group_sums[draws].sum(axis=1)/group_counts[draws].sum(axis=1)
            intervals[level]=[float(q) for q in np.quantile(vals,[.025,.975])]
        gains[name]={
            "heldout_multiclass_brier_reduction":float(dloss.mean()),
            "source_species_cluster_95CI":intervals["original_species"],
            "source_cell_cluster_95CI":intervals["original_photo_cell"],
            "positive_in_both_conditional_resamplings":bool(
                intervals["original_species"][0]>0 and intervals["original_photo_cell"][0]>0),
        }
    return {
        "status":"SOURCE_SPECIES_FIXED_WITH_NONLINEAR_SPATIAL_BASIS_ABIOTIC_EXPLORATION",
        "n_original_classified_abiotic_complete":len(frame),
        "n_original_outofcell_photo_species_heldout":n,
        "n_distinct_original_species_heldout":int(frame.loc[eligible,"inat_taxon_id"].nunique()),
        "n_heldout_geographic_cells":int(frame.loc[eligible,"cell_id"].nunique()),
        "spatial_basis":"up_to_degree_3_polynomial_of_true_photo_unit_sphere_xyz",
        "spatial_basis_dimensions":len(SPATIAL),
        "true_spatial_covariance_modelled":False,
        "all_model_source_test_photos_and_fold_assignments_identical":True,
        "source_photo_colour_unchanged":True,
        "five_source_region_blocked_folds":folds,
        "model_four_colour_brier_scores":{
            k:{"score":float(v.mean()),"n_same_test_photos":n}
            for k,v in errors.items()},
        "incremental_predictive_gains":gains,
        "estimation_limit":"Species means only from other source regions, spherical polynomial trends are not a spatial random field; environmental gain is predictive not adaptive or genetic",
    }


def run(taxon_cell:pd.DataFrame,original_sites:pd.DataFrame,*,strict=True,nboot=BOOT)->dict:
    merged,source=match_original(taxon_cell,original_sites,strict=strict)
    cohorts,coverage=populations(merged)
    climate=spatial_species_heldout(cohorts["CLIMATE_SOURCE_ONLY"],soil=False,nboot=nboot)
    full=spatial_species_heldout(cohorts["CLIMATE_AND_SOIL"],soil=True,nboot=nboot)
    return {
        "schema":SCHEMA,"date_jst":"2026-10-10",
        "status":"RETROSPECTIVE_ORIGINAL_PHOTO_SPECIES_AND_NONLINEAR_GEOGRAPHY_ENVIRONMENT_AUDIT",
        "original_source":source,
        "original_environment_complete_case_coverage":coverage,
        "source_7_environment_blocks":{k:list(v) for k,v in BLOCKS.items()},
        "spatial_polynomial_max_degree":DEGREE,
        "n_spatial_polynomial_columns":len(SPATIAL),
        "climate_without_soil_selection":climate,
        "climate_and_soil_complete":full,
        "all_original_photo_ids_and_flower_colour_labels_preserved":True,
        "source_phylogenetic_species_invariant_component_absorbed_by_species_intercepts":True,
        "genuine_spatial_autocorrelation_random_field_not_estimable_from_this_basis_alone":True,
        "old_full_original_42111_dated_phylogeny_claim_HOLD":True,
        "nonclaims":[
            "Spherical polynomial spatial basis is a finite global smooth trend, not flexible geodesic spatial covariance or a residual Moran test",
            "An original nominal species photographic intercept absorbs species-constant phylogeny but is not a genetic assay",
            "One photographed flower per species-cell does not estimate allelic polymorphism or phenotypic plasticity",
            "Source 162-cell heldout folds and genus-level photography selection can still bias any reported environmental gain",
            "A positive all-environment gain in source repeated photo data is not a causal ecological or evolutionary adaptation estimate",
        ],
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-taxon-cell",type=Path,required=True)
    p.add_argument("--all-original-expanded-photo-sites",type=Path,required=True)
    p.add_argument("--outdir",type=Path,required=True)
    a=p.parse_args()
    old=pd.read_csv(a.original_taxon_cell,low_memory=False)
    full=pd.read_csv(a.all_original_expanded_photo_sites,low_memory=False)
    result=run(old,full)
    a.outdir.mkdir(parents=True,exist_ok=True)
    (a.outdir/"result.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "schema":SCHEMA,"coverage":result["original_environment_complete_case_coverage"],
        "climate":result["climate_without_soil_selection"].get("incremental_predictive_gains"),
        "soil":result["climate_and_soil_complete"].get("incremental_predictive_gains"),
    },sort_keys=True),flush=True)

if __name__=="__main__":
    main()

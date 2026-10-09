#!/usr/bin/env python3
"""Photo classifiability selection-control for global FCP colour-climate ecology.

Outcome is ONLY whether the original September 2026 flower ROI produced a
classified four-state photo. Includes the unclassified original 23,654 taxa,
rather than silently dropping them. Spatially-heldout train-only genus
intercepts and simple fixed geographic/climate/soil predictor blocks.
This is photo-ascertainment ecology, not a test of biological flower colours.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd

PHOTO_STATES=("white","yellow_orange","red_pink","blue_purple")
GEO=("abs_latitude","lon_sin","lon_cos","wc_elevation_m")
TEMPERATURE=("wc_bio1","wc_bio5")
MOISTURE=("wc_bio12","wc_bio15")
SOIL=("soil_pH","soil_SOC","soil_N","soil_clay","soil_available_water_proxy")
FAMILIES={
    "GENUS_ONLY":(),
    "GENUS_GEO":GEO,
    "GENUS_GEO_TEMP":GEO+TEMPERATURE,
    "GENUS_GEO_MOISTURE":GEO+MOISTURE,
    "GENUS_GEO_ALL_CLIMATE":GEO+TEMPERATURE+MOISTURE,
    "GENUS_GEO_ALL_CLIMATE_SOIL":GEO+TEMPERATURE+MOISTURE+SOIL,
}
EXPECTED_SOURCE=42111
EXPECTED_CLASSIFIABLE=18457
EXPECTED_CLIMATE_COVER=42014
EXPECTED_SOIL_COVER=32231
FOLDS=5
RIDGE=10.0
BOOT=999
SEED=20261009


def source_sets(original:pd.DataFrame,*,strict:bool=True)->tuple[dict[str,pd.DataFrame],dict]:
    columns={"inat_taxon_id","species","morph","measurement_status","site_geo_status",
             "latitude","longitude","environment_climate_complete",
             "environment_soil_complete",*GEO[3:],*TEMPERATURE,*MOISTURE,*SOIL}
    if not columns.issubset(original):
        raise ValueError("Source original 42111 photo classification data missing "+str(sorted(columns-set(original))))
    if original.inat_taxon_id.duplicated().any():
        raise ValueError("Original nominal species replicated")
    if strict and len(original)!=EXPECTED_SOURCE:
        raise ValueError("Full source 42111 taxa missing")
    z=original.copy()
    z["classifiable"]=(z.measurement_status.eq("classified_four_state_morph")&z.morph.isin(PHOTO_STATES)).astype(int)
    if strict and int(z.classifiable.sum())!=EXPECTED_CLASSIFIABLE:
        raise ValueError("Original 18457 source photographic classifications altered")
    latitude=pd.to_numeric(z.latitude,errors="coerce")
    longitude=pd.to_numeric(z.longitude,errors="coerce")
    z["abs_latitude"]=latitude.abs()
    z["lon_sin"]=np.sin(np.deg2rad(longitude))
    z["lon_cos"]=np.cos(np.deg2rad(longitude))
    good=z.site_geo_status.eq("VALID_PUBLIC_ORIGINAL_PHOTO_POINT")&latitude.between(-90,90)&longitude.between(-180,180)
    # Standalone original photo geometry; no dependency on a pytest-named
    # analysis module that could be confused with the test harness.
    latv=latitude.to_numpy(float); lonv=longitude.to_numpy(float)
    source_cells=np.full(len(z),-1,dtype=int)
    finite=np.isfinite(latv)&np.isfinite(lonv)&(np.abs(latv)<=90)&(np.abs(lonv)<=180)
    y=np.clip(np.floor((np.sin(np.deg2rad(latv[finite]))+1)*4.5).astype(int),0,8)
    x=np.clip(np.floor((lonv[finite]+180)/20).astype(int),0,17)
    source_cells[finite]=18*y+x
    z["photo_cell_162"]=source_cells
    z["genus"]=z.species.astype(str).str.split().str[0]
    climate=(good&z.environment_climate_complete.fillna(False).astype(bool)&
             z[list(GEO+TEMPERATURE+MOISTURE)].notna().all(axis=1)&z.photo_cell_162.ge(0))
    soil=(climate&z.environment_soil_complete.fillna(False).astype(bool)&z[list(SOIL)].notna().all(axis=1))
    pools={"ALL_ORIGINAL_CLIMATE_PHOTO_OPPORTUNITY":z.loc[climate].copy(),
           "ALL_ORIGINAL_SOIL_PHOTO_OPPORTUNITY":z.loc[soil].copy()}
    stats={"original_source_taxa":len(z),"n_classifiable_original_photos":int(z.classifiable.sum()),
           "n_unclassifiable_original_photos":int((1-z.classifiable).sum()),
           "n_all_photos_with_original_climate":int(climate.sum()),
           "n_all_photos_with_original_climate_soil":int(soil.sum()),
           "n_original_unclassified_with_climate":int((climate&z.classifiable.eq(0)).sum()),
           "n_original_unclassified_with_climate_soil":int((soil&z.classifiable.eq(0)).sum())}
    if strict and (stats["n_all_photos_with_original_climate"]!=EXPECTED_CLIMATE_COVER or
                   stats["n_all_photos_with_original_climate_soil"]!=EXPECTED_SOIL_COVER):
        raise RuntimeError("Source 42014/32231 original observed climate/soil completeness drifted")
    for q,d in pools.items():
        if not d.classifiable.isin([0,1]).all() or d.genus.eq("").any() or d.inat_taxon_id.duplicated().any():
            raise ValueError("Missing source classifiability or genus")
    return pools,stats


def predict_train_only_genus(train:pd.DataFrame,test:pd.DataFrame,
                             feature_names:tuple[str,...],eligible:np.ndarray)->np.ndarray:
    # Same fixed original held-out photographed species in EVERY model.
    gtrain=train.genus.to_numpy(str)
    train_status=train.classifiable.to_numpy(float)
    yavg=pd.DataFrame({"genus":gtrain,"classifiable":train_status}).groupby("genus").classifiable.mean()
    te=test.iloc[eligible]
    yg=yavg.reindex(te.genus.to_numpy(str)).to_numpy(float)
    if len(feature_names)==0:return yg
    xt=train[list(feature_names)].to_numpy(float)
    xe=te[list(feature_names)].to_numpy(float)
    if not np.isfinite(xt).all() or not np.isfinite(xe).all():
        raise ValueError("Classifiability source soil/temperature photo values unavailable")
    mean=xt.mean(axis=0)
    sd=xt.std(axis=0)
    sd[sd<1e-8]=1
    xt=(xt-mean)/sd;xe=(xe-mean)/sd
    xm=pd.DataFrame(xt,columns=list(feature_names)).assign(genus=gtrain).groupby("genus")[list(feature_names)].mean()
    resid_x=xt-xm.reindex(gtrain).to_numpy(float)
    resid_y=train_status-yavg.reindex(gtrain).to_numpy(float)
    beta=np.linalg.solve(resid_x.T@resid_x+RIDGE*np.eye(len(feature_names)),resid_x.T@resid_y)
    return np.clip(yg+(xe-xm.reindex(te.genus.to_numpy(str)).to_numpy(float))@beta,0,1)


def blocked_classifiability(data:pd.DataFrame,*,with_soil:bool)->dict:
    from sklearn.model_selection import GroupKFold
    families={k:v for k,v in FAMILIES.items() if with_soil or k!="GENUS_GEO_ALL_CLIMATE_SOIL"}
    n=len(data)
    if n<500 or data.photo_cell_162.nunique()<FOLDS:
        raise ValueError("Not enough source geographic photographic opportunity")
    eligible=np.zeros(n,bool)
    prediction={key:np.full(n,np.nan) for key in families}
    folds=[]
    for i,(train,test) in enumerate(GroupKFold(n_splits=FOLDS).split(data,groups=data.photo_cell_162)):
        tr=data.iloc[train];te=data.iloc[test]
        freq=tr.genus.value_counts()
        elig=np.flatnonzero(te.genus.isin(freq[freq>=2].index).to_numpy())
        eligible[test[elig]]=True
        for name,cols in families.items():
            prediction[name][test[elig]]=predict_train_only_genus(tr,te,cols,elig)
        folds.append({"fold":i,"n_original_photos_heldout":len(test),
                      "n_train_supported_genus_original_photos":len(elig),
                      "n_missing_genus_training_support":int(len(test)-len(elig))})
    if int(eligible.sum())<100:
        raise RuntimeError("No within-genus photograph ascertainment support")
    true=data.loc[eligible,"classifiable"].to_numpy(float)
    errors={}
    summary={}
    for name,pred in prediction.items():
        score=(pred[eligible]-true)**2
        errors[name]=score
        summary[name]={"n_heldout_original_photos":int(eligible.sum()),
                       "mean_binary_classifiability_brier":float(score.mean())}
    groups=data.loc[eligible,"photo_cell_162"].to_numpy(int)
    cells,idx=np.unique(groups,return_inverse=True)
    counts=np.bincount(idx)
    effects={}
    comparisons={
        "geo_beyond_genus":("GENUS_ONLY","GENUS_GEO"),
        "climate_beyond_genus_geography":("GENUS_GEO","GENUS_GEO_ALL_CLIMATE"),
        "unique_temperature_beyond_genus_geography_moisture":("GENUS_GEO_MOISTURE","GENUS_GEO_ALL_CLIMATE"),
        "unique_moisture_beyond_genus_geography_temperature":("GENUS_GEO_TEMP","GENUS_GEO_ALL_CLIMATE"),
    }
    if with_soil:
        comparisons["soil_beyond_genus_geography_climate"]=("GENUS_GEO_ALL_CLIMATE","GENUS_GEO_ALL_CLIMATE_SOIL")
    for name,(base,added) in comparisons.items():
        delta=errors[base]-errors[added]
        by_cell=np.bincount(idx,weights=delta)
        rng=np.random.default_rng(SEED+len(name))
        s=rng.integers(0,len(cells),size=(BOOT,len(cells)))
        vals=by_cell[s].sum(axis=1)/counts[s].sum(axis=1)
        effects[name]={
            "heldout_classifiability_brier_improvement":float(delta.mean()),
            "original_geographic_cell_bootstrap_95CI":[float(v) for v in np.quantile(vals,[.025,.975])],
            "positive_gain_in_fixed_cell_resampling":bool(np.quantile(vals,.025)>0),
        }
    return {"n_source_original_photo_candidates":n,
            "n_classifiable_source_photos_in_candidate_pool":int(data.classifiable.sum()),
            "n_heldout_original_photos_genus_estimable":int(eligible.sum()),
            "n_source_genera_estimated":int(data.loc[eligible,"genus"].nunique()),
            "n_source_heldout_equal_area_cells":int(len(cells)),
            "folds":folds,
            "genus_means_estimated_on_training_photograph_opportunities_only":True,
            "same_original_source_photos_in_every_compared_model":True,
            "models":summary,"increments":effects,
            "outcome":"whether original photograph flower ROI can be classified, not colour of the flower"}


def run(source:pd.DataFrame,*,strict:bool=True)->dict:
    pools,coverage=source_sets(source,strict=strict)
    return {
        "schema":"fcp_global42111_colour_classifiability_moisture_negative_control_v1",
        "date_jst":"2026-10-09",
        "status":"PHOTOGRAPHIC_CLASSIFICATION_SELECTION_DIAGNOSTIC_NOT_FLOWER_COLOUR_ASSOCIATION",
        "original_source":coverage,
        "climate_opportunity":blocked_classifiability(pools["ALL_ORIGINAL_CLIMATE_PHOTO_OPPORTUNITY"],with_soil=False),
        "soil_opportunity":blocked_classifiability(pools["ALL_ORIGINAL_SOIL_PHOTO_OPPORTUNITY"],with_soil=True),
        "source_four_colour_labels_not_remeasured":True,
        "original_42111_taxa_and_23654_unclassifiable_denominators_preserved":True,
        "never_treat_unclassifiable_photographs_as_white_or_genetic_monomorphs":True,
        "hard_nonclaims":[
            "This predicts PHOTO classification success/failure, not whether an unknown natural plant lacks flower colour",
            "A precipitation signal in classifiability can expose selective photographic missingness but cannot quantify the induced bias by itself",
            "Genus and geographic blocking leave remaining phenology, photographer, flowering morphology and taxonomic confounding",
            "Environmental source variables are long-term climate, not conditions when the image was taken",
            "No genetic flower-colour variants, fitness, pigment biosynthesis or causal ecology is observed",
        ],
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-species-abiotic",type=Path,required=True)
    p.add_argument("--outdir",type=Path,required=True)
    a=p.parse_args()
    source=pd.read_csv(a.original_species_abiotic,low_memory=False)
    result=run(source)
    a.outdir.mkdir(parents=True,exist_ok=True)
    (a.outdir/"result.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,indent=2,sort_keys=True),flush=True)
if __name__=="__main__":
    main()

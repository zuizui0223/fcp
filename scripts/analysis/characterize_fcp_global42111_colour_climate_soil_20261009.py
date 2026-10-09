#!/usr/bin/env python3
"""Global 42111-species ecology: fixed-source flower-colour x climate-soil models.

Only source photograph outcomes; species-equal single breadth photo per taxon.
Three prespecified feature families, same complete-case species and identical
folds. Genus-held-out and 162 equal-area cell-held-out prediction checks.
NO causal climate, soil, pollinator or adaptation claims.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

COLOURS=("white","yellow_orange","red_pink","blue_purple")
GEO=["abs_latitude","lon_sin","lon_cos","wc_elevation_m"]
CLIM=["wc_bio1","wc_bio5","wc_bio12","wc_bio15"]
SOIL=["soil_pH","soil_SOC","soil_N","soil_clay","soil_available_water_proxy"]
FAMILIES={
    "GEO":GEO,
    "GEO_CLIMATE":GEO+CLIM,
    "GEO_CLIMATE_SOIL":GEO+CLIM+SOIL,
}
MIN_COMPLETE_SPECIES=1000
MIN_CLASS_PER_FOLD_TRAIN=4
SPLITS=5


def photo_equal_frame(data:pd.DataFrame)->tuple[pd.DataFrame,dict]:
    if len(data)!=42111 or data.inat_taxon_id.nunique()!=42111:
        raise ValueError("Original 42111 breadth photo-equal species denominator changed")
    req={"inat_taxon_id","species","morph","measurement_status","latitude","longitude",
         "site_geo_status","environment_climate_complete","environment_soil_complete",
         *GEO[3:],*CLIM,*SOIL}
    if not req.issubset(data):
        raise ValueError(f"Missing original source+environment columns {req-set(data)}")
    d=data.copy()
    classified=d.morph.isin(COLOURS)&d.measurement_status.eq("classified_four_state_morph")
    geolocated=d.site_geo_status.eq("VALID_PUBLIC_ORIGINAL_PHOTO_POINT")
    climate=geolocated & d.environment_climate_complete.astype(bool)
    soil=geolocated & d.environment_soil_complete.astype(bool)
    d["abs_latitude"]=pd.to_numeric(d.latitude,errors="coerce").abs()
    rad=np.deg2rad(pd.to_numeric(d.longitude,errors="coerce"))
    d["lon_sin"]=np.sin(rad)
    d["lon_cos"]=np.cos(rad)
    d["genus"]=d.species.fillna("").astype(str).str.split().str[0].fillna("")
    d["site_cell_162"]=np.where(geolocated,cell_id(d.latitude.to_numpy(float),d.longitude.to_numpy(float)),-1)
    support={
        "source_taxa":len(d),
        "original_colour_classifiable_taxa":int(classified.sum()),
        "original_colour_unclassifiable_taxa":int((~classified).sum()),
        "source_photo_geolocated":int(geolocated.sum()),
        "geo_plus_climate_eligible_classified":int((classified&climate).sum()),
        "geo_plus_soil_eligible_classified":int((classified&soil).sum()),
        "geo_climate_soil_eligible_classified":int((classified&climate&soil).sum()),
        "photo_classifiable_but_geolocated_missing":int((classified&~geolocated).sum()),
        "classified_complete_soil_fraction":float((classified&climate&soil).sum()/max(1,classified.sum())),
    }
    good=classified & climate & soil & d[GEO+CLIM+SOIL].notna().all(axis=1)
    f=d.loc[good].copy()
    if len(f)<MIN_COMPLETE_SPECIES:
        raise RuntimeError("Too few complete coloured species for honest global climate/soil model")
    if (f.genus=="").any() or f.site_cell_162.lt(0).any():
        raise RuntimeError("Missing genus or true point spatial cell in eligible rows")
    return f,support


def cell_id(latitude:np.ndarray,longitude:np.ndarray)->np.ndarray:
    lat=np.asarray(latitude,float);lon=np.asarray(longitude,float)
    finite=np.isfinite(lat)&np.isfinite(lon)
    result=np.full(len(lat),-1,dtype=int)
    row=np.clip(np.floor((np.sin(np.deg2rad(lat[finite]))+1.0)*4.5).astype(int),0,8)
    col=np.clip(np.floor((lon[finite]+180.0)/20.0).astype(int),0,17)
    result[finite]=18*row+col
    return result


def region_coverage(all_species:pd.DataFrame)->pd.DataFrame:
    """Include ALL original 42111 source taxa, even geolocation/classification missing."""
    d=all_species.copy()
    classified=(d.morph.isin(COLOURS)&d.measurement_status.eq("classified_four_state_morph"))
    loc=d.site_geo_status.eq("VALID_PUBLIC_ORIGINAL_PHOTO_POINT")
    abslat=pd.to_numeric(d.latitude,errors="coerce").abs()
    region=np.select(
        [loc & abslat.lt(30), loc & abslat.ge(30)&abslat.lt(60),
         loc & abslat.ge(60)&abslat.le(90)],
        ["0_30","30_60","60_90"],
        default="NO_EXACT_PUBLIC_GEO",
    )
    d["source_latitude_region"]=region
    d["classified"]=classified
    d["climate_complete"]=loc & d.environment_climate_complete.astype(bool)
    d["soil_complete"]=loc & d.environment_soil_complete.astype(bool)
    d["all_complete"]=d.climate_complete & d.soil_complete & d.wc_elevation_m.notna()
    return d.groupby("source_latitude_region",sort=True).agg(
        source_species=("inat_taxon_id","size"),
        classified_flower_photo=("classified","sum"),
        climate_complete_species=("climate_complete","sum"),
        soil_complete_species=("soil_complete","sum"),
        all_environment_complete_species=("all_complete","sum"),
    ).reset_index()


def split_loss(data:pd.DataFrame,group_by:str)->dict:
    from sklearn.model_selection import GroupKFold
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import make_pipeline
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import log_loss,accuracy_score
    if group_by not in ("genus","site_cell_162"):
        raise ValueError("Only prespecified taxonomic or spatial holdouts allowed")
    groups=data[group_by].to_numpy()
    if len(set(groups))<SPLITS:
        raise RuntimeError("Insufficient taxonomic/spatial groups for held-out prediction")
    y=pd.Categorical(data.morph,categories=COLOURS).codes
    if (y<0).any(): raise RuntimeError("Unclassified source photo leaked into model")
    foldplan=list(GroupKFold(n_splits=SPLITS).split(data,y,groups=groups))
    prediction={}
    metrics={}
    for name,features in FAMILIES.items():
        x=data[features].to_numpy(dtype=float)
        if not np.isfinite(x).all():
            raise ValueError("Feature family missing in supposedly complete species")
        predictions=np.zeros((len(x),len(COLOURS)),float)
        fold_metrics=[]
        for fold,(train,test) in enumerate(foldplan):
            if len(set(y[train]))!=MIN_CLASS_PER_FOLD_TRAIN:
                raise RuntimeError(f"Fold {fold} training missing flower-colour class")
            model=make_pipeline(StandardScaler(),LogisticRegression(
                max_iter=800,C=1.0,random_state=20261009))
            model.fit(x[train],y[train])
            predictions[test]=model.predict_proba(x[test])
            fold_metrics.append({
                "fold":fold,"n_test":int(len(test)),
                "n_heldout_groups":int(len(set(groups[test]))),
                "log_loss":float(log_loss(y[test],predictions[test],labels=np.arange(len(COLOURS)))),
            })
        if not np.allclose(predictions.sum(axis=1),1.0,atol=1e-6):
            raise RuntimeError("Incomplete predicted class probabilities")
        metrics[name]={
            "n_species_same_tested":int(len(data)),
            "pooled_heldout_log_loss":float(log_loss(y,predictions,labels=np.arange(len(COLOURS)))),
            "pooled_heldout_accuracy":float(accuracy_score(y,predictions.argmax(axis=1))),
            "heldout_fold_diagnostics":fold_metrics,
        }
    clim_gain=metrics["GEO"]["pooled_heldout_log_loss"]-metrics["GEO_CLIMATE"]["pooled_heldout_log_loss"]
    soil_gain=metrics["GEO_CLIMATE"]["pooled_heldout_log_loss"]-metrics["GEO_CLIMATE_SOIL"]["pooled_heldout_log_loss"]
    return {
        "heldout_group":group_by,
        "n_independent_heldout_groups":int(len(set(groups))),
        "heldout_species":len(data),
        "feature_families":metrics,
        "climate_predictive_gain_logloss_reduction":float(clim_gain),
        "soil_increment_predictive_gain_logloss_reduction":float(soil_gain),
        "soil_gain_positive_in_all_folds":all(
            metrics["GEO_CLIMATE"]["heldout_fold_diagnostics"][i]["log_loss"] >
            metrics["GEO_CLIMATE_SOIL"]["heldout_fold_diagnostics"][i]["log_loss"]
            for i in range(SPLITS)
        ),
        "interpretation":"Predictive association among SAME soil-complete classified photos, not an experimental ecological mechanism",
    }


def environment_profiles(data:pd.DataFrame)->pd.DataFrame:
    z=[]
    for morph,g in data.groupby("morph",observed=True):
        for key in GEO[3:]+CLIM+SOIL:
            p=g[key].to_numpy(float)
            z.append({
                "source_photo_colour":morph,
                "environmental_feature":key,
                "n_species":len(g),
                "median":float(np.median(p)),
                "q25":float(np.quantile(p,.25)),
                "q75":float(np.quantile(p,.75)),
            })
    return pd.DataFrame(z)


def run(source:pd.DataFrame)->tuple[dict,pd.DataFrame]:
    subset,support=photo_equal_frame(source)
    if len(subset.morph.unique())!=4:
        raise RuntimeError("Source complete cases miss entire photo-colour class")
    models=[split_loss(subset,"genus"),split_loss(subset,"site_cell_162")]
    report={
        "schema":"fcp_global42111_colour_climate_soil_predictive_characterization_v1",
        "date_jst":"2026-10-09",
        "status":"SOURCE_PHOTO_SPECIES_EQUAL_ENVIRONMENTAL_ASSOCIATION_ONLY",
        "global_original_species_denominator":42111,
        "global_one_photo_classified_species_denominator":18457,
        "all_original_photo_rows_preserved_in_source":True,
        "historical_source_completeness":support,
        "n_full_abiotic_colour_sample":len(subset),
        "n_genera_complete":int(subset.genus.nunique()),
        "n_source_photo_cells_complete":int(subset.site_cell_162.nunique()),
        "class_counts_complete_case":{c:int((subset.morph==c).sum()) for c in COLOURS},
        "identical_complete_case_rows_for_every_model":True,
        "photo_colour_outcomes_reclassified":False,
        "spatial_cell_centroid_as_plant_soil":False,
        "models":models,
        "interpretation_limits":[
            "One source photo of a species is not the species modal or allelic flower colour",
            "Out-of-genus and out-of-cell validation exposes species turnover or spatial portability; not local adaptation",
            "Environmental covariates are long-term means and 5km modelled topsoil, not experimental exposures",
            "Complete-case soil coverage is observationally selected and cannot represent all original 42111 taxa",
            "A positive predictive log-loss gain neither implies a true fitness effect nor causal soil pigment selection",
            "FCP old 1499 species and new 2000+730 prospective cohorts were not changed",
        ],
    }
    return report,environment_profiles(subset)


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--species-breadth-with-abiotic",required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    a=p.parse_args()
    d=pd.read_csv(a.species_breadth_with_abiotic,low_memory=False)
    report,profile=run(d)
    coverage=region_coverage(d)
    if int(coverage.source_species.sum())!=42111:
        raise RuntimeError('The source-wide latitude/environment missingness denominator drifted')
    a.outdir.mkdir(parents=True,exist_ok=True)
    coverage.to_csv(a.outdir/'original_42111_source_latitude_environment_coverage.csv',index=False)
    (a.outdir/"result.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    profile.to_csv(a.outdir/"species_equal_flower_colour_environmental_profiles.csv",index=False)
    print(json.dumps(report,indent=2,sort_keys=True),flush=True)
if __name__=="__main__":
    main()

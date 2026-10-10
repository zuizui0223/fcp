#!/usr/bin/env python3
"""Source-exposed FCP original 42111 species: solar, temperature, soil ranking.

Within the SAME exact original complete-case photographed species, fit
geography versus full 2.1 WorldClim BIO/temp/moisture/solar/elevation, plus
SoilGrids edaphic predictors. Block leave-one-out importance and named single
covariate sensitivity are out-of-genus and out-of-original-photo-cell scored.
Everything reported; no favourable feature selection, causal inference or
new independent taxa. 499 group-bootstrap percentile intervals conditional
on one realized fixed source cross-validation fit.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd

N_ORIGINAL=42111
N_CLASSIFIED=18457
CLASSES=("white","yellow_orange","red_pink","blue_purple")
N_FOLDS=5
BOOT=499
SEED=20261010

GEO=("abs_latitude","lon_sin","lon_cos")
BLOCKS={
    "elevation":("wc_elevation_m",),
    "temperature":("wc_bio1","wc_bio5","wc_bio6","wc_bio4","wc_bio2","wc_bio7"),
    "precipitation":("wc_bio12","wc_bio15","wc_bio14","wc_bio18"),
    "solar_radiation":("wc_srad_annual_kj_m2_day","wc_srad_monthly_cv"),
    "soil":("soil_pH","soil_SOC","soil_N","soil_clay","soil_available_water_proxy"),
}
ALL=GEO+tuple(f for v in BLOCKS.values() for f in v)
SCHEMA="fcp_original_42111_complete_case_multiabiotic_photo_colour_cv_v1"


def cell_ids(lat,lon)->np.ndarray:
    a=np.asarray(lat,float);b=np.asarray(lon,float)
    out=np.full(len(a),-1,int)
    yes=np.isfinite(a)&np.isfinite(b)&(abs(a)<=90)&(abs(b)<=180)
    r=np.clip(np.floor((np.sin(np.deg2rad(a[yes]))+1)*4.5).astype(int),0,8)
    c=np.clip(np.floor((b[yes]+180)/20).astype(int),0,17)
    out[yes]=18*r+c
    return out


def make_source(d:pd.DataFrame,strict:bool=True)->tuple[pd.DataFrame,dict]:
    req={"inat_taxon_id","species","morph","measurement_status","site_geo_status",
         "latitude","longitude",*ALL[3:]}
    if not req.issubset(d):
        raise ValueError("Original source photograph expanded-abiotic fields missing "+str(sorted(req-set(d))))
    if d.inat_taxon_id.duplicated().any() or (strict and len(d)!=N_ORIGINAL):
        raise ValueError("Original 42111 photo-equal species identities changed")
    x=d.copy()
    labeled=x.measurement_status.eq("classified_four_state_morph")&x.morph.isin(CLASSES)
    if strict and labeled.sum()!=N_CLASSIFIED:
        raise ValueError("Frozen original 18457 photo colour classification drifted")
    loc=x.site_geo_status.eq("VALID_PUBLIC_ORIGINAL_PHOTO_POINT")
    lat=pd.to_numeric(x.latitude,errors="coerce").to_numpy(float)
    lon=pd.to_numeric(x.longitude,errors="coerce").to_numpy(float)
    x["abs_latitude"]=abs(lat)
    x["lon_sin"]=np.sin(np.deg2rad(lon))
    x["lon_cos"]=np.cos(np.deg2rad(lon))
    x["source_cell_162"]=cell_ids(lat,lon)
    x["genus"]=x.species.fillna("").astype(str).str.split().str[0]
    complete=(labeled & loc & x[list(ALL)].notna().all(axis=1) & x.source_cell_162.ge(0))
    study=x.loc[complete].copy().reset_index(drop=True)
    if len(study)<1000 or set(study.morph)!=set(CLASSES) or study.genus.eq("").any():
        raise RuntimeError("Insufficient full solar temperature precip soil species-complete outcome sample")
    bycol=study.morph.value_counts().reindex(CLASSES,fill_value=0)
    coverage={
        "n_original_species_taxa":len(x),
        "n_original_classifiable_photo_species":int(labeled.sum()),
        "n_original_unclassifiable_photo_species":int((~labeled).sum()),
        "n_geolocated_source_photo_species":int(loc.sum()),
        "n_complete_all_environment_and_original_colour":len(study),
        "n_original_classified_excluded_from_common_complete_case":int(labeled.sum()-len(study)),
        "n_original_genera_in_common_complete_case":int(study.genus.nunique()),
        "n_original_geographic_cells_in_common_complete_case":int(study.source_cell_162.nunique()),
        "source_colour_classes_common_complete_case":{k:int(bycol[k]) for k in CLASSES},
        "n_original_classified_with_monthly_solar_complete":int((labeled&x[["wc_srad_annual_kj_m2_day","wc_srad_monthly_cv"]].notna().all(axis=1)).sum()),
    }
    return study,coverage


def features()->dict[str,tuple[str,...]]:
    families={
        "GEOGRAPHY_ONLY":GEO,
        "FULL_ALL_BLOCKS":ALL,
    }
    for name,fields in BLOCKS.items():
        families["FULL_MINUS_"+name.upper()]=tuple(f for f in ALL if f not in fields)
    for name,fields in BLOCKS.items():
        for f in fields:
            families["FULL_MINUS_SINGLE_"+f.upper()]=tuple(q for q in ALL if q!=f)
    # Follow source nomenclature, no data-driven candidate construction.
    if len(families)!=2+len(BLOCKS)+sum(len(x) for x in BLOCKS.values()):
        raise RuntimeError("Feature family count wrong")
    return families


def oof(d:pd.DataFrame,heldout:str,*,nboot:int=BOOT)->dict:
    from sklearn.model_selection import GroupKFold
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import make_pipeline
    from sklearn.metrics import log_loss
    if heldout not in ("genus","source_cell_162"):
        raise ValueError("Only source-defined genus or region independent holdout")
    y=pd.Categorical(d.morph,categories=CLASSES).codes
    group=d[heldout].to_numpy()
    if len(set(group))<N_FOLDS:
        raise ValueError("Original photo grouping cannot support 5 folds")
    plan=list(GroupKFold(n_splits=N_FOLDS).split(d,y,groups=group))
    families=features()
    loss={}
    metrics={}
    for k,cols in families.items():
        values=d[list(cols)].to_numpy(float)
        if not np.isfinite(values).all():
            raise RuntimeError("Compared solar/climate/soil families have different missingness")
        probabilities=np.full((len(d),4),np.nan)
        for tr,te in plan:
            if len(set(y[tr]))!=4:
                raise RuntimeError("Photo-colour training fold lost a class")
            model=make_pipeline(StandardScaler(),LogisticRegression(max_iter=550,C=1.0,random_state=SEED))
            model.fit(values[tr],y[tr])
            probabilities[te]=model.predict_proba(values[te])
        if not np.isfinite(probabilities).all() or not np.allclose(probabilities.sum(axis=1),1):
            raise RuntimeError("Incomplete 4-colour heldout source photo prediction")
        clipped=np.clip(probabilities[np.arange(len(d)),y],1e-12,1)
        errors=-np.log(clipped)
        loss[k]=errors
        metrics[k]={"heldout_mean_multiclass_logloss":float(np.mean(errors)),
                    "n_original_photo_taxa":len(d)}
    original,inv=np.unique(group,return_inverse=True)
    sizes=np.bincount(inv)
    rng=np.random.default_rng(SEED+(0 if heldout=="genus" else 1))
    draws=rng.integers(0,len(original),size=(nboot,len(original)))
    values={}
    comparisons={
        "full_minus_geography_only":("GEOGRAPHY_ONLY","FULL_ALL_BLOCKS"),
    }
    for block in BLOCKS:
        comparisons["conditional_"+block]=("FULL_MINUS_"+block.upper(),"FULL_ALL_BLOCKS")
    for fields in BLOCKS.values():
        for f in fields:
            comparisons["conditional_single_"+f]=("FULL_MINUS_SINGLE_"+f.upper(),"FULL_ALL_BLOCKS")
    for name,(base,extra) in comparisons.items():
        gain=loss[base]-loss[extra]
        total=np.bincount(inv,weights=gain)
        sampling=total[draws].sum(axis=1)/sizes[draws].sum(axis=1)
        ci=np.quantile(sampling,[.025,.975])
        values[name]={
            "heldout_mean_logloss_reduction":float(np.mean(gain)),
            "group_resampled_fixed_oof_95CI":[float(v) for v in ci],
            "nominal_conditional_interval_entirely_positive":bool(ci[0]>0),
            "n_test_species_unchanged":len(d),
        }
    return {
        "heldout_group":heldout,
        "original_source_n_independent_heldout_groups":len(original),
        "fivefold_grouping_kind":"original nominal genus" if heldout=="genus" else "original photo site's 162 equal-area cell",
        "all_models_same_original_photo_species_and_fixed_folds":True,
        "n_original_species_evaluated":len(d),
        "source_model_scores":metrics,
        "full_vs_geography_and_drop_one_predictor_gains":values,
        "confidence_intervals_conditional_on_realized_fixed_model_folds":True,
    }


def run(d:pd.DataFrame,*,strict=True,nboot=BOOT)->dict:
    sub,den=make_source(d,strict=strict)
    parts=[oof(sub,"genus",nboot=nboot),oof(sub,"source_cell_162",nboot=nboot)]
    return {
        "schema":SCHEMA,
        "date_jst":"2026-10-10",
        "status":"RETROSPECTIVE_EXPLORATORY_PHOTO_COLOUR_MULTIVARIABLE_IMPORTANCE",
        "original_source_denominators":den,
        "tested_named_blocks":{k:list(v) for k,v in BLOCKS.items()},
        "geography_features":list(GEO),
        "n_feature_families":len(features()),
        "n_block_importance_comparisons":len(BLOCKS),
        "n_single_predictor_importance_comparisons":sum(len(x) for x in BLOCKS.values()),
        "original_source_photo_pixels_and_colours_reclassified":False,
        "unclassified_original_photos_not_assigned_a_colour":True,
        "source_photo_sites_not_original_cell_centroids":True,
        "true_local_canopy_sunshine_not_observed":True,
        "full_source_872_1761_phylogenetic_claims_remain_HOLD":True,
        "models":parts,
        "multiple_comparison_guard":"All blocks and single features reported; 95pct intervals unadjusted exploratory sensitivity, not confirmatory tests; no winner-only p-value",
        "hard_nonclaims":[
            "An important covariate prediction does not demonstrate climate-driven floral adaptation",
            "Source one-photo-per-species four-state class is not a genetically fixed species flower colour",
            "WorldClim solar radiation 1970-2000 is grid-cell climatic incident radiation, not under-canopy exposure",
            "WorldClim 10 arc-minute values are coarse relative to the location error of some uploaded photos",
            "SoilGrids 5km topsoil predictors are modeled, not actual floral rhizosphere assays",
            "Block collinearity means conditional drop-one scores need not sum to total",
            "Same soil+solar complete-case sample has observational selection, never extrapolated to 42111 taxa",
        ],
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--expanded-source-breadth",type=Path,required=True)
    p.add_argument("--outdir",type=Path,required=True)
    a=p.parse_args()
    orig=pd.read_csv(a.expanded_source_breadth,low_memory=False)
    result=run(orig)
    a.outdir.mkdir(parents=True,exist_ok=True)
    (a.outdir/"result.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    short=[]
    for item in result["models"]:
        for name,value in item["full_vs_geography_and_drop_one_predictor_gains"].items():
            short.append({"heldout_group":item["heldout_group"],"comparison":name,
                "n_common_original_species":item["n_original_species_evaluated"],
                "gain_heldout_logloss":value["heldout_mean_logloss_reduction"],
                "ci_low":value["group_resampled_fixed_oof_95CI"][0],
                "ci_high":value["group_resampled_fixed_oof_95CI"][1],
                "nominal_interval_positive":value["nominal_conditional_interval_entirely_positive"]})
    pd.DataFrame(short).to_csv(a.outdir/"all_predictor_block_and_single_feature_comparisons.csv",index=False)
    print(json.dumps({
        "schema":SCHEMA,
        "support":result["original_source_denominators"],
        "n_families":result["n_feature_families"],
        "block_heldouts":[{
            "heldout_group":m["heldout_group"],
            "blocks":{key:val for key,val in m["full_vs_geography_and_drop_one_predictor_gains"].items()
                      if key.startswith("conditional_") and not key.startswith("conditional_single_")}
        } for m in result["models"]]},sort_keys=True),flush=True)


if __name__=="__main__":
    main()

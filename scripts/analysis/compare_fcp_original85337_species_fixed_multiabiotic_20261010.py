#!/usr/bin/env python3
"""Same-species cross-region global FCP original flower photograph abiotic CV.

85,337 ORIGINAL taxon×cell photographs (39,075 four-colour classifiable),
joined by exact photo+observation+taxon IDs to 100,543 original real-site
expanded BIO1-19, solar, wind, vapor, elevation and ten soil predictors.
Species intercepts come from DIFFERENT TRAINING CELLS of same nominal species.
Thus species fixed composition is explicitly removed; new single-photo
species cannot be used as a repeated species estimand and are kept missing.
Original geography-cell GroupKFold; each model fits the identical test rows
under climate-complete and separately soil-complete cohorts. Noncausal.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd

from compare_fcp_global42111_full_climate_sun_soil_blocks_20261010 import (
    BLOCKS,GEO,CLASSES
)
from extend_fcp_global42111_worldclim_solar_20261010 import NEW_FEATURES

N_CELL=85337
N_CLASSIFIED=39075
N_SOURCE_SITES=100543
FOLDS=5
SEED=20261010
BOOT=199
RIDGE=10.
MIN_TEST=200
OLD_KEY=("inat_taxon_id","observation_id","photo_id")
SCHEMA="fcp_original_85337_species_intercept_multiabiotic_photo_colour_cv_v1"


def match_original(old:pd.DataFrame,fullsites:pd.DataFrame,*,strict=True)->tuple[pd.DataFrame,dict]:
    if not set(OLD_KEY).issubset(old) or not set(OLD_KEY).issubset(fullsites):
        raise ValueError("Original species observation photo key missing")
    if strict and (len(old)!=N_CELL or len(fullsites)!=N_SOURCE_SITES):
        raise ValueError("Original 85337 taxon-cell / 100543 unique photo site population changed")
    if old[list(OLD_KEY)].duplicated().any() or fullsites.photo_id.duplicated().any():
        raise ValueError("Duplicate source taxon×cell or full-site photo identity")
    if not {"cell_id","morph","measurement_status","species"}.issubset(old):
        raise ValueError("Original repeated species photo colour ledger missing")
    classified=old.measurement_status.eq("classified_four_state_morph") & old.morph.isin(CLASSES)
    if strict and int(classified.sum())!=N_CLASSIFIED:
        raise ValueError("Original 39075 classifiable same-species source photos changed")
    new=list(NEW_FEATURES)
    if not set(new).issubset(fullsites):
        raise ValueError("Original 100543 genuine photo site expanded weather/soil missing")
    joined=old.merge(fullsites[list(OLD_KEY)+new],on=list(OLD_KEY),
                     how="left",validate="one_to_one",indicator=True)
    if len(joined)!=len(old) or not joined["_merge"].eq("both").all():
        raise RuntimeError("Original 85337 original source photo key lost while attaching environment")
    joined=joined.drop(columns=["_merge"])
    if strict and joined[["inat_taxon_id","cell_id"]].duplicated().any():
        raise ValueError("Duplicated historical species×geocell source photo")
    return joined,{
        "source_original_taxon_cell_photo_rows":len(old),
        "source_four_state_colour_classifiable":int(classified.sum()),
        "source_unclassifiable_photo_cell_records":int((~classified).sum()),
        "n_original_source_taxa_in_cell_ledger":int(old.inat_taxon_id.nunique()),
        "all_original_photo_ID_triples_matched_exactly":True,
    }


def populations(joined:pd.DataFrame)->tuple[dict[str,pd.DataFrame],dict]:
    d=joined.copy()
    req={"inat_taxon_id","cell_id","species","morph","measurement_status",
         "latitude","longitude","wc_elevation_m",*GEO[0:0],
         *(v for cols in BLOCKS.values() for v in cols)}
    if not req.issubset(d):
        raise ValueError("Missing full original photograph climate/soil status or covariates")
    labelled=d.measurement_status.eq("classified_four_state_morph")&d.morph.isin(CLASSES)
    real=d.site_geo_status.eq("VALID_PUBLIC_ORIGINAL_PHOTO_POINT")
    lat=pd.to_numeric(d.latitude,errors="coerce")
    lon=pd.to_numeric(d.longitude,errors="coerce")
    d["abs_latitude"]=lat.abs()
    d["lon_sin"]=np.sin(np.deg2rad(lon))
    d["lon_cos"]=np.cos(np.deg2rad(lon))
    if not pd.to_numeric(d.cell_id,errors="raise").between(0,161).all():
        raise ValueError("Old source equal-area region identity corrupt")
    climate_features=tuple(x for name,features in BLOCKS.items() if name!="soil" for x in features)
    base=labelled&real&d[list(GEO+climate_features)].notna().all(axis=1)
    full=base&d[list(BLOCKS["soil"])].notna().all(axis=1)
    pools={"CLIMATE_SOURCE_ONLY":d.loc[base].copy().reset_index(drop=True),
           "CLIMATE_AND_SOIL":d.loc[full].copy().reset_index(drop=True)}
    ledger={
        "n_original_source_photo_cell_records":len(d),
        "n_original_fourstate_classified_photo_cells":int(labelled.sum()),
        "n_original_fourstate_unclassified_photo_cells":int((~labelled).sum()),
        "n_original_photo_taxa_any_cell":int(d.inat_taxon_id.nunique()),
        "n_real_photo_site_coordinates":int(real.sum()),
        "n_original_classified_photo_cells_climate_complete":int(base.sum()),
        "n_original_classified_photo_cells_climate_soil_complete":int(full.sum()),
        "n_original_climate_complete_repeated_species":int(pools["CLIMATE_SOURCE_ONLY"].inat_taxon_id.value_counts().ge(2).sum()),
        "n_original_soil_complete_repeated_species":int(pools["CLIMATE_AND_SOIL"].inat_taxon_id.value_counts().ge(2).sum()),
        "original_photo_colour_states_reclassified":False,
        "original_unclassifiable_photos_never_white_or_monomorphic":True,
    }
    return pools,ledger


def source_features(*,with_soil:bool)->dict[str,tuple[str,...]]:
    blocks={k:v for k,v in BLOCKS.items() if with_soil or k!="soil"}
    allweather=tuple(x for vals in blocks.values() for x in vals)
    groups={
        "SPECIES_ONLY":(),
        "SPECIES_PLUS_GEOGRAPHY":GEO,
        "SPECIES_GEO_FULL_ABIOTIC":GEO+allweather,
    }
    for name,cols in blocks.items():
        groups["SPECIES_GEO_FULL_MINUS_"+name.upper()]=GEO+tuple(
            z for z in allweather if z not in cols)
    # Single features are exposed in global species-composition model; here
    # species control is the higher priority, with the same full block suite.
    return groups


def train_predict(tr:pd.DataFrame,te:pd.DataFrame,columns:tuple[str,...],eligible:np.ndarray)->np.ndarray:
    source=tr.inat_taxon_id.to_numpy(int)
    targets=te.iloc[eligible]
    if set(tr[["inat_taxon_id","cell_id"]].itertuples(index=False,name=None))&set(
        targets[["inat_taxon_id","cell_id"]].itertuples(index=False,name=None)):
        raise RuntimeError("Same original species×site photo appeared on both sides")
    y=pd.get_dummies(pd.Categorical(tr.morph,categories=CLASSES)).to_numpy(float)
    ym=pd.DataFrame(y,columns=list(CLASSES)).assign(species=source).groupby("species",sort=False)[list(CLASSES)].mean()
    baseline=ym.reindex(targets.inat_taxon_id.to_numpy(int)).to_numpy(float)
    if np.isnan(baseline).any():raise RuntimeError("Unknown original source species in heldout species intercept")
    if len(columns)==0:return baseline
    xtr=tr[list(columns)].to_numpy(float)
    xte=targets[list(columns)].to_numpy(float)
    if not np.isfinite(xtr).all() or not np.isfinite(xte).all():
        raise ValueError("One full feature was unavailable in climate/soil source")
    mean=xtr.mean(axis=0)
    sd=xtr.std(axis=0)
    sd[sd<1e-8]=1.0
    xtr=(xtr-mean)/sd;xte=(xte-mean)/sd
    xm=pd.DataFrame(xtr,columns=list(columns)).assign(species=source).groupby("species",sort=False)[list(columns)].mean()
    devx=xtr-xm.reindex(source).to_numpy(float)
    devy=y-ym.reindex(source).to_numpy(float)
    beta=np.linalg.solve(devx.T@devx+RIDGE*np.eye(len(columns)),devx.T@devy)
    predicted=baseline+(xte-xm.reindex(targets.inat_taxon_id.to_numpy(int)).to_numpy(float))@beta
    predicted=np.maximum(predicted,0)
    denom=predicted.sum(axis=1,keepdims=True)
    return np.divide(predicted,denom,out=np.full_like(predicted,.25),where=denom>0)


def species_out_of_cell(d:pd.DataFrame,*,with_soil:bool,nboot:int=BOOT)->dict:
    from sklearn.model_selection import GroupKFold
    if len(d)<MIN_TEST or d.cell_id.nunique()<FOLDS:
        return {"status":"HOLD_INSUFFICIENT_SOURCE_PHOTO_SAME_SPECIES_CROSSCELL_COVERAGE",
                "n_original_classified_environment_complete":len(d)}
    y=pd.Categorical(d.morph,categories=CLASSES).codes
    if np.any(y<0):raise RuntimeError("Unclassified photo label entered species fixed model")
    families=source_features(with_soil=with_soil)
    pred={k:np.full((len(d),4),np.nan,float) for k in families}
    eligible=np.zeros(len(d),bool)
    fit_support=[]
    for i,(train,test) in enumerate(GroupKFold(n_splits=FOLDS).split(d,y,d.cell_id)):
        tr=d.iloc[train]
        te=d.iloc[test]
        grouped=tr.inat_taxon_id.value_counts()
        ids=np.flatnonzero(te.inat_taxon_id.isin(grouped.index).to_numpy())
        eligible[test[ids]]=True
        n_multi=int(te.iloc[ids].inat_taxon_id.isin(grouped[grouped>=2].index).sum())
        fit_support.append({
            "fold":i,"n_test_region_source_photos":len(test),
            "n_test_photos_with_train_species_intercept":len(ids),
            "n_test_photos_with_at_least_two_training_same_species_other_cells":n_multi,
            "n_test_geographical_source_cells":int(te.cell_id.nunique())})
        for k,cols in families.items():
            if len(ids):
                pred[k][test[ids]]=train_predict(tr,te,cols,ids)
    n=int(eligible.sum())
    if n<MIN_TEST or len(np.unique(y[eligible]))<4:
        return {"status":"HOLD_INSUFFICIENT_SAME_SPECIES_HELDOUT_SOURCE_PHOTO_SUPPORT",
                "n_original_classified_environment_complete":len(d),
                "n_test_species_supported":n,"folds":fit_support}
    if not all(np.isfinite(v[eligible]).all() for v in pred.values()):
        raise RuntimeError("Compared models do not cover the same originally heldout species photos")
    truth=np.eye(4)[y[eligible]]
    losses={}
    scores={}
    for k,values in pred.items():
        err=np.square(values[eligible]-truth).sum(axis=1)
        losses[k]=err
        scores[k]={"heldout_multiclass_brier":float(err.mean()),"same_photo_site_species_holdout_n":n}
    groups=d.loc[eligible,"cell_id"].to_numpy(int)
    taxon=d.loc[eligible,"inat_taxon_id"].to_numpy(int)
    contrasts={}
    comps={"all_abiotic_beyond_species_geography":("SPECIES_PLUS_GEOGRAPHY","SPECIES_GEO_FULL_ABIOTIC")}
    blocks={k:v for k,v in BLOCKS.items() if with_soil or k!="soil"}
    for name in blocks:
        comps["conditional_"+name]=("SPECIES_GEO_FULL_MINUS_"+name.upper(),"SPECIES_GEO_FULL_ABIOTIC")
    for name,(base,full) in comps.items():
        diff=losses[base]-losses[full]
        cis={}
        for kind,g in (("source_cell_162",groups),("original_species",taxon)):
            unique,idx=np.unique(g,return_inverse=True)
            totals=np.bincount(idx,weights=diff)
            size=np.bincount(idx)
            rng=np.random.default_rng(SEED+len(name)+(1 if kind=="original_species" else 0))
            draw=rng.integers(0,len(unique),size=(nboot,len(unique)))
            boot=totals[draw].sum(axis=1)/size[draw].sum(axis=1)
            cis[kind]=[float(k) for k in np.quantile(boot,[.025,.975])]
        contrasts[name]={
            "heldout_same_species_photo_brier_improvement":float(diff.mean()),
            "original_geo_cell_bootstrap_95CI":cis["source_cell_162"],
            "original_species_bootstrap_95CI":cis["original_species"],
            "positive_both_conditional_CIs":bool(cis["source_cell_162"][0]>0 and cis["original_species"][0]>0),
            "same_heldout_source_photos_all_models":True,
        }
    return {
        "status":"SOURCE_ORIGINAL_SPECIES_INTERCEPT_OUTOFCELL_ABIOTIC_DIAGNOSTIC",
        "n_original_classified_environment_complete":len(d),
        "n_heldout_original_photos_with_species_training_intercept":n,
        "n_distinct_nominal_source_species_evaluated":int(d.loc[eligible,"inat_taxon_id"].nunique()),
        "n_heldout_photos_with_two_plus_training_same_species_photo_cells":sum(i["n_test_photos_with_at_least_two_training_same_species_other_cells"] for i in fit_support),
        "n_original_photo_cells_evaluated":int(d.loc[eligible,"cell_id"].nunique()),
        "n_discordant_original_photo_species_among_evaluated":int(d.loc[eligible].groupby("inat_taxon_id").morph.nunique().gt(1).sum()),
        "frozen_heldout_cell_plan":fit_support,
        "training_species_colour_means_only_from_other_cells":True,
        "species_fixed_intercept_eliminates_unidentified_lineage_static_baseline":True,
        "within_species_abiotic_slopes_are_original_photo_observational_not_genetic":True,
        "comparison_model_scores":scores,
        "source_block_incremental_prediction":contrasts,
        "confidence_intervals_conditional_on_fixed_model_folds":True,
    }


def run(old:pd.DataFrame,fullsites:pd.DataFrame,*,strict=True,nboot:int=BOOT)->dict:
    joined,source=match_original(old,fullsites,strict=strict)
    pop,coverage=populations(joined)
    return {
        "schema":SCHEMA,
        "date_jst":"2026-10-10",
        "status":"ORIGINAL_SAME_SPECIES_DIFFERENT_REGION_PHOTO_COLOUR_ENVIRONMENTAL_ASSOCIATION",
        "original_repeated_photo_source":source,
        "original_abiotic_coverage":coverage,
        "climate_without_soil_selection":species_out_of_cell(pop["CLIMATE_SOURCE_ONLY"],with_soil=False,nboot=nboot),
        "climate_soil_joint_complete":species_out_of_cell(pop["CLIMATE_AND_SOIL"],with_soil=True,nboot=nboot),
        "all_42111_one_photo_species_remain_distinct_from_85337_repeated_cell_source":True,
        "all_real_photo_site_environment_no_cell_centroid_imputation":True,
        "all_models_have_source_species_intercepts_not_nominal_genus_substitutes":True,
        "phylogenetic_main_effect_not_estimable_simultaneously_with_species_intercepts":True,
        "photographic_colour_cannot_prove_genetic_polymorphism":True,
        "old_prospective_2000_730_taxa_untouched":True,
        "nonclaims":[
            "Cross-cell original flower photos represent different observed plants, potentially varying photographer and taxonomic identity",
            "A fixed species colour mean does not identify genetic or fitness variation without repeated biological population samples",
            "Phylogeny's time-invariant species intercept is absorbed by the species factor and cannot be separately estimated",
            "Photo-cell regional folds do not remove finer-than-cell geographic confounding",
            "Soil and solar modelled climatic variables are spatial predictions, not in situ plant exposures",
            "Feature block intervals are multiple posthoc exploratory comparisons conditional on fixed trained model",
            "Direct-tip LCVP subset previous analysis failed full original species coverage and spatial robustness gates",
        ],
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--original-taxon-cell",required=True,type=Path)
    ap.add_argument("--expanded-all-original-photo-sites",required=True,type=Path)
    ap.add_argument("--outdir",required=True,type=Path)
    a=ap.parse_args()
    result=run(pd.read_csv(a.original_taxon_cell,low_memory=False),
               pd.read_csv(a.expanded_all_original_photo_sites,low_memory=False))
    a.outdir.mkdir(parents=True,exist_ok=True)
    (a.outdir/"result.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps({
        "schema":SCHEMA,
        "source":result["original_repeated_photo_source"],
        "coverage":result["original_abiotic_coverage"],
        "climate":result["climate_without_soil_selection"],
        "full":result["climate_soil_joint_complete"]},sort_keys=True),flush=True)


if __name__=="__main__":
    main()

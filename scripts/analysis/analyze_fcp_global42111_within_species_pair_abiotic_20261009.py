#!/usr/bin/env python3
"""Within-species (one fixed observer-disjoint pair per taxon) flower-colour ecology.

Join original 13,416 measured pair endpoints to ORIGINAL 100,543 photo identities
and real WorldClim/SoilGrids site covariates. Do not use equal-area centroids,
invent colour states, refit on future 2,000+730 taxa, or substitute missing soil.

Target: out-of-genus / out-of-real-midpoint-cell prediction of source-photo
four-colour DISCORDANCE, not genetic adaptation or local fitness selection.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

SOURCE_PAIR_SHA256="26755316b2cbd07e1e2241eb4ddc1cddc85219fd7e96424df8c14df32575a4d2"
N_PAIRS=13416
BOTH_CLASSIFIABLE=3196
DISCORDANT=797
MIN_COMPLETE_PAIRS=300
N_FOLDS=5
N_BOOT=999
SEED=20261009
KEY=("inat_taxon_id","observation_id","photo_id")
CLIMATE=("wc_bio1","wc_bio5","wc_bio12","wc_bio15")
SOIL=("soil_pH","soil_SOC","soil_N","soil_clay","soil_available_water_proxy")
ELEV="wc_elevation_m"
SOURCE_GEO="VALID_PUBLIC_ORIGINAL_PHOTO_POINT"
BASE=("log_geodesic_distance_km","delta_abs_latitude","delta_elevation_m")
CLIM_FEATURES=tuple("delta_"+v for v in CLIMATE)
SOIL_FEATURES=tuple("delta_"+v for v in SOIL)
FAMILIES={"GEOGRAPHY":BASE,
          "GEOGRAPHY_CLIMATE":BASE+CLIM_FEATURES,
          "GEOGRAPHY_CLIMATE_SOIL":BASE+CLIM_FEATURES+SOIL_FEATURES}


def sha256(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1<<20),b""):h.update(block)
    return h.hexdigest()


def read_original_pairs(path:Path)->pd.DataFrame:
    if sha256(path)!=SOURCE_PAIR_SHA256:
        raise ValueError("Original 2026-09 observer-disjoint pair SHA256 drift")
    return pd.read_csv(path,low_memory=False)


def pair_distance_and_midpoint(lat1,lon1,lat2,lon2)->tuple[np.ndarray,np.ndarray]:
    lat1,lon1,lat2,lon2=map(lambda a:np.deg2rad(np.asarray(a,float)),
                             (lat1,lon1,lat2,lon2))
    a=np.column_stack((np.cos(lat1)*np.cos(lon1),np.cos(lat1)*np.sin(lon1),np.sin(lat1)))
    b=np.column_stack((np.cos(lat2)*np.cos(lon2),np.cos(lat2)*np.sin(lon2),np.sin(lat2)))
    dot=np.clip((a*b).sum(axis=1),-1,1)
    d=np.arccos(dot)*6371.0088
    midway=a+b
    mag=np.linalg.norm(midway,axis=1)
    good=np.isfinite(mag)&(mag>1e-10)
    cell=np.full(len(d),-1,dtype=int)
    if good.any():
        v=midway[good]/mag[good,None]
        mid_lat=np.rad2deg(np.arcsin(np.clip(v[:,2],-1,1)))
        mid_lon=np.rad2deg(np.arctan2(v[:,1],v[:,0]))
        y=np.clip(np.floor((np.sin(np.deg2rad(mid_lat))+1)*4.5).astype(int),0,8)
        x=np.clip(np.floor((mid_lon+180)/20).astype(int),0,17)
        cell[np.flatnonzero(good)]=18*y+x
    return d,cell


def join_origins(pairs:pd.DataFrame,env:pd.DataFrame)->tuple[pd.DataFrame,dict]:
    req={"inat_taxon_id","species","photo_id_1","photo_id_2",
         "observation_id_1","observation_id_2","observer_id_1","observer_id_2",
         "cell_id_1","cell_id_2","pair_state","both_endpoints_classifiable"}
    if not req.issubset(pairs):raise ValueError(f"Missing source fixed-pair columns {sorted(req-set(pairs))}")
    if len(pairs)!=N_PAIRS or pairs.inat_taxon_id.duplicated().any():
        raise ValueError("One immutable FCP pair per 13416 different taxa required")
    if pairs.cell_id_1.eq(pairs.cell_id_2).any() or pairs.observer_id_1.astype(str).eq(pairs.observer_id_2.astype(str)).any():
        raise ValueError("Same cell/observer invalid original fixed cross-cell contrast")
    if pairs.photo_id_1.eq(pairs.photo_id_2).any() or pairs.observation_id_1.eq(pairs.observation_id_2).any():
        raise ValueError("Same original photo or observation appears in both endpoints")
    if not set(KEY).issubset(env):raise ValueError("Original source photo/observation identity unavailable")
    if env.photo_id.duplicated().any() or env[list(KEY)].duplicated().any():
        raise ValueError("Non-injective old photo environmental matching")
    needs={*CLIMATE,*SOIL,ELEV,"site_geo_status","latitude","longitude",
           "environment_climate_complete","environment_soil_complete","environment_all_complete"}
    if not needs.issubset(env):raise ValueError(f"Missing true-location original abiotic source columns {needs-set(env)}")
    if not np.isfinite(pd.to_numeric(env.loc[env.site_geo_status.eq(SOURCE_GEO),"latitude"],errors="coerce")).all():
        raise ValueError("Claimed valid original site has no latitude")
    current=pairs.copy()
    for side in (1,2):
        keymaps={k:f"{k}_{side}" for k in KEY}
        columns=list(KEY)+sorted(needs)
        root=env[columns].rename(columns={**keymaps,**{k:f"{k}_{side}" for k in needs}})
        current=current.merge(root,on=list(keymaps.values()),how="left",validate="one_to_one",indicator=f"_join_{side}")
        if not current[f"_join_{side}"].eq("both").all():
            raise RuntimeError(f"Source pair endpoint {side} missing exact photo identity in abiotic ledger")
        current=current.drop(columns=[f"_join_{side}"])
    both_classified=current.both_endpoints_classifiable.astype(str).str.strip().str.casefold().isin(("true","1","yes"))
    declared=current.pair_state.isin(("same","discordant"))
    if not (both_classified==declared).all():
        raise RuntimeError("Immutable source photo classification and pair-state disagree")
    if int(both_classified.sum())!=BOTH_CLASSIFIABLE or int(current.pair_state.eq("discordant").sum())!=DISCORDANT:
        raise RuntimeError("Frozen 3196 classifiable /797 discordant pair counts drifted")
    pos=current.site_geo_status_1.eq(SOURCE_GEO)&current.site_geo_status_2.eq(SOURCE_GEO)
    climate=pos&current.environment_climate_complete_1.eq(True)&current.environment_climate_complete_2.eq(True)
    soil=pos&current.environment_soil_complete_1.eq(True)&current.environment_soil_complete_2.eq(True)
    complete=both_classified&climate&soil&current.environment_all_complete_1.eq(True)&current.environment_all_complete_2.eq(True)
    for feature in (ELEV,*CLIMATE,*SOIL):
        current["delta_"+feature]=(pd.to_numeric(current[f"{feature}_1"],errors="coerce")-
                                    pd.to_numeric(current[f"{feature}_2"],errors="coerce")).abs()
    current["delta_abs_latitude"]=(pd.to_numeric(current.latitude_1,errors="coerce").abs()-
                                   pd.to_numeric(current.latitude_2,errors="coerce").abs()).abs()
    havegeo=pos & current[["latitude_1","longitude_1","latitude_2","longitude_2"]].notna().all(axis=1)
    current["geodesic_distance_km"]=np.nan
    current["pair_midpoint_cell_162"]=-1
    subset=current.loc[havegeo]
    if len(subset):
        d,c=pair_distance_and_midpoint(subset.latitude_1.to_numpy(float),subset.longitude_1.to_numpy(float),
                                      subset.latitude_2.to_numpy(float),subset.longitude_2.to_numpy(float))
        current.loc[havegeo,"geodesic_distance_km"]=d
        current.loc[havegeo,"pair_midpoint_cell_162"]=c
    current["log_geodesic_distance_km"]=np.log1p(current.geodesic_distance_km)
    current["genus"]=current.species.fillna("").astype(str).str.split().str[0]
    allfeatures=list(BASE+CLIM_FEATURES+SOIL_FEATURES)
    keep=complete & current[allfeatures].notna().all(axis=1)&current.pair_midpoint_cell_162.ge(0)
    summary={
        "n_original_species_unique_cross_cell_pairs":len(current),
        "n_original_pairs_both_four_state_classified":int(both_classified.sum()),
        "n_original_pairs_discordant":int(current.pair_state.eq("discordant").sum()),
        "n_both_source_photo_locations_available":int(pos.sum()),
        "n_classified_pairs_with_both_original_locations":int((both_classified&pos).sum()),
        "n_classified_pairs_two_sites_climate_complete":int((both_classified&climate).sum()),
        "n_classified_pairs_two_sites_soil_complete":int((both_classified&soil).sum()),
        "n_classified_pairs_all_geo_climate_soil":int(keep.sum()),
        "n_complete_cases_source_photo_discordant":int((keep&current.pair_state.eq("discordant")).sum()),
        "n_original_pairs_missing_classified_response":int((~both_classified).sum()),
        "no_measured_pair_or_missing_soil_imputed":True,
    }
    return current.loc[keep].copy().reset_index(drop=True),summary


def oof_compare(complete:pd.DataFrame,group_name:str)->dict:
    from sklearn.model_selection import GroupKFold
    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline
    from sklearn.metrics import log_loss,roc_auc_score
    if group_name not in ("genus","pair_midpoint_cell_162"):
        raise ValueError("Only source-defined taxa or spatial-fold groups supported")
    group=complete[group_name].to_numpy()
    if len(set(group))<N_FOLDS:
        return {"group":group_name,"status":"INSUFFICIENT_SPATIAL_OR_GENUS_GROUPS"}
    y=complete.pair_state.eq("discordant").to_numpy(dtype=int)
    folds=list(GroupKFold(n_splits=N_FOLDS).split(complete,y,groups=group))
    if any(len(set(y[train]))<2 or y[train].sum()<20 for train,test in folds):
        return {"group":group_name,"status":"INSUFFICIENT_OBSERVED_DISCORDANT_TRAIN_FOLDS"}
    prediction={}
    metrics={}
    for name,features in FAMILIES.items():
        x=complete[list(features)].to_numpy(float)
        if not np.isfinite(x).all():raise ValueError("Soil-complete source unexpectedly has missing gradients")
        pred=np.full(len(x),np.nan)
        foldloss=[]
        for i,(train,test) in enumerate(folds):
            model=make_pipeline(StandardScaler(),LogisticRegression(max_iter=750,C=1.0,random_state=SEED))
            model.fit(x[train],y[train])
            prob=np.clip(model.predict_proba(x[test])[:,1],1e-6,1-1e-6)
            pred[test]=prob
            foldloss.append({"fold":i,"n_test":len(test),"n_heldout_groups":len(set(group[test])),
                "heldout_log_loss":float(log_loss(y[test],prob,labels=[0,1]))})
        if np.isnan(pred).any():raise RuntimeError("A pair received no held-out probability")
        prediction[name]=pred
        metrics[name]={
            "heldout_log_loss":float(log_loss(y,pred,labels=[0,1])),
            "heldout_AUC":float(roc_auc_score(y,pred)) if len(set(y))==2 else None,
            "fold_diagnostics":foldloss,
        }
    individual={k:-(y*np.log(p)+(1-y)*np.log1p(-p)) for k,p in prediction.items()}
    rng=np.random.default_rng(SEED+ (1 if group_name=="genus" else 2))
    unique,groupind=np.unique(group,return_inverse=True)
    count=np.bincount(groupind)
    gain_clim=individual["GEOGRAPHY"]-individual["GEOGRAPHY_CLIMATE"]
    gain_soil=individual["GEOGRAPHY_CLIMATE"]-individual["GEOGRAPHY_CLIMATE_SOIL"]
    boot={}
    for key,arr in (("climate",gain_clim),("soil",gain_soil)):
        scores=np.bincount(groupind,weights=arr)
        draws=rng.integers(0,len(unique),size=(N_BOOT,len(unique)))
        vals=scores[draws].sum(axis=1)/count[draws].sum(axis=1)
        boot[key]=[float(v) for v in np.quantile(vals,[.025,.975])]
    return {
        "group":group_name,"status":"FIXED_ORIGINAL_PAIRS_BLOCKED_OOF_DIAGNOSTIC",
        "n_pairs":len(complete),"n_heldout_groups":len(unique),
        "n_discordant_source_photo_pairs":int(y.sum()),
        "n_coincident_source_colour_pairs":int((1-y).sum()),
        "models":metrics,
        "climate_delta_outofgroup_logloss":float(gain_clim.mean()),
        "soil_increment_outofgroup_logloss":float(gain_soil.mean()),
        "group_bootstrap_95CI_logloss_gain_climate":boot["climate"],
        "group_bootstrap_95CI_logloss_gain_soil":boot["soil"],
        "bootstrap_explanation":"Conditional group bootstrap of fixed out-of-fold prediction losses; no fold/model refitting, no confirmation or causal inference",
    }


def analyze(pairs:pd.DataFrame,env:pd.DataFrame)->dict:
    complete,cover=join_origins(pairs,env)
    result={
        "schema":"fcp_42111_same_species_crosscell_colour_discordance_abiotic_v1",
        "date_jst":"2026-10-09",
        "status":"RETROSPECTIVE_SOURCE_PHOTO_PAIR_ENVIRONMENTAL_COVERAGE",
        "fixed_photo_pair_denominator":N_PAIRS,
        "original_classifiable_pairs":BOTH_CLASSIFIABLE,
        "original_colour_discordant_pairs":DISCORDANT,
        "coverage":cover,
        "same_species_photo_pairs_not_genetic_or_populational_morphs":True,
        "original_pairs_unclassifiable_preserved_in_denominator":True,
        "no_species_cell_centroid_used_for_abiotic_location":True,
        "no_new_colour_pixel_or_prospective_species_opened":True,
        "original_fcp_manuscript_decisions_changed":False,
        "models":[],
    }
    if len(complete)<MIN_COMPLETE_PAIRS or complete.pair_state.nunique()<2:
        result["status"]="HOLD_SOURCE_PAIR_SOIL_COMPLETE_COVERAGE_OR_COLOUR_CLASS"
    else:
        result["models"]=[oof_compare(complete,v) for v in ("genus","pair_midpoint_cell_162")]
        result["status"]="SOURCE_PAIRS_GENUS_AND_SPATIAL_BLOCKED_ENVIRONMENT_DIAGNOSTIC"
    result["hard_nonclaims"]=[
        "Colour mismatch across two photographed places is not proved genetic polymorphism in the same breeding population",
        "Original 13416 species-level fixed pair denominator differs from 85337 taxon-cell rows and 42111 species breadth",
        "Photos without classified labels, public points or soil are not coded as matching colours",
        "Geographic distance may correlate with environment even after covariate adjustment; this is not causal adaptation",
        "Group-held-out cross-validation estimates prediction portability for one fixed source pair/species only",
        "SoilGrids 5km represents modelled soil, not root-zone measurements or direct pollinator fitness",
    ]
    return result


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--pairs",required=True,type=Path)
    p.add_argument("--photo-environment",required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    a=p.parse_args()
    original=read_original_pairs(a.pairs)
    env=pd.read_csv(a.photo_environment,low_memory=False)
    d=analyze(original,env)
    a.outdir.mkdir(parents=True,exist_ok=True)
    (a.outdir/"result.json").write_text(json.dumps(d,indent=2,sort_keys=True)+"\n")
    print(json.dumps(d,indent=2,sort_keys=True),flush=True)
if __name__=="__main__":
    main()

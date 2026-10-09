#!/usr/bin/env python3
"""Within-genus BETWEEN-SPECIES flower-photo colour/abiotic predictive audit.

Uses original source 42,111 one-photo-per-species measurements. In geographic
cell-held-out folds, estimate genus intercepts from TRAIN species only; predict
only species in genera with >=2 different training species. Genus centring
separates a within-genus BETWEEN-SPECIES environmental component from broad
between-genus floristic turnover. It is NOT a within-species or causal test.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

CLASSES=("white","yellow_orange","red_pink","blue_purple")
GEO=("abs_latitude","lon_sin","lon_cos","wc_elevation_m")
CLIMATE=("wc_bio1","wc_bio5","wc_bio12","wc_bio15")
SOIL=("soil_pH","soil_SOC","soil_N","soil_clay","soil_available_water_proxy")
FAMILIES={
    "GENUS_BASELINE":(),
    "GENUS_GEO":GEO,
    "GENUS_GEO_CLIMATE":GEO+CLIMATE,
    "GENUS_GEO_CLIMATE_SOIL":GEO+CLIMATE+SOIL,
}
N_SPECIES_ORIGINAL=42111
N_CLASSIFIED_SOURCE=18457
EXPECTED_CLIMATE_SOURCE=18413
EXPECTED_SOIL_SOURCE=14136
N_FOLDS=5
BOOT=999
RIDGE=10.
RNG_SEED=20261009


def source_populations(source:pd.DataFrame,*,strict:bool=True)->tuple[dict[str,pd.DataFrame],dict]:
    req={"inat_taxon_id","species","morph","measurement_status","site_geo_status",
         "latitude","longitude","wc_elevation_m","environment_climate_complete",
         "environment_soil_complete",*CLIMATE,*SOIL}
    if not req.issubset(source):
        raise ValueError("Original global photo/environment source schema missing "+str(sorted(req-set(source))))
    if source.inat_taxon_id.duplicated().any():
        raise ValueError("Repeated global source species/photo ID")
    if strict and len(source)!=N_SPECIES_ORIGINAL:
        raise ValueError("Original 42111 single-photo species denominator changed")
    s=source.copy()
    yes=s.morph.isin(CLASSES)&s.measurement_status.eq("classified_four_state_morph")
    if strict and int(yes.sum())!=N_CLASSIFIED_SOURCE:
        raise ValueError("Original source 18457 photo classifications drift")
    located=s.site_geo_status.eq("VALID_PUBLIC_ORIGINAL_PHOTO_POINT")
    lat=pd.to_numeric(s.latitude,errors="coerce")
    lon=pd.to_numeric(s.longitude,errors="coerce")
    s["abs_latitude"]=lat.abs()
    s["lon_sin"]=np.sin(np.deg2rad(lon))
    s["lon_cos"]=np.cos(np.deg2rad(lon))
    s["genus"]=s.species.astype(str).str.strip().str.split().str[0]
    s["photo_cell_162"]=cell_index(lat.to_numpy(float),lon.to_numpy(float))
    good_coord=located & lat.between(-90,90) & lon.between(-180,180) & s.photo_cell_162.ge(0)
    base=yes&good_coord&s.environment_climate_complete.fillna(False).astype(bool)&s[list(GEO+CLIMATE)].notna().all(axis=1)
    soil=base&s.environment_soil_complete.fillna(False).astype(bool)&s[list(SOIL)].notna().all(axis=1)
    pools={"CLIMATE_ALL":s.loc[base].copy(),"CLIMATE_SOIL_COMPLETE":s.loc[soil].copy()}
    stats={
        "source_single_photo_species":len(s),
        "original_photograph_colour_classifiable":int(yes.sum()),
        "source_photo_geolocated_classified":int((yes&good_coord).sum()),
        "source_classified_climate_complete":int(base.sum()),
        "source_classified_climate_soil_complete":int(soil.sum()),
        "unclassified_photo_species_preserved_outside_model":int((~yes).sum()),
        "climate_only_additional_species_not_available_in_soil_complete":int(base.sum()-soil.sum()),
    }
    if strict and (stats["source_classified_climate_complete"]!=EXPECTED_CLIMATE_SOURCE or
                   stats["source_classified_climate_soil_complete"]!=EXPECTED_SOIL_SOURCE):
        raise ValueError("Frozen 18413/14136 source climate/soil support drifted")
    for name,d in pools.items():
        if d.empty or (d.genus=="").any() or d.inat_taxon_id.duplicated().any():
            raise ValueError(f"{name}: malformed source genus / taxon ID")
        if len(set(d.morph))!=4:
            raise ValueError(f"{name}: missing photo-colour class")
    return pools,stats


def cell_index(latitude:np.ndarray,longitude:np.ndarray)->np.ndarray:
    lat=np.asarray(latitude,float);lon=np.asarray(longitude,float)
    out=np.full(len(lat),-1,dtype=int)
    finite=np.isfinite(lat)&np.isfinite(lon)&(np.abs(lat)<=90)&(np.abs(lon)<=180)
    x=np.minimum(17,np.maximum(0,np.floor((lon[finite]+180)/20).astype(int)))
    y=np.minimum(8,np.maximum(0,np.floor((np.sin(np.deg2rad(lat[finite]))+1)*4.5).astype(int)))
    out[finite]=18*y+x
    return out


def training_genus_predictors(train:pd.DataFrame,test:pd.DataFrame,
                              features:tuple[str,...],train_min_taxa:int=2)->tuple[np.ndarray,np.ndarray]:
    """Predict held-out original photo labels using TRAIN genus means only."""
    if train.inat_taxon_id.duplicated().any() or test.inat_taxon_id.duplicated().any():
        raise ValueError("Original source species cannot repeat in train or test")
    if set(train.inat_taxon_id)&set(test.inat_taxon_id):
        raise ValueError("Original source taxon ID entered both train and geographic holdout")
    gtrain=train.genus.to_numpy(str)
    counts=pd.Series(gtrain).value_counts()
    eligible=np.flatnonzero(test.genus.isin(counts[counts>=train_min_taxa].index).to_numpy())
    if len(eligible)==0:
        return eligible,np.zeros((0,len(CLASSES)),float)
    # A genus label/count is estimated using only the different TRAIN species;
    # never calculate a genus mean from the held-out photo's true colour.
    y=pd.get_dummies(pd.Categorical(train.morph,categories=CLASSES)).to_numpy(float)
    if y.shape[1]!=len(CLASSES):raise ValueError("Lost source photograph colour class")
    yg=pd.DataFrame(y,columns=CLASSES).assign(genus=gtrain).groupby("genus",sort=False)[list(CLASSES)].mean()
    observed=test.iloc[eligible]
    ymean=yg.reindex(observed.genus.to_numpy(str)).to_numpy(float)
    if len(features)==0:
        return eligible,ymean
    xt=train[list(features)].to_numpy(float)
    xe=observed[list(features)].to_numpy(float)
    if not np.isfinite(xt).all() or not np.isfinite(xe).all():
        raise ValueError("Incomplete source environmental values in study population")
    mu=xt.mean(axis=0)
    sigma=xt.std(axis=0)
    sigma[sigma<1e-8]=1.0
    xt=(xt-mu)/sigma; xe=(xe-mu)/sigma
    xm=pd.DataFrame(xt,columns=list(features)).assign(genus=gtrain).groupby("genus",sort=False)[list(features)].mean()
    xmeantrain=xm.reindex(gtrain).to_numpy(float)
    ymeantrain=yg.reindex(gtrain).to_numpy(float)
    xc=xt-xmeantrain
    yc=y-ymeantrain
    beta=np.linalg.solve(xc.T@xc+RIDGE*np.eye(len(features)),xc.T@yc)
    pred=ymean+(xe-xm.reindex(observed.genus.to_numpy(str)).to_numpy(float))@beta
    # All compared linear probability models use the same projection into the
    # 4-state simplex, avoiding false probability/negative Brier interpretations.
    pred=np.maximum(pred,0)
    denom=pred.sum(axis=1,keepdims=True)
    pred=np.divide(pred,denom,out=np.full_like(pred,1/len(CLASSES)),where=denom>0)
    return eligible,pred


def fivefold_within_genus_cv(source:pd.DataFrame,*,with_soil:bool)->dict:
    from sklearn.model_selection import GroupKFold
    if len(source)<100 or source.photo_cell_162.nunique()<N_FOLDS:
        raise ValueError("Too few original source species/spatial cells for fivefold")
    if source.inat_taxon_id.duplicated().any():raise ValueError("Repeated source taxa")
    features={k:v for k,v in FAMILIES.items() if with_soil or k!="GENUS_GEO_CLIMATE_SOIL"}
    n=len(source)
    y=pd.Categorical(source.morph,categories=CLASSES).codes
    if (y<0).any():raise ValueError("Non-classifiable photo entered ecological group")
    ys=np.eye(len(CLASSES))[y]
    pred={k:np.full((n,len(CLASSES)),np.nan,float) for k in features}
    eligible=np.zeros(n,bool)
    folds=[]
    for fold,(itr,ite) in enumerate(GroupKFold(n_splits=N_FOLDS).split(source,y,source.photo_cell_162)):
        tr=source.iloc[itr]
        te=source.iloc[ite]
        chosen=None
        for name,feat in features.items():
            ix,values=training_genus_predictors(tr,te,feat)
            if chosen is None:chosen=ix
            if not np.array_equal(chosen,ix):
                raise RuntimeError("Compared models switched included test species")
            pred[name][ite[ix]]=values
        eligible[ite[chosen]]=True
        folds.append({"fold":fold,"n_test_source_species":len(ite),
                      "n_genus_estimable_test_species":len(chosen),
                      "n_unseen_or_single_training_genus_test_species":len(ite)-len(chosen),
                      "n_test_source_geographic_cells":int(te.photo_cell_162.nunique())})
    if int(eligible.sum())<100:raise RuntimeError("Within-genus geographical prediction lacks test support")
    if len({k for k in pred if np.isfinite(pred[k][eligible]).all()})!=len(features):
        raise RuntimeError("Missing predictions among fixed within-genus evaluation rows")
    # Secondary equal-GENUS weight gives each test-supported lineage an equal
    # contribution, rather than allowing speciose genera to dominate the score.
    eval_genus=source.loc[eligible,"genus"].to_numpy(str)
    distinct_genus,genus_index=np.unique(eval_genus,return_inverse=True)
    genus_counts=np.bincount(genus_index)
    genus_weights=1.0/genus_counts[genus_index]
    scores={}
    loss={}
    for name,prob in pred.items():
        z=np.sum((prob[eligible]-ys[eligible])**2,axis=1)
        scores[name]={
            "n_same_source_test_species":int(eligible.sum()),
            "heldout_multiclass_brier":float(np.mean(z)),
            "heldout_genus_equal_multiclass_brier":float(np.average(z,weights=genus_weights)),
        }
        loss[name]=z
    keys=list(features)
    delta={}
    comparisons=[("GENUS_BASELINE","GENUS_GEO","geo_beyond_genus"),
                 ("GENUS_GEO","GENUS_GEO_CLIMATE","climate_beyond_genus_geography")]
    if with_soil:
        comparisons.append(("GENUS_GEO_CLIMATE","GENUS_GEO_CLIMATE_SOIL","soil_beyond_genus_geography_climate"))
    groups=source.loc[eligible,"photo_cell_162"].to_numpy(int)
    unique,gindex=np.unique(groups,return_inverse=True)
    count=np.bincount(gindex)
    for base,added,label in comparisons:
        gain=loss[base]-loss[added]
        v=float(gain.mean())
        sums=np.bincount(gindex,weights=gain)
        rng=np.random.default_rng(RNG_SEED+len(label))
        draw=rng.integers(0,len(unique),size=(BOOT,len(unique)))
        boot=sums[draw].sum(axis=1)/count[draw].sum(axis=1)
        # Two distinct sensitivity estimands, always on exactly the same OOF
        # species. Cell reweighting retains spatial grouping; genus resampling
        # tests the sensitivity to unequal genus richness in the photo sample.
        weighted_sums=np.bincount(gindex,weights=genus_weights*gain)
        weighted_denoms=np.bincount(gindex,weights=genus_weights)
        weighted_boot=weighted_sums[draw].sum(axis=1)/weighted_denoms[draw].sum(axis=1)
        genus_means=np.bincount(genus_index,weights=gain)/genus_counts
        genus_rng=np.random.default_rng(RNG_SEED+len(label)+100)
        gdraw=genus_rng.integers(0,len(distinct_genus),size=(BOOT,len(distinct_genus)))
        gboot=genus_means[gdraw].mean(axis=1)
        delta[label]={
            "mean_heldout_brier_reduction":v,
            "source_cell_block_bootstrap_95CI":[float(x) for x in np.quantile(boot,[.025,.975])],
            "positive_gain_supported_by_bootstrap":bool(np.quantile(boot,.025)>0),
            "genus_equal_heldout_brier_reduction":float(np.average(gain,weights=genus_weights)),
            "genus_equal_source_cell_bootstrap_95CI":[float(x) for x in np.quantile(weighted_boot,[.025,.975])],
            "genus_equal_genus_cluster_bootstrap_95CI":[float(x) for x in np.quantile(gboot,[.025,.975])],
            "fraction_evaluated_genera_with_positive_increment":float(np.mean(genus_means>0)),
            "n_genera_in_genus_equal_sensitivity":int(len(distinct_genus)),
        }
    return {
        "schema":"fcp_global42111_train_genus_only_within_genus_cell_holdout_v1",
        "source_species_in_pool":n,
        "n_original_source_genera":int(source.genus.nunique()),
        "n_original_source_162_geographic_cells":int(source.photo_cell_162.nunique()),
        "source_species_with_train_observed_genus_and_min_two_training_species":int(eligible.sum()),
        "n_distinct_source_genera_in_evaluation":int(source.loc[eligible,"genus"].nunique()),
        "n_distinct_original_cells_in_evaluation":int(len(unique)),
        "source_species_not_prediction_eligible_due_to_genus_support":int((~eligible).sum()),
        "same_exact_species_and_cell_folds_for_all_feature_families":True,
        "genus_means_estimated_only_from_training_source_species":True,
        "one_original_photo_per_taxon":True,
        "heldout_5fold_original_cell_status":folds,
        "source_colour_classes":list(CLASSES),
        "brier_score_models":scores,
        "fixed_fold_incremental_gains":delta,
        "block_bootstrap_conditional_on_one_realized_fivefold_cv":True,
        "equal_genus_weight_secondary_sensitivity":True,
        "genus_bootstrap_conditional_on_fixed_outofcell_predictions":True,
        "inference_boundary":"Between different species within the same nominal genus; does NOT estimate within-species evolution or a causal climate or soil selection effect",
    }


def run(source:pd.DataFrame,*,strict:bool=True)->dict:
    pops,coverage=source_populations(source,strict=strict)
    result={
        "schema":"fcp_global42111_within_genus_compositional_vs_climate_soil_v1",
        "date_jst":"2026-10-09",
        "status":"RETROSPECTIVE_OBSERVATIONAL_BETWEEN_SPECIES_WITHIN_GENUS",
        "original_global_species_denominator":len(source),
        "original_historic_colour_classifiable":coverage["original_photograph_colour_classifiable"],
        "original_source_photo_coverage":coverage,
        "climate_only_population":fivefold_within_genus_cv(pops["CLIMATE_ALL"],with_soil=False),
        "soil_complete_population":fivefold_within_genus_cv(pops["CLIMATE_SOIL_COMPLETE"],with_soil=True),
        "no_photo_pixels_or_new_colour_labels_opened":True,
        "no_environment_soil_imputation_or_cell_centroid_as_site":True,
        "old_fcp_highdepth_and_new_independent_species_unchanged":True,
        "hard_nonclaims":[
            "This is BETWEEN different species of the same named genus, not within-species photo-colour evolution",
            "Genus fixed effects do not control full evolutionary phylogeny, species turnover inside genera or photo misclassification",
            "Cell-heldout evaluation restricts to genera with two or more TRAIN species, so the support population is conditional",
            "Genus means always come from training species only; no test-label contamination",
            "A source photograph is not a species modal flower colour or genetically measured morph",
            "A small apparent predictive gain is not causation, selection or reproductive fitness",
            "Unequal soil-mask/classifiability across the original 42111 named species limits generalization",
        ],
    }
    return result


def main()->None:
    ap=argparse.ArgumentParser()
    ap.add_argument("--original-breadth-abiotic",required=True,type=Path)
    ap.add_argument("--outdir",required=True,type=Path)
    a=ap.parse_args()
    source=pd.read_csv(a.original_breadth_abiotic,low_memory=False)
    result=run(source)
    a.outdir.mkdir(parents=True,exist_ok=True)
    (a.outdir/"result.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2,sort_keys=True),flush=True)

if __name__=="__main__":
    main()

#!/usr/bin/env python3
"""Exploratory all-environment ecological prediction beyond TRUE LCVP patristic.

Only exact directly-tipped 342/649 original source species (250/500km parent
cohorts), and photo-site <=50/100km congeneric dyads. Reuse previously checked
dated tip tree, source photo colours and 100543-site expanded WorldClim/soil.
Always hold out whole genus in 5 folds, compare 7 blocks on SAME dyads, emit
HOLD below source support. This is no phylo-GLMM or genetic evolution claim.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo

from assess_fcp_direct_LCVP_phylogeny_climate_dyads_20261010 import original_dyads
from audit_fcp_global42111_direct_phylo_nearby_photo_pairs_20261010 import direct_taxa
from audit_fcp_global42111_local_congener_moisture_null_20261009 import fixed_source_cohort
from audit_fcp_global42111_local_congeners_20261009 import photo_populations
from compare_fcp_global42111_full_climate_sun_soil_blocks_20261010 import BLOCKS

SCHEMA="fcp_42111_direct_LCVP_photo_dyads_full_abiotic_source_cv_v1"
CAPS=(250,500)
MICRO=(50,100)
PHOTO_N={"250":342,"500":649}
MIN_DYADS=100
MIN_GEO_GENUS=20
FOLDS=5
BOOT=499
SEED=20261010

BASE=("log_km","delta_abs_latitude","log_LCVP_patristic")
PHOTO_VARIABLES=tuple(x for block in BLOCKS.values() for x in block)
DYAD_BLOCKS={name:tuple("delta_"+x for x in variables) for name,variables in BLOCKS.items()}
# Original local photo script uses delta_elevation, not delta_wc_elevation_m.
DYAD_BLOCKS["elevation"]=("delta_elevation",)


def build_source_dyads(direct:pd.DataFrame,tree,diameter:int)->tuple[pd.DataFrame,dict]:
    pair,coverage=original_dyads(direct,tree,diameter)
    if pair.empty:
        return pair,coverage
    source=direct.set_index("inat_taxon_id",verify_integrity=True)
    for name,cols in BLOCKS.items():
        for feature in cols:
            key="delta_elevation" if feature=="wc_elevation_m" else "delta_"+feature
            if key in pair.columns:
                continue
            if feature not in source:
                raise ValueError("Expanded full photo source missing "+feature)
            a=source.loc[pair.taxon_1.to_numpy(int),feature].to_numpy(float)
            b=source.loc[pair.taxon_2.to_numpy(int),feature].to_numpy(float)
            pair[key]=np.abs(a-b)
    return pair,coverage


def feature_families(*,soil:bool)->dict[str,tuple[str,...]]:
    blocks={k:v for k,v in DYAD_BLOCKS.items() if soil or k!="soil"}
    allfeatures=BASE+tuple(x for fs in blocks.values() for x in fs)
    out={"GEO_TRUE_LCVP_PHYLO_ONLY":BASE,
         "FULL_TRUE_LCVP_PHYLO_PLUS_ALL_ABIOTIC":allfeatures}
    for k,v in blocks.items():
        out["FULL_MINUS_"+k.upper()]=tuple(x for x in allfeatures if x not in v)
    return out


def blocked_models(d:pd.DataFrame,*,soil:bool,nboot:int=BOOT)->dict:
    from sklearn.model_selection import GroupKFold
    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline
    y=d.colour_photo_mismatch.to_numpy(int)
    genus=d.genus.to_numpy(str)
    n=len(d)
    if n<MIN_DYADS or len(set(genus))<MIN_GEO_GENUS or len(set(y))!=2:
        return {"status":"HOLD_DIRECT_PHYLOGENETIC_DYADS_OR_GENERA_INSUFFICIENT",
                "n_selected_original_photo_dyads":n}
    fold=list(GroupKFold(n_splits=FOLDS).split(d,y,groups=genus))
    if any(len(set(y[train]))!=2 for train,_ in fold):
        return {"status":"HOLD_DIRECT_PHYLOGENETIC_TRAINING_FOLD_CLASS",
                "n_selected_original_photo_dyads":n}
    families=feature_families(soil=soil)
    losses={}
    metrics={}
    for modelname,cols in families.items():
        X=d[list(cols)].to_numpy(float)
        if not np.isfinite(X).all():
            raise ValueError("Source complete-case dyadic abiotic model has missing feature")
        proba=np.full(n,np.nan)
        for train,test in fold:
            fit=make_pipeline(StandardScaler(),LogisticRegression(
                max_iter=1200,C=1.0,random_state=SEED))
            fit.fit(X[train],y[train])
            proba[test]=np.clip(fit.predict_proba(X[test])[:,1],1e-6,1-1e-6)
        if not np.isfinite(proba).all():
            raise RuntimeError("Original source dyad lacks heldout model prediction")
        l=-(y*np.log(proba)+(1-y)*np.log1p(-proba))
        losses[modelname]=l
        metrics[modelname]=float(np.mean(l))
    ref=losses["FULL_TRUE_LCVP_PHYLO_PLUS_ALL_ABIOTIC"]
    comparisons={}
    for name,base in [("all_abiotic_beyond_geo_phylogeny","GEO_TRUE_LCVP_PHYLO_ONLY")]+[
        ("conditional_"+k,"FULL_MINUS_"+k.upper()) for k in DYAD_BLOCKS if soil or k!="soil"]:
        gains=losses[base]-ref
        uniq,idx=np.unique(genus,return_inverse=True)
        weighted=np.bincount(idx,weights=gains)
        count=np.bincount(idx)
        rng=np.random.default_rng(SEED+len(name))
        draw=rng.integers(0,len(uniq),size=(nboot,len(uniq)))
        sampling=weighted[draw].sum(axis=1)/count[draw].sum(axis=1)
        ci=np.quantile(sampling,[.025,.975])
        comparisons[name]={
            "heldout_original_dyad_logloss_reduction":float(np.mean(gains)),
            "original_genus_bootstrap_95CI_conditional_fixed_models":[float(k) for k in ci],
            "ci_nominally_positive":bool(ci[0]>0),
        }
    return {
        "status":"EXPLORATORY_TRUE_PHYLO_CONDITIONAL_MULTIENVIRONMENT_DYAD_ASSOCIATION",
        "n_same_original_dyads_each_model":n,
        "n_unique_source_heldout_genera":len(set(genus)),
        "n_mismatch_original_photo_dyads":int(y.sum()),
        "all_models_same_original_dyads_and_genus_holdout_folds":True,
        "models_genus_holdout_binary_logloss":metrics,
        "all_block_conditional_incremental_gains":comparisons,
        "intervals_fix_one_realized_genus_holdout_fit_and_shared_dyad_dependencies_remain":True,
    }


def run(original:pd.DataFrame,ledger:pd.DataFrame,tree250:Path,tree500:Path,
        *,strict=True,nboot:int=BOOT)->dict:
    if strict and (len(original)!=42111 or original.inat_taxon_id.nunique()!=42111):
        raise ValueError("Original 42111 source photographed taxon identities changed")
    pools,support=photo_populations(original,strict=strict)
    outputs={}
    for cap,tree_file in ((250,tree250),(500,tree500)):
        full=fixed_source_cohort(pools["CLIMATE_ALL"],cap)
        direct,tree=direct_taxa(full,ledger,tree_file,cap,strict=strict)
        if strict and len(direct)!=PHOTO_N[str(cap)]:
            raise RuntimeError("Original direct LCVP taxon photo source drift")
        cohorts={}
        for diameter in MICRO:
            paired,cover=build_source_dyads(direct,tree,diameter)
            needed=tuple(x for k,v in DYAD_BLOCKS.items() if k!="soil" for x in v)
            extra=tuple(DYAD_BLOCKS["soil"])
            climate=paired[needed].notna().all(axis=1) if len(paired) else pd.Series([],dtype=bool)
            soil=climate&paired[list(extra)].notna().all(axis=1) if len(paired) else pd.Series([],dtype=bool)
            a=paired.loc[climate].copy()
            b=paired.loc[soil].copy()
            cohorts[str(diameter)]={
                "n_original_direct_LCVP_species":len(direct),
                "n_original_photo_pairs_all_original_4class":len(paired),
                "n_original_photo_pairs_full_climate_solar_wind_vapor":len(a),
                "n_original_photo_pairs_full_climate_and_soil":len(b),
                "source_original_pair_coverage":cover,
                "climate_only":blocked_models(a,soil=False,nboot=nboot),
                "all_environment_with_soil":blocked_models(b,soil=True,nboot=nboot),
            }
        outputs[str(cap)]={
            "original_parent_group_max_photo_diameter_km":cap,
            "n_original_LCVP_direct_tips":len(direct),
            "local_original_photo_dyads":cohorts,
        }
    return {
        "schema":SCHEMA,
        "date_jst":"2026-10-10",
        "status":"SOURCE_EXPOSED_DIRECT_LCVP_ABIOTIC_PHOTOGRAPHIC_DYAD_SENSITIVITY",
        "source_original_species_photo_denominator":len(original),
        "original_global_colour_classified_one_photo_species":support["source_classified_photo_taxa"],
        "n_environment_blocks":len(BLOCKS),
        "environment_blocks":{k:list(v) for k,v in BLOCKS.items()},
        "source_true_direct_LCVP_tips_only":True,
        "source_250_500_original_photo_dyads_nested_not_independent":True,
        "original_full_source_phylogenetic_50pct_tip_coverage_HOLD":True,
        "source_original_4class_photograph_labels_not_remeasured":True,
        "results":outputs,
        "hard_nonclaims":[
            "Genus-blocked logistic prediction on repeated-species dyads is not a phylogenetic generalized linear model",
            "Each source taxon can enter several dyads, so dyadic likelihood independence fails",
            "Direct LCVP tip subset covers only 37-39pct of original full local species cohorts and is climatically biased",
            "Original full source geographic 50/100km label-shuffle null does not support strong local precipitation-adaptation claim",
            "WorldClim incident solar, vapor and wind are climatic modeled exposures; soil not direct root-zone conditions",
            "Exploratory correlated 7-block comparisons yield unadjusted fixed-fit bootstrap intervals, not independent hypothesis confirmation",
        ],
    }


def main():
    a=argparse.ArgumentParser()
    a.add_argument("--expanded-breadth",type=Path,required=True)
    a.add_argument("--original-direct-tip-ledger",type=Path,required=True)
    a.add_argument("--tree-250",type=Path,required=True)
    a.add_argument("--tree-500",type=Path,required=True)
    a.add_argument("--outdir",type=Path,required=True)
    v=a.parse_args()
    report=run(pd.read_csv(v.expanded_breadth,low_memory=False),
               pd.read_csv(v.original_direct_tip_ledger,low_memory=False),v.tree_250,v.tree_500)
    v.outdir.mkdir(parents=True,exist_ok=True)
    (v.outdir/"result.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print(json.dumps(report,sort_keys=True),flush=True)


if __name__=="__main__":
    main()

#!/usr/bin/env python3
"""Exploratory direct-tip within-genus photo-colour mismatch prediction.

Only original dated LCVP species tips and original 42111-source measured
flower-photo colours + genuine photo-site geography and WorldClim. The
250/500km parent cohorts and <=50/100km source microgroups are fixed before
this run; the nested cohorts are not independent replications. Compare
within-microgroup DIFFERENT-species photo mismatch on geographic/thermal
distance, then dated-tree patristic distance, then rainfall difference.
Leave whole GENERA out of training, bootstrap genus-level held-out losses.
Do not infer phenotype genetics, climate causality or adaptation.
"""
from __future__ import annotations
import argparse
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo

from audit_fcp_global42111_direct_phylo_nearby_photo_pairs_20261010 import (
    canonical,direct_taxa,MIN_DIRECT_PHYLO_CONNECTED_TAXA,
)
from audit_fcp_global42111_local_congener_moisture_null_20261009 import (
    fixed_source_cohort,EXPECTED_ORIGINAL_N,
)
from audit_fcp_global42111_local_congeners_20261009 import photo_populations
from audit_fcp_global42111_microspatial_photo_shuffle_20261009 import (
    microgroups,great_circle_matrix_km
)

SOURCE_COHORTS=(250,500)
LOCAL_DIAMETERS=(50,100)
TRUE_DIRECT_TIP_N={"250":342,"500":649}
SOURCE_CLASSIFIABLE=18457
SOURCE_GLOBAL=42111
SOURCE_CLIMATE=18413
MIN_SPECIES=100
MIN_LOCAL_GROUPS=30
MIN_GENERA=20
MIN_PAIRS=100
FOLDS=5
BOOT=999
SEED=20261010
SCHEMA="fcp_42111_direct_LCVP_local_congener_colour_climate_prediction_v1"

BASE=("log_km","delta_abs_latitude","delta_elevation","delta_wc_bio1","delta_wc_bio5")
PHYLO=("log_LCVP_patristic",)
RAIN=("delta_wc_bio12","delta_wc_bio15")
MODELS={
    "GEO_THERMAL":BASE,
    "GEO_THERMAL_PHYLO":BASE+PHYLO,
    "GEO_THERMAL_RAIN":BASE+RAIN,
    "GEO_THERMAL_PHYLO_RAIN":BASE+PHYLO+RAIN,
}
COMPARE={
    "phylogeny_beyond_geo_thermal":("GEO_THERMAL","GEO_THERMAL_PHYLO"),
    "rain_beyond_geo_thermal_without_phylogeny":("GEO_THERMAL","GEO_THERMAL_RAIN"),
    "rain_beyond_geo_thermal_and_dated_phylogeny":("GEO_THERMAL_PHYLO","GEO_THERMAL_PHYLO_RAIN"),
}


def original_dyads(source:pd.DataFrame,tree,diameter:int)->tuple[pd.DataFrame,dict]:
    """Build dyads with both taxon IDs, colour states and real LCVP patristic
    paths, never use photo colour for inclusion or source group selection."""
    required={"inat_taxon_id","genus","genus_cell_id","morph","latitude","longitude",
              "wc_elevation_m","wc_bio1","wc_bio5","wc_bio12","wc_bio15",
              "original_direct_LCVP_tip"}
    if not required.issubset(source):raise ValueError("Missing old original photo/phylo inputs")
    if source.index.tolist()!=list(range(len(source))):raise ValueError("Original source IDs require zero index")
    if source.inat_taxon_id.duplicated().any():raise ValueError("Same old photographed species duplicated")
    micro=microgroups(source,float(diameter))
    records=[]
    represented=set()
    represented_genera=set()
    useful_groups=0
    variable_groups=0
    n_group_three=0
    varied_taxa=set()
    for ii,ix in enumerate(micro):
        if len(ix)<2:continue
        p=source.iloc[ix]
        if p.genus.nunique()!=1 or p.genus_cell_id.nunique()!=1:
            raise RuntimeError("Source local phylogenetic microgroup mixed original taxa")
        useful_groups+=1
        if len(ix)>=3:n_group_three+=1
        local_name=str(p.genus_cell_id.iloc[0])+"|micro_"+str(ii)
        positions=p[["latitude","longitude"]].to_numpy(float)
        distances=great_circle_matrix_km(positions[:,0],positions[:,1])
        branches=[]
        for a,b in itertools.combinations(range(len(ix)),2):
            q=p.iloc[a];r=p.iloc[b]
            km=float(distances[a,b])
            if km>diameter+1e-7:raise RuntimeError("Source microphoto pair exceeds capped great-circle distance")
            pat=float(tree.distance(q.original_direct_LCVP_tip,r.original_direct_LCVP_tip))
            if not np.isfinite(pat) or pat<0:
                raise RuntimeError("Invalid original direct LCVP dated path distance")
            branches.append(pat)
            records.append({
                "local_group":local_name,
                "genus":str(q.genus),
                "photo_cell":int(q.photo_cell_162),
                "taxon_1":int(q.inat_taxon_id),"taxon_2":int(r.inat_taxon_id),
                "colour_photo_mismatch":int(str(q.morph)!=str(r.morph)),
                "log_km":float(np.log1p(km)),
                "delta_abs_latitude":float(abs(abs(q.latitude)-abs(r.latitude))),
                "delta_elevation":float(abs(q.wc_elevation_m-r.wc_elevation_m)),
                "delta_wc_bio1":float(abs(q.wc_bio1-r.wc_bio1)),
                "delta_wc_bio5":float(abs(q.wc_bio5-r.wc_bio5)),
                "delta_wc_bio12":float(abs(q.wc_bio12-r.wc_bio12)),
                "delta_wc_bio15":float(abs(q.wc_bio15-r.wc_bio15)),
                "log_LCVP_patristic":float(np.log1p(pat)),
            })
        represented.update(p.inat_taxon_id.astype(int).tolist())
        represented_genera.add(str(p.genus.iloc[0]))
        if len(ix)>=3 and len(set(round(x,5) for x in branches))>=2:
            variable_groups+=1
            varied_taxa.update(p.inat_taxon_id.astype(int).tolist())
    table=pd.DataFrame.from_records(records,columns=[
        "local_group","genus","photo_cell","taxon_1","taxon_2",
        "colour_photo_mismatch",*MODELS["GEO_THERMAL_PHYLO_RAIN"]])
    if len(table)>0 and table[["taxon_1","taxon_2"]].duplicated().any():
        raise RuntimeError("The same pair was counted in multiple source microgroups")
    checks={
        "n_source_direct_tips":len(source),
        "n_original_microgroups_all_sizes":len(micro),
        "n_comparable_local_two_plus_tip_groups":useful_groups,
        "n_comparable_local_three_plus_tip_groups":n_group_three,
        "n_three_plus_local_groups_with_variable_dated_phylogeny":variable_groups,
        "n_original_species_in_local_direct_tip_comparisons":len(represented),
        "n_original_species_in_locally_variable_phylo_groups":len(varied_taxa),
        "n_original_genera_in_local_tip_comparisons":len(represented_genera),
        "n_original_congeneric_photo_pair_dyads":len(table),
        "n_photo_colour_mismatched_dyads":int(table.colour_photo_mismatch.sum()) if len(table) else 0,
        "all_pairs_direct_LCVP_tips":True,
        "all_pair_locations_from_original_photo_coordinates":True,
        "group_selection_did_not_use_photo_colour_or_rainfall":True,
    }
    return table,checks


def predictive_check(d:pd.DataFrame)->dict:
    from sklearn.model_selection import GroupKFold
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import make_pipeline
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import log_loss,roc_auc_score
    y=d.colour_photo_mismatch.to_numpy(int)
    group=d.genus.to_numpy(str)
    if len(d)<MIN_PAIRS or len(set(group))<MIN_GENERA or len(np.unique(y))<2:
        return {"status":"HOLD_INSUFFICIENT_ORIGINAL_DIRECT_PHYLO_DYAD_PREDICTION"}
    folds=list(GroupKFold(n_splits=FOLDS).split(d,y,groups=group))
    if any(len(set(y[train]))!=2 for train,_ in folds):
        return {"status":"HOLD_SOURCE_PHYLO_TRAIN_FOLD_MISSING_DISCORDANCE_CLASS"}
    pred={}
    loss={}
    metrics={}
    for name,features in MODELS.items():
        x=d[list(features)].to_numpy(float)
        if not np.isfinite(x).all():
            raise RuntimeError("Original photo phylogenetic or rainfall feature missing")
        p=np.full(len(d),np.nan,float)
        for train,test in folds:
            model=make_pipeline(StandardScaler(),LogisticRegression(
                C=1.0,max_iter=1000,random_state=SEED))
            model.fit(x[train],y[train])
            p[test]=np.clip(model.predict_proba(x[test])[:,1],1e-6,1-1e-6)
        if not np.isfinite(p).all():raise RuntimeError("Some original congeneric photo dyads not held out")
        pred[name]=p
        loss[name]=-(y*np.log(p)+(1-y)*np.log1p(-p))
        metrics[name]={
            "n_same_heldout_source_dyads":len(d),
            "genus_heldout_binary_logloss":float(log_loss(y,p,labels=[0,1])),
            "genus_heldout_AUC":float(roc_auc_score(y,p)),
        }
    distinct,idx=np.unique(group,return_inverse=True)
    count=np.bincount(idx)
    delta={}
    for nm,(a,b) in COMPARE.items():
        gain=loss[a]-loss[b]
        weight=np.bincount(idx,weights=gain)
        rng=np.random.default_rng(SEED+len(nm))
        draws=rng.integers(0,len(distinct),size=(BOOT,len(distinct)))
        samples=weight[draws].sum(axis=1)/count[draws].sum(axis=1)
        delta[nm]={
            "heldout_original_dyad_logloss_reduction":float(np.mean(gain)),
            "original_genus_block_bootstrap_95CI":[float(z) for z in np.quantile(samples,[.025,.975])],
            "positive_gain_supported_under_fixed_genus_bootstrap":bool(np.quantile(samples,.025)>0),
        }
    return {
        "status":"EXPLORATORY_DIRECT_DATED_PHYLO_CLIMATE_PAIRWISE_PREDICTION",
        "n_original_photo_pair_dyads":len(d),
        "n_photo_mismatched_dyads":int(y.sum()),
        "n_actual_original_genera_held_out":int(len(distinct)),
        "same_original_dyads_and_genus_cv_folds_all_models":True,
        "models":metrics,"incremental_prediction":delta,
        "warning":"Dyads share source species, genus-bootstrap of fixed predictions does not model all dyadic/phylogenetic dependence; no confirmatory p-value",
    }


def assess(original:pd.DataFrame,ledger:pd.DataFrame,trees:dict[int,Path],*,strict=True)->dict:
    pops,source=photo_populations(original,strict=strict)
    results={}
    for cap in (250,500):
        full=fixed_source_cohort(pops["CLIMATE_ALL"],cap)
        if strict and len(full)!=EXPECTED_ORIGINAL_N[str(cap)]:
            raise RuntimeError("Original photographed local genus source population altered")
        direct,tree=direct_taxa(full,ledger,trees[cap],cap,strict=strict)
        if strict and len(direct)!={"250":342,"500":649}[str(cap)]:
            raise RuntimeError("Original LCVP direct species tip support changed")
        geographic={}
        for radius in (50,100):
            pairs,cover=original_dyads(direct,tree,radius)
            ready=(cover["n_original_species_in_local_direct_tip_comparisons"]>=MIN_SPECIES
                   and cover["n_comparable_local_two_plus_tip_groups"]>=MIN_LOCAL_GROUPS
                   and cover["n_original_genera_in_local_tip_comparisons"]>=MIN_GENERA
                   and cover["n_original_congeneric_photo_pair_dyads"]>=MIN_PAIRS)
            geographic[str(radius)]={
                "diameter_km":radius,
                "coverage":cover,
                "support_gate":"PASS_SOURCE_LOCAL_DIRECT_PHYLOGENETIC_DYAD_OPPORTUNITY" if ready else "HOLD_LOCAL_DIRECT_PHYLOGENY_SAMPLE_COVERAGE",
                "models":predictive_check(pairs) if ready else {"status":"NOT_FIT_BECAUSE_COVERAGE_HOLD"},
            }
        results[str(cap)]={
            "source_original_local_group_diameter_km":cap,
            "n_original_250_or_500km_photo_species":len(full),
            "n_source_species_directly_tipped_in_LCVP":len(direct),
            "subcohorts":geographic,
        }
    return {
        "schema":SCHEMA,
        "date_jst":"2026-10-10",
        "status":"RETROSPECTIVE_EXPLORATORY_DIRECT_LCVP_DYADIC_PREDICTION_ONLY",
        "global_source_original_42111_photo_species":len(original),
        "global_original_photo_fourstate_classifiable":source["source_classified_photo_taxa"],
        "source_250_and_500km_cohorts_nested_not_independent":True,
        "source_same_species_3182_photo_pairs_not_reanalysed":True,
        "source_direct_LCVP_only_no_synthetic_genus_grafting":True,
        "all_original_photo_ids_and_colour_classes_unchanged":True,
        "no_new_photographic_phenotypes_or_observations_opened":True,
        "full_original_250_500km_phylogenetic_inference_remains_HOLD":True,
        "cohorts":results,
        "hard_nonclaims":[
            "Original photographed flower-colour mismatches are not confirmed heritable morph states or within-population polymorphisms",
            "A congeneric local original source pair may consist of two distantly related subclades and multiple habitats",
            "Old LCVP direct-tip subset is nonrandom and not a complete 1761 species phylogeny",
            "Group-heldout pairwise log-loss is not an independent pair likelihood due to taxa shared across dyads",
            "Conditional genus bootstrap intervals are not family tree uncertainty propagation",
            "Climate differences are observational and cannot establish selective or causal flower-pigment evolution",
        ],
    }


def main():
    p=argparse.ArgumentParser()
    for field in ("original-breadth-abiotic","original-LCVP-tip-ledger","tree-250","tree-500","outdir"):
        p.add_argument("--"+field,type=Path,required=True)
    x=p.parse_args()
    source=pd.read_csv(x.original_breadth_abiotic,low_memory=False)
    ledger=pd.read_csv(x.original_LCVP_tip_ledger,low_memory=False)
    result=assess(source,ledger,{250:x.tree_250,500:x.tree_500})
    x.outdir.mkdir(parents=True,exist_ok=True)
    (x.outdir/"result.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,indent=2,sort_keys=True),flush=True)


if __name__=="__main__":
    main()

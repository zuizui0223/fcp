#!/usr/bin/env python3
"""True spatial and true dated-LCVP covariance baseline vs abiotic flower colour.

Uses only original one-photo nominal taxa, a source-matched dated LCVP tip
subtree, original photo positions and frozen environment. Species-equal
four-colour kernel ridge predictions. Never label a genus fixed effect as
phylogenetic covariance or simple lat/lon as complete space adjustment.
Brownian shared-branch matrix and great-circle exponential spatial kernel
with fixed 50/250km ranges; 5fold whole-genus and whole-photo-cell holdouts.
Additional environment tested beyond JOINT spatial+phylogenetic kernels.
Only 342/649 direct tips; full 872/1761 ecological source remains HOLD.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from Bio import Phylo

from audit_fcp_global42111_direct_phylo_nearby_photo_pairs_20261010 import direct_taxa
from audit_fcp_global42111_local_congener_moisture_null_20261009 import fixed_source_cohort
from audit_fcp_global42111_local_congeners_20261009 import photo_populations

SCHEMA="fcp_42111_direct_LCVP_spatial_and_phylo_kernel_environment_v1"
CLASSES=("white","yellow_orange","red_pink","blue_purple")
GEO=("abs_lat","longitude_sin","longitude_cos")
BLOCKS={
    "elevation":("wc_elevation_m",),
    "temperature":("wc_bio1","wc_bio5"),
    "precipitation":("wc_bio12","wc_bio15"),
    "soil":("soil_pH","soil_SOC","soil_N","soil_clay","soil_available_water_proxy"),
}
N_FOLDS=5
MIN_SOURCE=150
BOOT=199
SEED=20261010
RIDGE=5.
SCALES_KM=(50.,250.)
COHORTS=(250,500)


def BM_shared_tree_cov(tree,tipnames:list[str])->np.ndarray:
    """Common-ancestral-edge Brownian covariance, not a genus proxy."""
    lookup={str(v):i for i,v in enumerate(tipnames)}
    if len(lookup)!=len(tipnames) or len(lookup)<2:
        raise ValueError("Duplicate or insufficient EXACT source LCVP tree tips")
    tips={str(t.name) for t in tree.get_terminals()}
    if tips!=set(lookup):
        raise ValueError("Original LCVP tree tips are not exactly the old source species")
    K=np.zeros((len(tipnames),len(tipnames)),float)
    def visit(node):
        if node.is_terminal():
            children=[lookup[str(node.name)]]
        else:
            children=[q for child in node.clades for q in visit(child)]
        length=0. if node.branch_length is None else float(node.branch_length)
        if not np.isfinite(length) or length<0:
            raise ValueError("Undated or negative LCVP tree edge")
        if length>0:
            K[np.ix_(children,children)]+=length
        return children
    visit(tree.root)
    denom=np.mean(np.diag(K))
    if not np.isfinite(denom) or denom<=0:
        raise ValueError("Dated phylogenetic tree did not yield positive path length")
    return K/denom


def geographic_kernel(d:pd.DataFrame,scale_km:float)->np.ndarray:
    latitude=np.deg2rad(d.latitude.to_numpy(float))
    longitude=np.deg2rad(d.longitude.to_numpy(float))
    if not np.isfinite(latitude).all() or not np.isfinite(longitude).all():
        raise ValueError("An original source photo has no real public geo position")
    sphere=np.stack((np.cos(latitude)*np.cos(longitude),
                     np.cos(latitude)*np.sin(longitude),np.sin(latitude)),axis=1)
    dot=np.clip(sphere@sphere.T,-1,1)
    distance=np.arccos(dot)*6371.0088
    return np.exp(-distance/scale_km)


def add_source_geo(d:pd.DataFrame)->pd.DataFrame:
    d=d.copy()
    lat=pd.to_numeric(d.latitude,errors="coerce").to_numpy(float)
    lon=pd.to_numeric(d.longitude,errors="coerce").to_numpy(float)
    d["abs_lat"]=np.abs(lat)
    d["longitude_sin"]=np.sin(np.deg2rad(lon))
    d["longitude_cos"]=np.cos(np.deg2rad(lon))
    d["true_photo_cell"]=np.clip(np.floor((np.sin(np.deg2rad(lat))+1)*4.5).astype(int),0,8)*18+np.clip(np.floor((lon+180)/20).astype(int),0,17)
    return d


def source_features(d:pd.DataFrame,soil:bool,blocks:dict|None=None)->tuple[pd.DataFrame,dict]:
    feats=tuple(f for k,v in (BLOCKS if blocks is None else blocks).items() if soil or k!="soil" for f in v)
    allcols=GEO+feats
    bad=d[list(allcols)].isna().any(axis=1)
    valid=d.loc[~bad].copy().reset_index(drop=True)
    receipt={
        "n_direct_original_photo_taxa_pre_missingness":len(d),
        "n_direct_original_photo_taxa_full_environment_common_case":len(valid),
        "n_direct_original_photo_taxa_excluded_by_environment_mask":int(bad.sum()),
        "n_direct_genus_represented":int(valid.genus.nunique()),
        "n_direct_original_photo_geo_cells":int(valid.true_photo_cell.nunique()),
        "n_original_taxa_each_tested_once":int(valid.inat_taxon_id.nunique()),
    }
    return valid,receipt


def environmental_kernel(X_train:np.ndarray,X_test:np.ndarray)->np.ndarray:
    mean=np.mean(X_train,axis=0)
    sd=np.std(X_train,axis=0)
    sd[sd<1e-9]=1.
    a=(X_train-mean)/sd
    b=(X_test-mean)/sd
    return (b@a.T)/len(mean)


def fixed_krr_scores(d:pd.DataFrame,Ksp:np.ndarray,Kphy:np.ndarray,
                     heldout:str,*,soil:bool,nboot:int=BOOT,blocks:dict|None=None)->dict:
    from sklearn.model_selection import GroupKFold
    if heldout not in ("genus","true_photo_cell"):
        raise ValueError("Cross-validation group is original species taxon genus or true geographic photo cell")
    groups=d[heldout].to_numpy()
    y=pd.Categorical(d.morph,categories=CLASSES).codes
    if len(d)<MIN_SOURCE or len(set(groups))<N_FOLDS or len(set(y))<4:
        return {"status":"HOLD_DIRECT_LCVP_MODEL_COVERAGE","n_source_complete_species":len(d)}
    candidates={k:v for k,v in (BLOCKS if blocks is None else blocks).items() if soil or k!="soil"}
    allfeatures=tuple(f for group in candidates.values() for f in group)
    methods={"GEO_SPATIAL_ONLY":(),"GEO_SPATIAL_PHYLOGENY":()}
    methods["GEO_SPATIAL_PHYLOGENY_ALL_ENV"]=allfeatures
    for k,cols in candidates.items():
        methods["GEO_SPATIAL_PHYLOGENY_MINUS_"+k.upper()]=tuple(
            f for f in allfeatures if f not in cols)
    predictions={m:np.full((len(d),4),np.nan) for m in methods}
    yy=np.eye(4)[y]
    plan=list(GroupKFold(n_splits=N_FOLDS).split(d,y,groups))
    geo=d[list(GEO)].to_numpy(float)
    for train,test in plan:
        center=yy[train].mean(axis=0)
        Kspace=Ksp[np.ix_(train,train)]
        Kspacecross=Ksp[np.ix_(test,train)]
        Ktree=Kphy[np.ix_(train,train)]
        Ktreecross=Kphy[np.ix_(test,train)]
        kg=environmental_kernel(geo[train],geo[train])
        kt=environmental_kernel(geo[train],geo[test])
        Kbase=Kspace+kg
        Kbase_test=Kspacecross+kt
        for name,cols in methods.items():
            havephy=name!="GEO_SPATIAL_ONLY"
            Ktrain=Kbase+(Ktree if havephy else 0)
            Kcross=Kbase_test+(Ktreecross if havephy else 0)
            if cols:
                X=d[list(cols)].to_numpy(float)
                Ktrain=Ktrain+environmental_kernel(X[train],X[train])
                Kcross=Kcross+environmental_kernel(X[train],X[test])
            fit=np.linalg.solve(Ktrain+RIDGE*np.eye(len(train)),yy[train]-center)
            p=np.maximum(center+Kcross@fit,0)
            total=p.sum(axis=1,keepdims=True)
            predictions[name][test]=np.divide(p,total,out=np.full_like(p,.25),where=total>0)
    if any(not np.isfinite(p).all() for p in predictions.values()):
        raise RuntimeError("Not every source direct-tip photo species got same-fold conditional predictions")
    loss={k:np.square(p-yy).sum(axis=1) for k,p in predictions.items()}
    metrics={k:{"mean_heldout_four_colour_brier":float(v.mean()),"n_same_species":len(d)}
             for k,v in loss.items()}
    comparisons={
        "spatial_plus_real_phylogeny_beyond_space_only":
            ("GEO_SPATIAL_ONLY","GEO_SPATIAL_PHYLOGENY"),
        "all_abiotic_beyond_spatial_plus_real_phylogeny":
            ("GEO_SPATIAL_PHYLOGENY","GEO_SPATIAL_PHYLOGENY_ALL_ENV"),
    }
    for k in candidates:
        comparisons["unique_"+k+"_beyond_space_real_phylogeny_and_other_abiotic"]=(
            "GEO_SPATIAL_PHYLOGENY_MINUS_"+k.upper(),"GEO_SPATIAL_PHYLOGENY_ALL_ENV")
    root,gi=np.unique(groups,return_inverse=True)
    nn=np.bincount(gi)
    result={}
    for name,(small,large) in comparisons.items():
        difference=loss[small]-loss[large]
        sums=np.bincount(gi,weights=difference)
        rng=np.random.default_rng(SEED+len(name))
        draw=rng.integers(0,len(root),size=(nboot,len(root)))
        bootstrap=sums[draw].sum(axis=1)/nn[draw].sum(axis=1)
        result[name]={
            "heldout_brier_reduction":float(difference.mean()),
            "original_group_cluster_95CI":[float(v) for v in np.quantile(bootstrap,[.025,.975])],
            "conditional_interval_entirely_positive":bool(np.quantile(bootstrap,.025)>0),
        }
    return {
        "status":"EXPLORATORY_EXACT_PHYLO_BROWNIAN_KERNEL_PLUS_SPATIAL_KERNEL_ENVIRONMENT",
        "n_direct_original_source_species":len(d),
        "n_heldout_groups":len(root),
        "same_original_species_and_fold_assignments":True,
        "modelled_true_phylogeny_BM_shared_branch_covariance":True,
        "modelled_original_photograph_spatial_exponential_covariance":True,
        "source_environment_blocks":list(candidates),
        "source_model_scores":metrics,
        "conditional_feature_gains":result,
        "not_a_phylogenetic_logistic_MCMC_or_genetic_assay":True,
    }


def run(source:pd.DataFrame,ledger:pd.DataFrame,trees:dict[int,Path],*,strict:bool=True,nboot:int=BOOT,blocks:dict|None=None)->dict:
    population,coverage=photo_populations(source,strict=strict)
    out={}
    for cap in COHORTS:
        subset=fixed_source_cohort(population["CLIMATE_ALL"],cap)
        direct,tree=direct_taxa(subset,ledger,trees[cap],cap,strict=strict)
        direct=add_source_geo(direct)
        for mode in ("climate","soil"):
            sample,receipt=source_features(direct,soil=mode=="soil",blocks=blocks)
            names=sample.original_direct_LCVP_tip.tolist()
            Kphy=BM_shared_tree_cov(tree,names) if len(sample)==len(direct) else BM_shared_tree_cov(
                Phylo.read(str(trees[cap]),"newick"),direct.original_direct_LCVP_tip.tolist())[np.ix_(
                    direct.inat_taxon_id.isin(sample.inat_taxon_id).to_numpy(),
                    direct.inat_taxon_id.isin(sample.inat_taxon_id).to_numpy())]
            tests={}
            for spatial_scale in SCALES_KM:
                Ks=geographic_kernel(sample,spatial_scale)
                for heldout in ("genus","true_photo_cell"):
                    tests[f"{int(spatial_scale)}km__{heldout}"]=fixed_krr_scores(
                        sample,Ks,Kphy,heldout,soil=mode=="soil",nboot=nboot,blocks=blocks)
            out[f"{cap}km_{mode}"]={"source":receipt,"spatial_kernel_scale_km_all_reported":list(SCALES_KM),"models":tests}
    return {
        "schema":SCHEMA,"date_jst":"2026-10-10",
        "status":"LIMITED_DIRECT_TIP_SOURCE_SPATIAL_AND_PHYLOGENETIC_COVARIANCE_PREDICTION",
        "original_42111_photo_species_denominator":len(source),
        "original_one_photo_classifiable":coverage["source_classified_photo_taxa"],
        "full_872_1761_direct_tip_phylogenetic_coverage_HOLD":True,
        "local_photo_source_250_500_nested_cohorts_not_independent":True,
        "source_model":"kernel ridge on four original photo-colour categories; Brownian shared branch, exponential geodesic spatial residual",
        "spatial_scales_km":[50,250],"fixed_ridge":RIDGE,
        "two_heldout_groups":["genus","true_photo_cell"],
        "no_phenotype_reclassification_or_prospective_source_opened":True,
        "cohorts":out,
        "hard_nonclaims":[
            "Brownian dated shared-branch kernel is one specific phylogenetic covariance hypothesis, not measured genotype effect",
            "Multiple overlapping 250/500km, 50/250km and genus/cell model comparisons are exploratory, not independent confirmatory tests",
            "Kernel ridge four-colour Brier improvement does not demonstrate environment-driven pigment adaptation",
            "The true direct LCVP matched 342/649 source species differ environmentally from untipped original 872/1761 species",
            "Spatially matched original photo-colour label exchange nulls previously did not support robust localized rain inference",
        ],
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--original-breadth-abiotic",required=True,type=Path)
    ap.add_argument("--original-LCVP-tip-ledger",required=True,type=Path)
    ap.add_argument("--tree-250",required=True,type=Path)
    ap.add_argument("--tree-500",required=True,type=Path)
    ap.add_argument("--outdir",required=True,type=Path)
    a=ap.parse_args()
    result=run(pd.read_csv(a.original_breadth_abiotic,low_memory=False),
               pd.read_csv(a.original_LCVP_tip_ledger,low_memory=False),
               {250:a.tree_250,500:a.tree_500})
    a.outdir.mkdir(parents=True,exist_ok=True)
    (a.outdir/"result.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"schema":SCHEMA,"cohorts":result["cohorts"]},sort_keys=True),flush=True)

if __name__=="__main__":
    main()

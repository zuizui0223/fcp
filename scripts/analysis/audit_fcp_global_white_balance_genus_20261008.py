#!/usr/bin/env python3
"""Genus-equal sensitivity to apparent ~50:50 white/chromatic FCP geography.

Source: previously measured 42,111-species breadth taxonomic identity and
85,337 taxon-cell photographed colour records, SHA256 pinned. A genus is an
equal-weight CLUSTER here; it is NOT a phylogeny or an independent causal test.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

import audit_fcp_global_white_balance_20261008 as white
import synthesize_fcp_all42111_measured_geography_20261008 as original

SOURCE_COMMIT=original.FROZEN_SOURCE_COMMIT
BREADTH_SHA=original.FILES_SHA["breadth"]
CELL_SHA=original.FILES_SHA["taxon_cell"]
N_BOOT=1999
SEED=20261008193
PAIR_THRESHOLDS=(1,3)


def sha(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):h.update(b)
    return h.hexdigest()


def source_genus(breadth:pd.DataFrame)->pd.DataFrame:
    expected={"inat_taxon_id","species"}
    if not expected.issubset(breadth):
        raise ValueError("Need original exact taxon species and ID")
    if len(breadth)!=original.N_WORLD_SPECIES or breadth.inat_taxon_id.duplicated().any():
        raise RuntimeError("Historical global species taxon registry is not one per species")
    result=breadth[["inat_taxon_id","species"]].copy()
    result["genus"]=result.species.astype(str).str.strip().str.split().str[0]
    if result.genus.isna().any() or result.genus.eq("").any():
        raise RuntimeError("Empty historical genus label")
    return result[["inat_taxon_id","genus"]]


def cluster_ci(x:np.ndarray,seed:int)->list[float]|None:
    x=np.asarray(x,float)
    if len(x)<2:return None
    r=np.random.default_rng(seed)
    b=np.empty(N_BOOT,float)
    for i in range(N_BOOT):
        b[i]=x[r.integers(0,len(x),size=len(x))].mean()
    return [float(v) for v in np.quantile(b,[.025,.975])]


def within_region_genus_stat(species_region:pd.DataFrame,
                             label:str,seed:int)->tuple[dict,pd.DataFrame]:
    region=species_region.loc[(species_region.region==label)&(
        species_region.n_classified>0)].copy()
    if region.inat_taxon_id.duplicated().any():
        raise RuntimeError("Duplicate source species in regional genus analysis")
    if not len(region):
        return {"region":label,"status":"HOLD_NO_CLASSIFIABLE_SPECIES"},pd.DataFrame()
    genus_means=region.groupby("genus",sort=True).agg(
        n_species=("inat_taxon_id","size"),
        mean_species_white_fraction=("fraction_white_when_classified","mean"),
        n_classified_taxon_cells=("n_classified","sum")
    ).reset_index()
    x=genus_means.mean_species_white_fraction.to_numpy(float)
    sizes=genus_means.n_species.to_numpy(int)
    return {
        "region":label,
        "status":"DESCRIPTIVE_GENUS_EQUAL_WHITE_FRACTION",
        "n_species_classifiable":int(len(region)),
        "n_genera_classifiable":int(len(genus_means)),
        "mean_species_equal_photo_white_fraction":float(region.fraction_white_when_classified.mean()),
        "mean_genus_equal_photo_white_fraction":float(x.mean()),
        "genus_cluster_bootstrap_95CI":cluster_ci(x,seed),
        "absolute_difference_genus_equal_minus_species_equal":float(
            x.mean()-region.fraction_white_when_classified.mean()),
        "max_species_per_genus":int(sizes.max()),
        "largest_genus_species_share":float(sizes.max()/len(region)),
        "n_genera_with_more_than_one_photo_species":int((sizes>1).sum()),
        "n_white_fraction_exactly_half_genus_means":int((x==.5).sum()),
        "biological_limit":"equal genus weight mitigates variable sampled genus richness, but not phylogenetic covariance, source missingness, image clipping or inherited morph frequencies"
    },genus_means


def genus_paired_region_delta(s:pd.DataFrame,label_a:str,label_b:str,
                              min_support:int)->dict:
    original_estimate,paired=white.species_pair_region(s,label_a,label_b,min_support)
    result={
        "region_a":label_a,"region_b":label_b,
        "minimum_classifiable_taxon_cells_each_region":min_support,
        "n_matched_species":len(paired),
        "status":original_estimate["status"]
    }
    if paired.empty:
        result["n_matched_genera"]=0
        result["genus_equal_paired_photo_white_delta"]=None
        result["bootstrap_CI_95_genus_cluster"]=None
        return result
    if paired.inat_taxon_id.duplicated().any():
        raise RuntimeError("A paired species duplicated by source join")
    tax= s[["inat_taxon_id","genus"]].drop_duplicates("inat_taxon_id")
    joined=paired.merge(tax,on="inat_taxon_id",how="left",validate="one_to_one")
    if joined.genus.isna().any():raise RuntimeError("Matched species with no genus identity")
    groups=joined.groupby("genus",sort=True).delta_white_region_b_minus_a.mean()
    x=groups.to_numpy(float)
    seed=SEED+list(original.REGIONS).index(label_a)*19+list(original.REGIONS).index(label_b)*5+min_support*100
    result.update({
        "n_matched_genera":int(len(groups)),
        "species_equal_paired_photo_white_delta":float(paired.delta_white_region_b_minus_a.mean()),
        "genus_equal_paired_photo_white_delta":float(x.mean()),
        "bootstrap_CI_95_genus_cluster":cluster_ci(x,seed),
        "fraction_genera_with_positive_paired_white_change":float((x>0).mean()),
        "fraction_genera_with_negative_paired_white_change":float((x<0).mean()),
        "no_assumption_of_independent_congeneric_species":True
    })
    return result


def analyze(breadth:pd.DataFrame,cell:pd.DataFrame)->tuple[dict,pd.DataFrame]:
    if len(cell)!=original.N_WORLD_CELL_ROWS:
        raise RuntimeError("Source geographic classifiable opportunity changed")
    taxonomy=source_genus(breadth)
    sr=white.source_species_regions(cell)
    joined=sr.merge(taxonomy,on="inat_taxon_id",how="left",validate="many_to_one")
    if len(joined)!=len(sr) or joined.genus.isna().any():
        raise RuntimeError("Global species-genus join failed")
    byregion=[];bygenus=[]
    for i,region in enumerate(original.REGIONS):
        stats,table=within_region_genus_stat(joined,region,SEED+i)
        byregion.append(stats)
        if len(table):
            table["region"]=region
            bygenus.append(table)
    contrasts=[]
    regs=list(original.REGIONS)
    for i,first in enumerate(regs):
        for other in regs[i+1:]:
            for k in PAIR_THRESHOLDS:
                contrasts.append(genus_paired_region_delta(joined,first,other,k))
    return {
        "schema":"fcp_global42111_genus_equal_white_geography_sensitivity_v1",
        "date_jst":"2026-10-08",
        "status":"RETROSPECTIVE_SOURCE_VERIFIED_GENUS_BALANCE_NOT_PHYLOGENETIC_MODEL",
        "original_source_commit":SOURCE_COMMIT,
        "source_breadth_SHA256":BREADTH_SHA,
        "source_taxon_cell_SHA256":CELL_SHA,
        "historical_species_n":original.N_WORLD_SPECIES,
        "historical_taxon_cell_n":original.N_WORLD_CELL_ROWS,
        "genus_equal_regional_white_balance":byregion,
        "genus_equal_paired_within_species_latitude_differences":contrasts,
        "source_terminology":"species-grouped by the original genus token; not a resolved phylogenetic tree",
        "limitations":[
            "photo white may be overexposure, not white anthocyanin-deficient alleles",
            "any-photo and classifiability denominator not random globally",
            "genus labels are not independent phylogenetic contrasts",
            "species regional sampling and within-genus sampling differ globally",
            "no direct genetic FCP, morph frequency, adaptation, net selection or opposing fitness costs inferred"
        ],
        "confirmatory_decisions_changed":False
    },pd.concat(bygenus,ignore_index=True)


def main()->None:
    p=argparse.ArgumentParser()
    p.add_argument("--breadth",required=True,type=Path)
    p.add_argument("--taxon-cell",required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    arg=p.parse_args()
    for path,expected in ((arg.breadth,BREADTH_SHA),(arg.taxon_cell,CELL_SHA)):
        if sha(path)!=expected:raise RuntimeError("Original photo colour class/identity source SHA drift")
    breadth=pd.read_csv(arg.breadth,low_memory=False)
    cell=pd.read_csv(arg.taxon_cell,low_memory=False)
    result,genus=analyze(breadth,cell)
    arg.outdir.mkdir(parents=True,exist_ok=True)
    path=arg.outdir/"genus_equal_region_white_photo_means.csv"
    genus.to_csv(path,index=False,lineterminator="\n")
    result["genus_region_row_ledger_SHA256"]=sha(path)
    (arg.outdir/"result.json").write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False))


if __name__=="__main__":
    main()

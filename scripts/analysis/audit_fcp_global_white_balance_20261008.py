#!/usr/bin/env python3
"""Is the roughly 50:50 white/nonwhite global FCP photo proportion biological?

Audit three separate estimands on ORIGINAL 85,337 taxon-cell photos:
- all taxon-cell observations (cell-weighted);
- each species equally weighted within a region (species turnover remains);
- species PAIRED across two latitude regions (species identity held fixed).

No assumption that a white *photo* is a white genotype or that global
white-flowered species are maintained in a stable 50:50 morph equilibrium.
All unclassified taxon cells remain in missingness bounds.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

import synthesize_fcp_all42111_measured_geography_20261008 as original

N_SOURCE=original.N_WORLD_CELL_ROWS
SOURCE_SHA=original.FILES_SHA["taxon_cell"]
REGION_IDS=tuple(original.REGIONS)
BOOTSTRAPS=1999
SEED=20261008153
PER_SPECIES_MIN_CELLS_SENSITIVITY=3
MIN_PAIR_SPECIES_WITH_BOTH_REGIONS=25


def file_sha(p:Path)->str:
    d=hashlib.sha256()
    with p.open("rb") as f:
        for v in iter(lambda:f.read(1<<20),b""):d.update(v)
    return d.hexdigest()


def source_species_regions(frame:pd.DataFrame)->pd.DataFrame:
    required={"inat_taxon_id","cell_id","morph","measurement_status"}
    if not required.issubset(frame):
        raise RuntimeError(f"Missing photograph species/cell/classification keys: {required-set(frame)}")
    if frame[["inat_taxon_id","cell_id"]].duplicated().any():
        raise RuntimeError("Repeated exact species/cell outcome must never count twice")
    df=frame[list(required)].copy()
    df["region"]=original.classify_region(
        original.row_cell_latitude(df.cell_id.to_numpy(int)))
    df["classifiable"]=(
        df.morph.astype(str).isin(original.MORPHS)&
        df.measurement_status.astype(str).eq("classified_four_state_morph"))
    df["white"]=df.classifiable & df.morph.astype(str).eq("white")
    s=df.groupby(["region","inat_taxon_id"],observed=True,sort=True).agg(
        n_cells=("cell_id","size"),
        n_classified=("classifiable","sum"),
        n_photo_white=("white","sum"),
    ).reset_index()
    if s.n_cells.sum()!=len(df):
        raise RuntimeError("Taxon-cell source rows lost in species/region aggregation")
    if (s.n_photo_white>s.n_classified).any():
        raise RuntimeError("White photos more numerous than photo classifications")
    s["fraction_white_when_classified"]=np.where(
        s.n_classified>0,
        s.n_photo_white/s.n_classified.replace(0,np.nan),np.nan)
    s["white_and_chromatic_images_in_same_region"]=(
        (s.n_photo_white>0)&(s.n_photo_white<s.n_classified))
    return s


def boot_mean(values:np.ndarray,seed:int)->list[float]|None:
    x=np.asarray(values,dtype=float)
    x=x[np.isfinite(x)]
    if len(x)<2:return None
    rng=np.random.default_rng(seed)
    out=np.empty(BOOTSTRAPS,float)
    for i in range(BOOTSTRAPS):
        out[i]=float(x[rng.integers(0,len(x),size=len(x))].mean())
    return [float(z) for z in np.quantile(out,[.025,.975])]


def region_coverage(s:pd.DataFrame)->list[dict]:
    out=[]
    for i,reg in enumerate(REGION_IDS):
        part=s.loc[s.region==reg]
        viewed=part.loc[part.n_classified>0]
        n=int(part.n_cells.sum())
        nclass=int(viewed.n_classified.sum())
        nwhite=int(viewed.n_photo_white.sum())
        unknown=n-nclass
        white_rate=nwhite/nclass if nclass else None
        species_equal=float(viewed.fraction_white_when_classified.mean()) if len(viewed) else None
        strict=viewed.loc[viewed.n_classified>=PER_SPECIES_MIN_CELLS_SENSITIVITY]
        out.append({
            "region":reg,
            "n_taxa_with_any_cell":len(part),
            "n_taxa_with_classifiable_cells":len(viewed),
            "n_taxa_with_ge3_classifiable_cells":len(strict),
            "n_taxon_cell_rows":n,
            "n_classifiable_cell":nclass,
            "n_unclassifiable_cell":unknown,
            "white_classified_cells":nwhite,
            "coloured_classified_cells":nclass-nwhite,
            "fraction_white_all_classified_cell_records":white_rate,
            "difference_from_0p5_all_classified_cell_records":(
                white_rate-.5 if white_rate is not None else None),
            "NO_assumption_all_cell_white_frequency_bound":[float(nwhite/n),float((nwhite+unknown)/n)] if n else None,
            "species_equal_mean_white_fraction_within_region":species_equal,
            "species_bootstrap_CI_95_species_equal_white_fraction":boot_mean(
                viewed.fraction_white_when_classified.to_numpy(float),SEED+i),
            "n_observed_species_with_white_plus_coloured_cells_within_region":int(
                viewed.white_and_chromatic_images_in_same_region.sum()),
            "n_strict_ge3_classifiable_cell_species_with_both_colors":int(
                strict.white_and_chromatic_images_in_same_region.sum()),
            "source_is_observer_photo_morph_not_genetic_fitness":True,
        })
    return out


def species_pair_region(s:pd.DataFrame,region_a:str,region_b:str,
                        min_classified_each:int)->tuple[dict,pd.DataFrame]:
    if region_a==region_b:
        raise ValueError("Independent regions required")
    r1=s.loc[(s.region==region_a)&(s.n_classified>=min_classified_each),
             ["inat_taxon_id","n_classified","fraction_white_when_classified","white_and_chromatic_images_in_same_region"]]
    r2=s.loc[(s.region==region_b)&(s.n_classified>=min_classified_each),
             ["inat_taxon_id","n_classified","fraction_white_when_classified","white_and_chromatic_images_in_same_region"]]
    paired=r1.merge(r2,on="inat_taxon_id",how="inner",validate="one_to_one",suffixes=("_a","_b"))
    paired["delta_white_region_b_minus_a"]=paired.fraction_white_when_classified_b-paired.fraction_white_when_classified_a
    n=len(paired)
    if not n:
        return {
            "n_matched_species":0,"status":"HOLD_ZERO_MATCHED_SPECIES",
            "region_a":region_a,"region_b":region_b,
            "min_classifiable_cell_images_each_region":min_classified_each,
        },paired
    x=paired.delta_white_region_b_minus_a.to_numpy(float)
    rngseed=SEED+int(list(REGION_IDS).index(region_a)*10+list(REGION_IDS).index(region_b))+min_classified_each*500
    # Do NOT count paired taxon-cell anchors as independent biological species.
    return {
        "region_a":region_a,"region_b":region_b,
        "min_classifiable_cell_images_each_region":min_classified_each,
        "n_matched_species":n,
        "n_observed_mixed_white_and_coloured_in_either_region":int(
            (paired.white_and_chromatic_images_in_same_region_a|
             paired.white_and_chromatic_images_in_same_region_b).sum()),
        "n_observed_mixed_white_and_coloured_in_both_regions":int(
            (paired.white_and_chromatic_images_in_same_region_a&
             paired.white_and_chromatic_images_in_same_region_b).sum()),
        "species_equal_white_difference_region_b_minus_region_a":float(x.mean()),
        "species_paired_bootstrap_95CI":boot_mean(x,rngseed),
        "species_fraction_positive_white_change":float(np.mean(x>0)),
        "species_fraction_negative_white_change":float(np.mean(x<0)),
        "species_fraction_equal_photo_white_fraction":float(np.mean(x==0)),
        "n_species_positive_photo_white_change":int(np.count_nonzero(x>0)),
        "n_species_negative_photo_white_change":int(np.count_nonzero(x<0)),
        "status":(
            "EXPLORATORY_SAME_SPECIES_CROSS_REGION_AVAILABLE"
            if n>=MIN_PAIR_SPECIES_WITH_BOTH_REGIONS
            else "HOLD_SPARSE_SAME_SPECIES_REGIONAL_PHOTOS"),
        "claim_limits":"paired species-controlled photographic contrast, not within-individual plasticity or adaptation; measurement origin, observer, source photo dates and clipped whites still confound",
    },paired


def analyze(frame:pd.DataFrame)->tuple[dict,pd.DataFrame,pd.DataFrame]:
    if len(frame)!=N_SOURCE:
        raise RuntimeError(f"Original region-wide measured source expected {N_SOURCE}, got {len(frame)}")
    s=source_species_regions(frame)
    regions=region_coverage(s)
    paired_rows=[]
    allpairs=[]
    for i,a in enumerate(REGION_IDS):
        for b in REGION_IDS[i+1:]:
            for threshold in (1,PER_SPECIES_MIN_CELLS_SENSITIVITY):
                payload,matches=species_pair_region(s,a,b,threshold)
                paired_rows.append(payload)
                if len(matches):
                    matches["region_a"]=a
                    matches["region_b"]=b
                    matches["threshold"]=threshold
                    allpairs.append(matches)
    result={
        "schema":"fcp_global42111_photo_white_balance_versus_species_turnover_audit_v1",
        "date_jst":"2026-10-08",
        "status":"SOURCE_VERIFIED_RETROSPECTIVE_WHITE_BALANCE_DESCRIPTION",
        "origin_global_measured_commit":original.FROZEN_SOURCE_COMMIT,
        "source_taxon_cell_SHA256":SOURCE_SHA,
        "original_species_universe":original.N_WORLD_SPECIES,
        "observed_taxon_cell_rows":N_SOURCE,
        "region_white_balance":regions,
        "same_species_paired_region_contrasts":paired_rows,
        "same_species_region_photo_diversity":"number of image colours observed in region, not genetic morph coexistence in one mating population",
        "what_half_white_means":"approx 50% is conditional on source photographs whose ROI was classifiable, among regionally sampled taxon-cell anchors; NOT 50:50 within a species, or selected white/pigmented genotypes in a population",
        "diagnostic_distinctions":[
            "taxon-cell-weighted share differs from species-equal within-region share when widespread species repeat across cells",
            "species-equal within-region still compares different species pools; paired crossregion contrasts constrain species identity",
            "rare-white photo non-detection does not prove regional absence of white morphs",
            "unclassified photo outcomes can be white or coloured; absolute no-assumption bounds are wide",
            "observer/light exposure coupling can falsely raise source photo-white classification",
            "different geographic cell photos may sample separate populations and even nonhomologous floral organs",
            "a photo-derived 50:50 ratio is not an ecological balancing-selection equilibrium"
        ],
        "confirmatory_decisions_changed":False,
    }
    detail=pd.concat(allpairs,ignore_index=True) if allpairs else pd.DataFrame()
    return result,s,detail


def main()->None:
    p=argparse.ArgumentParser()
    p.add_argument("--taxon-cell",required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    args=p.parse_args()
    if file_sha(args.taxon_cell)!=SOURCE_SHA:
        raise RuntimeError("Original entire taxon-cell source byte hash mismatch")
    d=pd.read_csv(args.taxon_cell,low_memory=False)
    outcome,s,matched=analyze(d)
    args.outdir.mkdir(parents=True,exist_ok=True)
    sm=args.outdir/"species_by_abs_latitude_region_photo_colour_summary.csv.gz"
    s.to_csv(sm,index=False,lineterminator="\n",
             compression={"method":"gzip","compresslevel":9,"mtime":0})
    pm=args.outdir/"same_species_paired_absolute_latitude_photo_colour.csv.gz"
    matched.to_csv(pm,index=False,lineterminator="\n",
                   compression={"method":"gzip","compresslevel":9,"mtime":0})
    outcome["species_region_ledger_SHA256"]=file_sha(sm)
    outcome["crossregion_matched_species_ledger_SHA256"]=file_sha(pm)
    (args.outdir/"result.json").write_text(
        json.dumps(outcome,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print(json.dumps(outcome,ensure_ascii=False,indent=2,allow_nan=False))


if __name__=="__main__":
    main()

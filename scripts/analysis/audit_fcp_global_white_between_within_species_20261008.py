#!/usr/bin/env python3
"""Decompose apparent global ~50% photographed white into BETWEEN- and WITHIN-species variation.

Original September 2026 FCP one-image-per-taxon×cell 85,337-row table.
Only source-classified four coarse photo-colour labels enter variance; all
unclassifiable photos remain in opportunity/coverage denominator.

A balanced global population of photos need NOT be balanced genetic morph
frequencies within any biological mating population.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

import synthesize_fcp_all42111_measured_geography_20261008 as source

INPUT_SHA=source.FILES_SHA["taxon_cell"]
N_INPUT=source.N_WORLD_CELL_ROWS
THRESHOLDS=(2,3,5,10)
BOOTSTRAPS=999
SEED=20261008143
MORPHS=set(source.MORPHS)


def sha(path:Path)->str:
    d=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1<<20),b""):d.update(block)
    return d.hexdigest()


def compress_species_photo_records(frame:pd.DataFrame)->pd.DataFrame:
    required={"inat_taxon_id","cell_id","measurement_status","morph"}
    if not required.issubset(frame):
        raise RuntimeError(f"original measured source keys absent: {required-set(frame)}")
    if frame[["inat_taxon_id","cell_id"]].duplicated().any():
        raise RuntimeError("Repeating same photo taxon×cell identity biases apparent FCP polymorphism")
    d=frame[["inat_taxon_id","cell_id","measurement_status","morph"]].copy()
    d["region"]=source.classify_region(source.row_cell_latitude(d.cell_id.to_numpy(int)))
    d["classifiable"]=d.morph.isin(MORPHS) & d.measurement_status.eq("classified_four_state_morph")
    d["white"]=d.classifiable & d.morph.eq("white")
    d["global_region"]="ALL_GLOBAL"
    parts=[]
    for label in ("global_region","region"):
        g=d.groupby([label,"inat_taxon_id"],sort=True,observed=True).agg(
            observed_cells=("cell_id","size"),
            n_classified=("classifiable","sum"),
            n_white=("white","sum"),
        ).reset_index().rename(columns={label:"region"})
        parts.append(g)
    x=pd.concat(parts,ignore_index=True)
    if (x.n_classified>x.observed_cells).any() or (x.n_white>x.n_classified).any():
        raise RuntimeError("Original image-classification opportunity counts inconsistent")
    return x


def decompose(eligible:pd.DataFrame, seed:int)->dict:
    n=len(eligible)
    if n<2:
        return {"n_species":n,"status":"HOLD_TOO_FEW_SPECIES_FOR_VAR_DECOMPOSITION"}
    ni=eligible.n_classified.to_numpy(float)
    w=eligible.n_white.to_numpy(float)
    pi=w/ni

    def metrics(ns:np.ndarray, ps:np.ndarray)->tuple[float,float,float,float,float,float]:
        ps_mean=float(ps.mean())
        equal_within=float(np.mean(ps*(1-ps)))
        equal_between=float(np.mean((ps-ps_mean)**2))
        p_weighted=float(np.dot(ns,ps)/ns.sum())
        within=float(np.dot(ns,ps*(1-ps))/ns.sum())
        between=float(np.dot(ns,(ps-p_weighted)**2)/ns.sum())
        return p_weighted,within,between,ps_mean,equal_within,equal_between

    pw,win,btw,peq,weq,beq=metrics(ni,pi)
    total_photo=pw*(1-pw)
    total_eq=peq*(1-peq)
    if abs((win+btw)-total_photo)>1e-10 or abs((weq+beq)-total_eq)>1e-10:
        raise RuntimeError("Expected exact observed binary photograph variance decomposition does not hold")
    rng=np.random.default_rng(seed)
    boot_photo=np.empty(BOOTSTRAPS,float)
    boot_equal=np.empty(BOOTSTRAPS,float)
    boot_mixed=np.empty(BOOTSTRAPS,float)
    mixed=(w>0)&(w<ni)
    for b in range(BOOTSTRAPS):
        idx=rng.integers(0,n,size=n)
        z,within,between,ze,we,be=metrics(ni[idx],pi[idx])
        boot_photo[b]=between/(z*(1-z)) if 0<z<1 else np.nan
        boot_equal[b]=be/(ze*(1-ze)) if 0<ze<1 else np.nan
        boot_mixed[b]=mixed[idx].mean()
    def CI(x:np.ndarray)->list[float]|None:
        v=x[np.isfinite(x)]
        return [float(t) for t in np.quantile(v,[.025,.975])] if len(v) else None

    n_mix=int(mixed.sum())
    num_balanced=int(np.count_nonzero((pi>=.25)&(pi<=.75)))
    n_fullwhite=int(np.count_nonzero(pi==1))
    n_no_white=int(np.count_nonzero(pi==0))
    return {
        "n_species":n,
        "n_classifiable_photo_taxon_cell_records":int(ni.sum()),
        "status":"OBSERVED_PHOTO_VARIANCE_DECOMPOSED_NOT_POPULATION_MORPH_FREQUENCY",
        "fraction_white_classified_photo_weighted":pw,
        "fraction_white_species_equal":peq,
        "photo_weighted_between_species_variance":btw,
        "photo_weighted_within_species_variance":win,
        "photo_weighted_between_species_share_of_observed_photo_variance":float(btw/total_photo) if total_photo else None,
        "photo_weighted_between_species_share_bootstrap_95CI":CI(boot_photo),
        "species_equal_between_species_variance":beq,
        "species_equal_within_species_variance":weq,
        "species_equal_between_species_share_of_observed_photo_variance":float(beq/total_eq) if total_eq else None,
        "species_equal_between_species_share_bootstrap_95CI":CI(boot_equal),
        "n_species_with_white_only_in_sampled_cells":n_fullwhite,
        "n_species_with_nonwhite_only_in_sampled_cells":n_no_white,
        "n_species_with_both_white_and_nonwhite_in_sampled_cells":n_mix,
        "fraction_species_with_both_observed_states":float(n_mix/n),
        "mixed_species_fraction_bootstrap_95CI":CI(boot_mixed),
        "n_species_with_observed_white_fraction_between_25_and_75_percent":num_balanced,
        "fraction_species_in_25_to_75_percent_range":float(num_balanced/n),
        "n_species_with_0_or_100_percent_sampled_white":int(n-n_mix),
        "estimator_boundary":"Descriptive variance of observed binary photo classifications. Finite images inflate apparent between-species variation; cells are not mating populations; colour categories not alleles."
    }


def analyze(frame:pd.DataFrame)->tuple[dict,pd.DataFrame]:
    if len(frame)!=N_INPUT:
        raise RuntimeError(f"expected 85337 historical source cells, got {len(frame)}")
    species=compress_species_photo_records(frame)
    regions=("ALL_GLOBAL",*source.REGIONS)
    summaries=[]
    for i,reg in enumerate(regions):
        x=species.loc[species.region==reg].copy()
        for t in THRESHOLDS:
            eligible=x.loc[x.n_classified>=t].copy()
            y=decompose(eligible,SEED+i*100+t)
            y.update({
                "region":reg,"minimum_classifiable_distinct_cells_per_species":t,
                "n_species_observed_any_cell":len(x),
                "n_species_with_at_least_one_classifiable_cell":int((x.n_classified>=1).sum()),
                "n_opportunity_taxon_cells":int(x.observed_cells.sum()),
                "n_classifiable_taxon_cells_any_depth":int(x.n_classified.sum()),
                "n_missing_photo_classifications":int((x.observed_cells-x.n_classified).sum())
            })
            summaries.append(y)
    for threshold in THRESHOLDS:
        if threshold==2:
            counts=[x for x in summaries if x["minimum_classifiable_distinct_cells_per_species"]==2 and x["region"]!="ALL_GLOBAL"]
            if sum(z["n_opportunity_taxon_cells"] for z in counts)!=N_INPUT:
                raise RuntimeError("3 latitude band opportunity strata not exhaustive")
    out={
        "schema":"fcp_global42111_white_near_half_between_within_species_photo_decomposition_v1",
        "date_jst":"2026-10-08",
        "status":"RETROSPECTIVE_SOURCE_FROZEN_PHOTO_VARIANCE_AND_DETECTION",
        "input_SHA256":INPUT_SHA,
        "original_species_frame":source.N_WORLD_SPECIES,
        "original_measured_taxon_cell_frame":N_INPUT,
        "regions_plus_world":list(regions),
        "minimum_classifiable_per_species_sensitivity":list(THRESHOLDS),
        "breakdown":summaries,
        "photo_white_state":"classified source photo white versus three other classified source colour categories pooled; NOT white versus pigmented genetic morph",
        "main_question":"Does globally ~half white among classified photos imply that sampled individual species have white and coloured flowers at ~half frequency? Answer requires species-level photo variance and detection sensitivity, not aggregate ratio alone.",
        "inference_limits":[
            "species repeatedly appearing in different cells are not replicate independent populations",
            "all-white or all-coloured after two or three photos is not verified natural species monomorphism",
            "observed within-species white+nonwhite in separated cells is not verified same-population sympatry",
            "finite per-species photo sampling itself biases between-species variance upward",
            "photo colour white/nonwhite depends on exposure, ROI, pigment chemistry and image classification",
            "observer/site/plant phenology and true species genotype are not separated in this atlas",
            "these already observed historical source photos cannot be reused as new prospective 2000+730 taxa colour confirmation",
            "this is not a null test of balancing selection, nor evidence that biological costs/benefits are equal"
        ],
        "confirmatory_decisions_changed":False
    }
    return out,species


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--taxon-cell",type=Path,required=True)
    p.add_argument("--outdir",type=Path,required=True)
    arg=p.parse_args()
    if sha(arg.taxon_cell)!=INPUT_SHA:
        raise RuntimeError("Original taxonomy×cell measured-photo source SHA drift")
    data=pd.read_csv(arg.taxon_cell,low_memory=False)
    j,x=analyze(data)
    arg.outdir.mkdir(parents=True,exist_ok=True)
    path=arg.outdir/"species_region_photo_white_counts.csv.gz"
    x.to_csv(path,index=False,lineterminator="\n",
             compression={"method":"gzip","compresslevel":9,"mtime":0})
    j["species_region_photo_summary_SHA256"]=sha(path)
    (arg.outdir/"result.json").write_text(json.dumps(j,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print(json.dumps(j,indent=2,ensure_ascii=False,allow_nan=False))


if __name__=="__main__":
    main()

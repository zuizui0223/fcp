#!/usr/bin/env python3
"""Label-blind photo-disjoint 50-km pair checks for geographic FCP observations.

Each photo can occur in ONE selected local pair within a species. All source
species, colour labels, exact photo geodesics and species-level colour frequencies
remain frozen. Greedy matching is MAXIMAL, NOT a globally maximum-cardinality match.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
ORIGINAL=ROOT/"scripts/analysis/run_polymorphism_distributed_polymorphism_posthoc_20261007.py"
_spec=importlib.util.spec_from_file_location("fcp_source_original",ORIGINAL)
base=importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(base)

RADIUS_KM=50.
MATCH_SCHEMES=("nearest","fixed_random")
POLICIES=("all_photos","different_observer")
MATCHING_THRESHOLDS=(5,10,15,20,30)
PRIMARY_THRESHOLD=10
N_BOOT=1999
MASTER_SEED=2026101049
EXPECTED_SOURCE={
 "discovery":(369,166,.020529254583812922),
 "validation":(363,181,.018672971642749295),
 "third":(377,204,.01468491968437185),
}

def seed(*terms):
    raw="|".join(map(str,(MASTER_SEED,*terms))).encode()
    return int.from_bytes(hashlib.sha256(raw).digest()[:8],"little")

def outcome_blind_matching(g:pd.DataFrame,cohort:str,scheme:str,policy:str):
    """Return label-agnostic vertex-disjoint photo pairs, never more than n//2."""
    if scheme not in MATCH_SCHEMES or policy not in POLICIES:
        raise ValueError("undefined pre-frozen photographic matching policy")
    n=len(g)
    dist=base.pairwise_geo_km(g.latitude.to_numpy(float),g.longitude.to_numpy(float))
    a,b=np.triu_indices(n,k=1)
    select=dist[a,b]<=RADIUS_KM
    if policy=="different_observer":
        if "observer_id" not in g.columns:
            raise RuntimeError("different observer pairing requires real observer IDs")
        obs=g.observer_id.fillna("").astype(str).to_numpy()
        select&=(obs[a]!="")&(obs[b]!="")&(obs[a]!=obs[b])
    a,b=a[select],b[select]
    if not len(a):
        return np.array([],int),np.array([],int),0
    rng=np.random.default_rng(seed(cohort,int(g.inat_taxon_id.iloc[0]),scheme,policy))
    if scheme=="nearest":
        # Lexicographically deterministic distance and source photo-order tie break.
        order=np.lexsort((b,a,dist[a,b]))
    else:
        order=rng.permutation(len(a))
    used=np.zeros(n,dtype=bool)
    chosen_a=[];chosen_b=[]
    for index in order:
        i=int(a[index]);j=int(b[index])
        if not used[i] and not used[j]:
            chosen_a.append(i);chosen_b.append(j)
            used[i]=True;used[j]=True
    if int(used.sum())!=2*len(chosen_a):
        raise RuntimeError("original photo used in more than one local pair")
    return np.array(chosen_a,dtype=int),np.array(chosen_b,dtype=int),int(len(a))

def species_test(g:pd.DataFrame,cohort:str)->tuple[dict,dict[tuple[str,str],np.ndarray]]:
    g=g.dropna(subset=["latitude","longitude"]).sort_values("photo_id",kind="stable")
    if len(g)<40: raise RuntimeError("Source high-depth geographical classifiability drift")
    original,null=base.analyse_species(g,cohort)
    colour=pd.Categorical(g.morph.astype(str),categories=base.MORPHS).codes.astype(np.int8)
    if (colour<0).any():raise RuntimeError("Source biological state drift")
    rng=np.random.default_rng(base.stable_seed(cohort,int(g.inat_taxon_id.iloc[0])))
    permutations=np.stack([rng.permutation(colour) for _ in range(base.PERMUTATIONS)])
    detail={"cohort":cohort,"inat_taxon_id":int(g.inat_taxon_id.iloc[0]),
            "species":str(g.species.iloc[0]),"n_source_photos":len(g),
            "source_50km_local_pair_edges":int(original["local_pairs_50km"]),
            "source_50km_local_depletion":original["depletion_50km"],
            "D_specieswide":float(original["D_pair"])}
    nulls={}
    for scheme in MATCH_SCHEMES:
        for policy in POLICIES:
            key=f"{scheme}__{policy}"
            a,b,n_edges=outcome_blind_matching(g,cohort,scheme,policy)
            m=len(a)
            detail[f"{key}__eligible_pair_edges"]=n_edges
            detail[f"{key}__disjoint_local_pairs"]=m
            detail[f"{key}__unique_photos"]=2*m
            if m:
                local=float(np.mean(colour[a]!=colour[b]))
                detail[f"{key}__mean_local_pair_discordance"]=local
                detail[f"{key}__depletion"]=float(original["D_pair"]-local)
                nulls[(scheme,policy)]=float(original["D_pair"])-np.mean(
                    permutations[:,a]!=permutations[:,b],axis=1)
            else:
                detail[f"{key}__mean_local_pair_discordance"]=np.nan
                detail[f"{key}__depletion"]=np.nan
                nulls[(scheme,policy)]=np.full(base.PERMUTATIONS,np.nan)
    return detail,nulls

def summarize(rows:pd.DataFrame,null_vectors:list[np.ndarray],scheme:str,policy:str,
              cohort:str,min_pairs:int,restrict_to_original_pair30:bool=False)->dict:
    name=f"{scheme}__{policy}"
    keep=rows[f"{name}__disjoint_local_pairs"].to_numpy(int)>=min_pairs
    if restrict_to_original_pair30:
        # Distinct supplementary estimand: preserve EXACT original >=30-edge
        # species universe before applying the new disjoint-pair threshold.
        keep &= rows["source_50km_local_pair_edges"].to_numpy(int)>=base.MIN_LOCAL_PAIRS
    n=int(keep.sum())
    result={"n_species":n,"minimum_photo_disjoint_pairs":min_pairs,
            "matching_scheme":scheme,"observer_policy":policy,
            "restricted_to_original_30_edge_species":restrict_to_original_pair30,
            "status":"HOLD_INSUFFICIENT_SPECIES" if n<30 else "EXPLORATORY_ESTIMABLE"}
    if n==0:return result
    used=rows.loc[keep]
    dep=used[f"{name}__depletion"].to_numpy(float)
    if not np.isfinite(dep).all():raise RuntimeError("Photo-disjoint source depletion missing")
    mat=np.vstack([x for x,yes in zip(null_vectors,keep) if yes])
    if mat.shape!=(n,base.PERMUTATIONS) or not np.isfinite(mat).all():
        raise RuntimeError("Matched original photographed label permutation matrix invalid")
    obs=float(dep.mean())
    null_means=mat.mean(axis=0)
    result.update({
        "mean_depletion":obs,
        "median_depletion":float(np.median(dep)),
        "mean_specieswide_discordance":float(used.D_specieswide.mean()),
        "mean_disjoint_local_discordance":float(used[f"{name}__mean_local_pair_discordance"].mean()),
        "positive_species_fraction":float(np.mean(dep>0)),
        "mean_unique_photos_in_selected_pairs":float(used[f"{name}__unique_photos"].mean()),
        "median_disjoint_pair_count":float(used[f"{name}__disjoint_local_pairs"].median()),
        "photo_reuse_across_selected_pairs":False,
        "permutation_null_mean":float(null_means.mean()),
        "permutation_p_upper":float((1+np.count_nonzero(null_means>=obs))/(base.PERMUTATIONS+1)),
    })
    if n>=30:
        rng=np.random.default_rng(seed(cohort,scheme,policy,min_pairs,"bootstrap"))
        boots=dep[rng.integers(n,size=(N_BOOT,n))].mean(axis=1)
        result["species_bootstrap_ci95"]=[float(z) for z in np.quantile(boots,[.025,.975])]
    else:
        result["species_bootstrap_ci95"]=None
    result["positive_bounded_evidence"]=bool(
        n>=30 and obs>0 and result["permutation_p_upper"]<.05
        and result["species_bootstrap_ci95"][0]>0)
    return result

def run_cohort(path:Path,cohort:str):
    d=base.load(path,cohort)
    records=[]
    nulls={pair:[] for pair in ((s,p) for s in MATCH_SCHEMES for p in POLICIES)}
    for _,g in d.groupby("inat_taxon_id",sort=True):
        row,nul=species_test(g,cohort)
        records.append(row)
        for pair in nulls:nulls[pair].append(nul[pair])
    df=pd.DataFrame(records)
    n_total,n_orig,dep_orig=EXPECTED_SOURCE[cohort]
    valid=df.source_50km_local_pair_edges>=base.MIN_LOCAL_PAIRS
    if len(df)!=n_total or int(valid.sum())!=n_orig or abs(
        float(df.loc[valid,"source_50km_local_depletion"].mean())-dep_orig)>1e-10:
        raise RuntimeError(f"Frozen original high-depth 50km source result mismatch: {cohort}")
    ladder={};restricted={}
    for scheme in MATCH_SCHEMES:
        for policy in POLICIES:
            name=f"{scheme}__{policy}"
            ladder[name]={
                str(t):summarize(df,nulls[(scheme,policy)],scheme,policy,cohort,t)
                for t in MATCHING_THRESHOLDS
            }
            restricted[name]={
                str(t):summarize(df,nulls[(scheme,policy)],scheme,policy,cohort,t,True)
                for t in MATCHING_THRESHOLDS
            }
    for name in restricted:
        for t in MATCHING_THRESHOLDS:
            if restricted[name][str(t)]["n_species"]>ladder[name][str(t)]["n_species"]:
                raise RuntimeError("restricted original source subset has more species")
    return {"high_depth_species":len(df),
            "original_50km_30edge_n_species":n_orig,
            "original_50km_30edge_mean_depletion":dep_orig,
            "ladder":ladder,
            "historical_original_30edge_species_conditioned_ladder":restricted},df

def main():
    p=argparse.ArgumentParser()
    for c in EXPECTED_SOURCE:p.add_argument("--"+c,type=Path,required=True)
    p.add_argument("--outdir",type=Path,required=True)
    a=p.parse_args()
    res={"schema":"fcp_photo_disjoint_50km_geography_posthoc_20261010_v1",
         "date_jst":"2026-10-10","role":"posthoc_photo_pseudoreplication_support_diagnostic",
         "status":"COMPLETE_SOURCE_FROZEN_POSTHOC",
         "source_sha256":base.SHA,"radius_km":RADIUS_KM,
         "shuffling":"original within-species whole-vertex composition-preserving 199-label null",
         "pair_selection":"outcome-blind greedy maximal vertex-disjoint photo matching",
         "matching_schemes":list(MATCH_SCHEMES),"observer_policies":list(POLICIES),
         "thresholds_disjoint_pair_count":list(MATCHING_THRESHOLDS),
         "primary_min_disjoint_photo_pairs":PRIMARY_THRESHOLD,
         "minimum_species_to_estimate":30,
         "secondary_source_populations":["all_high_depth","historical_original_30edge_species"],
         "permutations":base.PERMUTATIONS,
         "species_bootstraps":N_BOOT,
         "confirmatory_decisions_changed":False,
         "hard_nonclaims":[
            "Photo-disjoint pairs are not guaranteed independent plants, field sites or genotypes.",
            "Greedy matchings are maximal but not guaranteed maximum-cardinality independent pair sets.",
            "Matching nearest versus label-blind random order can select different geographic opportunities.",
            "Nested 5/10/15/20/30 threshold subsets select different species; p values not multiple-testing corrected.",
            "Post-outcome sensitivity cannot establish local adaptation, true genetic polymorphism or fitness.",
         ]}
    allrows=[]
    for cohort in EXPECTED_SOURCE:
        res[cohort],df=run_cohort(getattr(a,cohort),cohort)
        allrows.append(df)
    res["cross_cohort"]={}
    for cohort_population in ("ladder","historical_original_30edge_species_conditioned_ladder"):
        res["cross_cohort"][cohort_population]={
            f"{s}__{p}__min{t}":all(
                res[c][cohort_population][f"{s}__{p}"][str(t)]["positive_bounded_evidence"]
                for c in EXPECTED_SOURCE)
            for s in MATCH_SCHEMES for p in POLICIES for t in MATCHING_THRESHOLDS
        }
    out=a.outdir;out.mkdir(parents=True,exist_ok=True)
    pd.concat(allrows,ignore_index=True).to_csv(out/"photo_disjoint_species_support.csv",index=False)
    (out/"result.json").write_text(json.dumps(res,indent=2)+"\n",encoding="utf8")
    print(json.dumps({c:res[c]["ladder"] for c in EXPECTED_SOURCE},sort_keys=True))
if __name__=="__main__":main()

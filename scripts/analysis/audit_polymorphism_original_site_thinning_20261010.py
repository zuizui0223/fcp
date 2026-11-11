#!/usr/bin/env python3
"""Source-exact geographic hard-core thinning of FCP flower photographs.

Keep at most one chosen original photographed observation inside a radius
around any other selected observation of the SAME species, then replay a
composition-preserving photo-label geographic-null analysis. Geographical
selection is blind to colour outcome, flower metrics and observer identity.
The historical 50-km result and future selected source samples are untouched.
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
spec=importlib.util.spec_from_file_location("original_fcp_distributed",ORIGINAL)
orig=importlib.util.module_from_spec(spec)
spec.loader.exec_module(orig)

SPACINGS_KM=(1.,5.,10.)
N_REALIZATIONS=12
MIN_VALID_REALIZATIONS=8
MIN_RETAINED_PHOTOS=20
MIN_RETAINED_LOCAL_PAIRS=10
N_BOOT=1999
N_PERM=orig.PERMUTATIONS
BASE_SEED=20261010697
RADIUS_KM=50.
OBSERVER_POLICIES=("all_photos","different_observer")
EXPECTED={"discovery":(369,166,.020529254583812922),
          "validation":(363,181,.018672971642749295),
          "third":(377,204,.01468491968437185)}

def seed(*parts):
    payload="|".join(map(str,(BASE_SEED,*parts))).encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:8],"little")

def greedy_min_separation(dist:np.ndarray, spacing:float, order:np.ndarray)->np.ndarray:
    """Retain vertices whose spherical distances from ALL chosen > spacing."""
    dist=np.asarray(dist,float)
    n=dist.shape[0]
    if dist.shape!=(n,n) or not np.isfinite(dist).all():
        raise ValueError("source geographic matrix must be finite square")
    used=[]
    for k in order:
        i=int(k)
        if not used or bool(np.all(dist[i,np.asarray(used,int)]>spacing)):
            used.append(i)
    keep=np.sort(np.asarray(used,int))
    if len(keep)>1:
        a,b=np.triu_indices(len(keep),1)
        if (dist[keep[a],keep[b]]<=spacing).any():
            raise RuntimeError("source photo hard-core minimum separation violated")
    return keep

def study_one_realization(g:pd.DataFrame, distance_matrix:np.ndarray, labels:np.ndarray,
                          cohort:str, taxon:int, spacing:float, rep:int):
    n=len(g)
    order=np.random.default_rng(seed(cohort,taxon,spacing,rep,"photo_selection")).permutation(n)
    take=greedy_min_separation(distance_matrix,spacing,order)
    out={}
    n_retained=len(take)
    if n_retained<MIN_RETAINED_PHOTOS:
        return {p:None for p in OBSERVER_POLICIES},n_retained
    cc=labels[take]
    u,v=np.triu_indices(n_retained,1)
    near=distance_matrix[take[u],take[v]]<=RADIUS_KM
    D=float(np.mean(cc[u]!=cc[v]))
    observers=g.observer_id.fillna("").astype(str).to_numpy()[take]
    permutation_rng=np.random.default_rng(seed(cohort,taxon,spacing,rep,"photo_label_null"))
    nulls=np.stack([permutation_rng.permutation(cc) for _ in range(N_PERM)])
    for policy in OBSERVER_POLICIES:
        local_mask=near.copy()
        if policy=="different_observer":
            local_mask &= (observers[u]!="")&(observers[v]!="")&(observers[u]!=observers[v])
        n_edges=int(local_mask.sum())
        if n_edges<MIN_RETAINED_LOCAL_PAIRS:
            out[policy]=None
            continue
        x,y=u[local_mask],v[local_mask]
        loc=float(np.mean(cc[x]!=cc[y]))
        null_loc=np.mean(nulls[:,x]!=nulls[:,y],axis=1)
        out[policy]={
          "retained_photos":n_retained,
          "local_edges":n_edges,
          "mean_specieswide_discordance":D,
          "mean_local_discordance":loc,
          "depletion":D-loc,
          "null_depletion":D-null_loc,
        }
    return out,n_retained

def analyze_species(g:pd.DataFrame,cohort:str):
    g=g.dropna(subset=["latitude","longitude"]).sort_values("photo_id",kind="stable").reset_index(drop=True)
    if len(g)<40:
        raise RuntimeError("the original high-depth 40-photo classifiability gate drifted")
    row,original_nulls=orig.analyse_species(g,cohort)
    if row["local_pairs_50km"]<orig.MIN_LOCAL_PAIRS:
        return None
    if "observer_id" not in g.columns:
        raise RuntimeError("frozen source photo observer identities missing")
    lat=g.latitude.to_numpy(float)
    lon=g.longitude.to_numpy(float)
    dist=orig.pairwise_geo_km(lat,lon)
    labels=pd.Categorical(g.morph.astype(str),categories=orig.MORPHS).codes.astype(np.int8)
    if (labels<0).any():
        raise RuntimeError("four biological source photo classes changed")
    taxon=int(g.inat_taxon_id.iloc[0])
    record={"cohort":cohort,"inat_taxon_id":taxon,"species":str(g.species.iloc[0]),
            "source_n_classified_photos":len(g),
            "original_allphoto_local_pairs":int(row["local_pairs_50km"]),
            "original_allphoto_depletion":float(row["depletion_50km"])}
    source_result={}
    for spacing in SPACINGS_KM:
        series={p:[] for p in OBSERVER_POLICIES}
        all_n=[]
        for rep in range(N_REALIZATIONS):
            result,nret=study_one_realization(g,dist,labels,cohort,taxon,spacing,rep)
            all_n.append(nret)
            for policy in OBSERVER_POLICIES:
                if result[policy] is not None:
                    series[policy].append(result[policy])
        for policy in OBSERVER_POLICIES:
            key=f"{int(spacing)}km_{policy}"
            vals=series[policy]
            record[f"{key}_n_valid_repeats"]=len(vals)
            record[f"{key}_n_retained_median_all_repeats"]=float(np.median(all_n))
            record[f"{key}_n_photo_pairs_median_valid"]=float(np.median([x["local_edges"] for x in vals])) if vals else np.nan
            if len(vals)<MIN_VALID_REALIZATIONS:
                record[f"{key}_depletion"]=np.nan
                source_result[(int(spacing),policy)]=None
                continue
            dep=float(np.mean([x["depletion"] for x in vals]))
            originalD=float(np.mean([x["mean_specieswide_discordance"] for x in vals]))
            meanLoc=float(np.mean([x["mean_local_discordance"] for x in vals]))
            null_vec=np.mean(np.stack([x["null_depletion"] for x in vals]),axis=0)
            record[f"{key}_depletion"]=dep
            record[f"{key}_mean_specieswide_discordance"]=originalD
            record[f"{key}_mean_local_discordance"]=meanLoc
            source_result[(int(spacing),policy)]=null_vec
    return record,source_result

def summarize_source(records:pd.DataFrame,nulllist:list[dict])->dict:
    cohort=str(records.cohort.iloc[0])
    out={}
    for spacing in SPACINGS_KM:
        for policy in OBSERVER_POLICIES:
            key=f"{int(spacing)}km_{policy}"
            dep=records[f"{key}_depletion"].to_numpy(float)
            mask=np.isfinite(dep)
            n=int(mask.sum())
            item={"n_species":n,"n_original_50km_30edge_species":len(records),
                  "fraction_of_original_species_testable":n/len(records),
                  "radius_minimum_selected_site_separation_km":spacing,
                  "source_observer_policy":policy,
                  "minimum_valid_repeats":MIN_VALID_REALIZATIONS,
                  "status":"HOLD_INSUFFICIENT_GEOGRAPHIC_SUPPORT" if n<30 else "POSTHOC_ESTIMABLE"}
            if n:
                value=dep[mask]
                vecs=np.vstack([r[(int(spacing),policy)] for r,good in zip(nulllist,mask) if good])
                if vecs.shape!=(n,N_PERM) or not np.isfinite(vecs).all():
                    raise RuntimeError("source null arrays invalid")
                avg=vecs.mean(axis=0)
                effect=float(value.mean())
                item.update({
                    "species_equal_mean_depletion":effect,
                    "positive_species_fraction":float(np.mean(value>0)),
                    "median_spaced_photos_per_species":float(records.loc[mask,f"{key}_n_retained_median_all_repeats"].median()),
                    "median_source_50km_pairs_of_selected_photos":float(records.loc[mask,f"{key}_n_photo_pairs_median_valid"].median()),
                    "null_mean":float(avg.mean()),
                    "permutation_p_upper":float((1+np.count_nonzero(avg>=effect))/(N_PERM+1)),
                })
                if n>=30:
                    rng=np.random.default_rng(seed(cohort,spacing,policy,"species_bootstrap"))
                    boot=value[rng.integers(n,size=(N_BOOT,n))].mean(axis=1)
                    item["species_bootstrap_ci95"]=[float(x) for x in np.quantile(boot,[.025,.975])]
                else:
                    item["species_bootstrap_ci95"]=None
            item["positive_bounded_evidence"]=bool(n>=30 and item.get("species_equal_mean_depletion",0)>0
                and item.get("permutation_p_upper",1)<.05
                and item.get("species_bootstrap_ci95") is not None
                and item["species_bootstrap_ci95"][0]>0)
            out[key]=item
    return out

def run_cohort(path:Path,cohort:str):
    source=orig.load(path,cohort)
    rows=[];nulls=[]
    for _,g in source.groupby("inat_taxon_id",sort=True):
        result=analyze_species(g,cohort)
        if result is not None:
            row,nul=result
            rows.append(row);nulls.append(nul)
    sp=pd.DataFrame(rows)
    exp_high,exp_n,exp_dep=EXPECTED[cohort]
    if source.inat_taxon_id.nunique()!=exp_high or len(sp)!=exp_n:
        raise RuntimeError("source high-depth or original spatial species count drift")
    historical=float(sp.original_allphoto_depletion.mean())
    if abs(historical-exp_dep)>1e-10:
        raise RuntimeError("original source 50km colour-depletion baseline drift")
    return {"n_original_source_high_depth_species":exp_high,
            "n_original_50km_local_30pair_species":len(sp),
            "original_mean_depletion":historical,
            "thinning":summarize_source(sp,nulls)},sp

def main():
    parser=argparse.ArgumentParser()
    for c in EXPECTED:
        parser.add_argument("--"+c,required=True,type=Path)
    parser.add_argument("--outdir",required=True,type=Path)
    args=parser.parse_args()
    payload={
       "schema":"fcp_original_photo_geographic_site_thinning_posthoc_v1",
       "date_jst":"2026-10-10",
       "role":"posthoc spatial-pseudoreplication falsification of visible photographed four-state colour",
       "spatial_separation_km":list(SPACINGS_KM),
       "original_local_radius_km":RADIUS_KM,
       "realizations_per_original_species":N_REALIZATIONS,
       "minimum_valid_realizations":MIN_VALID_REALIZATIONS,
       "minimum_retained_photos":MIN_RETAINED_PHOTOS,
       "minimum_retained_local_pairs":MIN_RETAINED_LOCAL_PAIRS,
       "required_species_per_cohort":30,
       "outcome_blind_selection":True,
       "permutation_reps":N_PERM,"species_bootstraps":N_BOOT,
       "source_SHA256":orig.SHA,
       "confirmatory_decisions_changed":False,
       "hard_nonclaims":[
           "A photographed coordinate cluster is not a biological individual or genetically connected population.",
           "Observed coordinates have their original uncertainty and may not identify the plant precisely.",
           "Spaced photos can retain shared observer, site, phylogeny, phenotype plasticity and camera biases.",
           "Averaging repeated geographic thinnings does not create independent empirical cohorts.",
           "P values are uncorrected across 3 geographic spacings x 2 observer restrictions.",
           "This does not provide independent biological colour annotation, genotype or fitness.",
       ]
    }
    allsp=[]
    for cohort in EXPECTED:
        payload[cohort],sp=run_cohort(getattr(args,cohort),cohort)
        allsp.append(sp)
    payload["replication_all_three"]={
        k:all(payload[c]["thinning"][k]["positive_bounded_evidence"] for c in EXPECTED)
        for k in payload["discovery"]["thinning"]}
    out=args.outdir;out.mkdir(parents=True,exist_ok=True)
    pd.concat(allsp,ignore_index=True).to_csv(out/"species_geo_thinning_support.csv",index=False)
    (out/"result.json").write_text(json.dumps(payload,indent=2)+"\n")
    print(json.dumps({c:payload[c] for c in EXPECTED},sort_keys=True),flush=True)

if __name__=="__main__":
    main()

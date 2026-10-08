#!/usr/bin/env python3
"""Does the white–chromatic axis exhibit geographic sorting beyond generic hue sorting?

Outcome-exposed FCP screen; visible photo colours, NOT plant pigment chemistry.
Within each species, geographic latitude tertiles and calendar-month photograph
composition are held fixed. Identical geocoordinates are NEVER split by rank.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

import run_polymorphism_specieswide_space_season_20261008 as common

PHOTO_SHA256 = common.SOURCE_SHA256
NONWHITE = ("yellow_orange", "red_pink", "blue_purple")
N_PERM = 199
N_BOOT = 1999
N_SIGNFLIP = 9999
MIN_PER_LABEL = 5
MIN_DISTINCT_BINS = 2
MIN_BIN_PHOTOS = 3
MIN_LAT_RANGE_DEGREES = 1.0
MIN_GLOBAL_W_OBS = 75
MIN_GLOBAL_W_IDENT = 50
MIN_PAIRED_SPECIES = 25
MASTER_SEED = 2026100893


def rng_seed(*parts: object) -> int:
    msg = "|".join(map(str,(MASTER_SEED,*parts))).encode()
    return int.from_bytes(hashlib.sha256(msg).digest()[:8],"little")


def latitude_bins(g: pd.DataFrame) -> tuple[np.ndarray, dict]:
    x = np.abs(g.latitude.to_numpy(float))
    if np.any(~np.isfinite(x)) or np.ptp(x)<MIN_LAT_RANGE_DEGREES:
        return np.full(len(g),-1,int), {"lat_span":float(np.ptp(x)),
                                           "valid":False}
    edges=np.quantile(x,[1/3,2/3])
    # Equal latitudes get the SAME band. Using ranking would invent geographic
    # association among photos taken at exactly one site.
    labels=np.searchsorted(edges,x,side="right").astype(int)
    return labels, {"valid":True,"edges_degrees":edges.tolist(),
                    "lat_span":float(np.ptp(x)),
                    "band_photo_counts":np.bincount(labels,minlength=3).tolist()}


def bin_association(binary: np.ndarray, bins: np.ndarray) -> float:
    """Cramer's V^2 = normalized between-bin Bernoulli variance (0..1)."""
    binary=np.asarray(binary,dtype=np.int8)
    bins=np.asarray(bins,int)
    if len(binary)!=len(bins) or len(binary)<2:
        raise ValueError("Invalid latitude/morph vectors")
    p=float(binary.mean())
    if not 0<p<1:
        return 0.0
    total=0.0
    for k in np.unique(bins):
        b=binary[bins==k]
        total += len(b)*(float(b.mean())-p)**2
    return float(total/(len(binary)*p*(1-p)))


def measure_pair(g:pd.DataFrame,selected:tuple[str,...],target:str,
                 bins:np.ndarray,month:np.ndarray,seed:int,
                 n_permutations:int=N_PERM) -> tuple[dict,np.ndarray]:
    morph=g.morph.astype(str).to_numpy()
    mask=np.isin(morph,selected)
    n=int(mask.sum())
    targetmask=(morph[mask]==target).astype(np.int8)
    near=bins[mask]
    mm=month[mask]
    if n<2*MIN_PER_LABEL or int(targetmask.sum())<MIN_PER_LABEL or int((1-targetmask).sum())<MIN_PER_LABEL:
        return {"eligible":False,"reason":"insufficient_specific_state_photo_counts"},np.zeros(n_permutations)
    counts=np.bincount(near,minlength=3)
    if int(np.count_nonzero(counts>=MIN_BIN_PHOTOS))<MIN_DISTINCT_BINS:
        return {"eligible":False,"reason":"insufficient_latitude_bin_coverage"},np.zeros(n_permutations)
    obs=bin_association(targetmask,near)
    strata=[np.flatnonzero(mm==m) for m in np.unique(mm)]
    swappable=int(sum(len(idx) for idx in strata if len(idx)>1 and
                      len(np.unique(targetmask[idx]))>1))
    rng=np.random.default_rng(seed)
    null=np.empty(n_permutations,float)
    if swappable<10:
        null.fill(obs)
        identifiable=False
    else:
        for i in range(n_permutations):
            perm=targetmask.copy()
            for idx in strata:
                perm[idx]=rng.permutation(targetmask[idx])
            null[i]=bin_association(perm,near)
        identifiable=bool(np.ptp(null)>1e-12)
        if not identifiable: null.fill(obs)
    sd=float(null.std(ddof=1))
    delta=0.0 if not identifiable else float(obs-null.mean())
    standardized=0.0 if not identifiable or sd<=0 else float(delta/sd)
    return {
        "eligible":True,"identifiable":bool(identifiable),
        "n_photos":n,"n_target":int(targetmask.sum()),
        "n_comparison":int((1-targetmask).sum()),
        "n_month_morph_exchangeable_photos":swappable,
        "occupied_lat_bins":int(np.count_nonzero(counts>=MIN_BIN_PHOTOS)),
        "observed_V2":float(obs),"null_mean_V2":float(null.mean()),
        "excess_V2":delta,"null_sd_V2":sd,
        "excess_null_SD":standardized,
        "nonexchangeable_photo_month_convention":"retained with zero identified additional effect; not biological zero",
    },null


def species_analysis(g:pd.DataFrame,cohort:str) -> tuple[dict,dict]:
    g=g.loc[g.month.notna()].copy().reset_index(drop=True)
    tid=int(g.inat_taxon_id.iloc[0]) if len(g) else -1
    result={"cohort":cohort,"inat_taxon_id":tid,
            "species":str(g.species.iloc[0]) if len(g) else "",
            "genus":str(g.species.iloc[0]).split()[0] if len(g) else "",
            "n_date_valid":len(g)}
    if len(g)<40:
        return result,{}
    lat,meta=latitude_bins(g)
    result.update({"lat_range_degrees":meta["lat_span"],
                   "lat_bins_valid":meta["valid"]})
    if not meta["valid"]:return result,{}
    counts=g.morph.value_counts()
    # Fix the nonwhite control groups by observed counts *before* geographical
    # permutations; no post-hoc choice based on which hue has best z-score.
    hues=sorted(NONWHITE,key=lambda h:(-int(counts.get(h,0)),h))
    lead,next_hue=hues[:2]
    pairs={
        "W_all":(("white",*NONWHITE),"white"),
        "W_lead":(("white",lead),"white"),
        "hue_lead_next":((lead,next_hue),lead),
    }
    result["lead_nonwhite_hue"]=lead
    result["next_nonwhite_hue"]=next_hue
    result["n_white"]=int(counts.get("white",0))
    result["n_lead_hue"]=int(counts.get(lead,0))
    result["n_next_hue"]=int(counts.get(next_hue,0))
    month=g.month.to_numpy(int)
    nulls={}
    for key,(members,target) in pairs.items():
        feat,perm=measure_pair(
            g,members,target,lat,month,
            rng_seed("morphsort",cohort,tid,key))
        for col,val in feat.items():
            result[f"{key}__{col}"]=val
        if feat["eligible"]:
            nulls[key]=perm
    if (result.get("W_lead__eligible") and result.get("hue_lead_next__eligible")
        and result["W_lead__identifiable"] and result["hue_lead_next__identifiable"]):
        result["matched_z_diff_W_minus_hue"] = (
            result["W_lead__excess_null_SD"] -
            result["hue_lead_next__excess_null_SD"])
        result["matched_delta_V2_W_minus_hue"] = (
            result["W_lead__excess_V2"] -
            result["hue_lead_next__excess_V2"])
    return result,nulls


def bootstrap_means(vals:np.ndarray,key:object) -> list[float] | None:
    if len(vals)==0:return None
    rng=np.random.default_rng(rng_seed("bootstrap",key))
    samples=vals[rng.integers(0,len(vals),(N_BOOT,len(vals)))].mean(axis=1)
    return [float(x) for x in np.quantile(samples,[.025,.975])]


def summarise_axis(rows:pd.DataFrame,cohort:str, key:str)->dict:
    good=rows.loc[rows.get(f"{key}__eligible",pd.Series(False,index=rows.index)).fillna(False).astype(bool)].copy()
    n=len(good)
    if not n:
        return {"n_geographically_evaluable":0,"n_conditional_identifiable":0,
                "status":"HOLD_NO_LATITUDE_OPPORTUNITY"}
    good[f"{key}__identifiable"]=good[f"{key}__identifiable"].fillna(False).astype(bool)
    ident=good.loc[good[f"{key}__identifiable"]]
    x=good[f"{key}__excess_V2"].to_numpy(float)
    z=good[f"{key}__excess_null_SD"].to_numpy(float)
    genus=good.groupby("genus")[f"{key}__excess_V2"].mean()
    return {
        "n_geographically_evaluable":int(n),
        "n_conditional_identifiable":int(len(ident)),
        "n_genera":int(len(genus)),
        "n_nonidentifiable_kept_as_zero":int(n-len(ident)),
        "mean_excess_Cramer_V2":float(x.mean()),
        "mean_standardized_null_SD_excess":float(z.mean()),
        "genus_balanced_mean_excess_V2":float(genus.mean()),
        "fraction_species_positive":float(np.mean(x>0)),
        "species_bootstrap_CI_excess_V2":bootstrap_means(x,(cohort,key,"V2")),
        "status":"EXPLORATORY_POPULATION_COVERAGE" if n>=MIN_GLOBAL_W_OBS and len(ident)>=MIN_GLOBAL_W_IDENT
                 else "HOLD_SPARSE_GEOGRAPHY_OR_CALENDAR",
    }


def summary_match(rows:pd.DataFrame,cohort:str)->dict:
    x=rows.matched_z_diff_W_minus_hue.dropna().to_numpy(float) if "matched_z_diff_W_minus_hue" in rows else np.array([])
    d=rows.loc[rows.matched_z_diff_W_minus_hue.notna()].copy() if len(x) else pd.DataFrame()
    if not len(x):return{"n_matched_species":0,"status":"HOLD_NO_COMPARABLE_MORPH_AXIS"}
    rng=np.random.default_rng(rng_seed(cohort,"axis-paired-signflip"))
    null=(rng.choice((-1.,1.),size=(N_SIGNFLIP,len(x)))*x).mean(axis=1)
    observed=float(x.mean())
    p2=float((1+np.count_nonzero(abs(null)>=abs(observed)))/(N_SIGNFLIP+1))
    g=d.groupby("genus").matched_z_diff_W_minus_hue.mean()
    return{
        "n_matched_species":int(len(x)),"n_genera":int(len(g)),
        "mean_standardized_white_minus_nonwhite_axis":observed,
        "median_standardized_axis_difference":float(np.median(x)),
        "species_bootstrap_95CI":bootstrap_means(x,(cohort,"paired")),
        "paired_signflip_two_sided_p":p2,
        "genus_balanced_axis_difference":float(g.mean()),
        "positive_species_fraction":float(np.mean(x>0)),
        "n_minimum_25_matched":bool(len(x)>=MIN_PAIRED_SPECIES),
        "status":"EXPLORATORY_COMPARISON_ELIGIBLE" if len(x)>=MIN_PAIRED_SPECIES else "HOLD_SPARSE_PAIRED_AXIS",
        "null_model":"both pairs each have their own composition/month preserving geography shuffle; effect normalized by pair-specific null SD",
        "claim_boundary":"relative photographic latitude sorting, not white-pigment fitness cost; exposure and demography can confound",
    }


def main():
    parser=argparse.ArgumentParser()
    for cohort in PHOTO_SHA256:
        parser.add_argument("--"+cohort,required=True,type=Path)
    parser.add_argument("--outdir",required=True,type=Path)
    args=parser.parse_args()
    args.outdir.mkdir(parents=True,exist_ok=True)
    res={
        "schema":"fcp_white_vs_nonwhite_hue_unsigned_latitude_sorting_posthoc_v1",
        "date_jst":"2026-10-08",
        "status":"complete",
        "role":"retrospective_morph_axis_specific_geographical_sorting_not_benefit_cost_selection",
        "source_measured_SHA256":PHOTO_SHA256,
        "n_perm_fixed_per_species":N_PERM,
        "n_species_bootstrap":N_BOOT,
        "n_paired_signflip":N_SIGNFLIP,
        "latitude_bins":"within-species absolute-latitude tertiles, tied latitudes never rank-split",
        "minimum_abs_lat_span_degrees":MIN_LAT_RANGE_DEGREES,
        "minimum_photos_per_individual_morph":MIN_PER_LABEL,
        "white_chromatic_axis":"white vs three coarse nonwhite photographic states pooled",
        "paired_axis":"white vs most abundant coloured category; most abundant vs second most abundant coloured category, same species",
        "no_signed_global_cline_imposed":True,
        "cross_cohort_strict_paired_white_specific_gate":False,
        "hard_nonclaims":[
            "photographic white is light/exposure coupled and not anthocyanin loss validated",
            "no genetic polymorphism or metabolite cost/benefit measured",
            "latitude sorting under a month-preserving null is also compatible with neutral dispersal/history and observer/site bias",
            "white/colour and nonwhite/nonwhite have different photo opportunity, despite null SD standardization",
            "within-species 1-degree span/photo-terciles do not represent natural population boundaries",
            "selected >40-classifiable photo species not a random sampling frame for world flower diversity",
            "the prior photo H2 axis is conditional on classification and not evidence white mutations are cheaper",
            "this post-outcome ecological screen cannot upgrade historical manuscript H1/H2 claims",
        ],
        "cohorts":{},
    }
    dfs=[]
    for cohort in PHOTO_SHA256:
        d,cover=common.load_source(getattr(args,cohort),cohort)
        rows=[]
        for taxon,g in d.groupby("inat_taxon_id",sort=True):
            row,_=species_analysis(g,cohort)
            rows.append(row)
        df=pd.DataFrame(rows)
        dfs.append(df)
        stats={key:summarise_axis(df,cohort,key) for key in ("W_all","W_lead","hue_lead_next")}
        res["cohorts"][cohort]={
            "source_coverage":cover,
            "n_with_intraspecific_latitude_span":int(df.lat_bins_valid.fillna(False).sum()),
            "axis":stats,"paired_white_specific_comparison":summary_match(df,cohort),
        }
    res["cross_cohort_strict_paired_white_specific_gate"]=bool(all(
        res["cohorts"][c]["paired_white_specific_comparison"]["n_minimum_25_matched"]
        and res["cohorts"][c]["paired_white_specific_comparison"]["mean_standardized_white_minus_nonwhite_axis"]>0
        and res["cohorts"][c]["paired_white_specific_comparison"]["species_bootstrap_95CI"][0]>0
        and res["cohorts"][c]["paired_white_specific_comparison"]["paired_signflip_two_sided_p"]<.05
        for c in PHOTO_SHA256))
    pd.concat(dfs,ignore_index=True).to_csv(args.outdir/"species_latitude_sorting_axis_metrics.csv",index=False)
    (args.outdir/"result.json").write_text(json.dumps(res,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps(res,indent=2,ensure_ascii=False))


if __name__=="__main__":
    main()

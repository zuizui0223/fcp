#!/usr/bin/env python3
from __future__ import annotations

import argparse
from collections import Counter
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.discrete.conditional_models import ConditionalLogit

BIO=["white","yellow","orange","red","pink","magenta","purple","blue","bronze"]
FRACTIONS=[f"flower_fraction_{x}" for x in BIO]
MORPHS=["white","yellow_orange","red_pink","blue_purple"]
MIN_CLASSIFIABLE=40
THRESHOLD=0.10
N_NULL=999
SEED=20260915
EPS=1e-12
OR_LOW=0.80
OR_HIGH=1.25
REFERENCE_N=158

def bool_series(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False).astype(bool)
    return s.fillna("").astype(str).str.strip().str.lower().isin(["true","1","yes"])

def normalize_rows(x: np.ndarray) -> np.ndarray:
    x=np.asarray(x,dtype=float)
    mass=x.sum(axis=1)
    if np.any(~np.isfinite(x)) or np.any(mass<=0):
        raise RuntimeError("invalid palette rows")
    return x/mass[:,None]

def deterministic_two_means(x: np.ndarray) -> tuple[np.ndarray,np.ndarray]:
    x=np.asarray(x,float)
    grand=x.mean(axis=0)
    i0=int(np.argmax(np.sum((x-grand)**2,axis=1)))
    d0=np.sum((x-x[i0])**2,axis=1)
    i1=int(np.argmax(d0))
    if float(d0[i1])<=EPS:
        labels=np.zeros(len(x),dtype=int); labels[len(x)//2:]=1
    else:
        centers=np.vstack([x[i0],x[i1]])
        labels=np.full(len(x),-1,dtype=int)
        for _ in range(200):
            dist=np.sum((x[:,None,:]-centers[None,:,:])**2,axis=2)
            new=np.argmin(dist,axis=1).astype(int)
            if np.all(new==new[0]):
                only=int(new[0]); new[int(np.argmax(dist[:,only]))]=1-only
            if np.array_equal(new,labels):
                break
            labels=new
            for k in (0,1):
                members=x[labels==k]
                if len(members)==0: raise RuntimeError("empty two-means cluster")
                centers[k]=members.mean(axis=0)
        else:
            raise RuntimeError("two-means did not converge")
    return labels,np.bincount(labels,minlength=2)

def coarse_second_fraction(g: pd.DataFrame) -> float:
    c=Counter(g["morph"].astype(str))
    order=sorted(MORPHS,key=lambda m:(-c.get(m,0),m))
    return float(c.get(order[1],0)/len(g))

def vector_table(work: pd.DataFrame) -> tuple[pd.DataFrame,dict[str,int]]:
    rows=[]
    gates={"nclass_ge40":0,"coarse_second_pass":0,"continuous_minor_pass":0,"nonzero_delta_pass":0}
    for sp,g in work.groupby("species",sort=True):
        n=len(g)
        if n<MIN_CLASSIFIABLE: continue
        gates["nclass_ge40"]+=1
        if coarse_second_fraction(g)<THRESHOLD: continue
        gates["coarse_second_pass"]+=1
        p=g[FRACTIONS].to_numpy(float)
        labels,counts=deterministic_two_means(np.sqrt(p))
        if counts.min()/n<THRESHOLD: continue
        gates["continuous_minor_pass"]+=1
        if counts[0]>counts[1]: major,minor=0,1
        elif counts[1]>counts[0]: major,minor=1,0
        else:
            c0=p[labels==0].mean(axis=0); c1=p[labels==1].mean(axis=0)
            major,minor=(0,1) if tuple(c0.tolist())<=tuple(c1.tolist()) else (1,0)
        delta=p[labels==minor].mean(axis=0)-p[labels==major].mean(axis=0)
        if not np.isclose(delta.sum(),0.0,atol=1e-10):
            raise RuntimeError("Delta outside zero-sum subspace")
        norm=float(np.linalg.norm(delta))
        if not np.isfinite(norm) or norm<=EPS: continue
        gates["nonzero_delta_pass"]+=1
        rows.append({"species":str(sp),**{f"delta_{BIO[j]}":float(delta[j]) for j in range(9)}})
    return pd.DataFrame(rows),gates

def qwhite() -> np.ndarray:
    q=np.array([1.0]+[-1.0/8.0]*8)
    q/=np.linalg.norm(q)
    return q

def units_from_table(tab: pd.DataFrame) -> np.ndarray:
    x=tab[[f"delta_{c}" for c in BIO]].to_numpy(float)
    return x/np.linalg.norm(x,axis=1)[:,None]

def statistic(u: np.ndarray,q: np.ndarray) -> float:
    return float(np.mean((u@q)**2))

def prepare_selected(work: pd.DataFrame,species_labels:list[str]) -> dict:
    sub=work.loc[work["species"].astype(str).isin(set(species_labels))].copy().reset_index(drop=True)
    p=sub[FRACTIONS].to_numpy(float)
    morph=sub["morph"].astype(str).to_numpy()
    sp=sub["species"].astype(str).to_numpy()
    return {
        "p":p,
        "morph_indices":{m:np.flatnonzero(morph==m) for m in MORPHS},
        "species_indices":{s:np.flatnonzero(sp==s) for s in sorted(species_labels)},
    }

def null_units(p:np.ndarray,indices:dict[str,np.ndarray])->np.ndarray:
    out=[]
    for sp in sorted(indices):
        x=p[indices[sp]]
        labels,counts=deterministic_two_means(np.sqrt(x))
        if counts[0]>counts[1]: major,minor=0,1
        elif counts[1]>counts[0]: major,minor=1,0
        else:
            c0=x[labels==0].mean(axis=0); c1=x[labels==1].mean(axis=0)
            major,minor=(0,1) if tuple(c0.tolist())<=tuple(c1.tolist()) else (1,0)
        d=x[labels==minor].mean(axis=0)-x[labels==major].mean(axis=0)
        norm=float(np.linalg.norm(d))
        if norm<=EPS: raise RuntimeError(f"degenerate null Delta: {sp}")
        out.append(d/norm)
    return np.vstack(out)

def run_h2_sensitivity(measured:pd.DataFrame,high_photo_ids:set[int])->dict:
    df=measured.loc[~measured["photo_id"].astype(int).isin(high_photo_ids)].copy()
    keep=bool_series(df["global_classifiable"]) & df["morph"].astype(str).isin(MORPHS)
    work=df.loc[keep].copy()
    work.loc[:,FRACTIONS]=normalize_rows(work[FRACTIONS].to_numpy(float))
    tab,gates=vector_table(work)
    n=len(tab)
    if n<20:
        return {"evaluable":False,"vector_species":n,"gates":gates,"support":False}
    q=qwhite()
    observed=statistic(units_from_table(tab),q)
    prep=prepare_selected(work,tab["species"].astype(str).tolist())
    rng=np.random.default_rng(SEED)
    null=np.empty(N_NULL,float)
    for b in range(N_NULL):
        p=prep["p"].copy()
        for morph in MORPHS:
            idx=prep["morph_indices"][morph]
            if len(idx)>1:
                p[idx]=prep["p"][idx[rng.permutation(len(idx))]]
        null[b]=statistic(null_units(p,prep["species_indices"]),q)
    p_upper=float((1+np.sum(null>=observed))/(N_NULL+1))
    return {
        "evaluable":True,
        "vector_species":int(n),
        "vector_retention":float(n/REFERENCE_N),
        "gates":gates,
        "observed_W":observed,
        "structured_null_replicates":N_NULL,
        "structured_null_median":float(np.median(null)),
        "structured_null_q025":float(np.quantile(null,0.025)),
        "structured_null_q975":float(np.quantile(null,0.975)),
        "structured_null_upper_p":p_upper,
        "support":bool(p_upper<0.05),
    }

def find_one(root:Path,name:str)->Path:
    hits=list(root.rglob(name))
    if len(hits)!=1:
        raise RuntimeError(f"expected one {name}, found {len(hits)}")
    return hits[0]

def main()->None:
    p=argparse.ArgumentParser()
    p.add_argument("--technical-seal-dir",type=Path,required=True)
    p.add_argument("--biological-artifact-dir",type=Path,required=True)
    p.add_argument("--output-dir",type=Path,required=True)
    args=p.parse_args()

    tech=pd.read_csv(args.technical_seal_dir/"technical_table.csv.gz",compression="gzip",dtype={"measurement_id":str})
    join=pd.read_csv(args.technical_seal_dir/"sealed_join_key.csv",dtype={"measurement_id":str})
    high=pd.read_csv(args.technical_seal_dir/"high_clip_ids.csv",dtype={"measurement_id":str})
    tech_summary=json.loads((args.technical_seal_dir/"technical_summary.json").read_text())

    measured=pd.read_csv(find_one(args.biological_artifact_dir,"polymorphism_h2_third_cohort_measured_photos_v1.csv"))
    frozen_result=json.loads(find_one(args.biological_artifact_dir,"result.json").read_text())
    if len(measured)!=49_900 or measured["photo_id"].nunique()!=49_900:
        raise RuntimeError("biological artifact denominator drift")
    if frozen_result.get("decision",{}).get("verdict")!="H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED":
        raise RuntimeError("unexpected frozen H2 verdict")

    sealed=join.merge(tech,on="measurement_id",how="inner",validate="one_to_one")
    if len(sealed)!=49_900:
        raise RuntimeError("sealed technical join lost rows")
    combined=sealed.merge(measured,on="photo_id",how="inner",validate="one_to_one")
    if len(combined)!=49_900:
        raise RuntimeError("biological join lost rows")

    reacq=combined["reacquired_image_sha256"].fillna("").astype(str)
    old=combined["image_sha256"].fillna("").astype(str)
    comparable=reacq.ne("") & old.ne("")
    drift=comparable & reacq.ne(old)
    source_drift_count=int(drift.sum())

    near=pd.to_numeric(combined["near_clip_fraction"],errors="coerce")
    classifiable=bool_series(combined["global_classifiable"]) & combined["morph"].astype(str).isin(MORPHS)
    model=combined.loc[classifiable & near.notna()].copy()
    model["white_response"]=(model["morph"].astype(str)=="white").astype(int)

    eligible=[]
    standardized=[]
    for sp,g in model.groupby("species",sort=True):
        x=pd.to_numeric(g["near_clip_fraction"],errors="coerce").to_numpy(float)
        y=g["white_response"].to_numpy(int)
        if len(g)<2 or len(np.unique(y))<2 or not np.isfinite(x).all() or float(np.std(x,ddof=0))<=0:
            continue
        z=(x-float(np.mean(x)))/float(np.std(x,ddof=0))
        gg=g[["species","white_response"]].copy()
        gg["z_near_clip"]=z
        standardized.append(gg)
        eligible.append(str(sp))

    coupling={"estimable":False,"species":len(eligible),"rows":0}
    if standardized:
        dat=pd.concat(standardized,ignore_index=True)
        coupling["rows"]=int(len(dat))
        try:
            fit=ConditionalLogit(
                dat["white_response"].to_numpy(float),
                dat[["z_near_clip"]].to_numpy(float),
                groups=dat["species"].astype(str).to_numpy(),
            ).fit(disp=False)
            beta=float(fit.params[0]); se=float(fit.bse[0])
            lo=beta-1.959963984540054*se; hi=beta+1.959963984540054*se
            coupling.update({
                "estimable":True,
                "beta":beta,
                "se":se,
                "odds_ratio":float(math.exp(beta)),
                "or_ci_low":float(math.exp(lo)),
                "or_ci_high":float(math.exp(hi)),
                "within_species_standardization_ddof":0,
            })
        except Exception as exc:
            coupling["error"]=f"{type(exc).__name__}:{str(exc)[:400]}"

    high_ids=set(high["measurement_id"].astype(str))
    high_photo_ids=set(join.loc[join["measurement_id"].isin(high_ids),"photo_id"].astype(int))
    sensitivity=run_h2_sensitivity(measured,high_photo_ids)
    retention=float(sensitivity.get("vector_retention",0.0)) if sensitivity.get("evaluable") else 0.0

    state="INDETERMINATE"
    reasons=[]
    if source_drift_count:
        reasons.append("source_byte_drift")
    if not coupling.get("estimable"):
        reasons.append("coupling_model_not_estimable")
    if retention<0.90:
        reasons.append("h2_vector_retention_below_0.90")

    if not reasons:
        lo=float(coupling["or_ci_low"]); hi=float(coupling["or_ci_high"])
        outside=hi<OR_LOW or lo>OR_HIGH
        support_flip=not bool(sensitivity.get("support"))
        inside=lo>=OR_LOW and hi<=OR_HIGH
        if outside or (support_flip and retention>=0.90):
            state="FLAGGED"
            if outside: reasons.append("coupling_ci_entirely_outside_frozen_interval")
            if support_flip: reasons.append("h2_support_flipped_after_high_clip_exclusion")
        elif inside and bool(sensitivity.get("support")) and retention>=0.90:
            state="CLEAR"
            reasons.append("coupling_ci_inside_frozen_interval_and_h2_support_retained")
        else:
            reasons.append("coupling_ci_overlaps_frozen_interval_boundary")

    out={
        "schema":"third_cohort_highlight_validity_result_v1",
        "date_jst":"2026-09-22",
        "role":"one_shot_post_confirmatory_measurement_validity_control",
        "frozen_h2_verdict":"H2_PROSPECTIVE_WHITE_AXIS_CONFIRMED",
        "confirmatory_verdict_changed":False,
        "technical_seal":tech_summary,
        "source_identity":{
            "comparable_rows":int(comparable.sum()),
            "source_drift_rows":source_drift_count,
            "reacquisition_failed_rows":int((combined["acquisition_status"]!="acquired_and_decode_verified").sum()),
        },
        "coupling_model":coupling,
        "high_clip":{
            "rows":int(len(high)),
            "threshold":float(tech_summary["high_clip_threshold"]),
            "rule":tech_summary["high_clip_rule"],
        },
        "h2_high_clip_sensitivity":sensitivity,
        "decision":{
            "state":state,
            "reasons":reasons,
            "or_equivalence_interval":[OR_LOW,OR_HIGH],
            "minimum_vector_retention":0.90,
        },
        "claim_boundary":[
            "This is post-confirmatory validity control, not a replacement confirmatory test.",
            "No alternative clipping threshold, technical predictor, H2 axis, or rescue analysis is authorized after this result.",
        ],
    }
    args.output_dir.mkdir(parents=True,exist_ok=True)
    (args.output_dir/"result.json").write_text(
        json.dumps(out,indent=2,sort_keys=True,allow_nan=False)+"\n",encoding="utf-8"
    )
    print(json.dumps(out,indent=2,sort_keys=True,allow_nan=False))

if __name__=="__main__":
    main()

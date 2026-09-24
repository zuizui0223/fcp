#!/usr/bin/env python3
from __future__ import annotations

import argparse
from collections import Counter
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

BIO = ["white","yellow","orange","red","pink","magenta","purple","blue","bronze"]
MORPHS = ["white","yellow_orange","red_pink","blue_purple"]
FRACTION_SUFFIX = [f"fraction_{x}" for x in BIO]
ALLROW_CONDITIONS = ["base","fixed_ev_m1_0","fixed_ev_m0_5","fixed_ev_p0_0","fixed_ev_p0_5","fixed_ev_p1_0"]
PRIMARY_FIXED = ["fixed_ev_m1_0","fixed_ev_m0_5","fixed_ev_p0_5","fixed_ev_p1_0"]
HEAVY_CONDITIONS = [
    "full_ev_m1_0","full_ev_m0_5","full_ev_p0_5","full_ev_p1_0","neutral_bg",
    "jitter_0","jitter_1","jitter_2","jitter_3","jitter_4","jitter_5","jitter_6",
]
MIN_CLASSIFIABLE = 40
THRESHOLD = 0.10
EPS = 1e-12


def bool_series(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False).astype(bool)
    return s.fillna("").astype(str).str.strip().str.casefold().isin({"1","true","yes"})


def condition_classifiable(df: pd.DataFrame, prefix: str) -> pd.Series:
    return (
        df[f"{prefix}_status"].astype(str).eq("classified_four_state_morph")
        & df[f"{prefix}_morph"].astype(str).isin(MORPHS)
    )


def fraction_cols(prefix: str) -> list[str]:
    return [f"{prefix}_{suffix}" for suffix in FRACTION_SUFFIX]


def normalize_rows(x: np.ndarray) -> np.ndarray:
    x=np.asarray(x,float)
    mass=x.sum(axis=1)
    if np.any(~np.isfinite(x)) or np.any(mass<=0):
        raise ValueError("invalid palette row")
    return x/mass[:,None]


def hellinger_rows(a: np.ndarray,b: np.ndarray) -> np.ndarray:
    aa=normalize_rows(a)
    bb=normalize_rows(b)
    return np.sqrt(np.sum((np.sqrt(aa)-np.sqrt(bb))**2,axis=1))/math.sqrt(2.0)


def summarize_image_condition(df: pd.DataFrame, prefix: str) -> dict:
    base_ok=condition_classifiable(df,"base")
    cf_ok=condition_classifiable(df,prefix)
    paired=base_ok & cf_ok
    out={
        "rows":int(len(df)),
        "base_classifiable_rows":int(base_ok.sum()),
        "counterfactual_classifiable_rows":int(cf_ok.sum()),
        "paired_classifiable_rows":int(paired.sum()),
        "classifiable_to_nonclassifiable_rows":int((base_ok & ~cf_ok).sum()),
        "nonclassifiable_to_classifiable_rows":int((~base_ok & cf_ok).sum()),
    }
    if paired.any():
        b=df.loc[paired,"base_morph"].astype(str)
        c=df.loc[paired,f"{prefix}_morph"].astype(str)
        flip=(b.to_numpy()!=c.to_numpy())
        out["coarse_state_flip_probability_equal_image"]=float(np.mean(flip))
        trans={}
        for left in MORPHS:
            trans[left]={}
            for right in MORPHS:
                trans[left][right]=int(((b==left)&(c==right)).sum())
        out["transition_matrix"]=trans

        a=df.loc[paired,fraction_cols("base")].apply(pd.to_numeric,errors="coerce").to_numpy(float)
        z=df.loc[paired,fraction_cols(prefix)].apply(pd.to_numeric,errors="coerce").to_numpy(float)
        finite=np.isfinite(a).all(axis=1)&np.isfinite(z).all(axis=1)&(a.sum(axis=1)>0)&(z.sum(axis=1)>0)
        h=np.full(len(a),np.nan)
        if finite.any():
            h[finite]=hellinger_rows(a[finite],z[finite])
        out["hellinger_evaluable_rows"]=int(np.isfinite(h).sum())
        out["hellinger_mean_equal_image"]=float(np.nanmean(h)) if np.isfinite(h).any() else None
        out["hellinger_median_equal_image"]=float(np.nanmedian(h)) if np.isfinite(h).any() else None

        tmp=pd.DataFrame({
            "species":df.loc[paired,"species"].astype(str).to_numpy(),
            "flip":flip.astype(float),
            "hellinger":h,
        })
        species=tmp.groupby("species",sort=True).agg(
            flip=("flip","mean"),
            hellinger=("hellinger","mean"),
        )
        out["equal_species_n"]=int(len(species))
        out["coarse_state_flip_probability_equal_species"]=float(species["flip"].mean()) if len(species) else None
        hsp=species["hellinger"].dropna()
        out["hellinger_mean_equal_species"]=float(hsp.mean()) if len(hsp) else None
        out["hellinger_median_equal_species"]=float(hsp.median()) if len(hsp) else None
    else:
        out.update({
            "coarse_state_flip_probability_equal_image":None,
            "transition_matrix":{},
            "hellinger_evaluable_rows":0,
            "hellinger_mean_equal_image":None,
            "hellinger_median_equal_image":None,
            "equal_species_n":0,
            "coarse_state_flip_probability_equal_species":None,
            "hellinger_mean_equal_species":None,
            "hellinger_median_equal_species":None,
        })
    return out


def D_table(df: pd.DataFrame,prefix: str,min_n: int=1) -> pd.DataFrame:
    ok=condition_classifiable(df,prefix)
    rows=[]
    for (panel,sp),g in df.loc[ok].groupby(["panel","species"],sort=True):
        n=len(g)
        if n<min_n:
            continue
        counts=g[f"{prefix}_morph"].astype(str).value_counts().reindex(MORPHS,fill_value=0).to_numpy(float)
        p=counts/counts.sum()
        rows.append({"panel":str(panel),"species":str(sp),"n_classifiable":int(n),"D":float(1.0-np.sum(p*p))})
    return pd.DataFrame(rows)


def finite_or_none(value: float) -> float | None:
    v=float(value)
    return v if np.isfinite(v) else None


def describe_numeric(s: pd.Series) -> dict:
    if len(s)==0:
        return {}
    return {str(k): finite_or_none(v) for k,v in s.describe().to_dict().items()}


def spearman(x: np.ndarray,y: np.ndarray) -> float | None:
    rx=pd.Series(x).rank(method="average").to_numpy(float)
    ry=pd.Series(y).rank(method="average").to_numpy(float)
    return finite_or_none(np.corrcoef(rx,ry)[0,1]) if len(x)>1 else None


def ccc(x: np.ndarray,y: np.ndarray) -> float | None:
    x=np.asarray(x,float); y=np.asarray(y,float)
    if len(x)<2:
        return None
    vx=float(np.var(x,ddof=1)); vy=float(np.var(y,ddof=1))
    cov=float(np.cov(x,y,ddof=1)[0,1])
    denom=vx+vy+(float(np.mean(x))-float(np.mean(y)))**2
    return finite_or_none(2*cov/denom) if denom>0 else None


def D_comparison(df: pd.DataFrame,prefix: str) -> dict:
    base=D_table(df,"base",MIN_CLASSIFIABLE)
    cf=D_table(df,prefix,MIN_CLASSIFIABLE)
    joined=base.merge(cf,on=["panel","species"],suffixes=("_base","_cf"),validate="one_to_one")
    out={}
    for label,sub in [("P",joined.loc[joined["panel"]=="P"]),("N",joined.loc[joined["panel"]=="N"]),("pooled",joined)]:
        x=sub["D_base"].to_numpy(float); y=sub["D_cf"].to_numpy(float)
        delta=y-x
        out[label]={
            "paired_species":int(len(sub)),
            "spearman_rho":spearman(x,y),
            "ccc":ccc(x,y),
            "mean_absolute_D_difference":float(np.mean(np.abs(delta))) if len(sub) else None,
            "median_absolute_D_difference":float(np.median(np.abs(delta))) if len(sub) else None,
            "signed_mean_D_difference_cf_minus_base":float(np.mean(delta)) if len(sub) else None,
        }
    return out


def deterministic_two_means(x: np.ndarray) -> tuple[np.ndarray,np.ndarray]:
    x=np.asarray(x,float)
    grand=x.mean(axis=0)
    i0=int(np.argmax(np.sum((x-grand)**2,axis=1)))
    d0=np.sum((x-x[i0])**2,axis=1)
    i1=int(np.argmax(d0))
    if float(d0[i1])<=EPS:
        labels=np.zeros(len(x),dtype=int); labels[len(x)//2:]=1
        return labels,np.bincount(labels,minlength=2)
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
            if len(members)==0:
                raise RuntimeError("empty two-means cluster")
            centers[k]=members.mean(axis=0)
    else:
        raise RuntimeError("two-means did not converge")
    return labels,np.bincount(labels,minlength=2)


def vector_table(df: pd.DataFrame,prefix: str,min_classifiable: int=MIN_CLASSIFIABLE) -> tuple[pd.DataFrame,dict]:
    ok=condition_classifiable(df,prefix)
    work=df.loc[ok].copy()
    rows=[]
    gates={"nclass_pass":0,"coarse_second_pass":0,"continuous_minor_pass":0,"nonzero_delta_pass":0}
    for sp,g in work.groupby("species",sort=True):
        n=len(g)
        if n<min_classifiable:
            continue
        gates["nclass_pass"]+=1
        counts=Counter(g[f"{prefix}_morph"].astype(str))
        order=sorted(MORPHS,key=lambda m:(-counts.get(m,0),m))
        second=float(counts.get(order[1],0)/n)
        if second<THRESHOLD:
            continue
        gates["coarse_second_pass"]+=1
        p=normalize_rows(g[fraction_cols(prefix)].apply(pd.to_numeric,errors="coerce").to_numpy(float))
        labels,sizes=deterministic_two_means(np.sqrt(p))
        if sizes.min()/n<THRESHOLD:
            continue
        gates["continuous_minor_pass"]+=1
        if sizes[0]>sizes[1]: major,minor=0,1
        elif sizes[1]>sizes[0]: major,minor=1,0
        else:
            c0=p[labels==0].mean(axis=0); c1=p[labels==1].mean(axis=0)
            major,minor=(0,1) if tuple(c0.tolist())<=tuple(c1.tolist()) else (1,0)
        delta=p[labels==minor].mean(axis=0)-p[labels==major].mean(axis=0)
        norm=float(np.linalg.norm(delta))
        if not np.isfinite(norm) or norm<=EPS:
            continue
        gates["nonzero_delta_pass"]+=1
        row={"species":str(sp)}
        row.update({f"delta_{BIO[j]}":float(delta[j]) for j in range(9)})
        rows.append(row)
    return pd.DataFrame(rows),gates


def qwhite() -> np.ndarray:
    q=np.array([1.0]+[-1.0/8.0]*8,float)
    return q/np.linalg.norm(q)


def W_from_vectors(tab: pd.DataFrame) -> float | None:
    if len(tab)==0:
        return None
    x=tab[[f"delta_{b}" for b in BIO]].to_numpy(float)
    u=x/np.linalg.norm(x,axis=1)[:,None]
    return float(np.mean((u@qwhite())**2))


def main() -> None:
    p=argparse.ArgumentParser()
    p.add_argument("--biological-seal-dir",type=Path,required=True)
    p.add_argument("--technical-seal-dir",type=Path,required=True)
    p.add_argument("--output-dir",type=Path,required=True)
    args=p.parse_args()

    bsum=json.loads((args.biological_seal_dir/"biological_summary.json").read_text())
    if bsum.get("status")!="BIOLOGICAL_PASS_B_SEALED":
        raise RuntimeError("Pass B seal is not terminal")
    if int(bsum.get("rows",-1))!=40000 or int(bsum.get("partition_receipts",-1))!=256:
        raise RuntimeError("Pass B denominator drift")
    if bsum.get("source_identity",{}).get("replacement_rows")!=0:
        raise RuntimeError("Pass B replacement detected")

    df=pd.read_csv(args.biological_seal_dir/"biological_table.csv.gz",compression="gzip",dtype={"measurement_id":str})
    tech=pd.read_csv(args.technical_seal_dir/"technical_table.csv.gz",compression="gzip",dtype={"measurement_id":str})
    if len(df)!=40000 or df["measurement_id"].nunique()!=40000:
        raise RuntimeError("biological table census drift")
    if len(tech)!=40000 or tech["measurement_id"].nunique()!=40000:
        raise RuntimeError("technical table census drift")

    exact=df["source_identity_status"].astype(str).eq("exact_source_sha_match")
    work=df.loc[exact].copy()

    mv1={"all_row_fixed_mask":{},"heavy_subset":{}}
    for cond in PRIMARY_FIXED+["fixed_ev_p0_0"]:
        mv1["all_row_fixed_mask"][cond]=summarize_image_condition(work,cond)
        lab_col=f"{cond}_flower_lab_delta_e_to_base"
        if lab_col in tech.columns:
            z=pd.to_numeric(tech.loc[tech["measurement_id"].isin(work["measurement_id"]),lab_col],errors="coerce")
            mv1["all_row_fixed_mask"][cond]["pass_t_lab_delta_e_evaluable_rows"]=int(z.notna().sum())
            mv1["all_row_fixed_mask"][cond]["pass_t_lab_delta_e_mean"]=float(z.mean()) if z.notna().any() else None
            mv1["all_row_fixed_mask"][cond]["pass_t_lab_delta_e_median"]=float(z.median()) if z.notna().any() else None

    heavy=work.loc[bool_series(work["heavy_counterfactual"])].copy()
    for cond in HEAVY_CONDITIONS:
        if f"{cond}_status" in heavy.columns:
            mv1["heavy_subset"][cond]=summarize_image_condition(heavy,cond)

    mv2={"all_row_fixed_mask":{},"heavy_subset":{}}
    for cond in PRIMARY_FIXED+["fixed_ev_p0_0"]:
        mv2["all_row_fixed_mask"][cond]=D_comparison(work,cond)
    for cond in HEAVY_CONDITIONS:
        if f"{cond}_status" not in heavy.columns:
            continue
        base_h=D_table(heavy,"base",1)
        cf_h=D_table(heavy,cond,1)
        j=base_h.merge(cf_h,on=["panel","species"],suffixes=("_base","_cf"),validate="one_to_one")
        delta=j["D_cf"].to_numpy(float)-j["D_base"].to_numpy(float)
        mv2["heavy_subset"][cond]={
            "paired_species":int(len(j)),
            "mean_absolute_D_difference":float(np.mean(np.abs(delta))) if len(j) else None,
            "median_absolute_D_difference":float(np.median(np.abs(delta))) if len(j) else None,
            "signed_mean_D_difference_cf_minus_base":float(np.mean(delta)) if len(j) else None,
            "base_classifiable_n_distribution":describe_numeric(base_h["n_classifiable"]),
            "counterfactual_classifiable_n_distribution":describe_numeric(cf_h["n_classifiable"]),
        }

    mv3={"all_row":{},"technical_strata_exclusion":{}}
    base_tab,base_gates=vector_table(work,"base")
    base_W=W_from_vectors(base_tab)
    base_species=set(base_tab["species"].astype(str)) if len(base_tab) else set()
    mv3["all_row"]["base"]={"gates":base_gates,"vector_species":int(len(base_tab)),"W":base_W}
    for cond in PRIMARY_FIXED+["fixed_ev_p0_0"]:
        tab,gates=vector_table(work,cond)
        sp=set(tab["species"].astype(str)) if len(tab) else set()
        w=W_from_vectors(tab)
        mv3["all_row"][cond]={
            "gates":gates,
            "vector_species":int(len(tab)),
            "W":w,
            "W_minus_base":float(w-base_W) if w is not None and base_W is not None else None,
            "vector_retention_vs_base":float(len(sp)/len(base_species)) if base_species else None,
            "species_overlap_with_base":int(len(sp & base_species)),
        }

    strata=[
        "stratum_high_flower_near_clip",
        "stratum_high_background_near_clip",
        "stratum_low_flower_background_lab_separation",
        "stratum_roi_unstable_fixed",
    ]
    for col in strata:
        if col not in work.columns:
            continue
        sub=work.loc[~bool_series(work[col])].copy()
        tab,gates=vector_table(sub,"base")
        w=W_from_vectors(tab)
        mv3["technical_strata_exclusion"][col]={
            "remaining_rows":int(len(sub)),
            "excluded_rows":int(len(work)-len(sub)),
            "gates":gates,
            "vector_species":int(len(tab)),
            "W":w,
            "W_minus_full_base":float(w-base_W) if w is not None and base_W is not None else None,
        }

    result={
        "schema":"fcp_v2_mv1_mv3_result_v1",
        "status":"POST_PASS_B_MEASUREMENT_VALIDITY_SUMMARY",
        "rows":40000,
        "exact_source_sha_match_rows":int(exact.sum()),
        "source_byte_drift_rows":int(df["source_identity_status"].astype(str).eq("source_byte_drift").sum()),
        "pass_b_acquisition_failed_rows":int(df["source_identity_status"].astype(str).eq("pass_b_acquisition_failed").sum()),
        "MV1":mv1,
        "MV2":mv2,
        "MV3":mv3,
        "MV4":{"status":"DEFERRED_TO_EXISTING_SPATIAL_STATISTIC_REUSE"},
        "no_rescue":True,
    }
    args.output_dir.mkdir(parents=True,exist_ok=True)
    (args.output_dir/"result.json").write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2,sort_keys=True,allow_nan=False))


if __name__=="__main__":
    main()

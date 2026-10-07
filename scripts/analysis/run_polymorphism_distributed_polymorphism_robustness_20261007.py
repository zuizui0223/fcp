#!/usr/bin/env python3
"""Robustness tests for the FCP distributed-polymorphism result."""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd

EARTH_RADIUS_KM=6371.0088
MORPHS=["white","yellow_orange","red_pink","blue_purple"]
NONWHITE=["yellow_orange","red_pink","blue_purple"]
BIO9=["white","yellow","orange","red","pink","magenta","purple","blue","bronze"]
LEGACY_COLOUR=[f"palette_count_{x}" for x in BIO9]
THIRD_COLOUR=[f"flower_fraction_{x}" for x in BIO9]
PERMUTATIONS=199
MIN_LOCAL_PAIRS=30
MIN_NONWHITE_ROWS=30
RADIUS=50.0
MASTER_SEED=2026100743
SHA={
"discovery":"ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4",
"validation":"0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6",
"third":"57630fc9f281bce94a0c40a70aaf7bce879dde93d6154175adcd021e8f5c1186",
}

def sha256(p):
    h=hashlib.sha256()
    with Path(p).open("rb") as f:
        for ch in iter(lambda:f.read(1<<20),b""): h.update(ch)
    return h.hexdigest()

def as_bool(s):
    if pd.api.types.is_bool_dtype(s): return s.fillna(False).astype(bool)
    return s.fillna("").astype(str).str.strip().str.casefold().isin({"true","1","yes","y"})

def seed(*parts):
    b="|".join(map(str,(MASTER_SEED,*parts))).encode()
    return int.from_bytes(hashlib.sha256(b).digest()[:8],"little")

def geo(lat,lon):
    latr=np.deg2rad(np.asarray(lat,float)); lonr=np.deg2rad(np.asarray(lon,float))
    c=np.cos(latr); xyz=np.column_stack([c*np.cos(lonr),c*np.sin(lonr),np.sin(latr)])
    return np.arccos(np.clip(xyz@xyz.T,-1,1))*EARTH_RADIUS_KM

def load(path,cohort):
    if sha256(path)!=SHA[cohort]: raise RuntimeError(f"{cohort} sha mismatch")
    d=pd.read_csv(path,low_memory=False)
    colour_cols=THIRD_COLOUR if cohort=="third" else LEGACY_COLOUR
    req={"inat_taxon_id","species","photo_id","observer_id","latitude","longitude","observed_on","morph","global_classifiable",*colour_cols}
    miss=sorted(req-set(d.columns))
    if miss: raise RuntimeError(f"{cohort} missing {miss}")
    d=d.loc[as_bool(d.global_classifiable)&d.morph.astype(str).isin(MORPHS)].copy()
    n=d.groupby("inat_taxon_id").size()
    d=d.loc[d.inat_taxon_id.isin(n[n>=40].index)].copy()
    return d.sort_values(["inat_taxon_id","photo_id"],kind="stable").reset_index(drop=True)


def pairwise_jsd(prob):
    p=np.asarray(prob,float)
    mass=p.sum(axis=1)
    if p.ndim!=2 or np.any(~np.isfinite(p)) or np.any(p<0) or np.any(mass<=0):
        raise ValueError("invalid continuous colour rows")
    p=p/mass[:,None]
    a=p[:,None,:]; b=p[None,:,:]; m=.5*(a+b)
    with np.errstate(divide="ignore",invalid="ignore"):
        ka=np.where(a>0,a*np.log2(a/m),0).sum(axis=2)
        kb=np.where(b>0,b*np.log2(b/m),0).sum(axis=2)
    return np.clip(.5*(ka+kb),0,1)

def continuous_species_test(g,cohort):
    cols=THIRD_COLOUR if cohort=="third" else LEGACY_COLOUR
    gg=g.dropna(subset=["latitude","longitude",*cols]).copy()
    if len(gg)<40: return None,None
    gg=gg.sort_values("photo_id",kind="stable").reset_index(drop=True)
    n=len(gg); u,v=np.triu_indices(n,k=1)
    dist=geo(gg.latitude.to_numpy(float),gg.longitude.to_numpy(float))[u,v]
    mask=dist<=RADIUS
    if int(mask.sum())<MIN_LOCAL_PAIRS: return None,None
    jsd=pairwise_jsd(gg[cols].to_numpy(float))
    allv=jsd[u,v]
    localv=allv[mask]
    overall=float(np.mean(allv)); local=float(np.mean(localv)); dep=overall-local
    # Permuting complete rows only reorders the fixed pairwise-distance matrix.
    rng=np.random.default_rng(seed(cohort,int(gg.inat_taxon_id.iloc[0]),"continuous"))
    local_null=np.empty(PERMUTATIONS,float)
    uu=u[mask]; vv=v[mask]
    for i in range(PERMUTATIONS):
        p=rng.permutation(n)
        local_null[i]=float(np.mean(jsd[p[uu],p[vv]]))
    dep_null=overall-local_null
    row={"cohort":cohort,"mode":"continuous_nine_colour","inat_taxon_id":int(gg.inat_taxon_id.iloc[0]),
         "species":str(gg.species.iloc[0]),"n_rows":int(n),"n_local_pairs":int(mask.sum()),
         "overall_mean_jsd":overall,"local_mean_jsd":local,"depletion":dep}
    return row,dep_null

def summarize_continuous(d,cohort):
    rows=[]; nulls=[]
    for _,g in d.groupby("inat_taxon_id",sort=True):
        row,nul=continuous_species_test(g,cohort)
        if row is not None:
            rows.append(row); nulls.append(nul)
    if not rows: return {"evaluable":False,"n_species":0},pd.DataFrame()
    df=pd.DataFrame(rows); mat=np.vstack(nulls); nd=mat.mean(axis=0)
    obs=float(df.depletion.mean())
    p=float((1+np.count_nonzero(nd>=obs))/(PERMUTATIONS+1))
    return {"evaluable":True,"n_species":int(len(df)),"mean_depletion":obs,
            "median_depletion":float(df.depletion.median()),
            "positive_species_fraction":float(np.mean(df.depletion>0)),
            "null_mean":float(nd.mean()),"null_q025":float(np.quantile(nd,.025)),
            "null_q975":float(np.quantile(nd,.975)),"p_upper":p,
            "supported":bool(obs>0 and p<.05)},df


def quarter_stratified_species_test(g,cohort):
    gg=g.dropna(subset=["latitude","longitude","observed_on"]).copy()
    if len(gg)<40: return None,None
    dates=pd.to_datetime(gg.observed_on,errors="coerce")
    gg=gg.loc[dates.notna()].copy()
    dates=pd.to_datetime(gg.observed_on,errors="raise")
    if len(gg)<40: return None,None
    gg=gg.assign(_quarter=dates.dt.quarter.to_numpy())
    gg=gg.sort_values("photo_id",kind="stable").reset_index(drop=True)
    labels=pd.Categorical(gg.morph.astype(str),categories=MORPHS).codes.astype(np.int8)
    quarters=gg._quarter.to_numpy(int)
    n=len(gg); u,v=np.triu_indices(n,k=1)
    dist=geo(gg.latitude.to_numpy(float),gg.longitude.to_numpy(float))[u,v]
    mask=dist<=RADIUS
    if int(mask.sum())<MIN_LOCAL_PAIRS: return None,None
    uu=u[mask]; vv=v[mask]
    D_pair=float(np.mean(labels[u]!=labels[v]))
    D_local=float(np.mean(labels[uu]!=labels[vv]))
    dep=D_pair-D_local
    rng=np.random.default_rng(seed(cohort,int(gg.inat_taxon_id.iloc[0]),"quarter_stratified"))
    null=np.empty(PERMUTATIONS,float)
    groups=[np.flatnonzero(quarters==q) for q in sorted(np.unique(quarters))]
    for i in range(PERMUTATIONS):
        p=labels.copy()
        for idx in groups:
            p[idx]=rng.permutation(labels[idx])
        null[i]=D_pair-float(np.mean(p[uu]!=p[vv]))
    if np.ptp(null)<=1e-15:
        return None,None
    return {
        "cohort":cohort,"mode":"quarter_stratified","inat_taxon_id":int(gg.inat_taxon_id.iloc[0]),
        "species":str(gg.species.iloc[0]),"n_rows":int(n),"n_local_pairs":int(mask.sum()),
        "D_pair":D_pair,"D_local":D_local,"depletion":dep
    },null

def crossyear_species_test(g,cohort):
    gg=g.dropna(subset=["latitude","longitude","observed_on"]).copy()
    if len(gg)<40: return None,None
    dates=pd.to_datetime(gg.observed_on,errors="coerce")
    gg=gg.loc[dates.notna()].copy()
    dates=pd.to_datetime(gg.observed_on,errors="raise")
    if len(gg)<40: return None,None
    gg=gg.assign(_year=dates.dt.year.to_numpy())
    gg=gg.sort_values("photo_id",kind="stable").reset_index(drop=True)
    labels=pd.Categorical(gg.morph.astype(str),categories=MORPHS).codes.astype(np.int8)
    years=gg._year.to_numpy(int)
    n=len(gg); u,v=np.triu_indices(n,k=1)
    dist=geo(gg.latitude.to_numpy(float),gg.longitude.to_numpy(float))[u,v]
    mask=(dist<=RADIUS)&(years[u]!=years[v])
    if int(mask.sum())<MIN_LOCAL_PAIRS: return None,None
    uu=u[mask]; vv=v[mask]
    D_pair=float(np.mean(labels[u]!=labels[v]))
    D_local=float(np.mean(labels[uu]!=labels[vv]))
    dep=D_pair-D_local
    rng=np.random.default_rng(seed(cohort,int(gg.inat_taxon_id.iloc[0]),"crossyear"))
    perms=np.stack([rng.permutation(labels) for _ in range(PERMUTATIONS)])
    null=D_pair-np.mean(perms[:,uu]!=perms[:,vv],axis=1)
    return {
        "cohort":cohort,"mode":"crossyear","inat_taxon_id":int(gg.inat_taxon_id.iloc[0]),
        "species":str(gg.species.iloc[0]),"n_rows":int(n),"n_local_pairs":int(mask.sum()),
        "D_pair":D_pair,"D_local":D_local,"depletion":dep
    },null

def summarize_custom(d,cohort,fn):
    rows=[]; nulls=[]
    for _,g in d.groupby("inat_taxon_id",sort=True):
        row,nul=fn(g,cohort)
        if row is not None:
            rows.append(row); nulls.append(nul)
    if not rows: return {"evaluable":False,"n_species":0},pd.DataFrame()
    df=pd.DataFrame(rows); mat=np.vstack(nulls); nd=mat.mean(axis=0)
    obs=float(df.depletion.mean())
    p=float((1+np.count_nonzero(nd>=obs))/(PERMUTATIONS+1))
    return {"evaluable":True,"n_species":int(len(df)),"mean_depletion":obs,
            "median_depletion":float(df.depletion.median()),
            "positive_species_fraction":float(np.mean(df.depletion>0)),
            "null_mean":float(nd.mean()),"null_q025":float(np.quantile(nd,.025)),
            "null_q975":float(np.quantile(nd,.975)),"p_upper":p,
            "supported":bool(obs>0 and p<.05)},df

def species_test(g,cohort,mode):
    if mode=="different_observer":
        gg=g.copy()
        allowed=MORPHS
        minrows=40
    else:
        gg=g.loc[g.morph.astype(str).isin(NONWHITE)].copy()
        allowed=NONWHITE
        minrows=MIN_NONWHITE_ROWS
    gg=gg.dropna(subset=["latitude","longitude"]).copy()
    if len(gg)<minrows: return None,None
    gg=gg.sort_values("photo_id",kind="stable").reset_index(drop=True)
    labels=pd.Categorical(gg.morph.astype(str),categories=allowed).codes.astype(np.int8)
    n=len(gg); u,v=np.triu_indices(n,k=1)
    dist=geo(gg.latitude.to_numpy(float),gg.longitude.to_numpy(float))[u,v]
    mask=dist<=RADIUS
    if mode=="different_observer":
        obs=gg.observer_id.fillna("").astype(str).to_numpy()
        mask &= (obs[u]!="") & (obs[v]!="") & (obs[u]!=obs[v])
    if int(mask.sum())<MIN_LOCAL_PAIRS: return None,None
    uu=u[mask]; vv=v[mask]
    all_diff=(labels[u]!=labels[v])
    D_pair=float(np.mean(all_diff))
    D_local=float(np.mean(labels[uu]!=labels[vv]))
    depletion=D_pair-D_local
    rng=np.random.default_rng(seed(cohort,int(gg.inat_taxon_id.iloc[0]),mode))
    perms=np.stack([rng.permutation(labels) for _ in range(PERMUTATIONS)])
    local_null=np.mean(perms[:,uu]!=perms[:,vv],axis=1)
    dep_null=D_pair-local_null
    row={
      "cohort":cohort,"mode":mode,"inat_taxon_id":int(gg.inat_taxon_id.iloc[0]),
      "species":str(gg.species.iloc[0]),"n_rows":int(n),"n_local_pairs":int(mask.sum()),
      "D_pair":D_pair,"D_local":D_local,"depletion":depletion
    }
    return row,dep_null

def summarize(d,cohort,mode):
    rows=[]; nulls=[]
    for _,g in d.groupby("inat_taxon_id",sort=True):
        row,nul=species_test(g,cohort,mode)
        if row is not None:
            rows.append(row); nulls.append(nul)
    if not rows: return {"evaluable":False,"n_species":0},pd.DataFrame()
    df=pd.DataFrame(rows)
    mat=np.vstack(nulls)
    obs=float(df.depletion.mean())
    nd=mat.mean(axis=0)
    p=float((1+np.count_nonzero(nd>=obs))/(PERMUTATIONS+1))
    return {
      "evaluable":True,"n_species":int(len(df)),
      "mean_depletion":obs,"median_depletion":float(df.depletion.median()),
      "positive_species_fraction":float(np.mean(df.depletion>0)),
      "null_mean":float(nd.mean()),"null_q025":float(np.quantile(nd,.025)),
      "null_q975":float(np.quantile(nd,.975)),"p_upper":p,
      "supported":bool(obs>0 and p<.05)
    },df

def run(path,cohort):
    d=load(path,cohort)
    r1,d1=summarize(d,cohort,"different_observer")
    r2,d2=summarize(d,cohort,"nonwhite_only")
    r3,d3=summarize_continuous(d,cohort)
    r4,d4=summarize_custom(d,cohort,quarter_stratified_species_test)
    r5,d5=summarize_custom(d,cohort,crossyear_species_test)
    return {"R1_different_observer":r1,"R2_nonwhite_only":r2,"R3_continuous_nine_colour":r3,
            "R4_quarter_stratified_null":r4,"R5_crossyear_local_pairs":r5},pd.concat([d1,d2,d3,d4,d5],ignore_index=True)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--discovery",required=True); ap.add_argument("--validation",required=True)
    ap.add_argument("--third",required=True); ap.add_argument("--outdir",required=True)
    a=ap.parse_args()
    dr,dd=run(Path(a.discovery),"discovery")
    vr,vd=run(Path(a.validation),"validation")
    tr,td=run(Path(a.third),"third")
    result={
      "schema":"fcp_distributed_polymorphism_robustness_v1","date_jst":"2026-10-07",
      "status":"complete_posthoc_robustness_audit","confirmatory_decisions_changed":False,
      "radius_km":RADIUS,"permutations":PERMUTATIONS,
      "discovery":dr,"validation":vr,"third":tr,
      "cross_cohort":{
        "different_observer_supported_all_three":bool(all(x["R1_different_observer"].get("supported",False) for x in (dr,vr,tr))),
        "nonwhite_only_supported_all_three":bool(all(x["R2_nonwhite_only"].get("supported",False) for x in (dr,vr,tr))),
        "continuous_nine_colour_supported_all_three":bool(all(x["R3_continuous_nine_colour"].get("supported",False) for x in (dr,vr,tr))),
        "quarter_stratified_supported_all_three":bool(all(x["R4_quarter_stratified_null"].get("supported",False) for x in (dr,vr,tr))),
        "crossyear_supported_all_three":bool(all(x["R5_crossyear_local_pairs"].get("supported",False) for x in (dr,vr,tr)))
      },
      "hard_nonclaims":["does not establish adaptation","does not distinguish genetic from plastic differentiation","post hoc robustness audit"]
    }
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    pd.concat([dd,vd,td],ignore_index=True).to_csv(out/"species_robustness_metrics.csv",index=False)
    (out/"result.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
if __name__=="__main__": main()

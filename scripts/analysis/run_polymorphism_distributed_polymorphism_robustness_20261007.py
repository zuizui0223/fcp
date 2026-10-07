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
    req={"inat_taxon_id","species","photo_id","observer_id","latitude","longitude","morph","global_classifiable"}
    miss=sorted(req-set(d.columns))
    if miss: raise RuntimeError(f"{cohort} missing {miss}")
    d=d.loc[as_bool(d.global_classifiable)&d.morph.astype(str).isin(MORPHS)].copy()
    n=d.groupby("inat_taxon_id").size()
    d=d.loc[d.inat_taxon_id.isin(n[n>=40].index)].copy()
    return d.sort_values(["inat_taxon_id","photo_id"],kind="stable").reset_index(drop=True)

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
    return {"R1_different_observer":r1,"R2_nonwhite_only":r2},pd.concat([d1,d2],ignore_index=True)

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
        "nonwhite_only_supported_all_three":bool(all(x["R2_nonwhite_only"].get("supported",False) for x in (dr,vr,tr)))
      },
      "hard_nonclaims":["does not establish adaptation","does not distinguish genetic from plastic differentiation","post hoc robustness audit"]
    }
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    pd.concat([dd,vd,td],ignore_index=True).to_csv(out/"species_robustness_metrics.csv",index=False)
    (out/"result.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
if __name__=="__main__": main()

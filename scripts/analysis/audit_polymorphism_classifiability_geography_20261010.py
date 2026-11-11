#!/usr/bin/env python3
"""Audit spatial organization of photographic classifiability, not biological morphs."""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

SHA256={
 "discovery":"ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4",
 "validation":"0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6",
 "third":"57630fc9f281bce94a0c40a70aaf7bce879dde93d6154175adcd021e8f5c1186"}
COLOURS=("white","yellow_orange","red_pink","blue_purple")
RADIUS_KM=50.
N_PERM=199
N_BOOT=999
MIN_LOCAL=30
SEED=2026101059

def file_sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""): h.update(chunk)
    return h.hexdigest()

def as_bool(col):
    if pd.api.types.is_bool_dtype(col): return col.fillna(False).astype(bool)
    return col.fillna("").astype(str).str.strip().str.lower().isin(("true","yes","1","y"))

def seed(*args):
    s="|".join(map(str,(SEED,*args))).encode("utf8")
    return int.from_bytes(hashlib.sha256(s).digest()[:8],"little")

def distance_matrix(lat,lon):
    lat=np.deg2rad(np.asarray(lat,float))
    lon=np.deg2rad(np.asarray(lon,float))
    xyz=np.column_stack((np.cos(lat)*np.cos(lon),np.cos(lat)*np.sin(lon),np.sin(lat)))
    return 6371.0088*np.arccos(np.clip(xyz@xyz.T,-1,1))

def load(path,cohort):
    if file_sha(path)!=SHA256[cohort]:
        raise RuntimeError(f"{cohort}: frozen measured photo SHA mismatch")
    d=pd.read_csv(path,low_memory=False)
    required={"inat_taxon_id","species","photo_id","latitude","longitude","observed_on",
              "global_classifiable","morph"}
    if required-set(d): raise RuntimeError(f"Missing: {sorted(required-set(d))}")
    d["_classifiable"]=(as_bool(d.global_classifiable)&d.morph.astype(str).isin(COLOURS)).astype(int)
    counts=d.groupby("inat_taxon_id")["_classifiable"].sum()
    d=d[d.inat_taxon_id.isin(counts[counts>=40].index)].copy()
    if d.photo_id.duplicated().any(): raise RuntimeError("Duplicate immutable photo ID")
    return d.sort_values(["inat_taxon_id","photo_id"],kind="stable")

def measure(g,cohort):
    """Use all photos for technical status; retain original classified positions for colour."""
    g=g.dropna(subset=["latitude","longitude"]).sort_values("photo_id",kind="stable")
    y=g._classifiable.to_numpy(np.int8)
    n=len(g)
    if n<40: raise RuntimeError("Unexpected geography missingness")
    u,v=np.triu_indices(n,1)
    distances=distance_matrix(g.latitude.to_numpy(),g.longitude.to_numpy())
    near=distances[u,v]<=RADIUS_KM
    mu,mv=u[near],v[near]
    row={"cohort":cohort,"inat_taxon_id":int(g.inat_taxon_id.iloc[0]),
         "species":str(g.species.iloc[0]),"n_photos":n,
         "n_classifiable":int(y.sum()),"n_technical_local_pairs":len(mu)}
    nulls={}
    if len(mu)>=MIN_LOCAL:
        glob=float(np.mean(y[u]!=y[v]))
        loc=float(np.mean(y[mu]!=y[mv]))
        row.update(technical_global_discordance=glob,technical_local_discordance=loc,
                   technical_depletion=glob-loc)
        rng=np.random.default_rng(seed(cohort,row["inat_taxon_id"],"technical"))
        perms=np.stack([rng.permutation(y) for _ in range(N_PERM)])
        nulls["unconditional"]=glob-np.mean(perms[:,mu]!=perms[:,mv],axis=1)
        dates=pd.to_datetime(g.observed_on,errors="coerce")
        months=dates.dt.month.fillna(-1).to_numpy(int)
        groups=[np.flatnonzero(months==m) for m in np.unique(months)]
        restricted=np.repeat(y[None,:],N_PERM,axis=0)
        for ix in range(N_PERM):
            for grp in groups: restricted[ix,grp]=rng.permutation(y[grp])
        nulls["month"]=glob-np.mean(restricted[:,mu]!=restricted[:,mv],axis=1)
        row["technical_month_exchangeable_photos"]=int(sum(len(grp) for grp in groups
                                     if len(grp)>1 and len(np.unique(y[grp]))>1))
    else:
        row["technical_depletion"]=np.nan
    colour=g[g._classifiable.eq(1)]
    if len(colour)>=40:
        u2,v2=np.triu_indices(len(colour),1)
        near2=distance_matrix(colour.latitude.to_numpy(),colour.longitude.to_numpy())[u2,v2]<=RADIUS_KM
        row["n_colour_local_pairs"]=int(near2.sum())
        if near2.sum()>=MIN_LOCAL:
            codes=pd.Categorical(colour.morph,categories=COLOURS).codes
            row["colour_depletion"]=float(np.mean(codes[u2]!=codes[v2])-
                                      np.mean(codes[u2[near2]]!=codes[v2[near2]]))
        else:
            row["colour_depletion"]=np.nan
    else:
        row["n_colour_local_pairs"]=0
        row["colour_depletion"]=np.nan
    return row,nulls

def summarize(rows,nulls,cohort):
    valid=np.isfinite(rows.technical_depletion.to_numpy(float))
    out={"n_high_depth_eligible_species":len(rows),"n_technical_evaluable_species":int(valid.sum()),
         "n_colour_evaluable_species":int(np.isfinite(rows.colour_depletion).sum())}
    if not valid.any(): return {**out,"status":"HOLD_NO_TECHNICAL_LOCAL_SUPPORT"}
    subset=rows.loc[valid]
    obs=subset.technical_depletion.to_numpy(float)
    out["mean_technical_depletion"]=float(obs.mean())
    out["mean_technical_global_discordance"]=float(subset.technical_global_discordance.mean())
    out["mean_technical_local_discordance"]=float(subset.technical_local_discordance.mean())
    rng=np.random.default_rng(seed(cohort,"bootstrap"))
    boots=np.mean(obs[rng.integers(len(obs),size=(N_BOOT,len(obs)))],axis=1)
    out["species_bootstrap_ci_95"]=[float(x) for x in np.quantile(boots,[.025,.975])]
    for mode in ("unconditional","month"):
        mat=np.vstack([null[mode] for null,ok in zip(nulls,valid) if ok])
        avg=mat.mean(axis=0)
        out[f"{mode}_null_mean"]=float(avg.mean())
        out[f"{mode}_p_upper"]=float((1+np.count_nonzero(avg>=obs.mean()))/(N_PERM+1))
    matched=rows.loc[valid & np.isfinite(rows.colour_depletion)].copy()
    out["n_matched_technical_colour_species"]=len(matched)
    rho=None
    if len(matched)>=30 and matched.technical_depletion.nunique()>1 and matched.colour_depletion.nunique()>1:
        rho=float(spearmanr(matched.technical_depletion,matched.colour_depletion).statistic)
    out["technical_vs_colour_depletion_spearman_rho_descriptive"]=rho
    out["status"]="POSTHOC_OBSERVATION_PROCESS_DIAGNOSTIC"
    return out

def execute(paths,outdir):
    result={"schema":"fcp_classifiability_geography_posthoc_v1","date_jst":"2026-10-10",
            "radius_km":RADIUS_KM,"permutations":N_PERM,"source_sha256":SHA256,
            "confirmatory_decisions_changed":False,
            "nonclaims":["Technical classifiability is not a validated biological negative control.",
                         "Spatial classifiability alone does not explain colour-specific clustering.",
                         "Unmeasured morph-dependent missingness, site-specific photo conditions and phenotype plasticity remain.",
                         "No genetic morph, fitness, selection or local adaptation identification."]}
    allrows=[]
    for cohort,path in paths.items():
        frame=load(path,cohort)
        rows=[];nulls=[]
        for _,g in frame.groupby("inat_taxon_id",sort=True):
            row,null=measure(g,cohort)
            rows.append(row);nulls.append(null)
        df=pd.DataFrame(rows)
        result[cohort]=summarize(df,nulls,cohort)
        allrows.append(df)
    outdir=Path(outdir)
    outdir.mkdir(parents=True,exist_ok=True)
    pd.concat(allrows,ignore_index=True).to_csv(outdir/"species_classifiability_spatial_audit.csv",index=False)
    (outdir/"result.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf8")
    return result

def main():
    p=argparse.ArgumentParser()
    for c in SHA256: p.add_argument("--"+c,type=Path,required=True)
    p.add_argument("--outdir",type=Path,required=True)
    args=p.parse_args()
    result=execute({c:getattr(args,c) for c in SHA256},args.outdir)
    print(json.dumps({c:result[c] for c in SHA256},indent=2))

if __name__=="__main__": main()

#!/usr/bin/env python3
"""Post hoc test of geographically distributed flower-colour polymorphism."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

EARTH_RADIUS_KM = 6371.0088
MORPHS = ["white", "yellow_orange", "red_pink", "blue_purple"]
MORPH_CODE = {m:i for i,m in enumerate(MORPHS)}
RADII_KM = [25.0, 50.0, 100.0, 250.0]
PRIMARY_RADIUS = 50.0
MIN_LOCAL_PAIRS = 30
PERMUTATIONS = 199
MASTER_SEED = 2026100731
SHA = {
    "discovery":"ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4",
    "validation":"0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6",
    "third":"57630fc9f281bce94a0c40a70aaf7bce879dde93d6154175adcd021e8f5c1186",
}


def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda:fh.read(1<<20), b""):
            h.update(chunk)
    return h.hexdigest()


def as_bool(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False).astype(bool)
    return s.fillna("").astype(str).str.strip().str.casefold().isin({"true","1","yes","y"})


def stable_seed(cohort: str, taxon: int) -> int:
    b=f"{MASTER_SEED}|{cohort}|{taxon}".encode()
    return int.from_bytes(hashlib.sha256(b).digest()[:8],"little")


def pairwise_geo_km(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    latr=np.deg2rad(np.asarray(lat,float))
    lonr=np.deg2rad(np.asarray(lon,float))
    c=np.cos(latr)
    xyz=np.column_stack([c*np.cos(lonr),c*np.sin(lonr),np.sin(latr)])
    dot=np.clip(xyz@xyz.T,-1.0,1.0)
    return np.arccos(dot)*EARTH_RADIUS_KM


def load(path: Path, cohort: str) -> pd.DataFrame:
    if sha256(path)!=SHA[cohort]:
        raise RuntimeError(f"{cohort} SHA256 mismatch")
    d=pd.read_csv(path,low_memory=False)
    required={"inat_taxon_id","species","photo_id","latitude","longitude","morph","global_classifiable"}
    miss=sorted(required-set(d.columns))
    if miss:
        raise RuntimeError(f"{cohort} missing columns {miss}")
    keep=as_bool(d["global_classifiable"]) & d["morph"].astype(str).isin(MORPHS)
    d=d.loc[keep].copy()
    n=d.groupby("inat_taxon_id").size()
    d=d.loc[d["inat_taxon_id"].isin(n[n>=40].index)].copy()
    return d.sort_values(["inat_taxon_id","photo_id"],kind="stable").reset_index(drop=True)


def analyse_species(g: pd.DataFrame, cohort: str) -> tuple[dict,dict[float,np.ndarray]]:
    g=g.dropna(subset=["latitude","longitude"]).copy()
    labels=np.array([MORPH_CODE[str(x)] for x in g["morph"]],dtype=np.int8)
    n=len(g)
    if n<40:
        raise RuntimeError("eligible species lost below 40 rows after coordinate check")
    geo=pairwise_geo_km(g["latitude"].to_numpy(float),g["longitude"].to_numpy(float))
    u,v=np.triu_indices(n,k=1)
    dist=geo[u,v]
    different=(labels[u]!=labels[v])
    D_pair=float(np.mean(different))

    rng=np.random.default_rng(stable_seed(cohort,int(g["inat_taxon_id"].iloc[0])))
    perms=np.stack([rng.permutation(labels) for _ in range(PERMUTATIONS)])

    row={
        "cohort":cohort,
        "inat_taxon_id":int(g["inat_taxon_id"].iloc[0]),
        "species":str(g["species"].iloc[0]),
        "n_classifiable":int(n),
        "D_pair":D_pair,
        "sampled_span_km":float(np.max(dist)),
    }
    nulls={}
    for radius in RADII_KM:
        mask=dist<=radius
        m=int(mask.sum())
        key=str(int(radius))
        row[f"local_pairs_{key}km"]=m
        if m<MIN_LOCAL_PAIRS:
            row[f"D_local_{key}km"]=np.nan
            row[f"depletion_{key}km"]=np.nan
            nulls[radius]=np.full(PERMUTATIONS,np.nan)
            continue
        uu=u[mask]; vv=v[mask]
        obs_local=float(np.mean(labels[uu]!=labels[vv]))
        dep=D_pair-obs_local
        # Each row in perms is one complete vertex permutation.
        local_null=np.mean(perms[:,uu]!=perms[:,vv],axis=1)
        dep_null=D_pair-local_null
        row[f"D_local_{key}km"]=obs_local
        row[f"depletion_{key}km"]=dep
        row[f"null_mean_depletion_{key}km"]=float(np.mean(dep_null))
        nulls[radius]=dep_null
    return row,nulls


def cohort_summary(species: pd.DataFrame, null_by_radius: dict[float,list[np.ndarray]]) -> dict:
    out={}
    for radius in RADII_KM:
        key=str(int(radius))
        vals=species[f"depletion_{key}km"].to_numpy(float)
        good=np.isfinite(vals)
        arrays=[a for a,ok in zip(null_by_radius[radius],good) if ok and np.isfinite(a).all()]
        if int(good.sum())==0 or not arrays:
            out[key]={"evaluable":False,"n_species":0}
            continue
        mat=np.vstack(arrays)
        mean_null=mat.mean(axis=0)
        obs_mean=float(np.mean(vals[good]))
        p_mean=float((1+np.count_nonzero(mean_null>=obs_mean))/(PERMUTATIONS+1))

        D=species.loc[good,"D_pair"].to_numpy(float)
        dep=vals[good]
        rho_obs=float(spearmanr(D,dep).statistic) if np.ptp(D)>1e-15 and np.ptp(dep)>1e-15 else 0.0
        rho_null=np.empty(PERMUTATIONS,float)
        for i in range(PERMUTATIONS):
            z=mat[:,i]
            rho_null[i]=float(spearmanr(D,z).statistic) if np.ptp(z)>1e-15 else 0.0
        p_rho=float((1+np.count_nonzero(rho_null>=rho_obs))/(PERMUTATIONS+1))

        out[key]={
            "evaluable":True,
            "radius_km":radius,
            "n_species":int(good.sum()),
            "mean_D_pair":float(np.mean(D)),
            "mean_D_local":float(np.mean(D-dep)),
            "mean_depletion":obs_mean,
            "positive_depletion_species_fraction":float(np.mean(dep>0)),
            "mean_depletion_null_mean":float(np.mean(mean_null)),
            "mean_depletion_null_q025":float(np.quantile(mean_null,.025)),
            "mean_depletion_null_q975":float(np.quantile(mean_null,.975)),
            "mean_depletion_p_upper":p_mean,
            "mean_depletion_supported":bool(obs_mean>0 and p_mean<.05),
            "rho_D_vs_depletion":rho_obs,
            "rho_D_vs_depletion_null_mean":float(np.mean(rho_null)),
            "rho_D_vs_depletion_null_q025":float(np.quantile(rho_null,.025)),
            "rho_D_vs_depletion_null_q975":float(np.quantile(rho_null,.975)),
            "rho_D_vs_depletion_p_upper":p_rho,
            "rho_D_vs_depletion_supported":bool(rho_obs>0 and p_rho<.05),
        }
    return out


def run(path: Path, cohort: str) -> tuple[dict,pd.DataFrame]:
    d=load(path,cohort)
    rows=[]
    null_by_radius={r:[] for r in RADII_KM}
    for _,g in d.groupby("inat_taxon_id",sort=True):
        row,nulls=analyse_species(g,cohort)
        rows.append(row)
        for r in RADII_KM:
            null_by_radius[r].append(nulls[r])
    sp=pd.DataFrame(rows).sort_values("inat_taxon_id").reset_index(drop=True)
    return {
        "eligible_species":int(len(sp)),
        "radii":cohort_summary(sp,null_by_radius),
    },sp


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--discovery",type=Path,required=True)
    ap.add_argument("--validation",type=Path,required=True)
    ap.add_argument("--third",type=Path,required=True)
    ap.add_argument("--outdir",type=Path,required=True)
    a=ap.parse_args()

    disc,ds=run(a.discovery,"discovery")
    val,vs=run(a.validation,"validation")
    third,ts=run(a.third,"third")
    pk=str(int(PRIMARY_RADIUS))

    result={
        "schema":"fcp_distributed_polymorphism_posthoc_v1",
        "date_jst":"2026-10-07",
        "status":"complete_posthoc_distributed_polymorphism_test",
        "confirmatory_decisions_changed":False,
        "primary_radius_km":PRIMARY_RADIUS,
        "minimum_local_pairs":MIN_LOCAL_PAIRS,
        "permutations":PERMUTATIONS,
        "discovery":disc,
        "validation":val,
        "third":third,
        "replication":{
            "local_depletion_500_plus_500":bool(
                disc["radii"][pk].get("mean_depletion_supported",False)
                and val["radii"][pk].get("mean_depletion_supported",False)
            ),
            "D_depletion_coupling_500_plus_500":bool(
                disc["radii"][pk].get("rho_D_vs_depletion_supported",False)
                and val["radii"][pk].get("rho_D_vs_depletion_supported",False)
            ),
            "local_depletion_all_three":bool(
                all(x["radii"][pk].get("mean_depletion_supported",False) for x in (disc,val,third))
            ),
            "D_depletion_coupling_all_three":bool(
                all(x["radii"][pk].get("rho_D_vs_depletion_supported",False) for x in (disc,val,third))
            ),
        },
        "interpretation":{
            "if_local_depletion":"Nearby conspecific observations contain less coarse colour-state diversity than expected from the species-wide colour composition, consistent with geographic partitioning of ITV.",
            "if_D_coupling":"Species with greater species-wide colour diversity show stronger local depletion even relative to null worlds that preserve each species' D exactly.",
            "hard_nonclaims":[
                "does not distinguish selection from drift or dispersal limitation",
                "does not identify genetic versus plastic differentiation",
                "does not establish local adaptation or fitness differences",
                "post hoc analysis cannot alter frozen confirmatory decisions"
            ]
        }
    }
    out=a.outdir; out.mkdir(parents=True,exist_ok=True)
    pd.concat([ds,vs,ts],ignore_index=True).to_csv(out/"species_distributed_polymorphism_metrics.csv",index=False)
    (out/"result.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2))
    return 0


if __name__=="__main__":
    raise SystemExit(main())

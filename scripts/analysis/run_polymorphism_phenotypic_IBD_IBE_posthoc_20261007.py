#!/usr/bin/env python3
"""Post hoc phenotypic IBD versus IBE-like decomposition for FCP."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from scipy.stats import rankdata, spearmanr

EARTH_RADIUS_KM=6371.0088
MORPHS=["white","yellow_orange","red_pink","blue_purple"]
BIO9=["white","yellow","orange","red","pink","magenta","purple","blue","bronze"]
LEGACY_COLS=[f"palette_count_{x}" for x in BIO9]
THIRD_COLS=[f"flower_fraction_{x}" for x in BIO9]
PERMUTATIONS=199
SIGNFLIP=9999
MASTER_SEED=2026100737
SHA={
    "discovery":"ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4",
    "validation":"0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6",
    "third":"57630fc9f281bce94a0c40a70aaf7bce879dde93d6154175adcd021e8f5c1186",
}


def sha256(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def stable_seed(*parts:object)->int:
    b="|".join(map(str,(MASTER_SEED,*parts))).encode()
    return int.from_bytes(hashlib.sha256(b).digest()[:8],"little")


def as_bool(s:pd.Series)->pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False).astype(bool)
    return s.fillna("").astype(str).str.strip().str.casefold().isin({"true","1","yes","y"})


def sample_raster(path:Path,lon:np.ndarray,lat:np.ndarray)->np.ndarray:
    out=np.full(len(lon),np.nan,float)
    good=[i for i,(x,y) in enumerate(zip(lon,lat)) if np.isfinite(x) and np.isfinite(y)]
    if not good:
        return out
    with rasterio.open(path) as src:
        vals=list(src.sample([(float(lon[i]),float(lat[i])) for i in good]))
        for i,v in zip(good,vals):
            x=float(v[0])
            if src.nodata is not None and np.isclose(x,src.nodata):
                x=np.nan
            out[i]=x
    return out


def pairwise_geo(lat:np.ndarray,lon:np.ndarray)->np.ndarray:
    latr=np.deg2rad(np.asarray(lat,float)); lonr=np.deg2rad(np.asarray(lon,float))
    c=np.cos(latr)
    xyz=np.column_stack([c*np.cos(lonr),c*np.sin(lonr),np.sin(latr)])
    dot=np.clip(xyz@xyz.T,-1,1)
    return np.arccos(dot)*EARTH_RADIUS_KM


def pairwise_jsd(prob:np.ndarray)->np.ndarray:
    p=np.asarray(prob,float)
    mass=p.sum(axis=1)
    if p.ndim!=2 or np.any(~np.isfinite(p)) or np.any(p<0) or np.any(mass<=0):
        raise ValueError("invalid colour composition")
    p=p/mass[:,None]
    a=p[:,None,:]; b=p[None,:,:]; m=.5*(a+b)
    with np.errstate(divide="ignore",invalid="ignore"):
        ka=np.where(a>0,a*np.log2(a/m),0).sum(axis=2)
        kb=np.where(b>0,b*np.log2(b/m),0).sum(axis=2)
    return np.clip(.5*(ka+kb),0,1)


def partial_from_centered(yc:np.ndarray,xc:np.ndarray,zc:np.ndarray)->float:
    z2=float(np.dot(zc,zc))
    yr=yc-(0.0 if z2<=1e-15 else float(np.dot(zc,yc)/z2))*zc
    xr=xc-(0.0 if z2<=1e-15 else float(np.dot(zc,xc)/z2))*zc
    den=float(np.linalg.norm(yr)*np.linalg.norm(xr))
    return 0.0 if den<=1e-14 else float(np.dot(yr,xr)/den)


def r2(y:np.ndarray,X:np.ndarray)->float:
    y=np.asarray(y,float)
    X=np.asarray(X,float)
    A=np.column_stack([np.ones(len(y)),X])
    b,*_=np.linalg.lstsq(A,y,rcond=None)
    fit=A@b
    ss=float(np.dot(y-y.mean(),y-y.mean()))
    if ss<=1e-15:
        return 0.0
    return max(0.0,min(1.0,1-float(np.dot(y-fit,y-fit))/ss))


def prepare(path:Path,cohort:str,bio5:Path)->pd.DataFrame:
    if sha256(path)!=SHA[cohort]:
        raise RuntimeError(f"{cohort} SHA256 mismatch")
    d=pd.read_csv(path,low_memory=False)
    cols=THIRD_COLS if cohort=="third" else LEGACY_COLS
    req={"inat_taxon_id","species","photo_id","latitude","longitude","morph","global_classifiable",*cols}
    miss=sorted(req-set(d.columns))
    if miss:
        raise RuntimeError(f"{cohort} missing {miss}")
    keep=as_bool(d["global_classifiable"]) & d["morph"].astype(str).isin(MORPHS)
    d=d.loc[keep].copy()
    counts=d.groupby("inat_taxon_id").size()
    d=d.loc[d["inat_taxon_id"].isin(counts[counts>=40].index)].copy()
    d=d.sort_values(["inat_taxon_id","photo_id"],kind="stable").reset_index(drop=True)
    d["bio5"]=sample_raster(
        bio5,
        pd.to_numeric(d["longitude"],errors="coerce").to_numpy(float),
        pd.to_numeric(d["latitude"],errors="coerce").to_numpy(float),
    )
    return d


def species_stats(g:pd.DataFrame,cohort:str)->tuple[dict,np.ndarray,np.ndarray]:
    cols=THIRD_COLS if cohort=="third" else LEGACY_COLS
    g=g.dropna(subset=["latitude","longitude","bio5",*cols]).copy()
    if len(g)<40:
        raise RuntimeError("species fell below eligibility after BIO5 sampling")
    g=g.sort_values("photo_id",kind="stable").reset_index(drop=True)
    n=len(g)
    u,v=np.triu_indices(n,k=1)
    geo=pairwise_geo(g["latitude"].to_numpy(float),g["longitude"].to_numpy(float))
    gp=geo[u,v]
    bio=g["bio5"].to_numpy(float)
    ep=np.abs(bio[u]-bio[v])
    colour=pairwise_jsd(g[cols].to_numpy(float))
    yp=colour[u,v]

    gr=rankdata(gp,method="average"); er=rankdata(ep,method="average"); yr=rankdata(yp,method="average")
    gc=gr-gr.mean(); ec=er-er.mean(); yc=yr-yr.mean()

    ibe=partial_from_centered(yc,ec,gc)
    ibd=partial_from_centered(yc,gc,ec)

    Rg=r2(yr,gr[:,None])
    Re=r2(yr,er[:,None])
    Rge=r2(yr,np.column_stack([gr,er]))
    unique_g=Rge-Re
    unique_e=Rge-Rg
    shared=Rg+Re-Rge

    rank_matrix=np.zeros((n,n),float)
    rank_matrix[u,v]=yr
    rank_matrix[v,u]=yr
    rng=np.random.default_rng(stable_seed(cohort,int(g["inat_taxon_id"].iloc[0])))
    ibe_null=np.empty(PERMUTATIONS,float)
    ibd_null=np.empty(PERMUTATIONS,float)
    batch=32
    for start in range(0,PERMUTATIONS,batch):
        stop=min(PERMUTATIONS,start+batch)
        pp=np.stack([rng.permutation(n) for _ in range(stop-start)])
        vals=rank_matrix[pp[:,u],pp[:,v]]
        vals=vals-vals.mean(axis=1,keepdims=True)

        # IBE: y residual after geography, against env residual after geography.
        g2=float(np.dot(gc,gc))
        e_res=ec-(0.0 if g2<=1e-15 else float(np.dot(gc,ec)/g2))*gc
        en=float(np.linalg.norm(e_res))
        beta_g=np.zeros(stop-start) if g2<=1e-15 else (vals@gc)/g2
        y_rg=vals-beta_g[:,None]*gc
        yn=np.linalg.norm(y_rg,axis=1)
        ibe_null[start:stop]=np.where((yn<=1e-14)|(en<=1e-14),0.0,(y_rg@e_res)/(yn*en))

        # IBD: y residual after environment, against geography residual after environment.
        e2=float(np.dot(ec,ec))
        g_res=gc-(0.0 if e2<=1e-15 else float(np.dot(ec,gc)/e2))*ec
        gn=float(np.linalg.norm(g_res))
        beta_e=np.zeros(stop-start) if e2<=1e-15 else (vals@ec)/e2
        y_re=vals-beta_e[:,None]*ec
        yn2=np.linalg.norm(y_re,axis=1)
        ibd_null[start:stop]=np.where((yn2<=1e-14)|(gn<=1e-14),0.0,(y_re@g_res)/(yn2*gn))

    row={
        "cohort":cohort,
        "inat_taxon_id":int(g["inat_taxon_id"].iloc[0]),
        "species":str(g["species"].iloc[0]),
        "n":int(n),
        "rho_IBE":ibe,
        "rho_IBD":ibd,
        "delta_IBE_minus_IBD":ibe-ibd,
        "R2_geo":Rg,
        "R2_env":Re,
        "R2_both":Rge,
        "unique_geo":unique_g,
        "unique_env":unique_e,
        "shared_geo_env":shared,
        "rho_geo_env":float(spearmanr(gp,ep).statistic),
    }
    return row,ibe_null,ibd_null


def summarize(rows:pd.DataFrame,ibe_nulls:list[np.ndarray],ibd_nulls:list[np.ndarray],cohort:str)->dict:
    ibe=rows["rho_IBE"].to_numpy(float)
    ibd=rows["rho_IBD"].to_numpy(float)
    mat_e=np.vstack(ibe_nulls)
    mat_g=np.vstack(ibd_nulls)
    null_e=mat_e.mean(axis=0)
    null_g=mat_g.mean(axis=0)
    obs_e=float(np.mean(ibe)); obs_g=float(np.mean(ibd))
    p_e=float((1+np.count_nonzero(null_e>=obs_e))/(PERMUTATIONS+1))
    p_g=float((1+np.count_nonzero(null_g>=obs_g))/(PERMUTATIONS+1))

    delta=ibe-ibd
    obs_delta=float(np.mean(delta))
    null_delta=(mat_e-mat_g).mean(axis=0)
    p_delta_null=float((1+np.count_nonzero(null_delta>=obs_delta))/(PERMUTATIONS+1))

    rng=np.random.default_rng(stable_seed("signflip",cohort))
    sf=np.empty(SIGNFLIP,float)
    for i in range(SIGNFLIP):
        signs=rng.choice(np.array([-1.0,1.0]),size=len(delta))
        sf[i]=float(np.mean(delta*signs))
    p_delta_sign=float((1+np.count_nonzero(sf>=obs_delta))/(SIGNFLIP+1))

    return {
        "n_species":int(len(rows)),
        "IBE_like":{
            "mean_partial_rho":obs_e,
            "median_partial_rho":float(np.median(ibe)),
            "positive_species_fraction":float(np.mean(ibe>0)),
            "matched_vertex_p_upper":p_e,
            "null_mean":float(np.mean(null_e)),
            "null_q025":float(np.quantile(null_e,.025)),
            "null_q975":float(np.quantile(null_e,.975)),
            "supported":bool(obs_e>0 and p_e<.05),
        },
        "IBD_like":{
            "mean_partial_rho":obs_g,
            "median_partial_rho":float(np.median(ibd)),
            "positive_species_fraction":float(np.mean(ibd>0)),
            "matched_vertex_p_upper":p_g,
            "null_mean":float(np.mean(null_g)),
            "null_q025":float(np.quantile(null_g,.025)),
            "null_q975":float(np.quantile(null_g,.975)),
            "supported":bool(obs_g>0 and p_g<.05),
        },
        "relative_balance":{
            "mean_IBE_minus_IBD":obs_delta,
            "median_IBE_minus_IBD":float(np.median(delta)),
            "fraction_IBE_gt_IBD":float(np.mean(delta>0)),
            "matched_null_p_upper_for_positive_delta":p_delta_null,
            "species_signflip_p_upper_for_positive_delta":p_delta_sign,
            "IBE_stronger_supported":bool(obs_delta>0 and p_delta_null<.05 and p_delta_sign<.05),
        },
        "commonality_descriptive":{
            "mean_unique_geo":float(rows["unique_geo"].mean()),
            "mean_unique_env":float(rows["unique_env"].mean()),
            "mean_shared_geo_env":float(rows["shared_geo_env"].mean()),
            "median_unique_geo":float(rows["unique_geo"].median()),
            "median_unique_env":float(rows["unique_env"].median()),
            "mean_pairwise_geo_env_rho":float(rows["rho_geo_env"].mean()),
        },
    }


def run(path:Path,cohort:str,bio5:Path)->tuple[dict,pd.DataFrame]:
    d=prepare(path,cohort,bio5)
    rows=[]; en=[]; gn=[]
    for _,g in d.groupby("inat_taxon_id",sort=True):
        row,a,b=species_stats(g,cohort)
        rows.append(row); en.append(a); gn.append(b)
    df=pd.DataFrame(rows).sort_values("inat_taxon_id").reset_index(drop=True)
    return summarize(df,en,gn,cohort),df


def main()->int:
    ap=argparse.ArgumentParser()
    ap.add_argument("--discovery",type=Path,required=True)
    ap.add_argument("--validation",type=Path,required=True)
    ap.add_argument("--third",type=Path,required=True)
    ap.add_argument("--bio5",type=Path,required=True)
    ap.add_argument("--outdir",type=Path,required=True)
    a=ap.parse_args()
    disc,dd=run(a.discovery,"discovery",a.bio5)
    val,vd=run(a.validation,"validation",a.bio5)
    third,td=run(a.third,"third",a.bio5)

    result={
        "schema":"fcp_phenotypic_IBD_IBE_posthoc_v1",
        "date_jst":"2026-10-07",
        "status":"complete_posthoc_phenotypic_IBD_IBE_decomposition",
        "confirmatory_decisions_changed":False,
        "environment":"WorldClim 2.1 BIO5",
        "colour_representation":"continuous nine-colour biological palette JSD",
        "permutations":PERMUTATIONS,
        "discovery":disc,
        "validation":val,
        "third":third,
        "replication":{
            "IBE_like_500_plus_500":bool(disc["IBE_like"]["supported"] and val["IBE_like"]["supported"]),
            "IBE_like_all_three":bool(all(x["IBE_like"]["supported"] for x in (disc,val,third))),
            "IBD_like_500_plus_500":bool(disc["IBD_like"]["supported"] and val["IBD_like"]["supported"]),
            "IBD_like_all_three":bool(all(x["IBD_like"]["supported"] for x in (disc,val,third))),
            "IBE_stronger_than_IBD_500_plus_500":bool(
                disc["relative_balance"]["IBE_stronger_supported"] and val["relative_balance"]["IBE_stronger_supported"]
            ),
            "IBE_stronger_than_IBD_all_three":bool(
                all(x["relative_balance"]["IBE_stronger_supported"] for x in (disc,val,third))
            ),
        },
        "interpretation":{
            "IBE_like":"Positive BIO5-associated colour turnover after geographic distance is removed is consistent with phenotypic environmental sorting.",
            "IBD_like":"Positive geographic colour turnover after BIO5 difference is removed is consistent with nonredundant spatial/dispersal/history structure.",
            "hard_nonclaims":[
                "not genetic isolation by environment",
                "does not establish causal temperature selection",
                "does not measure fitness or local adaptation",
                "cannot distinguish genetic differentiation from plasticity",
                "post hoc analysis cannot alter frozen decisions"
            ]
        }
    }
    out=a.outdir; out.mkdir(parents=True,exist_ok=True)
    pd.concat([dd,vd,td],ignore_index=True).to_csv(out/"species_IBD_IBE_metrics.csv",index=False)
    (out/"result.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2))
    return 0


if __name__=="__main__":
    raise SystemExit(main())

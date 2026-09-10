#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import rasterio
from numpy.linalg import pinv, slogdet
from scipy.spatial import ConvexHull, QhullError
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[2]
DISC = ROOT / "data/derived/global_monte_carlo_measured_photos_v1.csv"
RES = ROOT / "data/derived/rgfca_reserve_replication_measured_photos_v1.csv"
STATES = ROOT / "results/rgfca_local_cooccurrence_segregation_step7a_20260910/photo_organization_states.csv"
OUT = ROOT / "results/rgfca_worldclim_old34_transfer_step7d_20260910"
OUT.mkdir(parents=True, exist_ok=True)
SELECT = [1, 4, 5, 6, 7, 12, 14, 15, 17]
CLIMATE = [f"bio{x}" for x in SELECT]
N_PERM = 20_000
HV_DRAWS = 99
SEED = 20260910
R_EARTH_KM = 6371.0088
PRED = {
    "H1a_mean_bio4": ("mean_bio4", "less"),
    "H1b_mean_bio15": ("mean_bio15", "less"),
    "H3a_spatial_niche_turnover": ("spatial_niche_turnover", "greater"),
    "H3b_regional_centroid_sep": ("regional_centroid_sep", "greater"),
    "H3c_regional_gaussian_overlap": ("regional_gaussian_overlap", "less"),
}


def haversine(lat1, lon1, lat2, lon2):
    p1 = np.radians(lat1); p2 = np.radians(lat2)
    dp = np.radians(lat2-lat1); dl = np.radians(lon2-lon1)
    a = np.sin(dp/2)**2 + np.cos(p1)*np.cos(p2)*np.sin(dl/2)**2
    return 2*R_EARTH_KM*np.arcsin(np.sqrt(np.clip(a,0,1)))


def spherical_xyz(lat, lon):
    lat=np.radians(lat); lon=np.radians(lon)
    return np.c_[np.cos(lat)*np.cos(lon), np.cos(lat)*np.sin(lon), np.sin(lat)]


def bhattacharyya_overlap(x1: np.ndarray, x2: np.ndarray, eps: float=1e-4) -> float:
    d=x1.shape[1]; mu1,mu2=x1.mean(0),x2.mean(0)
    s1=np.cov(x1,rowvar=False)+eps*np.eye(d); s2=np.cov(x2,rowvar=False)+eps*np.eye(d); s=(s1+s2)/2
    delta=mu1-mu2; sg,ld=slogdet(s); sg1,ld1=slogdet(s1); sg2,ld2=slogdet(s2)
    if min(sg,sg1,sg2)<=0: return np.nan
    db=.125*delta@pinv(s)@delta + .5*(ld-.5*(ld1+ld2))
    return float(np.exp(-db))


def extract_climate(path: Path, tranche: str, wcdir: Path, admitted: set[str]) -> tuple[pd.DataFrame, dict[str, Any]]:
    d=pd.read_csv(path, usecols=["species","latitude","longitude"])
    d=d[d.species.isin(admitted)].dropna(subset=["latitude","longitude"]).copy()
    coords=list(zip(d.longitude.astype(float), d.latitude.astype(float)))
    for bio in SELECT:
        candidates=list(wcdir.glob(f"*bio_{bio}.tif"))+list(wcdir.glob(f"*bio{bio}.tif"))
        if not candidates: raise FileNotFoundError(f"BIO{bio} raster missing")
        with rasterio.open(candidates[0]) as src:
            vals=[]
            for v in src.sample(coords):
                x=v[0] if len(v) else np.nan
                vals.append(np.nan if src.nodata is not None and x==src.nodata else x)
            d[f"bio{bio}"]=vals
    before=len(d); d=d.dropna(subset=CLIMATE).copy(); after_climate=len(d)
    d=d.drop_duplicates(["species",*CLIMATE]).copy(); dedup=len(d)
    scaler=StandardScaler(); Z=scaler.fit_transform(d[CLIMATE]); pca=PCA(n_components=3, random_state=SEED).fit(Z); pcs=pca.transform(Z)
    for j in range(3): d[f"pc{j+1}"]=pcs[:,j]
    d["tranche"]=tranche
    qc={"raw_coordinate_rows":int(before),"rows_with_climate":int(after_climate),"unique_species_climate_cells":int(dedup),"species":int(d.species.nunique()),"pca_variance_explained":pca.explained_variance_ratio_.tolist()}
    return d,qc


def species_metrics(cells: pd.DataFrame, states: pd.DataFrame, tranche: str, seed: int) -> pd.DataFrame:
    st=states[states.tranche.eq(tranche)].set_index("species")
    rng_pairs=np.random.default_rng(seed); rng_hv=np.random.default_rng(seed+100)
    rows=[]
    for sp,g in cells.groupby("species",sort=True):
        state=st.loc[sp]
        n=len(g)
        row={"tranche":tranche,"species":sp,"organization_state":state.organization_state,"C_star":bool(state.C_star),"S_star":bool(state.S_star),"D":float(state.D),"n_climate_cells":int(n)}
        if n<20:
            row["metric_status"]="insufficient_cells"; rows.append(row); continue
        lat=g.latitude.to_numpy(float); lon=g.longitude.to_numpy(float); x=g[["pc1","pc2","pc3"]].to_numpy(float)
        row["mean_bio4"]=float(g.bio4.mean()); row["mean_bio15"]=float(g.bio15.mean())
        total=n*(n-1)//2
        if total<=20000: ii,jj=np.triu_indices(n,1)
        else:
            ii=rng_pairs.integers(0,n,size=25000); jj=rng_pairs.integers(0,n,size=25000); ok=ii!=jj; ii=ii[ok]; jj=jj[ok]
        gd=haversine(lat[ii],lon[ii],lat[jj],lon[jj]); ed=np.linalg.norm(x[ii]-x[jj],axis=1)
        q25,q75=np.quantile(gd,[.25,.75]); short=np.median(ed[gd<=q25]); long=np.median(ed[gd>=q75]); med=np.median(ed)
        row["spatial_niche_turnover"]=float((long-short)/(med+1e-9))
        vols=[]
        for _ in range(HV_DRAWS):
            idx=rng_hv.choice(n,20,replace=False)
            try: vols.append(float(ConvexHull(x[idx]).volume))
            except QhullError: pass
        row["hv3d_rarefied20"]=float(np.median(vols)) if vols else np.nan
        row["regional_centroid_sep"]=np.nan; row["regional_gaussian_overlap"]=np.nan
        if n>=40:
            lab=KMeans(n_clusters=2,random_state=seed,n_init=20).fit_predict(spherical_xyz(lat,lon)); counts=np.bincount(lab)
            if len(counts)==2 and counts.min()>=15:
                x1,x2=x[lab==0],x[lab==1]; centroid=np.linalg.norm(x1.mean(0)-x2.mean(0)); rms=np.sqrt(np.mean(np.sum((x-x.mean(0))**2,axis=1)))
                row["regional_centroid_sep"]=float(centroid/(rms+1e-9)); row["regional_gaussian_overlap"]=bhattacharyya_overlap(x1,x2)
        row["metric_status"]="complete"; rows.append(row)
    return pd.DataFrame(rows)


def perm_contrast(values: np.ndarray, labels_s: np.ndarray, alt: str, seed: int) -> dict[str,Any]:
    v=np.asarray(values,float); s=np.asarray(labels_s,bool); keep=np.isfinite(v); v=v[keep]; s=s[keep]
    if s.sum()<1 or (~s).sum()<1: return {"n":int(len(v)),"n_S":int(s.sum()),"n_C":int((~s).sum()),"delta_S_minus_C":None,"p":1.0}
    obs=float(np.median(v[s])-np.median(v[~s])); rng=np.random.default_rng(seed); null=np.empty(N_PERM,float)
    for b in range(N_PERM):
        q=s[rng.permutation(len(s))]; null[b]=float(np.median(v[q])-np.median(v[~q]))
    if alt=="greater": p=float((1+np.sum(null>=obs-1e-15))/(N_PERM+1))
    elif alt=="less": p=float((1+np.sum(null<=obs+1e-15))/(N_PERM+1))
    else: p=float((1+np.sum(np.abs(null)>=abs(obs)-1e-15))/(N_PERM+1))
    return {"n":int(len(v)),"n_S":int(s.sum()),"n_C":int((~s).sum()),"delta_S_minus_C":obs,"p":p,"null_q025":float(np.quantile(null,.025)),"null_q975":float(np.quantile(null,.975))}


def holm(pvals: dict[str,float]) -> dict[str,float]:
    keys=list(pvals); p=np.array([pvals[k] for k in keys],float); order=np.argsort(p,kind="mergesort"); m=len(p)
    vals=np.array([(m-i)*p[order[i]] for i in range(m)],float); vals=np.minimum(1,np.maximum.accumulate(vals)); out=np.empty(m); out[order]=vals
    return {k:float(out[i]) for i,k in enumerate(keys)}


def test_tranche(m: pd.DataFrame,tranche: str,seed:int) -> dict[str,Any]:
    pure=m[m.organization_state.isin(["local_cooccurrence_only","spatial_segregation_only"])].copy(); pure=pure[pure.metric_status.eq("complete")]
    s=pure.organization_state.eq("spatial_segregation_only").to_numpy(bool)
    tests={}
    for j,(hid,(metric,alt)) in enumerate(PRED.items()): tests[hid]=perm_contrast(pure[metric].to_numpy(float),s,alt,seed+j)
    adj=holm({k:v["p"] for k,v in tests.items()})
    for k in tests:
        tests[k]["holm_p"]=adj[k]; metric,alt=PRED[k]; d=tests[k]["delta_S_minus_C"]
        tests[k]["direction_ok"]=bool(d is not None and ((alt=="greater" and d>0) or (alt=="less" and d<0)))
        tests[k]["pass"]=bool(tests[k]["direction_ok"] and adj[k]<=.05)
    h0=perm_contrast(pure.hv3d_rarefied20.to_numpy(float),s,"two-sided",seed+50)
    return {"tranche":tranche,"pure_species_complete":int(len(pure)),"pure_C":int((~s).sum()),"pure_S":int(s.sum()),"tests":tests,"H0_hypervolume_negative_control":h0}


def main():
    import argparse
    ap=argparse.ArgumentParser(); ap.add_argument("--worldclim-dir",required=True); a=ap.parse_args(); wcdir=Path(a.worldclim_dir)
    states=pd.read_csv(STATES); states["C_star"]=states.C_star.astype(str).str.lower().eq("true"); states["S_star"]=states.S_star.astype(str).str.lower().eq("true")
    admitted=set(states.species.astype(str))
    dc,dqc=extract_climate(DISC,"discovery",wcdir,admitted); rc,rqc=extract_climate(RES,"reserve",wcdir,admitted)
    dc.to_csv(OUT/"discovery_worldclim_cells.csv.gz",index=False,compression="gzip"); rc.to_csv(OUT/"reserve_worldclim_cells.csv.gz",index=False,compression="gzip")
    dm=species_metrics(dc,states,"discovery",SEED+1); rm=species_metrics(rc,states,"reserve",SEED+2); metrics=pd.concat([dm,rm],ignore_index=True); metrics.to_csv(OUT/"species_environment_metrics.csv",index=False)
    ds=test_tranche(dm,"discovery",SEED+101); rs=test_tranche(rm,"reserve",SEED+201)
    recurrent={hid:bool(ds["tests"][hid]["pass"] and rs["tests"][hid]["pass"]) for hid in PRED}
    result={"analysis":"rgfca_worldclim_old34_transfer_step7d","source":{"product":"WorldClim 2.1","resolution":"10 arc-minute","bio_variables":SELECT,"url":"https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_10m_bio.zip"},"discovery_qc":dqc,"reserve_qc":rqc,"discovery":ds,"reserve":rs,"two_tranche_recurrent":recurrent,"claim_boundary":"Transferred historical climate estimands on photo-derived pure C*/S* states; not confirmatory replication, causal selection, or dynamic temporal climate."}
    (OUT/"result.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    lines=["# RGFCA Step 7D — historical WorldClim transfer","",f"- discovery pure C/S complete: **{ds['pure_C']} C* / {ds['pure_S']} S***",f"- reserve pure C/S complete: **{rs['pure_C']} C* / {rs['pure_S']} S***","","## Directional hypotheses",""]
    for hid in PRED:
        d=ds['tests'][hid]; r=rs['tests'][hid]
        lines += [f"### {hid}",f"- discovery delta(S-C) **{d['delta_S_minus_C']}**, raw p **{d['p']:.6g}**, Holm **{d['holm_p']:.6g}**, pass **{d['pass']}**",f"- reserve delta(S-C) **{r['delta_S_minus_C']}**, raw p **{r['p']:.6g}**, Holm **{r['holm_p']:.6g}**, pass **{r['pass']}**",f"- two-tranche recurrent: **{recurrent[hid]}**",""]
    lines += ["## H0 total niche-size negative control","",f"- discovery delta(S-C) **{ds['H0_hypervolume_negative_control']['delta_S_minus_C']}**, p **{ds['H0_hypervolume_negative_control']['p']:.6g}**",f"- reserve delta(S-C) **{rs['H0_hypervolume_negative_control']['delta_S_minus_C']}**, p **{rs['H0_hypervolume_negative_control']['p']:.6g}**","","BIO4/BIO15 are climatological seasonality, not dynamic year-specific climate. H8 remains unopened."]
    (OUT/"RESULT.md").write_text("\n".join(lines)+"\n",encoding="utf-8")

if __name__=="__main__": main()

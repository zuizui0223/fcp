#!/usr/bin/env python3
from __future__ import annotations
import argparse,json,warnings
from pathlib import Path
import numpy as np
import pandas as pd
import rasterio
from scipy.stats import wilcoxon
from statsmodels.discrete.conditional_models import ConditionalLogit

MORPHS={"white","yellow_orange","red_pink","blue_purple"}
BG=[f"background_palette_count_{x}" for x in ["white","yellow","orange","red","pink","magenta","purple","blue","bronze","green","brown","black"]]

def truthy(s):
    return s.astype(str).str.strip().str.lower().isin({"true","1","yes","y"})

def wz(df,col,group="species"):
    g=df.groupby(group)[col]
    mu=g.transform("mean"); sd=g.transform(lambda x:x.std(ddof=0)).replace(0,np.nan)
    return (df[col]-mu)/sd

def sample(path,lon,lat):
    with rasterio.open(path) as src:
        v=np.array([z[0] for z in src.sample(list(zip(lon.astype(float),lat.astype(float))))],float)
        if src.nodata is not None: v[np.isclose(v,float(src.nodata),equal_nan=False)]=np.nan
        return v

def cells(path,lon,lat):
    with rasterio.open(path) as src:
        rc=[src.index(float(x),float(y)) for x,y in zip(lon,lat)]
    return np.array([r for r,c in rc]),np.array([c for r,c in rc])

def clogit(df,predictors,group):
    q=df.dropna(subset=["white",group]+predictors).copy()
    st=q.groupby(group).white.agg(["sum","count"])
    good=set(st.index[(st["sum"]>0)&(st["sum"]<st["count"])])
    q=q[q[group].isin(good)].copy()
    vv=q.groupby(group)[predictors[0]].std(ddof=0).fillna(0)
    q=q[q[group].isin(set(vv.index[vv>0]))].copy()
    if q.empty: return {"estimable":False}
    groups=q[group].astype("category").cat.codes.to_numpy()
    X=q[predictors].to_numpy(float)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fit=ConditionalLogit(q.white.to_numpy(float),X,groups=groups).fit(method="bfgs",maxiter=500,disp=False)
    return {
      "estimable":True,"n_rows":int(len(q)),"n_strata":int(q[group].nunique()),"n_species":int(q.species.nunique()),
      "beta":{p:float(fit.params[i]) for i,p in enumerate(predictors)},
      "p":{p:float(fit.pvalues[i]) for i,p in enumerate(predictors)},
      "OR":{p:float(np.exp(fit.params[i])) for i,p in enumerate(predictors)}
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--reserve",required=True); ap.add_argument("--bio5",required=True); ap.add_argument("--tmax-dir",required=True); ap.add_argument("--outdir",required=True)
    a=ap.parse_args(); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    use=["species","morph","global_classifiable","latitude","longitude","observed_on","nuisance_pixel_fraction"]+BG
    d=pd.read_csv(a.reserve,usecols=use,low_memory=False)
    d=d[truthy(d.global_classifiable)&d.morph.isin(MORPHS)].copy()
    d["white"]=(d.morph=="white").astype(int)
    d["month"]=pd.to_datetime(d.observed_on,errors="coerce").dt.month
    for c in BG+["nuisance_pixel_fraction"]: d[c]=pd.to_numeric(d[c],errors="coerce")
    denom=d[BG].sum(axis=1,min_count=1)
    d["background_white_fraction"]=np.where(denom>0,d["background_palette_count_white"]/denom,np.nan)
    d=d[d.month.notna()].copy()

    d["bio5"]=sample(Path(a.bio5),d.longitude,d.latitude)
    mats=[]
    for m in range(1,13):
        p=Path(a.tmax_dir)/f"wc2.1_10m_tmax_{m:02d}.tif"
        if not p.exists(): raise SystemExit(f"missing {p}")
        mats.append(sample(p,d.longitude,d.latitude))
    tm=np.vstack(mats).T
    d["tmax_mean12"]=np.nanmean(tm,axis=1)
    mm=d.month.astype(int).to_numpy()
    d["tmax_month"]=tm[np.arange(len(d)),mm-1]
    d["seasonal_heat_anomaly"]=d.tmax_month-d.tmax_mean12
    rr,cc=cells(Path(a.bio5),d.longitude,d.latitude)
    d["species_cell"]=d.species.astype(str)+"|"+rr.astype(str)+"|"+cc.astype(str)

    for c in ["seasonal_heat_anomaly","bio5","background_white_fraction","nuisance_pixel_fraction"]:
        d["wz_"+c]=wz(d,c)

    counts=d.groupby("species").white.agg(["sum","count"]); counts["nonwhite"]=counts["count"]-counts["sum"]
    eligible=set(counts.index[(counts["sum"]>=5)&(counts.nonwhite>=5)])
    q=d[d.species.isin(eligible)].copy()
    if len(eligible)<100: raise SystemExit(f"reserve species gate failed: {len(eligible)}")

    rows=[]
    for sp,g in q.dropna(subset=["wz_seasonal_heat_anomaly"]).groupby("species"):
        if g.white.sum()<5 or len(g)-g.white.sum()<5: continue
        delta=float(g.loc[g.white.eq(1),"wz_seasonal_heat_anomaly"].mean()-g.loc[g.white.eq(0),"wz_seasonal_heat_anomaly"].mean())
        rows.append({"species":sp,"delta_white_minus_nonwhite_SD":delta,"n":len(g),"n_white":int(g.white.sum())})
    dd=pd.DataFrame(rows)
    if len(dd)<100: raise SystemExit(f"reserve delta gate failed: {len(dd)}")
    vals=dd.delta_white_minus_nonwhite_SD.to_numpy(float)
    w=wilcoxon(vals,zero_method="wilcox",alternative="two-sided")

    spmodel=clogit(q,["wz_seasonal_heat_anomaly","wz_bio5","wz_background_white_fraction","wz_nuisance_pixel_fraction"],"species")
    sd=float(q.seasonal_heat_anomaly.std(ddof=0)); q["z_seasonal_heat_global"]=(q.seasonal_heat_anomaly-float(q.seasonal_heat_anomaly.mean()))/sd
    local=clogit(q,["z_seasonal_heat_global","wz_background_white_fraction","wz_nuisance_pixel_fraction"],"species_cell")
    local_estimable=bool(local.get("estimable") and local.get("n_strata",0)>=50 and local.get("n_species",0)>=30)

    species_ok=bool(np.median(vals)>0 and w.pvalue<0.05)
    model_ok=bool(spmodel.get("estimable") and spmodel["beta"]["wz_seasonal_heat_anomaly"]>0 and spmodel["p"]["wz_seasonal_heat_anomaly"]<0.05)
    local_ok=bool(local_estimable and local["beta"]["z_seasonal_heat_global"]>0 and local["p"]["z_seasonal_heat_global"]<0.05)
    if not local_estimable: verdict="RESERVE_LOCAL_GATE_NOT_ESTIMABLE"
    elif species_ok and model_ok and local_ok: verdict="RESERVE_SEASONAL_HEAT_CONFIRMED"
    else: verdict="RESERVE_SEASONAL_HEAT_NOT_CONFIRMED"

    result={
      "schema":"fcp_white_seasonal_heat_reserve_v1","status":verdict,
      "reserve_species":int(len(dd)),"reserve_rows":int(len(q)),
      "species_level":{"median_delta":float(np.median(vals)),"mean_delta":float(np.mean(vals)),"fraction_gt0":float(np.mean(vals>0)),"wilcoxon_p":float(w.pvalue)},
      "species_model":spmodel,"same_cell_model":local,"same_cell_gate_estimable":local_estimable,
      "decision_components":{"species_level":species_ok,"species_model":model_ok,"same_cell":local_ok},
      "source_commit":"db2514f622468747a7d5896c2696fcbf13729ef6",
      "reserve_file_sha256_expected":"0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6",
      "hard_nonclaims":["same iNaturalist/measurement programme, not independent-source replication","does not establish temperature-induced plasticity","does not establish selection or pigment mechanism","does not establish evolutionary transition direction"]
    }
    dd.to_csv(out/"reserve_species_deltas.csv",index=False)
    (out/"result.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
if __name__=="__main__": main()

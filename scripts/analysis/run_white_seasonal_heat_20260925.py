#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import rasterio
from scipy.stats import wilcoxon
from statsmodels.discrete.conditional_models import ConditionalLogit

MORPHS={"white","yellow_orange","red_pink","blue_purple"}

def truthy(s):
    if s.dtype == bool:
        return s.fillna(False)
    return s.astype(str).str.strip().str.lower().isin({"true","1","yes","y"})

def wz(df,col,group="species"):
    g=df.groupby(group)[col]
    mu=g.transform("mean")
    sd=g.transform(lambda x:x.std(ddof=0)).replace(0,np.nan)
    return (df[col]-mu)/sd

def sample_raster(path,lon,lat):
    with rasterio.open(path) as src:
        vals=np.array([x[0] for x in src.sample(list(zip(lon.astype(float),lat.astype(float))))],dtype=float)
        if src.nodata is not None:
            vals[np.isclose(vals,float(src.nodata),equal_nan=False)]=np.nan
        return vals

def raster_cells(path,lon,lat):
    with rasterio.open(path) as src:
        rc=[src.index(float(x),float(y)) for x,y in zip(lon,lat)]
    return np.array([x[0] for x in rc],int),np.array([x[1] for x in rc],int)

def clogit(df,predictors,group_col):
    q=df.dropna(subset=["white",group_col]+predictors).copy()
    st=q.groupby(group_col).white.agg(["sum","count"])
    good=set(st.index[(st["sum"]>0)&(st["sum"]<st["count"])])
    q=q[q[group_col].isin(good)].copy()
    # Groups with no predictor variation carry no slope information and can cause numerical singularity.
    var_ok=q.groupby(group_col)[predictors[0]].std(ddof=0)
    good2=set(var_ok.index[var_ok.fillna(0)>0])
    q=q[q[group_col].isin(good2)].copy()
    if q.empty:
        return {"estimable":False,"reason":"no informative strata"}
    groups=q[group_col].astype("category").cat.codes.to_numpy()
    X=q[predictors].to_numpy(float)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fit=ConditionalLogit(q.white.to_numpy(float),X,groups=groups).fit(method="bfgs",maxiter=500,disp=False)
    out={
      "estimable":True,
      "n_rows":int(len(q)),
      "n_strata":int(q[group_col].nunique()),
      "n_species":int(q.species.nunique()),
      "predictors":predictors,
      "beta":{p:float(fit.params[i]) for i,p in enumerate(predictors)},
      "se":{p:float(fit.bse[i]) for i,p in enumerate(predictors)},
      "p":{p:float(fit.pvalues[i]) for i,p in enumerate(predictors)},
      "OR":{p:float(np.exp(fit.params[i])) for i,p in enumerate(predictors)}
    }
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--measured",required=True)
    ap.add_argument("--technical-table",required=True)
    ap.add_argument("--join-key",required=True)
    ap.add_argument("--high-clip-ids",required=True)
    ap.add_argument("--bio5",required=True)
    ap.add_argument("--tmax-dir",required=True)
    ap.add_argument("--outdir",required=True)
    a=ap.parse_args()
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)

    cols=["photo_id","species","morph","global_classifiable","latitude","longitude","observed_on"]
    d=pd.read_csv(a.measured,usecols=cols,low_memory=False)
    d=d[truthy(d.global_classifiable)&d.morph.isin(MORPHS)].copy()
    d["white"]=(d.morph=="white").astype(int)
    d["month"]=pd.to_datetime(d.observed_on,errors="coerce").dt.month

    tech=pd.read_csv(a.technical_table,compression="gzip",dtype={"measurement_id":str})
    join=pd.read_csv(a.join_key,dtype={"measurement_id":str})
    high=pd.read_csv(a.high_clip_ids,dtype={"measurement_id":str})
    t=tech.merge(join,on="measurement_id",how="left",validate="one_to_one")
    highset=set(high.measurement_id.astype(str))
    t["high_clip"]=t.measurement_id.astype(str).isin(highset)
    d=d.merge(t[["photo_id","near_clip_fraction","high_clip"]],on="photo_id",how="left",validate="one_to_one")
    d=d[d.near_clip_fraction.notna() & ~d.high_clip.fillna(False) & d.month.notna()].copy()

    d["bio5"]=sample_raster(Path(a.bio5),d.longitude,d.latitude)
    tmax=[]
    for m in range(1,13):
        p=Path(a.tmax_dir)/f"wc2.1_10m_tmax_{m:02d}.tif"
        if not p.exists():
            raise SystemExit(f"missing {p}")
        tmax.append(sample_raster(p,d.longitude,d.latitude))
    tm=np.vstack(tmax).T
    d["tmax_annual_month_mean"]=np.nanmean(tm,axis=1)
    months=d.month.astype(int).to_numpy()
    d["tmax_observation_month"]=tm[np.arange(len(d)),months-1]
    d["seasonal_heat_anomaly"]=d.tmax_observation_month-d.tmax_annual_month_mean

    rr,cc=raster_cells(Path(a.bio5),d.longitude,d.latitude)
    d["climate_cell_row"]=rr; d["climate_cell_col"]=cc
    d["species_cell"]=d.species.astype(str)+"|"+d.climate_cell_row.astype(str)+"|"+d.climate_cell_col.astype(str)

    for c in ["seasonal_heat_anomaly","bio5","near_clip_fraction"]:
        d["wz_"+c]=wz(d,c)

    counts=d.groupby("species").white.agg(["sum","count"])
    counts["nonwhite"]=counts["count"]-counts["sum"]
    eligible=set(counts.index[(counts["sum"]>=5)&(counts["nonwhite"]>=5)])
    primary=d[d.species.isin(eligible)].copy()
    if len(eligible)<100:
        raise SystemExit(f"primary species gate failed: {len(eligible)}")

    per=[]
    for sp,g in primary.dropna(subset=["wz_seasonal_heat_anomaly"]).groupby("species"):
        if g.white.sum()<5 or (len(g)-g.white.sum())<5:
            continue
        delta=float(g.loc[g.white.eq(1),"wz_seasonal_heat_anomaly"].mean()-g.loc[g.white.eq(0),"wz_seasonal_heat_anomaly"].mean())
        per.append({"species":sp,"delta_white_minus_nonwhite_seasonal_heat_SD":delta,"n":len(g),"n_white":int(g.white.sum())})
    perdf=pd.DataFrame(per)
    if len(perdf)<100:
        raise SystemExit(f"species delta gate failed: {len(perdf)}")
    vals=perdf.delta_white_minus_nonwhite_seasonal_heat_SD.to_numpy(float)
    w=wilcoxon(vals,zero_method="wilcox",alternative="two-sided")

    species_model=clogit(
      primary,
      ["wz_seasonal_heat_anomaly","wz_bio5","wz_near_clip_fraction"],
      "species"
    )

    # Local same-cell model: use globally standardized seasonal anomaly to retain a common unit;
    # the species×cell stratum absorbs all time-invariant local climate.
    sd=float(primary.seasonal_heat_anomaly.std(ddof=0))
    primary["z_seasonal_heat_global"]=(primary.seasonal_heat_anomaly-float(primary.seasonal_heat_anomaly.mean()))/sd
    local_model=clogit(primary,["z_seasonal_heat_global","wz_near_clip_fraction"],"species_cell")
    local_gate_estimable=bool(local_model.get("estimable") and local_model.get("n_strata",0)>=50 and local_model.get("n_species",0)>=30)

    primary_positive=(float(np.median(vals))>0 and float(w.pvalue)<0.05)
    species_positive=bool(species_model.get("estimable") and species_model["beta"]["wz_seasonal_heat_anomaly"]>0 and species_model["p"]["wz_seasonal_heat_anomaly"]<0.05)
    local_positive=bool(local_gate_estimable and local_model["beta"]["z_seasonal_heat_global"]>0 and local_model["p"]["z_seasonal_heat_global"]<0.05)

    if not local_gate_estimable:
        verdict="SEASONAL_HEAT_LOCAL_GATE_NOT_ESTIMABLE"
    elif primary_positive and species_positive and local_positive:
        verdict="SEASONAL_HEAT_ASSOCIATION_SUPPORTED"
    else:
        verdict="SEASONAL_HEAT_ASSOCIATION_NOT_SUPPORTED_UNDER_THIS_TEST"

    result={
      "schema":"fcp_white_seasonal_heat_v1",
      "status":verdict,
      "role":"prospectively specified seasonal-heat discrimination after long-term BIO5 association was opened",
      "primary_species":int(len(perdf)),
      "primary_rows":int(len(primary)),
      "species_level":{
        "median_delta_white_minus_nonwhite_SD":float(np.median(vals)),
        "mean_delta_white_minus_nonwhite_SD":float(np.mean(vals)),
        "fraction_delta_gt_0":float(np.mean(vals>0)),
        "wilcoxon_two_sided_p":float(w.pvalue)
      },
      "species_stratified_model":species_model,
      "same_cell_local_model":local_model,
      "same_cell_local_gate_estimable":local_gate_estimable,
      "decision_components":{
        "species_level_positive":primary_positive,
        "species_stratified_positive":species_positive,
        "same_cell_positive":local_positive
      },
      "predictor_definition":"WorldClim monthly tmax at observed calendar month minus mean of 12 monthly tmax values at the same coordinate",
      "hard_nonclaims":[
        "WorldClim monthly values are climatology, not weather in the observation year",
        "does not establish temperature-induced plastic colour change",
        "does not establish selection on pigment chemistry",
        "does not establish pigmented-to-white evolutionary direction",
        "does not establish pollinator causation"
      ]
    }
    perdf.to_csv(out/"species_seasonal_heat_deltas.csv",index=False)
    (out/"result.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))

if __name__=="__main__":
    main()

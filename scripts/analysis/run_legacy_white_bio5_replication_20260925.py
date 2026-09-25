#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from scipy.stats import wilcoxon
from statsmodels.discrete.conditional_models import ConditionalLogit

MORPHS = ["white","yellow_orange","red_pink","blue_purple"]

def as_bool(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False)
    return s.fillna("").astype(str).str.strip().str.lower().isin({"true","1","yes","y"})

def coordinate_columns(df: pd.DataFrame) -> tuple[str,str]:
    pairs=[
        ("latitude","longitude"),
        ("decimalLatitude","decimalLongitude"),
        ("lat","lon"),
    ]
    for a,b in pairs:
        if a in df.columns and b in df.columns:
            return a,b
    raise RuntimeError("No recognized latitude/longitude columns")

def sample_raster(path: Path, lon: pd.Series, lat: pd.Series) -> np.ndarray:
    with rasterio.open(path) as src:
        pts=list(zip(pd.to_numeric(lon,errors="coerce"),pd.to_numeric(lat,errors="coerce")))
        vals=np.full(len(pts),np.nan,dtype=float)
        good=[i for i,(x,y) in enumerate(pts) if np.isfinite(x) and np.isfinite(y)]
        if good:
            samp=list(src.sample([pts[i] for i in good]))
            for i,v in zip(good,samp):
                x=float(v[0])
                if src.nodata is not None and np.isclose(x,src.nodata):
                    x=np.nan
                vals[i]=x
        return vals

def conditional_fit(df: pd.DataFrame, group_col: str) -> dict:
    x=df.dropna(subset=["white","bio5_z",group_col]).copy()
    stats=x.groupby(group_col).white.agg(["sum","count"])
    good=set(stats.index[(stats["sum"]>0)&(stats["sum"]<stats["count"])])
    x=x[x[group_col].isin(good)].copy()
    if x.empty or x[group_col].nunique()<5:
        return {"estimable":False,"n_rows":int(len(x)),"n_strata":int(x[group_col].nunique()),"n_species":int(x.species.nunique())}
    groups=x[group_col].astype("category").cat.codes.to_numpy()
    X=x[["bio5_z"]].to_numpy(float)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fit=ConditionalLogit(x.white.to_numpy(float),X,groups=groups).fit(method="bfgs",maxiter=500,disp=False)
    b=float(fit.params[0]); se=float(fit.bse[0]); p=float(fit.pvalues[0])
    return {
        "estimable":True,
        "n_rows":int(len(x)),
        "n_strata":int(x[group_col].nunique()),
        "n_species":int(x.species.nunique()),
        "beta":b,
        "se":se,
        "OR_per_within_species_SD":float(np.exp(b)),
        "ci_low":float(np.exp(b-1.96*se)),
        "ci_high":float(np.exp(b+1.96*se)),
        "p":p,
    }

def analyze(path: Path, cohort: str, bio5_path: Path) -> tuple[dict,pd.DataFrame]:
    d=pd.read_csv(path,low_memory=False)
    required={"species","morph","global_classifiable","observer_id"}
    missing=sorted(required-set(d.columns))
    if missing:
        raise RuntimeError(f"{cohort}: missing required columns {missing}")
    latc,lonc=coordinate_columns(d)
    d=d.copy()
    d["species"]=d.species.fillna("").astype(str).str.strip()
    d["observer_id"]=d.observer_id.fillna("").astype(str).str.strip()
    keep=(as_bool(d.global_classifiable)&d.morph.isin(MORPHS)&d.species.ne(""))
    d=d.loc[keep,["species","observer_id","morph",latc,lonc]].copy()
    d["white"]=(d.morph=="white").astype(int)
    d["bio5"]=sample_raster(bio5_path,d[lonc],d[latc])
    d=d.dropna(subset=["bio5"]).copy()
    mu=d.groupby("species").bio5.transform("mean")
    sd=d.groupby("species").bio5.transform(lambda s:s.std(ddof=0)).replace(0,np.nan)
    d["bio5_z"]=(d.bio5-mu)/sd
    d=d.dropna(subset=["bio5_z"]).copy()

    counts=d.groupby("species").white.agg(["sum","count"])
    counts["nonwhite"]=counts["count"]-counts["sum"]
    eligible=set(counts.index[(counts["sum"]>=5)&(counts["nonwhite"]>=5)])
    p=d[d.species.isin(eligible)].copy()

    per=[]
    for sp,g in p.groupby("species"):
        delta=float(g.loc[g.white.eq(1),"bio5_z"].mean()-g.loc[g.white.eq(0),"bio5_z"].mean())
        per.append({
            "cohort":cohort,
            "species":sp,
            "delta_white_minus_nonwhite_SD":delta,
            "n":int(len(g)),
            "n_white":int(g.white.sum()),
            "n_nonwhite":int(len(g)-g.white.sum()),
        })
    per=pd.DataFrame(per)
    estimable=len(per)>=100
    if estimable:
        w=wilcoxon(per.delta_white_minus_nonwhite_SD.to_numpy(float),zero_method="wilcox",alternative="two-sided")
        med=float(per.delta_white_minus_nonwhite_SD.median())
        mean=float(per.delta_white_minus_nonwhite_SD.mean())
        frac=float((per.delta_white_minus_nonwhite_SD>0).mean())
        wp=float(w.pvalue)
    else:
        med=mean=frac=wp=np.nan

    species_fit=conditional_fit(p.assign(species_group=p.species),"species_group") if estimable else {"estimable":False}
    p["species_observer"]=p.species.astype(str)+"|||"+p.observer_id.astype(str)
    obs_fit=conditional_fit(p.loc[p.observer_id.ne("")],"species_observer")
    obs_fit["coverage_gate_30_species"]=bool(obs_fit.get("n_species",0)>=30)

    support=bool(
        estimable and med>0 and wp<0.05 and
        species_fit.get("estimable",False) and
        species_fit.get("beta",np.nan)>0 and
        species_fit.get("p",1)>=0 and species_fit.get("p",1)<0.05
    )

    out={
        "cohort":cohort,
        "coordinate_columns":[latc,lonc],
        "classifiable_rows_with_bio5":int(len(d)),
        "eligible_species":int(len(per)),
        "estimable":bool(estimable),
        "species_level":{
            "median_delta_white_minus_nonwhite_SD":med,
            "mean_delta_white_minus_nonwhite_SD":mean,
            "fraction_delta_gt_0":frac,
            "wilcoxon_two_sided_p":wp,
        },
        "species_stratified_model":species_fit,
        "species_observer_sensitivity":obs_fit,
        "primary_support":support,
    }
    return out,per

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--discovery",required=True)
    ap.add_argument("--reserve",required=True)
    ap.add_argument("--bio5",required=True)
    ap.add_argument("--outdir",required=True)
    a=ap.parse_args()
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    disc,dp=analyze(Path(a.discovery),"discovery",Path(a.bio5))
    res,rp=analyze(Path(a.reserve),"reserve",Path(a.bio5))
    both=bool(disc["primary_support"] and res["primary_support"])
    verdict="LEGACY_BIO5_WHITE_REPLICATION_SUPPORTED" if both else "LEGACY_BIO5_WHITE_REPLICATION_NOT_SUPPORTED_UNDER_THIS_TEST"
    result={
        "schema":"fcp_legacy_white_bio5_replication_v1",
        "status":"complete",
        "verdict":verdict,
        "source_commit":"5142f7951af0dde5364bb047a566d67e8c479e51",
        "predictor":"WorldClim 2.1 10 arc-minute BIO5",
        "prediction":"white records occupy higher BIO5 within species",
        "discovery":disc,
        "reserve":res,
        "hard_nonclaims":[
            "does not establish causal heat selection",
            "does not distinguish genetic from plastic white states",
            "does not remove known image-exposure coupling of the coarse white classifier",
            "same iNaturalist and measurement system, so this is species-disjoint transport rather than independent-source causation"
        ]
    }
    pd.concat([dp,rp],ignore_index=True).to_csv(out/"species_bio5_deltas.csv",index=False)
    (out/"result.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))

if __name__=="__main__":
    main()

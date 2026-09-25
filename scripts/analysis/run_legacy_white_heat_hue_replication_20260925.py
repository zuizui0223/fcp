#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import rasterio
from scipy.stats import wilcoxon
from statsmodels.discrete.conditional_models import ConditionalLogit

MORPHS=["white","yellow_orange","red_pink","blue_purple"]
ANTHO={"red_pink","blue_purple"}
MIN_A=100
MIN_Y=50
MIN_S=30

def as_bool(s):
    if pd.api.types.is_bool_dtype(s): return s.fillna(False)
    return s.fillna("").astype(str).str.strip().str.lower().isin({"true","1","yes","y"})

def coords(df):
    for a,b in [("latitude","longitude"),("decimalLatitude","decimalLongitude"),("lat","lon")]:
        if a in df.columns and b in df.columns: return a,b
    raise RuntimeError("no coordinate columns")

def sample(path,lon,lat):
    with rasterio.open(path) as src:
        xy=list(zip(pd.to_numeric(lon,errors="coerce"),pd.to_numeric(lat,errors="coerce")))
        out=np.full(len(xy),np.nan,float)
        good=[i for i,(x,y) in enumerate(xy) if np.isfinite(x) and np.isfinite(y)]
        vv=list(src.sample([xy[i] for i in good])) if good else []
        for i,v in zip(good,vv):
            z=float(v[0])
            if src.nodata is not None and np.isclose(z,src.nodata): z=np.nan
            out[i]=z
        return out

def wz(df,col):
    mu=df.groupby("species")[col].transform("mean")
    sd=df.groupby("species")[col].transform(lambda x:x.std(ddof=0)).replace(0,np.nan)
    return (df[col]-mu)/sd

def signed(vals,min_n):
    x=np.asarray(vals,float); x=x[np.isfinite(x)]
    r={"n_species":int(len(x)),"minimum_required":min_n,"estimable":bool(len(x)>=min_n)}
    if not len(x):
        r.update(median_delta=None,mean_delta=None,fraction_positive=None,wilcoxon_one_sided_greater_p=None,wilcoxon_two_sided_p=None)
        return r
    r.update(
        median_delta=float(np.median(x)),
        mean_delta=float(np.mean(x)),
        fraction_positive=float(np.mean(x>0)),
        wilcoxon_one_sided_greater_p=float(wilcoxon(x,zero_method="wilcox",alternative="greater",method="auto").pvalue),
        wilcoxon_two_sided_p=float(wilcoxon(x,zero_method="wilcox",alternative="two-sided",method="auto").pvalue),
    )
    return r

def clogit(df):
    x=df.dropna(subset=["white","bio5_z","species"]).copy()
    st=x.groupby("species").white.agg(["sum","count"])
    good=set(st.index[(st["sum"]>0)&(st["sum"]<st["count"])])
    x=x[x.species.isin(good)].copy()
    if x.species.nunique()<20:
        return {"estimable":False,"n_species":int(x.species.nunique()),"n_rows":int(len(x))}
    groups=x.species.astype("category").cat.codes.to_numpy()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fit=ConditionalLogit(x.white.to_numpy(float),x[["bio5_z"]].to_numpy(float),groups=groups).fit(method="bfgs",maxiter=500,disp=False)
    b=float(fit.params[0]); se=float(fit.bse[0])
    return {"estimable":True,"n_species":int(x.species.nunique()),"n_rows":int(len(x)),
            "beta_BIO5":b,"OR_BIO5_per_within_species_SD":float(np.exp(b)),
            "ci_low":float(np.exp(b-1.96*se)),"ci_high":float(np.exp(b+1.96*se)),
            "p_BIO5":float(fit.pvalues[0])}

def analyze(path,cohort,bio5):
    d=pd.read_csv(path,low_memory=False)
    need={"species","morph","global_classifiable","observer_id"}
    miss=sorted(need-set(d.columns))
    if miss: raise RuntimeError(f"{cohort}: missing {miss}")
    latc,lonc=coords(d)
    d=d.copy()
    d["species"]=d.species.fillna("").astype(str).str.strip()
    d["observer_id"]=d.observer_id.fillna("").astype(str).str.strip()
    d=d[as_bool(d.global_classifiable)&d.morph.isin(MORPHS)&d.species.ne("")].copy()
    d["BIO5"]=sample(bio5,d[lonc],d[latc])
    d=d.dropna(subset=["BIO5"]).copy()
    d["bio5_z"]=wz(d,"BIO5")
    d=d.dropna(subset=["bio5_z"]).copy()

    c=(d.assign(w=d.morph.eq("white"),a=d.morph.isin(ANTHO),y=d.morph.eq("yellow_orange"))
         .groupby("species").agg(n_white=("w","sum"),n_antho=("a","sum"),n_yellow=("y","sum")).reset_index())
    eA=set(c.loc[c.n_white.ge(5)&c.n_antho.ge(5),"species"])
    eY=set(c.loc[c.n_white.ge(5)&c.n_yellow.ge(5),"species"])
    eS=set(c.loc[c.n_white.ge(5)&c.n_antho.ge(5)&c.n_yellow.ge(5),"species"])

    da={}; dy={}; rows=[]
    for sp,g in d.groupby("species"):
        w=g.loc[g.morph.eq("white"),"bio5_z"]
        aa=g.loc[g.morph.isin(ANTHO),"bio5_z"]
        yy=g.loc[g.morph.eq("yellow_orange"),"bio5_z"]
        if sp in eA:
            z=float(w.mean()-aa.mean()); da[sp]=z
            rows.append({"cohort":cohort,"species":sp,"contrast":"white_minus_anthocyanin_proxy","delta":z})
        if sp in eY:
            z=float(w.mean()-yy.mean()); dy[sp]=z
            rows.append({"cohort":cohort,"species":sp,"contrast":"white_minus_yellow_orange","delta":z})
    common=sorted(eS & set(da) & set(dy))
    sval=[da[s]-dy[s] for s in common]

    A=signed(list(da.values()),MIN_A)
    Y=signed(list(dy.values()),MIN_Y)
    S=signed(sval,MIN_S)

    Arow=d[d.species.isin(eA)&(d.morph.eq("white")|d.morph.isin(ANTHO))].copy()
    Arow["white"]=Arow.morph.eq("white").astype(int)
    Yrow=d[d.species.isin(eY)&(d.morph.eq("white")|d.morph.eq("yellow_orange"))].copy()
    Yrow["white"]=Yrow.morph.eq("white").astype(int)
    Am=clogit(Arow); Ym=clogit(Yrow)

    A_support=bool(A["estimable"] and A["median_delta"]>0 and A["wilcoxon_one_sided_greater_p"]<0.05 and Am.get("estimable") and Am["beta_BIO5"]>0 and Am["p_BIO5"]<0.05)
    S_support=bool(S["estimable"] and S["median_delta"]>0 and S["wilcoxon_one_sided_greater_p"]<0.05)

    obs={}
    # Same-observer sensitivity for A only.
    q=Arow[Arow.observer_id.ne("")].copy()
    q["species_observer"]=q.species+"|||"+q.observer_id
    st=q.groupby("species_observer").white.agg(["sum","count"])
    good=set(st.index[(st["sum"]>0)&(st["sum"]<st["count"])])
    q=q[q.species_observer.isin(good)].copy()
    if q.species.nunique()>=30:
        groups=q.species_observer.astype("category").cat.codes.to_numpy()
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                fit=ConditionalLogit(q.white.to_numpy(float),q[["bio5_z"]].to_numpy(float),groups=groups).fit(method="bfgs",maxiter=500,disp=False)
            b=float(fit.params[0]); se=float(fit.bse[0])
            obs={"estimable":True,"n_species":int(q.species.nunique()),"n_strata":int(q.species_observer.nunique()),
                 "OR_BIO5":float(np.exp(b)),"ci_low":float(np.exp(b-1.96*se)),"ci_high":float(np.exp(b+1.96*se)),"p_BIO5":float(fit.pvalues[0])}
        except Exception as e:
            obs={"estimable":False,"reason":str(e)[:300]}
    else:
        obs={"estimable":False,"n_species":int(q.species.nunique()),"n_strata":int(q.species_observer.nunique())}

    return {"cohort":cohort,"coordinate_columns":[latc,lonc],
            "anthocyanin_proxy":A,"yellow_orange_comparator":Y,"direct_specificity":S,
            "row_level_anthocyanin_proxy":Am,"row_level_yellow_orange":Ym,
            "same_observer_anthocyanin_proxy_sensitivity":obs,
            "anthocyanin_proxy_support":A_support,"specificity_support":S_support}, pd.DataFrame(rows)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--discovery",required=True); p.add_argument("--reserve",required=True)
    p.add_argument("--bio5",required=True); p.add_argument("--outdir",required=True)
    a=p.parse_args(); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    D,dr=analyze(Path(a.discovery),"discovery",Path(a.bio5))
    R,rr=analyze(Path(a.reserve),"reserve",Path(a.bio5))
    at=bool(D["anthocyanin_proxy_support"] and R["anthocyanin_proxy_support"])
    sp=bool(D["specificity_support"] and R["specificity_support"])
    result={"schema":"fcp_legacy_white_heat_hue_replication_v1","status":"complete",
            "source_commit":"5142f7951af0dde5364bb047a566d67e8c479e51",
            "opened_after_third_cohort_hue_result":True,
            "anthocyanin_proxy_transport_verdict":"ANTHOCYANIN_PROXY_TRANSPORT_SUPPORTED" if at else "ANTHOCYANIN_PROXY_TRANSPORT_NOT_SUPPORTED_UNDER_THIS_TEST",
            "hue_specificity_transport_verdict":"HUE_SPECIFICITY_TRANSPORT_SUPPORTED" if sp else "HUE_SPECIFICITY_TRANSPORT_NOT_SUPPORTED_UNDER_THIS_TEST",
            "discovery":D,"reserve":R,
            "hard_nonclaims":["post-result species-disjoint replication, not untouched confirmation",
                              "coarse hue groups are pigment proxies, not chemical measurements",
                              "legacy cohorts lack direct highlight technical control",
                              "does not change the failed overall BIO5 replication verdict",
                              "does not establish causal heat selection"]}
    pd.concat([dr,rr],ignore_index=True).to_csv(out/"species_hue_bio5_deltas.csv",index=False)
    (out/"result.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
if __name__=="__main__": main()

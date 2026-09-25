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
MIN_SPEC=30

def bool_series(s):
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False)
    return s.fillna("").astype(str).str.strip().str.lower().isin({"true","1","yes","y"})

def sample_raster(path, lon, lat):
    with rasterio.open(path) as src:
        pts=list(zip(pd.to_numeric(lon,errors="coerce"),pd.to_numeric(lat,errors="coerce")))
        vals=np.full(len(pts),np.nan,float)
        good=[i for i,(x,y) in enumerate(pts) if np.isfinite(x) and np.isfinite(y)]
        if good:
            vv=list(src.sample([pts[i] for i in good]))
            for i,v in zip(good,vv):
                x=float(v[0])
                if src.nodata is not None and np.isclose(x,src.nodata):
                    x=np.nan
                vals[i]=x
        return vals

def wz(df,col):
    mu=df.groupby("species")[col].transform("mean")
    sd=df.groupby("species")[col].transform(lambda x:x.std(ddof=0)).replace(0,np.nan)
    return (df[col]-mu)/sd

def signed_summary(vals, min_n):
    vals=np.asarray(vals,float)
    vals=vals[np.isfinite(vals)]
    out={
        "n_species":int(len(vals)),
        "minimum_required":int(min_n),
        "estimable":bool(len(vals)>=min_n),
        "median_delta":float(np.median(vals)) if len(vals) else None,
        "mean_delta":float(np.mean(vals)) if len(vals) else None,
        "fraction_positive":float(np.mean(vals>0)) if len(vals) else None,
    }
    if len(vals):
        try:
            out["wilcoxon_one_sided_greater_p"]=float(wilcoxon(vals,zero_method="wilcox",alternative="greater",method="auto").pvalue)
            out["wilcoxon_two_sided_p"]=float(wilcoxon(vals,zero_method="wilcox",alternative="two-sided",method="auto").pvalue)
        except ValueError:
            out["wilcoxon_one_sided_greater_p"]=None
            out["wilcoxon_two_sided_p"]=None
    else:
        out["wilcoxon_one_sided_greater_p"]=None
        out["wilcoxon_two_sided_p"]=None
    return out

def conditional_fit(df, comparator_name):
    x=df.copy()
    x["near_z"]=wz(x,"near_clip_fraction")
    x=x.dropna(subset=["bio5_z","near_z","white","species"])
    st=x.groupby("species").white.agg(["sum","count"])
    good=set(st.index[(st["sum"]>0)&(st["sum"]<st["count"])])
    x=x[x.species.isin(good)].copy()
    if len(x)==0 or x.species.nunique()<20:
        return {"estimable":False,"contrast":comparator_name,"n_species":int(x.species.nunique()),"n_rows":int(len(x))}
    groups=x.species.astype("category").cat.codes.to_numpy()
    X=x[["bio5_z","near_z"]].to_numpy(float)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fit=ConditionalLogit(x.white.to_numpy(float),X,groups=groups).fit(method="bfgs",maxiter=500,disp=False)
    b=float(fit.params[0]); se=float(fit.bse[0])
    return {
        "estimable":True,
        "contrast":comparator_name,
        "n_species":int(x.species.nunique()),
        "n_rows":int(len(x)),
        "beta_BIO5":b,
        "OR_BIO5_per_within_species_SD":float(np.exp(b)),
        "ci_low":float(np.exp(b-1.96*se)),
        "ci_high":float(np.exp(b+1.96*se)),
        "p_BIO5":float(fit.pvalues[0]),
        "beta_near_clip":float(fit.params[1]),
        "p_near_clip":float(fit.pvalues[1]),
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--measured",required=True)
    p.add_argument("--technical-table",required=True)
    p.add_argument("--join-key",required=True)
    p.add_argument("--high-clip-ids",required=True)
    p.add_argument("--bio5",required=True)
    p.add_argument("--outdir",required=True)
    a=p.parse_args()

    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)

    use=["photo_id","species","morph","global_classifiable","latitude","longitude"]
    d=pd.read_csv(a.measured,usecols=use,low_memory=False)
    d=d[bool_series(d.global_classifiable)&d.morph.isin(MORPHS)].copy()

    tech=pd.read_csv(a.technical_table,compression="gzip",dtype={"measurement_id":str})
    join=pd.read_csv(a.join_key,dtype={"measurement_id":str})
    high=pd.read_csv(a.high_clip_ids,dtype={"measurement_id":str})
    tech=tech.merge(join,on="measurement_id",how="left",validate="one_to_one")
    highset=set(high.measurement_id.astype(str))
    tech["high_clip"]=tech.measurement_id.astype(str).isin(highset)

    d=d.merge(tech[["photo_id","near_clip_fraction","high_clip"]],on="photo_id",how="left",validate="one_to_one")
    d=d[d.near_clip_fraction.notna() & ~d.high_clip.fillna(False)].copy()
    d["BIO5"]=sample_raster(Path(a.bio5),d.longitude,d.latitude)
    d=d.dropna(subset=["BIO5"]).copy()
    d["bio5_z"]=wz(d,"BIO5")
    d=d.dropna(subset=["bio5_z"]).copy()

    # species supports
    counts=(d.assign(is_white=d.morph.eq("white"),
                     is_antho=d.morph.isin(ANTHO),
                     is_yellow=d.morph.eq("yellow_orange"))
              .groupby("species")
              .agg(n_white=("is_white","sum"),
                   n_antho=("is_antho","sum"),
                   n_yellow=("is_yellow","sum"))
              .reset_index())

    elig_A=set(counts.loc[counts.n_white.ge(5)&counts.n_antho.ge(5),"species"])
    elig_Y=set(counts.loc[counts.n_white.ge(5)&counts.n_yellow.ge(5),"species"])
    elig_S=set(counts.loc[counts.n_white.ge(5)&counts.n_antho.ge(5)&counts.n_yellow.ge(5),"species"])

    rows=[]
    deltaA={}
    deltaY={}
    for sp,g in d.groupby("species"):
        wzv=g.bio5_z
        w=g.loc[g.morph.eq("white"),"bio5_z"]
        aa=g.loc[g.morph.isin(ANTHO),"bio5_z"]
        yy=g.loc[g.morph.eq("yellow_orange"),"bio5_z"]
        if sp in elig_A:
            da=float(w.mean()-aa.mean()); deltaA[sp]=da
            rows.append({"species":sp,"contrast":"white_minus_anthocyanin_proxy","delta":da,"n_white":len(w),"n_comparator":len(aa)})
        if sp in elig_Y:
            dy=float(w.mean()-yy.mean()); deltaY[sp]=dy
            rows.append({"species":sp,"contrast":"white_minus_yellow_orange","delta":dy,"n_white":len(w),"n_comparator":len(yy)})

    common=sorted(elig_S & set(deltaA) & set(deltaY))
    spec_vals=[deltaA[sp]-deltaY[sp] for sp in common]
    specificity_rows=[{"species":sp,"delta_A":deltaA[sp],"delta_Y":deltaY[sp],"specificity_delta_A_minus_Y":deltaA[sp]-deltaY[sp]} for sp in common]

    A_summary=signed_summary(list(deltaA.values()),MIN_A)
    Y_summary=signed_summary(list(deltaY.values()),MIN_Y)
    S_summary=signed_summary(spec_vals,MIN_SPEC)

    # row-level binary contrast panels
    Arow=d[d.species.isin(elig_A) & (d.morph.eq("white")|d.morph.isin(ANTHO))].copy()
    Arow["white"]=d.loc[Arow.index,"morph"].eq("white").astype(int)
    Yrow=d[d.species.isin(elig_Y) & (d.morph.eq("white")|d.morph.eq("yellow_orange"))].copy()
    Yrow["white"]=d.loc[Yrow.index,"morph"].eq("white").astype(int)

    result={
        "schema":"fcp_white_heat_hue_mechanism_discriminator_v1",
        "status":"complete",
        "role":"post_confirmatory_mechanism_discrimination",
        "hue_definitions":{
            "white":["white"],
            "anthocyanin_associated_proxy":["red_pink","blue_purple"],
            "yellow_orange_comparator":["yellow_orange"]
        },
        "technical_control":"frozen response-blind high-clip exclusion plus continuous near-clip adjustment in row-level models",
        "anthocyanin_proxy":A_summary,
        "yellow_orange_comparator":Y_summary,
        "direct_specificity_delta_A_minus_Y":S_summary,
        "row_level_anthocyanin_proxy":conditional_fit(Arow,"white_vs_redpink_bluepurple"),
        "row_level_yellow_orange":conditional_fit(Yrow,"white_vs_yellow_orange"),
        "interpretation_gate":{
            "anthocyanin_proxy_positive":bool(A_summary["estimable"] and A_summary["median_delta"] is not None and A_summary["median_delta"]>0 and (A_summary["wilcoxon_one_sided_greater_p"] or 1)<0.05),
            "specificity_positive":bool(S_summary["estimable"] and S_summary["median_delta"] is not None and S_summary["median_delta"]>0 and (S_summary["wilcoxon_one_sided_greater_p"] or 1)<0.05)
        },
        "hard_nonclaims":[
            "coarse hue groups are pigment-system proxies rather than chemical measurements",
            "does not establish causal heat selection",
            "does not rescue the failed legacy BIO5 replication",
            "does not establish universal anthocyanin causation or evolutionary transition direction"
        ]
    }

    pd.DataFrame(rows).to_csv(out/"species_hue_specific_bio5_deltas.csv",index=False)
    pd.DataFrame(specificity_rows).to_csv(out/"species_direct_specificity.csv",index=False)
    counts.to_csv(out/"species_support_counts.csv",index=False)
    (out/"result.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))

if __name__=="__main__":
    main()

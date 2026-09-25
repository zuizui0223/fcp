#!/usr/bin/env python3
"""Exploratory white-state/pollinator screen.

Post-outcome sensitivity only. This script cannot produce a confirmatory verdict.
Inputs are a frozen third-cohort flower measurement table, the frozen highlight
technical seal, and an outcome-blind standardized pollinator-record export.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm

MORPHS={"white","yellow_orange","red_pink","blue_purple"}

def b(s):
    if s.dtype==bool: return s.fillna(False)
    return s.astype(str).str.strip().str.lower().isin({"true","1","yes","y"})

def z(s):
    x=pd.to_numeric(s,errors="coerce").astype(float)
    sd=float(x.std(ddof=0))
    return (x-float(x.mean()))/sd if np.isfinite(sd) and sd>0 else pd.Series(np.nan,index=x.index)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--measured",required=True)
    ap.add_argument("--technical-table",required=True)
    ap.add_argument("--join-key",required=True)
    ap.add_argument("--high-clip-ids",required=True)
    ap.add_argument("--pollinator-records",required=True)
    ap.add_argument("--coverage",required=True)
    ap.add_argument("--outdir",required=True)
    a=ap.parse_args()
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)

    d=pd.read_csv(a.measured,usecols=["photo_id","species","morph","global_classifiable"])
    d=d.loc[b(d.global_classifiable)&d.morph.isin(MORPHS)].copy()
    d["white"]=(d.morph=="white").astype(int)

    tech=pd.read_csv(a.technical_table,compression="gzip",dtype={"measurement_id":str})
    join=pd.read_csv(a.join_key,dtype={"measurement_id":str})
    high=set(pd.read_csv(a.high_clip_ids,dtype={"measurement_id":str}).measurement_id.astype(str))
    tech=tech.merge(join,on="measurement_id",validate="one_to_one")
    tech["high_clip"]=tech.measurement_id.astype(str).isin(high)
    d=d.merge(tech[["photo_id","near_clip_fraction","high_clip"]],on="photo_id",how="left",validate="one_to_one")
    d=d.loc[d.near_clip_fraction.notna() & ~d.high_clip.fillna(False)].copy()

    outcome=d.groupby("species").agg(
        n_img=("white","size"), n_white=("white","sum"), mean_near=("near_clip_fraction","mean")
    ).reset_index()
    outcome["white_prop"]=outcome.n_white/outcome.n_img

    p=pd.read_csv(a.pollinator_records)
    u=p.dropna(subset=["plant_species","pollinator_name","guild"]).drop_duplicates(
        ["plant_species","pollinator_name","guild"]
    )
    cnt=u.groupby(["plant_species","guild"]).size().unstack(fill_value=0)
    for guild in ["bee","lepidoptera"]:
        if guild not in cnt: cnt[guild]=0
    cnt["total_taxa"]=cnt.sum(axis=1)
    cnt["bee_frac"]=cnt.bee/cnt.total_taxa.replace(0,np.nan)
    cnt["lep_frac"]=cnt.lepidoptera/cnt.total_taxa.replace(0,np.nan)
    sph=p.assign(sph=p.pollinator_path.fillna("").str.contains("Sphingidae",case=False)).groupby("plant_species").sph.any()
    cov=pd.read_csv(a.coverage)

    q=cnt.join(sph.rename("has_sphingidae"),how="outer").reset_index().rename(columns={"plant_species":"species"})
    q=q.merge(cov[["species","all_four_queries_valid","any_query_truncated","n_distinct_pollinator_taxa"]],
              on="species",how="outer")
    m=outcome.merge(q,on="species",how="inner")
    m=m.loc[
        m.all_four_queries_valid.eq(True)
        & m.any_query_truncated.eq(False)
        & (m.n_distinct_pollinator_taxa>=10)
        & (m.n_img>=20)
    ].copy()

    for c in ["bee_frac","lep_frac","mean_near","n_distinct_pollinator_taxa"]:
        m["z_"+c]=z(m[c])

    rows=[]
    for pred,direction in [("has_sphingidae",1),("z_bee_frac",-1),("z_lep_frac",1)]:
        q=m.dropna(subset=[pred,"z_mean_near","z_n_distinct_pollinator_taxa"]).copy()
        X=pd.DataFrame({"const":1.0,"pred":q[pred].astype(float),
                        "near":q.z_mean_near,"effort":q.z_n_distinct_pollinator_taxa},index=q.index)
        fit=sm.GLM(q.white_prop,X,family=sm.families.Binomial(),freq_weights=q.n_img).fit(
            cov_type="cluster",cov_kwds={"groups":q.species})
        bb=float(fit.params["pred"]); se=float(fit.bse["pred"])
        q["logit_white"]=np.log((q.n_white+0.5)/(q.n_img-q.n_white+0.5))
        X2=sm.add_constant(pd.DataFrame({"pred":q[pred].astype(float),
                                        "near":q.z_mean_near,"effort":q.z_n_distinct_pollinator_taxa}))
        fit2=sm.OLS(q.logit_white,X2).fit(cov_type="HC3")
        rows.append({
            "predictor":pred,"prediction_direction":direction,"n_species":int(len(q)),
            "beta":bb,"OR":float(np.exp(bb)),
            "ci_low":float(np.exp(bb-1.96*se)),"ci_high":float(np.exp(bb+1.96*se)),
            "clustered_p":float(fit.pvalues["pred"]),
            "species_equal_beta":float(fit2.params["pred"]),
            "species_equal_p":float(fit2.pvalues["pred"]),
        })

    result={
        "schema":"fcp_white_pollinator_exploratory_v1",
        "status":"complete",
        "role":"post_outcome_exploratory_pollinator_screen",
        "strict_species":int(len(m)),
        "hard_nonclaims":[
            "not confirmatory because flower-colour outcomes were already open",
            "live GloBI query output is not a stable final inference source",
            "non-detection is not pollinator absence",
            "does not establish pigmented-to-white transition direction",
            "does not establish pollinator causation"
        ],
        "results":rows,
    }
    pd.DataFrame(rows).to_csv(out/"models.csv",index=False)
    (out/"result.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))

if __name__=="__main__":
    main()

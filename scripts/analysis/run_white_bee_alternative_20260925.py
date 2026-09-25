#!/usr/bin/env python3
from __future__ import annotations
import argparse, hashlib, json, math, warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
import statsmodels.api as sm
from scipy.stats import norm

MORPHS={"white","yellow_orange","red_pink","blue_purple"}

def truthy(s):
    if s.dtype == bool:
        return s.fillna(False)
    return s.astype(str).str.strip().str.lower().isin({"true","1","yes","y"})

def z(s):
    x=pd.to_numeric(s,errors="coerce").astype(float)
    sd=float(x.std(ddof=0))
    if not np.isfinite(sd) or sd<=0:
        return pd.Series(np.nan,index=x.index)
    return (x-float(x.mean()))/sd

def sample_raster(path,lon,lat):
    with rasterio.open(path) as src:
        vals=np.array([v[0] for v in src.sample(list(zip(lon.astype(float),lat.astype(float))))],dtype=float)
        if src.nodata is not None:
            vals[np.isclose(vals,float(src.nodata),equal_nan=False)]=np.nan
        return vals

def derive_fcp(measured,technical,join_key,high_clip,bio5):
    use=["photo_id","species","morph","global_classifiable","latitude","longitude"]
    d=pd.read_csv(measured,usecols=use,low_memory=False)
    d=d[truthy(d.global_classifiable)&d.morph.isin(MORPHS)].copy()
    tech=pd.read_csv(technical,compression="gzip",dtype={"measurement_id":str})
    key=pd.read_csv(join_key,dtype={"measurement_id":str})
    high=pd.read_csv(high_clip,dtype={"measurement_id":str})
    t=tech.merge(key,on="measurement_id",how="left",validate="one_to_one")
    highset=set(high.measurement_id.astype(str))
    t["high_clip"]=t.measurement_id.astype(str).isin(highset)
    d=d.merge(t[["photo_id","near_clip_fraction","high_clip"]],on="photo_id",how="left",validate="one_to_one")
    d=d[d.near_clip_fraction.notna() & ~d.high_clip.fillna(False)].copy()
    d["BIO5"]=sample_raster(Path(bio5),d.longitude,d.latitude)
    d=d.dropna(subset=["BIO5"])
    d["white"]=(d.morph=="white").astype(int)
    s=d.groupby("species").agg(
        n_classifiable=("white","size"),
        n_white=("white","sum"),
        mean_BIO5=("BIO5","mean")
    ).reset_index()
    s=s[s.n_classifiable>=40].copy()
    s["n_nonwhite"]=s.n_classifiable-s.n_white
    s["white_fraction"]=s.n_white/s.n_classifiable
    s["empirical_logit_white"]=np.log((s.n_white+0.5)/(s.n_nonwhite+0.5))
    return s

def scan_globi(path,focal,chunksize=150000):
    header=pd.read_csv(path,nrows=0)
    required={"plant_species","plant_family"}
    missing=required-set(header.columns)
    if missing:
        raise SystemExit(f"GloBI missing required columns: {sorted(missing)}")
    # Zenodo v3.1 renamed the bee-name field relative to the authors' internal
    # analysis scripts. Resolve only an explicitly bee/species semantic column;
    # this is schema adaptation, not biological model selection.
    bee_exact=[
      "scientificName","bee_species","bee_scientific_name","bee_species_name",
      "bee_name","bee_taxon_name","targetTaxonName"
    ]
    bee_col=next((x for x in bee_exact if x in header.columns),None)
    if bee_col is None:
        semantic=[
          x for x in header.columns
          if "bee" in x.lower() and ("species" in x.lower() or "scientific" in x.lower() or "taxon" in x.lower())
        ]
        if len(semantic)==1:
            bee_col=semantic[0]
    if bee_col is None:
        raise SystemExit(f"Could not identify unique bee species column. GloBI columns: {list(header.columns)}")
    n_records=defaultdict(int)
    bees=defaultdict(set)
    families=defaultdict(set)
    matched_rows=0
    for chunk in pd.read_csv(path,usecols=["plant_species",bee_col,"plant_family"],chunksize=chunksize,low_memory=False):
        chunk["plant_species"]=chunk.plant_species.astype(str).str.strip()
        q=chunk[chunk.plant_species.isin(focal)].copy()
        if q.empty:
            continue
        q=q[q[bee_col].notna()]
        matched_rows += len(q)
        for sp,g in q.groupby("plant_species",sort=False):
            n_records[sp]+=len(g)
            bees[sp].update(x for x in g[bee_col].astype(str).str.strip() if x and x.lower()!="nan")
            families[sp].update(x for x in g.plant_family.dropna().astype(str).str.strip() if x and x.lower()!="nan")
    rows=[]
    for sp in sorted(focal):
        if sp not in n_records:
            continue
        fam=sorted(families.get(sp,set()))
        rows.append({
            "species":sp,
            "n_interaction_records":int(n_records[sp]),
            "bee_species_richness":int(len(bees.get(sp,set()))),
            "n_plant_families_reported":int(len(fam)),
            "plant_family":fam[0] if len(fam)==1 else np.nan,
        })
    return pd.DataFrame(rows),matched_rows

def cluster_ols(df,formula_cols):
    y=df.empirical_logit_white.astype(float)
    X=sm.add_constant(df[formula_cols].astype(float),has_constant="add")
    fit=sm.OLS(y,X).fit(cov_type="cluster",cov_kwds={"groups":df.plant_family.astype(str),"use_correction":True})
    out={"n_species":int(len(df)),"n_families":int(df.plant_family.nunique()),"r2":float(fit.rsquared)}
    for col in formula_cols:
        b=float(fit.params[col]); se=float(fit.bse[col]); p2=float(fit.pvalues[col]); zz=b/se
        out[col]={
          "beta":b,"se":se,"z":zz,"two_sided_p":p2,
          "ci_low":float(b-1.96*se),"ci_high":float(b+1.96*se)
        }
    return out

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--measured",required=True)
    ap.add_argument("--technical-table",required=True)
    ap.add_argument("--join-key",required=True)
    ap.add_argument("--high-clip-ids",required=True)
    ap.add_argument("--bio5",required=True)
    ap.add_argument("--globi",required=True)
    ap.add_argument("--outdir",required=True)
    a=ap.parse_args()
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)

    fcp=derive_fcp(a.measured,a.technical_table,a.join_key,a.high_clip_ids,a.bio5)
    globi,matched_rows=scan_globi(a.globi,set(fcp.species.astype(str)))
    joined=fcp.merge(globi,on="species",how="left",validate="one_to_one")
    joined["globi_covered"]=joined.n_interaction_records.notna()
    eligible=joined[
      joined.n_interaction_records.ge(5) &
      joined.bee_species_richness.ge(2) &
      joined.plant_family.notna() &
      joined.n_plant_families_reported.eq(1)
    ].copy()

    endpoint_low=int((eligible.white_fraction<=0.10).sum())
    endpoint_high=int((eligible.white_fraction>=0.90).sum())
    coverage={
      "fcp_species_after_technical_filters":int(len(fcp)),
      "exact_globi_species_matches":int(joined.globi_covered.sum()),
      "analysis_species":int(len(eligible)),
      "analysis_families":int(eligible.plant_family.nunique()),
      "white_fraction_le_0_10":endpoint_low,
      "white_fraction_ge_0_90":endpoint_high,
      "globi_focal_matched_rows":int(matched_rows)
    }
    gate=bool(len(eligible)>=100 and eligible.plant_family.nunique()>=20 and endpoint_low>=20 and endpoint_high>=20)
    joined.to_csv(out/"white_bee_coverage_all_fcp_species.csv",index=False)
    eligible.to_csv(out/"white_bee_analysis_species.csv",index=False)
    if not gate:
        result={
          "schema":"fcp_white_bee_alternative_v1","status":"not_estimable",
          "coverage":coverage,
          "reason":"predeclared coverage gate failed",
          "missing_is_zero":False
        }
        (out/"result.json").write_text(json.dumps(result,indent=2)+"\n")
        print(json.dumps(result,indent=2)); return

    eligible["z_log_bee_richness"]=z(np.log1p(eligible.bee_species_richness))
    eligible["z_log_records"]=z(np.log1p(eligible.n_interaction_records))
    eligible["z_BIO5"]=z(eligible.mean_BIO5)
    primary=cluster_ols(eligible,["z_log_bee_richness","z_log_records","z_BIO5"])
    bee=primary["z_log_bee_richness"]
    bee["one_sided_negative_p"]=float(norm.cdf(bee["z"]))
    support=bool(bee["beta"]<0 and bee["one_sided_negative_p"]<0.05 and bee["ci_high"]<0)

    # Non-rescuing effort-residual sensitivity.
    Xeff=sm.add_constant(eligible[["z_log_records"]],has_constant="add")
    eff=sm.OLS(eligible.z_log_bee_richness,Xeff).fit()
    eligible["bee_breadth_effort_residual"]=eff.resid
    residual=cluster_ols(eligible,["bee_breadth_effort_residual","z_BIO5"])
    residual_bee=residual["bee_breadth_effort_residual"]
    residual_bee["one_sided_negative_p"]=float(norm.cdf(residual_bee["z"]))
    eligible.to_csv(out/"white_bee_analysis_species.csv",index=False)

    result={
      "schema":"fcp_white_bee_alternative_v1",
      "status":"BEE_BREADTH_SUPPORTED" if support else "BEE_BREADTH_NOT_SUPPORTED_UNDER_THIS_TEST",
      "coverage":coverage,
      "coverage_gate_pass":gate,
      "primary_model":"empirical_logit_white ~ z(log1p(bee_species_richness)) + z(log1p(interaction_records)) + z(mean_BIO5), family-clustered SE",
      "primary":primary,
      "effort_residual_sensitivity":residual,
      "bee_support_gate_pass":support,
      "globi_missing_semantics":"unmatched species are missing coverage, never zero bee interaction",
      "hard_nonclaims":[
        "does not test local pollinator preference",
        "does not establish bee dominance in pollinator community",
        "does not treat every GloBI interaction as successful pollination",
        "does not test non-bee pollinators",
        "does not establish causal selection",
        "does not establish evolutionary transition direction"
      ]
    }
    (out/"result.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))

if __name__=="__main__":
    main()

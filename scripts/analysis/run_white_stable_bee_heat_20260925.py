#!/usr/bin/env python3
from __future__ import annotations
import argparse, json
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from scipy.stats import wilcoxon
from sklearn.neighbors import BallTree
from statsmodels.discrete.conditional_models import ConditionalLogit

EARTH_KM=6371.0088
RADII=[50,100,250,500]
MORPHS={"white","yellow_orange","red_pink","blue_purple"}

def bool_series(s):
    if s.dtype==bool:
        return s.fillna(False)
    return s.astype(str).str.strip().str.lower().isin({"true","1","yes","y"})

def wz(df,col):
    g=df.groupby("species")[col]
    mu=g.transform("mean")
    sd=g.transform(lambda x:x.std(ddof=0)).replace(0,np.nan)
    return (df[col]-mu)/sd

def holm(ps):
    vals=[(i,float(p)) for i,p in enumerate(ps) if np.isfinite(p)]
    vals.sort(key=lambda z:z[1]); m=len(vals); out=[np.nan]*len(ps); prev=0
    for rank,(i,p) in enumerate(vals):
        v=min(1.0,(m-rank)*p); v=max(v,prev); out[i]=v; prev=v
    return out

def sample_raster(path, lon, lat):
    with rasterio.open(path) as src:
        vals=np.array([x[0] for x in src.sample(list(zip(lon.astype(float),lat.astype(float))))],float)
        if src.nodata is not None:
            vals[np.isclose(vals,src.nodata)]=np.nan
        return vals

def load_bees(path, plant_names):
    cols=["bee_species","bee_genus","plant_species","x","y","coordinated","plant_level"]
    pieces=[]
    for chunk in pd.read_csv(path,usecols=cols,chunksize=150000,low_memory=False,encoding="utf-8",encoding_errors="replace"):
        keep=chunk.plant_species.astype(str).isin(plant_names)
        keep &= chunk.plant_level.astype(str).eq("species_level")
        keep &= bool_series(chunk.coordinated)
        q=chunk.loc[keep].copy()
        q["x"]=pd.to_numeric(q.x,errors="coerce")
        q["y"]=pd.to_numeric(q.y,errors="coerce")
        q=q.dropna(subset=["bee_species","plant_species","x","y"])
        q=q.loc[q.x.between(-180,180)&q.y.between(-90,90)]
        if len(q):
            pieces.append(q)
    if not pieces:
        return pd.DataFrame(columns=cols)
    d=pd.concat(pieces,ignore_index=True)
    return d.drop_duplicates(["plant_species","bee_species","x","y"])

def build_trees(bees):
    out={}
    for sp,g in bees.groupby("plant_species"):
        coords=np.deg2rad(g[["y","x"]].to_numpy(float))
        out[sp]=(BallTree(coords,metric="haversine"),g.reset_index(drop=True))
    return out

def coverage_preflight(flowers,trees):
    rows=[]
    for sp,g in flowers.groupby("species"):
        if len(g)<20 or sp not in trees:
            continue
        tree,_=trees[sp]
        xy=np.deg2rad(g[["latitude","longitude"]].to_numpy(float))
        for radius in RADII:
            counts=tree.query_radius(xy,r=radius/EARTH_KM,count_only=True)
            rows.append({
                "species":sp,"radius_km":radius,"n_flower_rows":int(len(g)),
                "fraction_rows_ge3_bee_records":float((counts>=3).mean()),
                "median_local_bee_records":float(np.median(counts)),
                "qualified":bool((counts>=3).mean()>=0.5),
            })
    cov=pd.DataFrame(rows)
    summary=[]
    chosen=None
    for r in RADII:
        q=cov.loc[cov.radius_km.eq(r)]
        n=int(q.qualified.sum()) if len(q) else 0
        summary.append({"radius_km":r,"qualified_species":n,"species_evaluated":int(len(q))})
        if chosen is None and n>=100:
            chosen=r
    return cov,summary,chosen

def local_predictors(flowers,trees,radius):
    parts=[]
    for sp,g in flowers.groupby("species"):
        h=g.copy()
        h["local_bee_records"]=np.nan
        h["local_bee_richness"]=np.nan
        h["local_bombus_fraction"]=np.nan
        if sp in trees:
            tree,b=trees[sp]
            xy=np.deg2rad(h[["latitude","longitude"]].to_numpy(float))
            idxs=tree.query_radius(xy,r=radius/EARTH_KM,return_distance=False)
            rec=[]; rich=[]; bomb=[]
            for ids in idxs:
                if len(ids)<3:
                    rec.append(np.nan); rich.append(np.nan); bomb.append(np.nan); continue
                sub=b.iloc[np.asarray(ids,dtype=int)]
                taxa=sub[["bee_species","bee_genus"]].drop_duplicates("bee_species")
                rec.append(float(len(ids)))
                rich.append(float(len(taxa)))
                bomb.append(float(taxa.bee_genus.astype(str).eq("Bombus").mean()) if len(taxa) else np.nan)
            h["local_bee_records"]=rec
            h["local_bee_richness"]=rich
            h["local_bombus_fraction"]=bomb
        parts.append(h)
    return pd.concat(parts,ignore_index=True)

def cond_model(d,predictors,min_species=60):
    x=d.dropna(subset=["white","species"]+predictors).copy()
    st=x.groupby("species").white.agg(["sum","count"])
    good=set(st.index[(st["sum"]>0)&(st["sum"]<st["count"])])
    x=x.loc[x.species.isin(good)].copy()
    nsp=x.species.nunique()
    if nsp<min_species:
        return {"status":"not_estimable","n_species":int(nsp),"n_rows":int(len(x))}
    groups=x.species.astype("category").cat.codes.to_numpy()
    X=x[predictors].to_numpy(float)
    fit=ConditionalLogit(x.white.to_numpy(float),X,groups=groups).fit(method="bfgs",maxiter=500,disp=False)
    rr={"status":"complete","n_species":int(nsp),"n_rows":int(len(x))}
    for j,p in enumerate(predictors):
        beta=float(fit.params[j]); se=float(fit.bse[j])
        rr[p]={"beta":beta,"OR":float(np.exp(beta)),"ci_low":float(np.exp(beta-1.96*se)),
               "ci_high":float(np.exp(beta+1.96*se)),"p":float(fit.pvalues[j])}
    return rr

def species_delta(d,col):
    rows=[]
    for sp,g in d.dropna(subset=[col]).groupby("species"):
        nw=int(g.white.sum()); nn=int(len(g)-nw)
        if nw<5 or nn<5: continue
        delta=float(g.loc[g.white.eq(1),col].mean()-g.loc[g.white.eq(0),col].mean())
        rows.append((sp,delta,nw,nn))
    vals=np.array([r[1] for r in rows],float)
    if len(vals)<20:
        return {"n_species":len(vals),"status":"not_estimable"},rows
    p=float(wilcoxon(vals,zero_method="wilcox").pvalue)
    return {"n_species":len(vals),"median_delta":float(np.median(vals)),
            "mean_delta":float(np.mean(vals)),"positive_fraction":float((vals>0).mean()),
            "wilcoxon_two_sided_p":p,"status":"complete"},rows

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--measured",required=True)
    ap.add_argument("--technical-table",required=True)
    ap.add_argument("--join-key",required=True)
    ap.add_argument("--high-clip-ids",required=True)
    ap.add_argument("--bee-curated",required=True)
    ap.add_argument("--bio5",required=True)
    ap.add_argument("--outdir",required=True)
    a=ap.parse_args()
    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)

    # Outcome-blind preflight: never read morph/classifiability here.
    f0=pd.read_csv(a.measured,usecols=["photo_id","species","latitude","longitude"])
    for c in ["latitude","longitude"]:
        f0[c]=pd.to_numeric(f0[c],errors="coerce")
    f0=f0.dropna(subset=["species","latitude","longitude"])
    plants=set(f0.species.astype(str).unique())
    bees=load_bees(a.bee_curated,plants)
    trees=build_trees(bees)
    cov,summary,radius=coverage_preflight(f0,trees)
    cov.to_csv(out/"stable_bee_spatial_preflight.csv",index=False)
    preflight={
        "schema":"fcp_white_stable_bee_preflight_v1",
        "flower_species":int(f0.species.nunique()),
        "bee_matched_species":int(bees.plant_species.nunique()),
        "filtered_bee_records":int(len(bees)),
        "candidate_radii":RADII,
        "radius_summary":summary,
        "chosen_radius_km":radius,
        "rule":"smallest radius with >=100 species having >=20 flower coordinates and >=50% flower coordinates with >=3 same-plant local bee interaction records",
        "morph_opened":False,
    }
    (out/"preflight.json").write_text(json.dumps(preflight,indent=2)+"\n")
    if radius is None:
        (out/"result.json").write_text(json.dumps({
            "schema":"fcp_white_stable_bee_heat_v1","status":"NOT_ESTIMABLE_STABLE_BEE_COVERAGE",
            "preflight":preflight},indent=2)+"\n")
        print(json.dumps(preflight,indent=2)); return

    local=local_predictors(f0,trees,radius)

    # Only after radius and local predictors are frozen do we open colour outcomes.
    y=pd.read_csv(a.measured,usecols=["photo_id","morph","global_classifiable"])
    d=local.merge(y,on="photo_id",how="left",validate="one_to_one")
    d=d.loc[bool_series(d.global_classifiable)&d.morph.isin(MORPHS)].copy()
    d["white"]=(d.morph=="white").astype(int)

    tech=pd.read_csv(a.technical_table,compression="gzip",dtype={"measurement_id":str})
    join=pd.read_csv(a.join_key,dtype={"measurement_id":str})
    high=set(pd.read_csv(a.high_clip_ids,dtype={"measurement_id":str}).measurement_id.astype(str))
    tech=tech.merge(join,on="measurement_id",how="left",validate="one_to_one")
    tech["high_clip"]=tech.measurement_id.astype(str).isin(high)
    d=d.merge(tech[["photo_id","near_clip_fraction","high_clip"]],on="photo_id",how="left",validate="one_to_one")
    d["bio5"]=sample_raster(a.bio5,d.longitude,d.latitude)
    d=d.loc[d.near_clip_fraction.notna() & ~d.high_clip.fillna(False)
            & d.local_bee_richness.notna() & d.bio5.notna()].copy()
    d["log_local_bee_richness"]=np.log1p(d.local_bee_richness)
    for c in ["bio5","log_local_bee_richness","local_bombus_fraction","near_clip_fraction"]:
        d["wz_"+c]=wz(d,c)

    m1=cond_model(d,["wz_bio5","wz_log_local_bee_richness","wz_near_clip_fraction"])
    m2=cond_model(d,["wz_bio5","wz_local_bombus_fraction","wz_log_local_bee_richness","wz_near_clip_fraction"])

    audit_rich,rows_rich=species_delta(d,"wz_log_local_bee_richness")
    audit_bomb,rows_bomb=species_delta(d,"wz_local_bombus_fraction")
    poll_ps=[
        m1.get("wz_log_local_bee_richness",{}).get("p",np.nan),
        m2.get("wz_local_bombus_fraction",{}).get("p",np.nan),
    ]
    adj=holm(poll_ps)
    if m1.get("status")=="complete": m1["wz_log_local_bee_richness"]["holm_pollinator_p"]=adj[0]
    if m2.get("status")=="complete": m2["wz_local_bombus_fraction"]["holm_pollinator_p"]=adj[1]

    pd.DataFrame(
        [{"species":sp,"predictor":"wz_log_local_bee_richness","delta":de,"n_white":nw,"n_nonwhite":nn}
         for sp,de,nw,nn in rows_rich]+
        [{"species":sp,"predictor":"wz_local_bombus_fraction","delta":de,"n_white":nw,"n_nonwhite":nn}
         for sp,de,nw,nn in rows_bomb]
    ).to_csv(out/"species_local_bee_deltas.csv",index=False)

    result={
        "schema":"fcp_white_stable_bee_heat_v1",
        "status":"complete",
        "role":"pre_specified_joint_heat_bee_discriminator_after_white_outcome_already_known",
        "preflight":preflight,
        "biological_rows":int(len(d)),
        "biological_species":int(d.species.nunique()),
        "model_B1_heat_plus_local_bee_richness":m1,
        "model_B2_heat_plus_bombus_fraction_plus_richness":m2,
        "species_audit_local_bee_richness":audit_rich,
        "species_audit_local_bombus_fraction":audit_bomb,
        "hard_nonclaims":[
            "not untouched confirmation because white outcomes and prior BIO5 result already existed",
            "interaction record density is not pollination effectiveness",
            "does not measure non-bee pollinator guilds",
            "does not establish pigmented-to-white transition direction",
            "does not establish causal selection"
        ]
    }
    (out/"result.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))

if __name__=="__main__":
    main()

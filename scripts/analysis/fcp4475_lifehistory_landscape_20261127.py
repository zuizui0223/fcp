#!/usr/bin/env python3
"""Full-photo macroecology: life history, woodiness and self-fertilization.

Use authenticated existing GIFT source join, not an independently selected
photo panel. Compare amount of photographic colour variation and spatial
allocation. Treat missing/variable traits as unobserved for a fixed contrast.
No observer splitting or new molecular inference.
"""
import argparse
import hashlib
import io
import json
import zipfile
from pathlib import Path
import numpy as np
import pandas as pd

ART_SHA="a3398a827f9d8b5a3478c3096f94d15752d858cfd95c62be2fa4e48c0b45adf4"
GIFT_SHA="087755d43dbf87de57049c9a065d7019df2bf3f3ddd36fcbb3e795d90eebed9f"
TRAITS={"lifecycle":("annual","perennial"),"woodiness":("woody","non-woody"),"self_fertilization":("present","absent")}
BOOT=3999
SEED=20261127
# 3 traits x 6 primary contrasts: diversity, white fraction, geographic
# allocation (50,250km), and thermal-landscape moderation (50,250km).
FAMILY=18
BASE=["log_env_range_q90_km","log_colour_range_q90_km","env_n","colour_valid_n",
      "mean_temp","mean_precip","abs_latitude","log_thermal_breadth","log_moisture_breadth",
      "temp_seasonality","precip_seasonality"]

def verify(path,sha):
    b=path.read_bytes()
    assert hashlib.sha256(b).hexdigest()==sha, f"source digest mismatch: {path}"
    return b

def load(artifact,gift):
    with zipfile.ZipFile(io.BytesIO(verify(artifact,ART_SHA))) as z:
        names=[n for n in z.namelist() if n.rsplit("/",1)[-1]=="species_climate_colour_landscape.csv.gz"]
        assert len(names)==1
        a=pd.read_csv(io.BytesIO(z.read(names[0])),compression="gzip")
    b=pd.read_csv(io.BytesIO(verify(gift,GIFT_SHA)))
    assert len(a)==len(b)==4475
    assert a.inat_taxon_id.is_unique and b.inat_taxon_id.is_unique
    assert a.source_n.sum()==447500 and a.classifiable_n.sum()==227107
    d=a.merge(b[["inat_taxon_id","species","source_cohort","name_match_status",*TRAITS]],
              on="inat_taxon_id",validate="one_to_one",suffixes=("","_gift"))
    assert len(d)==4475 and d.species.eq(d.species_gift).all()
    assert d.cohort.eq(d.source_cohort).all()
    assert d.name_match_status.eq("unique_exact_name").sum()==4334
    for t in TRAITS:
        assert d.loc[~d.name_match_status.eq("unique_exact_name"),t].isna().all()
    for k in ("env_range_q90_km","colour_range_q90_km","thermal_breadth","moisture_breadth"):
        d["log_"+k]=np.log1p(d[k])
    return d

def fit(d,outcome,trait,level,base,interaction=None,within=False,seed=SEED):
    cols=[outcome,trait,*base]
    x=d.loc[d[trait].isin(level)].dropna(subset=cols).copy()
    if within:
        # Only genera that have both trait levels provide identification.
        x=x.loc[x.groupby("genus")[trait].transform("nunique")==2].copy()
    counts=x[trait].value_counts()
    genera=x.groupby(trait).genus.nunique()
    n=len(x);G=x.genus.nunique()
    res={"n":n,"genera":int(G),"level_counts":{s:int(counts.get(s,0)) for s in level},
         "level_genera":{s:int(genera.get(s,0)) for s in level},
         "within_genus":within,"contrast":f"{level[0]} minus {level[1]}"}
    if min(res["level_counts"].values())<30 or min(res["level_genera"].values())<15 or G<30:
        return dict(res,status="insufficient_support")
    A=x[base].to_numpy(float)
    sd=A.std(axis=0)
    if not np.isfinite(A).all() or (sd<=1e-12).any():
        return dict(res,status="constant_or_invalid_predictor")
    A=(A-A.mean(axis=0))/sd
    T=x[trait].eq(level[0]).to_numpy(float)
    C=pd.get_dummies(x.cohort,drop_first=True,dtype=float).to_numpy()
    X=np.column_stack([T,A,C])
    target=[0]
    labels=["trait_difference"]
    if interaction:
        # Trait contrast evaluated at mean thermal partition, and difference
        # in thermal-partition slope per one pooled panel SD.
        X=np.column_stack([X,T*A[:,base.index(interaction)]])
        target.append(X.shape[1]-1)
        labels.append("trait_x_thermal")
    y=x[outcome].to_numpy(float)
    codes,_=pd.factorize(x.genus,sort=True)
    if within:
        # Frisch-Waugh-Lovell transformation for genus intercepts.
        cnt=np.bincount(codes,minlength=G)
        sums=np.zeros((G,X.shape[1]));sy=np.zeros(G)
        np.add.at(sums,codes,X);np.add.at(sy,codes,y)
        X=X-sums[codes]/cnt[codes,None]
        y=y-sy[codes]/cnt[codes]
        keep=np.std(X,axis=0)>1e-10
        if not keep[target].all():
            return dict(res,status="within_genus_target_constant")
        target=[int(np.sum(keep[:i])) for i in target]
        X=X[:,keep]
    else:
        X=np.column_stack([np.ones(n),X])
        target=[i+1 for i in target]
    p=X.shape[1]
    if np.linalg.matrix_rank(X)!=p:
        return dict(res,status="rank_deficient")
    beta=np.linalg.lstsq(X,y,rcond=None)[0]
    gram=np.zeros((G,p,p));rhs=np.zeros((G,p))
    np.add.at(gram,codes,np.einsum("ni,nj->nij",X,X))
    np.add.at(rhs,codes,X*y[:,None])
    rng=np.random.default_rng(seed)
    draws=[]
    for _ in range(BOOT):
        w=np.bincount(rng.integers(G,size=G),minlength=G)
        a=np.einsum("g,gij->ij",w,gram)
        if np.linalg.matrix_rank(a)==p:
            draws.append(np.linalg.solve(a,w@rhs)[target])
    if len(draws)<.9*BOOT:
        return dict(res,status="unstable_bootstrap",valid_boot=len(draws))
    bs=np.asarray(draws)
    ci=np.quantile(bs,[.025,.975],axis=0)
    sim=np.quantile(bs,[.025/FAMILY,1-.025/FAMILY],axis=0)
    return dict(res,status="estimable",valid_boot=len(bs),
                effects={name:{"coefficient":float(beta[k]),"marginal95":ci[:,i].tolist(),
                               "bonferroni18_approx95":sim[:,i].tolist()}
                         for i,(name,k) in enumerate(zip(labels,target))})

def analyse(d):
    eligible=d.loc[(d.env_n>=20)&(d.env_cells>=5)&(d.colour_valid_n>=20)].copy()
    spatial=eligible.loc[eligible.repeated5.eq(True)].dropna(subset=[
        f"{k}_{r}" for r in (50,250) for k in ("thermal_partition","moisture_partition","colour_partition")
    ]).copy()
    assert len(eligible)==4333 and len(spatial)==1258
    out={}
    for i,(trait,levels) in enumerate(TRAITS.items()):
        q={}
        for j,outcome in enumerate(("colour_D","white_fraction")):
            q[outcome]={mode:fit(eligible,outcome,trait,levels,BASE,within=(mode=="within_genus"),seed=SEED+i*1000+j*10+int(mode=="within_genus"))
                        for mode in ("pooled","within_genus")}
        for j,r in enumerate((50,250)):
            thermal=f"thermal_partition_{r}"
            base=BASE+["colour_D",f"env_near_fraction_{r}",f"colour_near_fraction_{r}",thermal,f"moisture_partition_{r}"]
            q[f"colour_partition_{r}"]={mode:fit(spatial,f"colour_partition_{r}",trait,levels,base,interaction=thermal,
                         within=(mode=="within_genus"),seed=SEED+i*1000+100+j*10+int(mode=="within_genus"))
                         for mode in ("pooled","within_genus")}
        out[trait]=q
    return {"source_known_traits":{t:int(d[t].notna().sum()) for t in TRAITS},
            "eligible_diversity":len(eligible),"matched_spatial":len(spatial),"tests":out}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--landscape",required=True,type=Path)
    p.add_argument("--gift",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    d=load(a.landscape,a.gift)
    j={"source_landscape_artifact":15184882647,"source_landscape_sha256":ART_SHA,
       "gift_commit":"f005220272f84c2d6d696f984cfb284d60485375","gift_sha256":GIFT_SHA,
       "source_taxa":len(d),"bootstrap":BOOT,"primary_family":FAMILY,
       "scope":"Post hoc macroecological trait comparisons. Pooled regressions are primary; within-genus same-direction checks are sensitivity, not full phylogenetic correction. Species-equal, no observer splitting. GIFT self-fertilization is a recorded trait not a realized selfing rate; lifecycle and woodiness are broad traits. Photo-colour variation not inherited morphs. No causal or pigment-cost inference. Family intervals approximate Bonferroni for 18 pooled coefficients; within-genus sensitivity separately descriptive.",
       **analyse(d)}
    a.out.mkdir(parents=True,exist_ok=True)
    d.to_csv(a.out/"fcp4475_lifehistory_landscape.csv.gz",index=False,compression={"method":"gzip","mtime":0},float_format="%.12g")
    (a.out/"results.json").write_text(json.dumps(j,indent=2,allow_nan=False)+"\n")
    print("FCP_LIFEHISTORY_LANDSCAPE_COMPLETE",json.dumps(j,allow_nan=False),flush=True)

if __name__=="__main__":
    main()

#!/usr/bin/env python3
"""Macroecology of within-species photographic colour and sampled climate landscapes.

Use every original photo location (including colour-unclassifiable records)
for the environmental template. Use all classifiable photographs for colour.
No observer division, no inferred genetic morphs or causal fitness effects.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd

from fcp_full4475_hypothesis_census_20261021 import (
    SOURCE_SHA, BIO_ARCHIVE_SHA, TOTAL_ROWS, TOTAL_TAXA, CLASSIFIABLE_ROWS,
    PHOTO_LABELS, load_original, sample_climate, all_coords_ok, distance_km,
)

RADII=(50,250)
MIN_PAIRS=15
MIN_COLOUR=20
MIN_ENV=20
MIN_CELLS=5
BOOT=1999
SEED=20261127
BIO=("BIO1","BIO5","BIO12","BIO14","BIO4","BIO15")
DOMAINS={"thermal":(0,1),"moisture":(2,3)}

def prepare_climate(raw, archive):
    # load_original authenticates the exact full original photo table.
    # Climate support deliberately includes the unclassifiable photographs.
    a=raw.copy()
    for k in ("latitude","longitude","positional_accuracy_m"):
        a[k]=pd.to_numeric(a[k],errors="coerce")
    a=sample_climate(a,archive)
    v=a[list(BIO)].to_numpy(float)
    good=all_coords_ok(a)&np.isfinite(v).all(axis=1)&(v[:,2:]>=0).all(axis=1)
    a=a.loc[good].copy()
    v=v[good]
    v[:,2:]=np.log1p(v[:,2:])
    mu=v.mean(axis=0);sd=v.std(axis=0)
    assert np.isfinite(mu).all() and (sd>1e-10).all()
    z=(v-mu)/sd
    for j,k in enumerate(BIO):
        a["z"+k]=z[:,j]
    return a,dict(columns=list(BIO),mean=mu.tolist(),sd=sd.tolist(),
                  transformation="BIO12,14,4,15 log1p; all axes standardized on all climate-complete original sampled photographs, not global raster background",
                  complete_photos=len(a))

def pairs(g):
    u,v=np.triu_indices(len(g),1)
    d=distance_km(g.latitude.to_numpy(float),g.longitude.to_numpy(float),u,v)
    return u,v,d

def climate_metrics(g):
    out={"env_n":len(g),"env_cells":0}
    if len(g)<2:
        return out
    cells=np.column_stack([np.floor((90-g.latitude.to_numpy(float))*6),
                           np.floor((g.longitude.to_numpy(float)+180)*6)])
    out["env_cells"]=len(np.unique(cells,axis=0))
    out["abs_latitude"]=float(np.mean(np.abs(g.latitude)))
    out["mean_temp"]=float(g[["zBIO1","zBIO5"]].to_numpy().mean())
    out["mean_precip"]=float(g[["zBIO12","zBIO14"]].to_numpy().mean())
    out["temp_seasonality"]=float(g.zBIO4.mean())
    out["precip_seasonality"]=float(g.zBIO15.mean())
    u,v,d=pairs(g)
    out["env_range_q90_km"]=float(np.quantile(d,.9))
    z=g[["z"+b for b in BIO]].to_numpy(float)
    for r in RADII:
        near=d<=r
        out[f"env_near_n_{r}"]=int(near.sum())
        out[f"env_far_n_{r}"]=int((~near).sum())
        out[f"env_near_fraction_{r}"]=float(near.mean())
    for name,ix in DOMAINS.items():
        delta=np.mean((z[u][:,ix]-z[v][:,ix])**2,axis=1)
        global_d=float(delta.mean())
        out[name+"_breadth"]=float(np.sqrt(global_d))
        for r in RADII:
            near=d<=r
            if near.sum()>=MIN_PAIRS and (~near).sum()>=MIN_PAIRS and global_d>1e-10:
                out[f"{name}_partition_{r}"]=1-float(delta[near].mean())/global_d
    return out

def colour_metrics(g):
    s=g.loc[all_coords_ok(g)].reset_index(drop=True)
    out={"classifiable_n":len(g),"colour_valid_n":len(s),
         "repeated5":bool((s.morph.value_counts()>=5).sum()>=2)}
    if len(s)<2:
        return out
    u,v,d=pairs(s)
    c=pd.Categorical(s.morph,categories=PHOTO_LABELS).codes
    assert (c>=0).all()
    diff=(c[u]!=c[v]).astype(float)
    pal=np.stack(s.palette_norm9.to_numpy()).astype(float)
    assert (pal>=0).all() and np.allclose(pal.sum(axis=1),1)
    # Squared Hellinger dissimilarity, continuous bounded nine-hue measure.
    sq=np.sqrt(pal)
    hell=.5*np.sum((sq[u]-sq[v])**2,axis=1)
    out.update(colour_D=float(diff.mean()),hellinger_D=float(hell.mean()),
               white_fraction=float(pal[:,0].mean()),
               colour_range_q90_km=float(np.quantile(d,.9)))
    for r in RADII:
        near=d<=r
        out[f"colour_near_n_{r}"]=int(near.sum())
        out[f"colour_far_n_{r}"]=int((~near).sum())
        out[f"colour_near_fraction_{r}"]=float(near.mean())
        if near.sum()>=MIN_PAIRS and (~near).sum()>=MIN_PAIRS:
            out[f"colour_partition_{r}"]=float(diff.mean()-diff[near].mean())
            out[f"hellinger_partition_{r}"]=float(hell.mean()-hell[near].mean())
    return out

def ledger(raw,cl,env):
    cgroups={int(t):g for t,g in cl.groupby("inat_taxon_id",sort=False)}
    egroups={int(t):g for t,g in env.groupby("inat_taxon_id",sort=False)}
    rows=[]
    for tid,g in raw.groupby("inat_taxon_id",sort=True):
        row=dict(inat_taxon_id=int(tid),species=str(g.species.iloc[0]),
                 genus=str(g.species.iloc[0]).split()[0],
                 cohort=str(g.fcp_source_cohort.iloc[0]),source_n=len(g))
        row.update(climate_metrics(egroups.get(int(tid),env.iloc[:0])))
        row.update(colour_metrics(cgroups.get(int(tid),cl.iloc[:0])))
        rows.append(row)
    d=pd.DataFrame(rows)
    assert len(d)==TOTAL_TAXA and d.source_n.sum()==TOTAL_ROWS
    assert d.classifiable_n.sum()==CLASSIFIABLE_ROWS
    return d

def fit(d,outcome,base,extra,seed):
    """Species equal OLS, all outcomes in original units, predictors per panel SD.
    Report within-sample incremental R2, not out-of-sample prediction.
    """
    x=d.dropna(subset=[outcome,*base,*extra]).copy()
    n=len(x)
    G=int(x.genus.nunique())
    result={"n":n,"genera":G,"outcome":outcome}
    if n<100 or G<50:
        return dict(result,status="insufficient_support")
    names=base+extra
    A=x[names].to_numpy(float)
    sd=A.std(axis=0)
    if not (np.isfinite(A).all() and (sd>1e-10).all()):
        return dict(result,status="constant_or_invalid_predictor")
    A=(A-A.mean(axis=0))/sd
    C=pd.get_dummies(x.cohort,drop_first=True,dtype=float).to_numpy()
    X0=np.column_stack([np.ones(n),A[:,:len(base)],C])
    X=np.column_stack([X0,A[:,len(base):]])
    y=x[outcome].to_numpy(float)
    p=X.shape[1]
    if np.linalg.matrix_rank(X)!=p:
        return dict(result,status="rank_deficient")
    beta=np.linalg.lstsq(X,y,rcond=None)[0]
    beta0=np.linalg.lstsq(X0,y,rcond=None)[0]
    sst=float(((y-y.mean())**2).sum())
    if sst<1e-12:
        return dict(result,status="constant_outcome")
    r2=1-float(((y-X@beta)**2).sum())/sst
    r20=1-float(((y-X0@beta0)**2).sum())/sst
    codes,_=pd.factorize(x.genus,sort=True)
    gram=np.zeros((G,p,p));rhs=np.zeros((G,p))
    np.add.at(gram,codes,np.einsum("ni,nj->nij",X,X))
    np.add.at(rhs,codes,X*y[:,None])
    rng=np.random.default_rng(seed)
    draws=[]
    for _ in range(BOOT):
        w=np.bincount(rng.integers(G,size=G),minlength=G)
        g=np.einsum("g,gij->ij",w,gram)
        if np.linalg.matrix_rank(g)==p:
            draws.append(np.linalg.solve(g,w@rhs))
    if len(draws)<.9*BOOT:
        return dict(result,status="unstable_bootstrap",valid_boot=len(draws))
    ci=np.quantile(draws,[.025,.975],axis=0).T
    indices=list(range(1,1+len(base)))+list(range(X0.shape[1],p))
    return dict(result,status="estimable",valid_boot=len(draws),r2_base=r20,
                r2_full=r2,incremental_r2=r2-r20,
                coefficient={k:{"per_predictor_sd":float(beta[i]),"genus_ci95":ci[i].tolist()}
                             for k,i in zip(names,indices)})

def analysis(d):
    d=d.copy()
    for k in ("env_range_q90_km","colour_range_q90_km","thermal_breadth","moisture_breadth"):
        d["log_"+k]=np.log1p(d[k])
    common_base=["log_env_range_q90_km","log_colour_range_q90_km","env_n",
                 "colour_valid_n","mean_temp","mean_precip","abs_latitude"]
    breadth=["log_thermal_breadth","log_moisture_breadth","temp_seasonality","precip_seasonality"]
    eligible=d.loc[(d.env_n>=MIN_ENV)&(d.env_cells>=MIN_CELLS)&
                   (d.colour_valid_n>=MIN_COLOUR)].copy()
    diversity=fit(eligible,"colour_D",common_base,breadth,SEED)
    p={}
    panels={}
    for i,r in enumerate(RADII):
        part=[f"thermal_partition_{r}",f"moisture_partition_{r}"]
        base=common_base+breadth+["colour_D",f"env_near_fraction_{r}",f"colour_near_fraction_{r}"]
        # Main population: repeated photographic categories, not prior
        # literature-defined morphs. All 4475 taxa remain in the ledger.
        q=eligible.loc[eligible.repeated5].dropna(subset=[f"colour_partition_{r}",*base,*part]).copy()
        panels[r]=q
        p[str(r)]={"categorical":fit(q,f"colour_partition_{r}",base,part,SEED+10+i),
                   "continuous":fit(q,f"hellinger_partition_{r}",
                     [k if k!="colour_D" else "hellinger_D" for k in base],part,SEED+20+i)}
    ids=set(panels[50].inat_taxon_id)&set(panels[250].inat_taxon_id)
    for i,r in enumerate(RADII):
        q=panels[r].loc[panels[r].inat_taxon_id.isin(ids)]
        base=common_base+breadth+["colour_D",f"env_near_fraction_{r}",f"colour_near_fraction_{r}"]
        p[str(r)]["matched_scales"]=fit(q,f"colour_partition_{r}",base,
                 [f"thermal_partition_{r}",f"moisture_partition_{r}"],SEED+30+i)
    return {"eligible_diversity_n":len(eligible),"repeated5_source_n":int(d.repeated5.sum()),
            "matched_scale_n":len(ids),"diversity":diversity,"spatial_partition":p}

def selftest():
    g=pd.DataFrame(dict(latitude=[0.,0.,0.,0.],longitude=[0.,.01,10.,10.01],
                       positional_accuracy_m=[10.]*4,morph=["white","white","red_pink","red_pink"]))
    g["palette_norm9"]=[np.eye(9)[i] for i in [0,0,3,3]]
    u,v,d=pairs(g)
    assert (d<50).sum()==2 and np.all(d[d>=50]>1000)
    g2=pd.concat([g]*10,ignore_index=True)
    cm=colour_metrics(g2)
    assert cm["colour_partition_50"]>0 and cm["hellinger_partition_50"]>0
    for k in BIO:
        g2["z"+k]=np.tile([0.,0.,1.,1.],10)
    em=climate_metrics(g2)
    assert em["thermal_partition_50"]==1 and em["moisture_partition_50"]==1
    print("FCP_CLIMATE_LANDSCAPE_SELFTEST_PASS")

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--photos",type=Path)
    p.add_argument("--worldclim-bio",type=Path)
    p.add_argument("--out",type=Path)
    p.add_argument("--self-test",action="store_true")
    a=p.parse_args()
    selftest()
    if a.self_test:
        return
    assert a.photos and a.worldclim_bio and a.out
    raw,cl=load_original(a.photos)
    env,scales=prepare_climate(raw,a.worldclim_bio)
    d=ledger(raw,cl,env)
    result={"source_photo_sha256":SOURCE_SHA,"source_bio_sha256":BIO_ARCHIVE_SHA,
            "source_photos":TOTAL_ROWS,"source_taxa":TOTAL_TAXA,"source_classified":CLASSIFIABLE_ROWS,
            "climate_scaling":scales,"radii_km":RADII,"genus_bootstrap":BOOT,
            "scope":"Cross-species exploratory associations among sampled climate landscapes and photo-colour diversity/spatial structure. All original photo locations for climate, all classifiable photos for colour. No observer split. Spatial partition is a scale-specific ratio of near vs range-wide environmental dissimilarity, not a niche preference or causal IBE. Seasonality is a climate normal, not observed interannual change. No true polymorphism prevalence, genetic inheritance, fitness or pigment costs.",
            **analysis(d)}
    a.out.mkdir(parents=True,exist_ok=True)
    d.to_csv(a.out/"species_climate_colour_landscape.csv.gz",index=False,
             compression={"method":"gzip","mtime":0},float_format="%.12g")
    (a.out/"results.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("FCP_CLIMATE_LANDSCAPE_COMPLETE",json.dumps(result,allow_nan=False),flush=True)

if __name__=="__main__":
    main()

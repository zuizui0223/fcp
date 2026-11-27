#!/usr/bin/env python3
"""Does thermal/moisture landscape organization predict white-specific colour sorting?

Post hoc macroecological analysis, no observer splitting. Exact decomposition of
four-label mismatch and nine-palette squared Hellinger dissimilarity.
Absolute component differences are not normalized effects; compare each
component relative to its own range-wide dissimilarity on identical taxa.
"""
import argparse
import hashlib
import io
import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

from fcp_full4475_hypothesis_census_20261021 import (
    SOURCE_SHA, TOTAL_ROWS, TOTAL_TAXA, CLASSIFIABLE_ROWS, load_original,
    all_coords_ok, PHOTO_LABELS,
)
from fcp4475_climate_landscape_organization_20261127 import pairs

ARTIFACT_SHA="a3398a827f9d8b5a3478c3096f94d15752d858cfd95c62be2fa4e48c0b45adf4"
RADII=(50,250)
MIN_PAIRS=15
MIN_BASE=.01
SEED=20261127
BOOT=1999

def load_landscape(path):
    raw=path.read_bytes()
    assert hashlib.sha256(raw).hexdigest()==ARTIFACT_SHA, "landscape ZIP hash mismatch"
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        f=[n for n in z.namelist() if n.rsplit("/",1)[-1]=="species_climate_colour_landscape.csv.gz"]
        assert len(f)==1
        d=pd.read_csv(io.BytesIO(z.read(f[0])),compression="gzip")
    assert len(d)==TOTAL_TAXA and d.inat_taxon_id.nunique()==TOTAL_TAXA
    assert d.source_n.sum()==TOTAL_ROWS and d.classifiable_n.sum()==CLASSIFIABLE_ROWS
    return d

def decompose(labels, p, u, v):
    # Coarse four-state categories: mismatch exactly partitions into
    # white/nonwhite pairs and two-different-nonwhite-category pairs.
    w=labels==0
    mismatch=(labels[u]!=labels[v]).astype(float)
    cat_white=(w[u]!=w[v]).astype(float)
    cat_chrom=mismatch-cat_white
    # Squared Hellinger identity:
    # H^2(p,q)=H^2(binary white/rest)+sqrt(nonwhite_p*nonwhite_q)*H^2(nonwhite conditional palettes).
    # This is defined even when either vector has zero nonwhite mass.
    sw=np.sqrt(p[:,0])
    sc=np.sqrt(np.clip(1-p[:,0],0,1))
    total=.5*np.sum((np.sqrt(p[u])-np.sqrt(p[v]))**2,axis=1)
    hw=.5*((sw[u]-sw[v])**2+(sc[u]-sc[v])**2)
    hc=total-hw
    assert np.min(hc)>=-1e-10
    hc=np.maximum(hc,0)
    assert np.allclose(mismatch,cat_white+cat_chrom,rtol=0,atol=1e-10)
    assert np.allclose(total,hw+hc,rtol=0,atol=1e-10)
    return {"categorical":(cat_white,cat_chrom),"continuous":(hw,hc)}

def measure(g):
    s=g.loc[all_coords_ok(g)].reset_index(drop=True)
    if len(s)<2:
        return {}
    u,v,d=pairs(s)
    labels=pd.Categorical(s.morph,categories=PHOTO_LABELS).codes
    p=np.stack(s.palette_norm9.to_numpy()).astype(float)
    assert (labels>=0).all() and np.allclose(p.sum(axis=1),1)
    out={}
    for typ,(white,chrom) in decompose(labels,p,u,v).items():
        for comp,a in (("white",white),("chromatic",chrom)):
            b=float(a.mean())
            out[f"{typ}_{comp}_baseline"]=b
            for r in RADII:
                near=d<=r
                if near.sum()>=MIN_PAIRS and (~near).sum()>=MIN_PAIRS:
                    delta=b-float(a[near].mean())
                    out[f"{typ}_{comp}_absolute_{r}"]=delta
                    if b>=MIN_BASE:
                        out[f"{typ}_{comp}_relative_{r}"]=delta/b
    return out

def build_ledger(raw,cl,land):
    rows=[]
    groups={int(t):g for t,g in cl.groupby("inat_taxon_id",sort=False)}
    for tid in sorted(raw.inat_taxon_id.unique()):
        rows.append({"inat_taxon_id":int(tid),**measure(groups.get(int(tid),cl.iloc[:0]))})
    d=land.merge(pd.DataFrame(rows),on="inat_taxon_id",validate="one_to_one")
    assert len(d)==TOTAL_TAXA
    for typ,old in (("categorical","colour"),("continuous","hellinger")):
        for r in RADII:
            a=d[f"{typ}_white_absolute_{r}"]+d[f"{typ}_chromatic_absolute_{r}"]
            b=d[f"{old}_partition_{r}"]
            valid=a.notna()&b.notna()
            assert valid.any() and np.allclose(a[valid],b[valid],atol=2e-10,rtol=1e-8)
    return d

def estimate(d,typ,r,metric,seed):
    w=f"{typ}_white_{metric}_{r}"
    c=f"{typ}_chromatic_{metric}_{r}"
    extra=[f"thermal_partition_{r}",f"moisture_partition_{r}"]
    base=["log_env_range_q90_km","log_colour_range_q90_km","env_n","colour_valid_n",
          "mean_temp","mean_precip","abs_latitude","log_thermal_breadth","log_moisture_breadth",
          "temp_seasonality","precip_seasonality",f"env_near_fraction_{r}",f"colour_near_fraction_{r}",
          f"{typ}_white_baseline",f"{typ}_chromatic_baseline"]
    names=base+extra
    x=d.dropna(subset=names+[w,c]).copy()
    n=len(x)
    G=int(x.genus.nunique())
    result={"n":n,"genera":G}
    if n<100 or G<50:
        return dict(result,status="insufficient_support")
    A=x[names].to_numpy(float)
    sd=A.std(axis=0)
    if (sd<=1e-10).any():
        return dict(result,status="constant_predictor")
    A=(A-A.mean(axis=0))/sd
    C=pd.get_dummies(x.cohort,drop_first=True,dtype=float).to_numpy()
    X=np.column_stack([np.ones(n),A,C])
    Y=x[[w,c]].to_numpy(float)
    Y=np.column_stack([Y,Y[:,0]-Y[:,1]])
    p=X.shape[1]
    if np.linalg.matrix_rank(X)!=p:
        return dict(result,status="rank_deficient")
    beta=np.linalg.lstsq(X,Y,rcond=None)[0]
    codes,_=pd.factorize(x.genus,sort=True)
    gram=np.zeros((G,p,p));rhs=np.zeros((G,p,3))
    np.add.at(gram,codes,np.einsum("ni,nj->nij",X,X))
    np.add.at(rhs,codes,X[:,:,None]*Y[:,None,:])
    rng=np.random.default_rng(seed)
    draws=[]
    for _ in range(BOOT):
        weights=np.bincount(rng.integers(G,size=G),minlength=G)
        g=np.einsum("g,gij->ij",weights,gram)
        if np.linalg.matrix_rank(g)==p:
            draws.append(np.linalg.solve(g,np.einsum("g,gij->ij",weights,rhs)))
    if len(draws)<.9*BOOT:
        return dict(result,status="unstable_bootstrap",valid_boot=len(draws))
    bs=np.asarray(draws)
    ci=np.quantile(bs,[.025,.975],axis=0)
    # Same-species paired differences: both component outcomes fitted
    # with precisely the same X and genus bootstrap draw.
    output={}
    for predictor in extra:
        j=1+names.index(predictor)
        output[predictor]={k:{"coefficient_per_predictor_sd":float(beta[j,i]),
                             "genus_ci95":ci[:,j,i].tolist()}
                          for i,k in enumerate(("white","chromatic","white_minus_chromatic"))}
    return dict(result,status="estimable",valid_boot=len(draws),mean_white=float(Y[:,0].mean()),
                mean_chromatic=float(Y[:,1].mean()),effects=output)

def analyse(d):
    for k in ("env_range_q90_km","colour_range_q90_km","thermal_breadth","moisture_breadth"):
        d["log_"+k]=np.log1p(d[k])
    eligible=d.loc[(d.env_n>=20)&(d.env_cells>=5)&(d.colour_valid_n>=20)&(d.repeated5==True)].copy()
    # The common set must support both components at both spatial scales;
    # it need not be the same across categorical and continuous definitions.
    out={}
    for i,typ in enumerate(("categorical","continuous")):
        required=[f"{typ}_{c}_relative_{r}" for c in ("white","chromatic") for r in RADII]
        required += [f"{c}_partition_{r}" for c in ("thermal","moisture") for r in RADII]
        q=eligible.dropna(subset=required)
        panels={}
        for j,r in enumerate(RADII):
            panels[str(r)]={m:estimate(q,typ,r,m,SEED+i*100+j*10+k)
                            for k,m in enumerate(("relative","absolute"))}
        out[typ]={"common_component_scale_taxa":len(q),"panels":panels}
    return out

def selftest():
    labels=np.array([0,1,2,0])
    p=np.array([[1,0,0],[0,1,0],[0,0,1],[.5,.5,0]])
    u,v=np.triu_indices(4,1)
    z=decompose(labels,p,u,v)
    assert np.array_equal(z["categorical"][0]+z["categorical"][1],(labels[u]!=labels[v]).astype(float))
    assert np.all(z["continuous"][1]>=0)
    assert z["continuous"][1][0]==0
    print("FCP_WHITE_CHROMATIC_LANDSCAPE_SELFTEST_PASS")

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--photos",type=Path)
    p.add_argument("--landscape",type=Path)
    p.add_argument("--out",type=Path)
    p.add_argument("--self-test",action="store_true")
    a=p.parse_args()
    selftest()
    if a.self_test:
        return
    assert a.photos and a.landscape and a.out
    raw,cl=load_original(a.photos)
    d=build_ledger(raw,cl,load_landscape(a.landscape))
    j={"source_photo_sha256":SOURCE_SHA,"landscape_artifact_id":15184882647,"landscape_zip_sha256":ARTIFACT_SHA,
       "source_photos":TOTAL_ROWS,"source_classified":CLASSIFIABLE_ROWS,"source_taxa":TOTAL_TAXA,
       "min_component_baseline":MIN_BASE,"radii_km":RADII,"genus_bootstrap":BOOT,
       "scope":"Exploratory macroecological comparative association using full source photos; no observer splitting. White/rest and chromatic component decomposition describes photographed palette geometry, not pigment biochemistry. Relative contrasts conditioned on positive baseline support; absolute contrasts evaluated on the same taxa. Marginal genus intervals, not independent tests or full phylogenetic model.",
       "results":analyse(d)}
    a.out.mkdir(parents=True,exist_ok=True)
    d.to_csv(a.out/"species_white_chromatic_landscape.csv.gz",index=False,compression={"method":"gzip","mtime":0},float_format="%.12g")
    (a.out/"results.json").write_text(json.dumps(j,indent=2,allow_nan=False)+"\n")
    print("FCP_WHITE_CHROMATIC_LANDSCAPE_COMPLETE",json.dumps(j,allow_nan=False),flush=True)

if __name__=="__main__":
    main()

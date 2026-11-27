#!/usr/bin/env python3
"""Separate growth form from life cycle in the full-photo macroecology.

Compare perennial woody vs perennial nonwoody, and annual vs perennial
nonwoody. The second contrast does not conflate annuality with woodiness.
Two photo-colour measures on identical species, and two spatial measures
at 250 km on identical species. No observer splitting. Conditional genus
wild-score uncertainty avoids singular random genus resampling of rare
nuisance covariates; fixed-design sensitivity is not a phylogenetic model.
"""
import argparse
import hashlib
import io
import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

ARTIFACT=15189521587
DIGEST="19e4daf3c7c1a73bf1ae5f5a98e0eb0c689408ff0cb4a737fc37329a2aa34239"
SEED=20261128
BOOT=9999
FAMILY=12
BASE=["log_env_range_q90_km","log_colour_range_q90_km","env_n","colour_valid_n",
      "mean_temp","mean_precip","abs_latitude","log_thermal_breadth","log_moisture_breadth",
      "temp_seasonality","precip_seasonality"]
R=250
THERMAL="thermal_partition_250"
GROUPS={
    "woodiness_within_perennials":("perennial_woody","perennial_nonwoody"),
    "annuality_within_nonwoody":("annual_nonwoody","perennial_nonwoody"),
}

def load(path):
    data=path.read_bytes()
    assert hashlib.sha256(data).hexdigest()==DIGEST, "life-history source artifact hash mismatch"
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        f=[n for n in z.namelist() if n.rsplit("/",1)[-1]=="fcp4475_lifehistory_landscape.csv.gz"]
        assert len(f)==1
        d=pd.read_csv(io.BytesIO(z.read(f[0])),compression="gzip")
    assert len(d)==4475 and d.inat_taxon_id.is_unique
    assert d.source_n.sum()==447500 and d.classifiable_n.sum()==227107
    d["strategy"]="other_or_unknown"
    for life,wood,label in [
        ("perennial","woody","perennial_woody"),
        ("perennial","non-woody","perennial_nonwoody"),
        ("annual","non-woody","annual_nonwoody")]:
        d.loc[d.lifecycle.eq(life)&d.woodiness.eq(wood),"strategy"]=label
    return d

def demean(A,codes,G):
    a=np.asarray(A,float)
    if a.ndim==1:
        return demean(a[:,None],codes,G)[:,0]
    sums=np.zeros((G,a.shape[1]))
    np.add.at(sums,codes,a)
    return a-sums[codes]/np.bincount(codes,minlength=G)[codes,None]

def orth(A):
    U,S,_=np.linalg.svd(A,full_matrices=False)
    if not len(S) or S[0]<1e-12:
        return np.empty((len(A),0))
    return U[:,S>S[0]*max(A.shape)*np.finfo(float).eps*10]

def wild_fit(x,levels,outcome,base,interaction,within,seed):
    if within:
        x=x.loc[x.groupby("genus").strategy.transform("nunique")==2].copy()
    n=len(x)
    G=x.genus.nunique()
    counts=x.strategy.value_counts()
    level_counts={k:int(counts.get(k,0)) for k in levels}
    result={"n":n,"genera":int(G),"counts":level_counts,"levels":list(levels),"within_genus":within}
    if min(level_counts.values())<30 or G<20:
        return dict(result,status="insufficient_support")
    A=x[base].to_numpy(float)
    sd=A.std(axis=0)
    assert np.isfinite(A).all()
    A=(A-A.mean(axis=0))/np.where(sd>1e-10,sd,1.)
    C=pd.get_dummies(x.cohort,drop_first=True,dtype=float).to_numpy()
    Z=np.column_stack([np.ones(n),A,C])
    t=x.strategy.eq(levels[0]).to_numpy(float)
    T=t[:,None]
    names=["strategy_difference"]
    if interaction:
        if sd[base.index(THERMAL)]<=1e-10:
            return dict(result,status="constant_thermal_partition")
        T=np.column_stack([t,t*A[:,base.index(THERMAL)]])
        names.append("strategy_x_thermal")
        result["thermal_partition_sd"]=float(sd[base.index(THERMAL)])
    Y=x[outcome].to_numpy(float)
    codes,_=pd.factorize(x.genus,sort=True)
    fe_h=np.zeros(n)
    if within:
        # Exact fixed-effect elimination; only mixed-strategy genera included.
        Z=demean(Z,codes,G);T=demean(T,codes,G);Y=demean(Y,codes,G)
        fe_h=1/np.bincount(codes,minlength=G)[codes]
    Q=orth(Z)
    W=T-Q@(Q.T@T)
    y=Y-Q@(Q.T@Y)
    if np.linalg.matrix_rank(W)<W.shape[1] or np.linalg.cond(W)>1e8:
        return dict(result,status="target_not_identifiable")
    inv=np.linalg.inv(W.T@W)
    beta=inv@(W.T@y)
    e=y-W@beta
    leverage=fe_h+np.sum(Q*Q,axis=1)+np.einsum("ni,ij,nj->n",W,inv,W)
    # HC2 leverage adjustment, clustered Rademacher multiplier on fixed design.
    # Clusters can span both strategy levels; their scores stay together.
    if np.max(leverage)>1+1e-7:
        raise ValueError("invalid projection leverage")
    e=e/np.sqrt(np.maximum(1-leverage,1e-8))
    score=np.zeros((G,W.shape[1]))
    np.add.at(score,codes,W*e[:,None])
    influence=score@inv
    influence*=np.sqrt(G/(G-1))
    rng=np.random.default_rng(seed)
    draws=[]
    for start in range(0,BOOT,250):
        m=min(250,BOOT-start)
        signs=2*rng.integers(0,2,size=(m,G))-1
        draws.append(beta+signs@influence)
    bs=np.concatenate(draws,axis=0)
    ci=np.quantile(bs,[.025,.975],axis=0)
    sim=np.quantile(bs,[.025/FAMILY,1-.025/FAMILY],axis=0)
    result.update(status="estimable",nuisance_rank=Q.shape[1],target_rank=W.shape[1],
                  bootstrap=BOOT,mean_by_level={k:float(x.loc[x.strategy.eq(k),outcome].mean()) for k in levels},
                  effects={name:{"coefficient":float(beta[i]),"marginal95":ci[:,i].tolist(),
                                 "bonferroni12_approx95":sim[:,i].tolist()} for i,name in enumerate(names)})
    return result

def analysis(d):
    eligible=d.loc[(d.env_n>=20)&(d.env_cells>=5)&(d.colour_valid_n>=20)].copy()
    spatial=eligible.loc[eligible.repeated5.eq(True)].dropna(subset=[
        f"{k}_{r}" for r in (50,250) for k in ("thermal_partition","moisture_partition","colour_partition")
    ]).copy()
    assert len(eligible)==4333 and len(spatial)==1258
    out={}
    for i,(name,levels) in enumerate(GROUPS.items()):
        # Exact common response support for categorical vs continuous colour.
        a=eligible.loc[eligible.strategy.isin(levels)].dropna(subset=BASE+["colour_D","hellinger_D"])
        extra=["env_near_fraction_250","colour_near_fraction_250","thermal_partition_250","moisture_partition_250"]
        b=spatial.loc[spatial.strategy.isin(levels)].dropna(subset=BASE+extra+[
            "colour_D","hellinger_D","colour_partition_250","hellinger_partition_250"])
        results={}
        for j,(outcome,x,base,interaction) in enumerate([
            ("colour_D",a,BASE,False),
            ("hellinger_D",a,BASE,False),
            ("colour_partition_250",b,BASE+extra+["colour_D"],True),
            ("hellinger_partition_250",b,BASE+extra+["hellinger_D"],True)]):
            results[outcome]={mode:wild_fit(x,levels,outcome,base,interaction,mode=="within_genus",SEED+i*100+j*10+int(mode=="within_genus"))
                              for mode in ("pooled","within_genus")}
        out[name]={"diversity_common_n":len(a),"spatial_common_n":len(b),"results":results}
    return {"strategies_all4475":d.strategy.value_counts().to_dict(),
            "eligible_diversity":len(eligible),"matched_spatial":len(spatial),"tests":out}

def selftest():
    # Known paired within-genus effect even when a nuisance variable is
    # constant inside every genus and therefore rank deficient after demeaning.
    rows=[]
    for i in range(40):
        for t,level in enumerate(("perennial_nonwoody","perennial_woody")):
            rows.append(dict(genus=str(i),strategy=level,cohort="x",z=float(i),
                             y=float(i+2*t)))
    x=pd.DataFrame(rows)
    q=wild_fit(x,("perennial_woody","perennial_nonwoody"),"y",["z"],False,True,SEED)
    assert q["status"]=="estimable"
    assert abs(q["effects"]["strategy_difference"]["coefficient"]-2)<1e-10
    print("FCP_GROWTH_FORM_LIFECYCLE_SELFTEST_PASS")

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--artifact",type=Path)
    p.add_argument("--out",type=Path)
    p.add_argument("--self-test",action="store_true")
    a=p.parse_args()
    selftest()
    if a.self_test:return
    assert a.artifact and a.out
    d=load(a.artifact)
    j={"source_artifact":ARTIFACT,"source_zip_sha256":DIGEST,"source_taxa":len(d),
       "source_photos":447500,"source_classified":227107,"bootstrap":BOOT,"primary_family":FAMILY,
       "scope":"Post hoc full-photo macroecological trait comparisons on fixed sampled landscapes. Growth form comparison restricted to perennials, annuality comparison restricted to nonwoody plants. Pooled 12 coefficients with approximate Bonferroni intervals; within-genus mixed-strategy contrasts are supplementary. Fixed-design genus wild-score bootstrap is conditional on observed design; not a full phylogeny, not a measurement of generation time, heritability, or selection.",
       **analysis(d)}
    a.out.mkdir(parents=True,exist_ok=True)
    d[["inat_taxon_id","species","genus","cohort","strategy","lifecycle","woodiness"]].to_csv(a.out/"strategy_ledger.csv",index=False)
    (a.out/"results.json").write_text(json.dumps(j,indent=2,allow_nan=False)+"\n")
    print("FCP_GROWTH_FORM_LIFECYCLE_COMPLETE",json.dumps(j,allow_nan=False),flush=True)

if __name__=="__main__":
    main()

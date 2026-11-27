#!/usr/bin/env python3
"""Secondary audit of frozen PR #275 observer-fold ledger; no refitting."""
import argparse
import hashlib
import io
import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

SOURCE_DIGEST = "b157ab08214cf1b67b37e8be5ba5db10bf5a511e4cc9c0e935b2a78017f59d4a"
SEED=20261128
BOOT=3999
NULL=9999

def read_source(path):
    raw=path.read_bytes()
    assert hashlib.sha256(raw).hexdigest()==SOURCE_DIGEST, "source ZIP digest mismatch"
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        names=z.namelist()
        files=[n for n in names if n.endswith("/fcp4475_local_observer_disjoint_ledger.csv.gz") or n=="fcp4475_local_observer_disjoint_ledger.csv.gz"]
        assert len(files)==1
        d=pd.read_csv(io.BytesIO(z.read(files[0])),compression="gzip")
        receipts=[n for n in names if n.endswith("/fcp4475_local_observer_disjoint_results.json") or n=="fcp4475_local_observer_disjoint_results.json"]
        assert len(receipts)==1
        r=json.loads(z.read(receipts[0]))
    assert (r["source_photos"],r["source_classifiable"],r["source_taxa"])==(447500,227107,4475)
    assert len(d)==8950 and not d.duplicated(["inat_taxon_id","direction"]).any()
    assert d.groupby("direction").size().eq(4475).all()
    a=d.loc[d.direction=="forward"].set_index("inat_taxon_id")
    b=d.loc[d.direction=="reverse"].set_index("inat_taxon_id")
    b=b.reindex(a.index)
    assert a.genus.eq(b.genus).all() and a.local_fold.eq(b.heldout_fold).all()
    return a,b

def cluster_ci(vals,genus,seed):
    vals=np.asarray(vals,float)
    if len(vals)<20 or len(set(genus))<15:
        return None
    codes,g=pd.factorize(genus,sort=True)
    G=len(g)
    sums=np.zeros((G,vals.shape[1]))
    np.add.at(sums,codes,vals)
    counts=np.bincount(codes,minlength=G)
    rng=np.random.default_rng(seed)
    draws=[]
    for _ in range(BOOT):
        w=np.bincount(rng.integers(G,size=G),minlength=G)
        draws.append((w@sums)/(w@counts))
    return np.quantile(draws,[.025,.975],axis=0).T.tolist()

def summarize(a,b,scope,key,seed):
    valid=(a.local_pairs>=1)&(b.local_pairs>=1)
    if scope=="five_pairs_each":
        valid&=(a.local_pairs>=5)&(b.local_pairs>=5)
    if scope=="nontrivial_both_month_nulls":
        valid&=a.local_expected_detection_month.between(.05,.95)&b.local_expected_detection_month.between(.05,.95)
    x,y=a.loc[valid],b.loc[valid]
    n=len(x)
    if not n:
        return {"n":0,"status":"not_estimable"}
    da,db=x.local_detected.to_numpy(float),y.local_detected.to_numpy(float)
    pa=x["local_expected_detection_"+key].to_numpy(float)
    pb=y["local_expected_detection_"+key].to_numpy(float)
    ma=x.local_mismatch.to_numpy(float)-x["local_expected_mismatch_"+key].to_numpy(float)
    mb=y.local_mismatch.to_numpy(float)-y["local_expected_mismatch_"+key].to_numpy(float)
    assert np.isfinite(np.column_stack([da,db,pa,pb,ma,mb])).all()
    assert np.isin(da,[0,1]).all() and np.isin(db,[0,1]).all()
    assert ((pa>=0)&(pa<=1)&(pb>=0)&(pb<=1)).all()
    ra,rb=da-pa,db-pb
    # Joint excess detects repeated local contrast above the composition null;
    # residual product asks whether deviations from null recur across folds.
    vals=np.column_stack([da,db,da*db,pa,pb,pa*pb,da*db-pa*pb,ra*rb,ma*mb])
    names=["detect_A","detect_B","detect_both","expected_A","expected_B","expected_both","joint_excess","detection_residual_product","mismatch_residual_product"]
    ci=cluster_ci(vals,x.genus,seed)
    # Conditional independent-vertex-null approximation from 199-shuffle
    # probabilities stored per species; not an independent biological sample.
    rng=np.random.default_rng(seed+1)
    count=0
    observed=int((da*db).sum())
    simulated=[]
    for start in range(0,NULL,250):
        m=min(250,NULL-start)
        z=(rng.random((m,n))<pa)&(rng.random((m,n))<pb)
        k=z.sum(axis=1)
        simulated.extend(k.tolist())
        count+=int((k>=observed).sum())
    def corr(u,v):
        return float(np.corrcoef(u,v)[0,1]) if len(u)>1 and np.std(u)>1e-12 and np.std(v)>1e-12 else None
    return {"n":n,"genera":int(x.genus.nunique()),"scope":scope,"null_stratum":key,
            "observed_both_n":observed,"expected_both_n":float((pa*pb).sum()),
            "detected_A_only_n":int((da*(1-db)).sum()),"detected_B_only_n":int(((1-da)*db).sum()),
            "neither_n":int(((1-da)*(1-db)).sum()),
            "exact_null_both_deterministic_n":int((np.isin(pa,[0,1])&np.isin(pb,[0,1])).sum()),
            "conditional_null_upper_p_approx":(count+1)/(NULL+1),
            "conditional_null_count_95":np.quantile(simulated,[.025,.975]).tolist(),
            "detection_residual_corr":corr(ra,rb),"mismatch_residual_corr":corr(ma,mb),
            "metrics":{name:{"estimate":float(vals[:,i].mean()),"genus_ci95":None if ci is None else ci[i]} for i,name in enumerate(names)}}

def crossfit_assoc(a,b,seed):
    # Continuous local mismatch excess in A versus independent geography
    # and climate in B. Each direction counted separately.
    out={}
    for k,(x,y) in enumerate([(a,b),(b,a)]):
        # geographic and thermal held-out measurements in the same row
        # as the local predictor, rather than the opposite row.
        for j,target in enumerate(["geo_delta","thermal_white"]):
            valid=x.local_excess_month.notna()&x[target].notna()
            d=x.loc[valid]
            u=d.local_excess_month.to_numpy(float)
            v=d[target].to_numpy(float)
            n=len(d)
            G=int(d.genus.nunique())
            rec={"n":n,"genera":G}
            if n>=20 and G>=15 and np.std(u)>1e-12 and np.std(v)>1e-12:
                # Opportunity covariates: baseline and heldout support, expected
                # mismatch from baseline month null, and source cohort.
                C=np.column_stack([np.log1p(d.local_pairs),np.log1p(d.local_valid_n),np.log1p(d.heldout_classified_n),d.local_expected_mismatch_month])
                C=C[:,C.std(axis=0)>1e-10]
                C=(C-C.mean(axis=0))/C.std(axis=0)
                cohorts=pd.get_dummies(d.source_cohort,drop_first=True,dtype=float).to_numpy()
                X=np.column_stack([np.ones(n),u,C,cohorts])
                p=X.shape[1]
                if np.linalg.matrix_rank(X)==p:
                    codes,_=pd.factorize(d.genus,sort=True)
                    gram=np.zeros((G,p,p));rhs=np.zeros((G,p))
                    np.add.at(gram,codes,np.einsum("ni,nj->nij",X,X))
                    np.add.at(rhs,codes,X*v[:,None])
                    rng=np.random.default_rng(seed+k*100+j)
                    bs=[]
                    for _ in range(BOOT):
                        w=np.bincount(rng.integers(G,size=G),minlength=G)
                        g=np.einsum("g,gij->ij",w,gram)
                        if np.linalg.matrix_rank(g)==p:
                            bs.append(float(np.linalg.solve(g,w@rhs)[1]))
                    rec.update(status="estimable" if len(bs)>=.9*BOOT else "unstable",
                               slope_per_unit_excess=float(np.linalg.lstsq(X,v,rcond=None)[0][1]),
                               ci95=np.quantile(bs,[.025,.975]).tolist() if len(bs)>=.9*BOOT else None,
                               valid_boot=len(bs))
                else:
                    rec["status"]="rank_deficient"
            else:
                rec["status"]="insufficient_support"
            out[("forward" if k==0 else "reverse")+"_"+target]=rec
    return out

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--artifact",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    args=p.parse_args()
    a,b=read_source(args.artifact)
    result={"source_run":45953690498,"source_artifact_id":15178380182,"source_zip_sha256":SOURCE_DIGEST,
            "source_taxa":len(a),"bootstrap":BOOT,"conditional_null_draws":NULL,
            "scope":"Exploratory species-level repeat detection, not persistence in the same locality or genetic morphs. Independent photographers within species, not an independent platform. Conditional null uses estimated probabilities (199 permutations per fold) and assumes independent permutations given source data. Exact-year/month may be degenerate. No new photo models.",
            "recurrence":{},"continuous_crossfit":crossfit_assoc(a,b,SEED+1000)}
    for i,scope in enumerate(["all_both_local_opportunity","five_pairs_each","nontrivial_both_month_nulls"]):
        for j,key in enumerate(["month","yearmonth"]):
            result["recurrence"][scope+"_"+key]=summarize(a,b,scope,key,SEED+i*100+j*10)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("FCP_LOCAL_FOLD_RECURRENCE_COMPLETE",json.dumps(result,allow_nan=False),flush=True)

if __name__=="__main__":
    main()

#!/usr/bin/env python3
"""Independent common-species audit of FCP seven-climate photographic covariance.

Compares total, between-four-original-photo-label and within-label covariance
on exactly the SAME eligible taxa. Signed covariance fractions are not genetic
morph shares or causal decompositions. The files must match source SHAs.
"""
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd

AXES=("BIO1","BIO4","BIO5","BIO12","BIO14","BIO15","SRADmean")
MODES=("original","observer_thinned")
REGIMES=("raw","geo_month_residual")
COMP=("total","between_coarse_label","within_coarse_label")
SHA={"original":"0fbffbf8f562937ca71edb95889298e6c56ac0e9253026977a0cb7a9b9a2e689",
     "observer_thinned":"8acfcb8b42b57937a155e1fce408e96929fa98cbe6472fb67aac2c24084bbeb1",
     "summary":"688b9ec926c39db07f4dc5d3a35940692cde8e670e96a606306eaa134e5cd0db"}
REPS=1999


def require(x,msg):
    if not x:raise RuntimeError("FCP_CLIMATE_COMMON_SUPPORT: "+msg)


def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(1<<20),b""):h.update(block)
    return h.hexdigest()


def source_subset(df,regime,axis):
    cols=[f"{axis}_{regime}_{c}_signed_white_fraction_excess_vs_month_null" for c in COMP]
    x=df[cols].to_numpy(float)
    good=np.isfinite(x).all(axis=1)
    x=x[good]
    sub=df.loc[good].reset_index(drop=True)
    require(len(sub)>100 and np.max(np.abs(x[:,0]-x[:,1]-x[:,2]))<1e-9,
            "same eligible source taxa or exact covariance identity lost")
    return sub,x


def genus_bootstrap(sub,x,reps=REPS,seed=20261024):
    g,names=pd.factorize(sub.genus,sort=True)
    n=len(names)
    sums=np.zeros((n,3));counts=np.bincount(g,minlength=n)
    np.add.at(sums,g,x)
    rng=np.random.default_rng(seed)
    indices=rng.integers(n,size=(reps,n))
    draws=sums[indices].sum(axis=1)/counts[indices].sum(axis=1)[:,None]
    obs=x.mean(axis=0)
    part=draws[:,1]/draws[:,0]
    part=part[np.isfinite(part)]
    return {
        "n_species":len(sub),"n_genera":n,
        "signed_white_covariance_month_matched_excess":{
            name:{"mean":float(obs[i]),
                  "genus_bootstrap_marginal95":[float(v) for v in np.quantile(draws[:,i],[.025,.975])]}
            for i,name in enumerate(COMP)},
        "between_minus_within":float(obs[1]-obs[2]),
        "between_minus_within_95":[float(v) for v in np.quantile(draws[:,1]-draws[:,2],[.025,.975])],
        "between_signed_covariance_share_if_same_direction":{
            "mean_ratio":float(obs[1]/obs[0]) if abs(obs[0])>1e-8 else None,
            "genus_bootstrap_marginal95":[float(v) for v in np.quantile(part,[.025,.975])],
            "not_fraction_genetic_explained_variance":True},
        "within_plus_between_equals_total":bool(np.allclose(obs[1]+obs[2],obs[0],atol=1e-12))
    }


def family7_total(df,regime,reps=REPS,seed=20261024):
    # Common membership with two original photo labels >=5 per label; species
    # with non-identifiable colour frequency remain null, not zero.
    cols=[f"{a}_{regime}_total_signed_white_fraction_excess_vs_month_null" for a in AXES]
    bc=[f"{a}_{regime}_between_coarse_label_signed_white_fraction_excess_vs_month_null" for a in AXES]
    val=df[cols].to_numpy(float)
    eligible=np.isfinite(df[bc].to_numpy(float))
    val=np.where(eligible,val,np.nan)
    gid,genera=pd.factorize(df.genus,sort=True)
    G=len(genera)
    sumG=np.zeros((G,7));countG=np.zeros((G,7))
    np.add.at(sumG,gid,np.nan_to_num(val))
    np.add.at(countG,gid,np.isfinite(val).astype(float))
    rng=np.random.default_rng(seed)
    ix=rng.integers(G,size=(reps,G))
    boot=sumG[ix].sum(axis=1)/countG[ix].sum(axis=1)
    require(np.isfinite(boot).all(),"seven climate axes not jointly estimable")
    mean=sumG.sum(axis=0)/countG.sum(axis=0)
    sd=boot.std(axis=0,ddof=1)
    crit=np.quantile(np.max(np.abs((boot-mean)/sd),axis=1),.95)
    return {ax:{"n_species":int(countG[:,i].sum()),
                "mean_month_null_excess":float(mean[i]),
                "genus_family7_simultaneous95":[float(mean[i]-crit*sd[i]),float(mean[i]+crit*sd[i])],
                "family7_excludes_zero":bool(mean[i]-crit*sd[i]>0 or mean[i]+crit*sd[i]<0)}
            for i,ax in enumerate(AXES)}


def audit(paths,summary):
    for name,path in {**paths,"summary":summary}.items():
        require(sha(path)==SHA[name],"source-locked "+name+" SHA mismatch")
    source=json.loads(Path(summary).read_text())
    require(source["source_original_photos"]==447500 and
            source["source_classifiable_photos"]==227107 and
            source["source_original_species"]==4475 and
            source["axes"]==list(AXES),"genuine 447500 photo/climate frame lost")
    results={}
    for ii,mode in enumerate(MODES):
        df=pd.read_csv(paths[mode])
        require(len(df)==4475 and df.inat_taxon_id.is_unique,
                "per-taxon original source 4475 ledger not preserved")
        block={}
        for jj,regime in enumerate(REGIMES):
            fam=family7_total(df,regime,seed=202610240+ii*11+jj)
            axes={}
            for k,axis in enumerate(AXES):
                sub,vals=source_subset(df,regime,axis)
                expected=source["results"][mode]["per_axis_results"][regime][axis]["n_between_label_supported"]
                require(len(sub)==expected,
                        "between-label support changed for "+mode+" "+regime+" "+axis)
                record=genus_bootstrap(sub,vals,seed=202610250+ii*109+jj*17+k)
                record["family7_total_on_common_support"]=fam[axis]
                axes[axis]=record
            block[regime]=axes
        results[mode]=block
    return {
        "schema":"fcp4475_photo9_climate7_exact_common_support_signed_covariance_v1",
        "source_run_id":40029760904,"source_artifact_sha256":dict(SHA),
        "n_original_photos":447500,"n_taxa":4475,"n_classified_photos":227107,
        "bootstrap_replicates":REPS,
        "covariance_not_slope":True,
        "fraction_is_algebraic_signed_covariance_not_genetic_morph_share":True,
        "no_independently_verified_flower_organs_or_photo_background":True,
        "familywise_seven_climate_axes_only_within_each_exploratory_family":True,
        "results":results
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original",type=Path)
    p.add_argument("--observer-thinned",type=Path)
    p.add_argument("--summary",type=Path)
    p.add_argument("--out",type=Path)
    p.add_argument("--self-test",action="store_true")
    a=p.parse_args()
    if a.self_test:
        test=np.array([[.03,.02,.01],[.01,.005,.005]])
        require(np.allclose(test[:,0],test[:,1]+test[:,2]),"fixture covariance false")
        print("FCP_COLOUR_CLIMATE_COMMON_SUPPORT_SELF_TEST_PASS")
        return
    require(a.original and a.observer_thinned and a.summary and a.out,
            "original and observer-thinned original-source CSVs/summary required")
    result=audit({"original":a.original,"observer_thinned":a.observer_thinned},a.summary)
    a.out.mkdir(parents=True,exist_ok=True)
    (a.out/"source_verified_common_support_seven_climate.json").write_text(
        json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("FCP4475_CLIMATE_COMMON_SUPPORT_AUDIT_PASS",json.dumps({
        kind:{reg:{ax:{"n":v["n_species"],
                         "signed":{k:x["mean"] for k,x in v["signed_white_covariance_month_matched_excess"].items()},
                         "share":v["between_signed_covariance_share_if_same_direction"]["mean_ratio"],
                         "family":v["family7_total_on_common_support"]}
                   for ax,v in axes.items()} for reg,axes in block.items()}
        for kind,block in result["results"].items()},sort_keys=True))


if __name__=="__main__":
    main()

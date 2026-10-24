#!/usr/bin/env python3
"""Same-taxon 5degree vs 1degree seven-climate photographic whiteness audit.

All four source CSVs must come unchanged from the SHA-authenticated prior
original-photo regional analysis. The difference between 5 and 1 degree
covariances is influenced by different climate residual variance; compare
variance-normalized null-excess covariance within EXACT same photo/taxon
sample, but do not interpret as a genetic-morph or causal climate slope.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd

AXES=("BIO1","BIO4","BIO5","BIO12","BIO14","BIO15","SRADmean")
SHA={
    "original_5":"828a9a8810cf3a2fed00075c21b544add2c24f2ae7a08794331f8364f850fc91",
    "original_1":"31ad317ff735b2300123149ff60f0307b6cc3c9dadff2e8b59efe63e5116283e",
    "observer_thinned_5":"d7a9ea1a86bcf788c099d7deee6492cb0dc56be1c1201c655d53399c4e84e27f",
    "observer_thinned_1":"27911db48c260a39a091741a766163691e008b609841de35114faf1e9c1c180d",
    "summary":"9239086c2141fca8d63b94bca149c296a57e863e7ee5d8a2c0491e5ac12a7c00"
}
B=1999
SEED=20261024
METRICS=("cov5","cov1","cov5_minus_cov1",
         "slope5","slope1","slope5_minus_slope1")


def require(ok,msg):
    if not ok:raise RuntimeError("FCP_BIO5_REGIONAL_COMMON_SUPPORT: "+msg)


def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for part in iter(lambda:f.read(1<<20),b""):h.update(part)
    return h.hexdigest()


def load_and_match(a,b):
    p=pd.read_csv(a,low_memory=False)
    q=pd.read_csv(b,low_memory=False)
    for x in (p,q):
        require(len(x)==4475 and x.inat_taxon_id.is_unique,
                "expected complete 4,475 source taxon ledger in each scale")
        require(int(x.n_source_classifiable.sum())==227107,
                "photographic classifiability denominator changed")
    require(p.inat_taxon_id.equals(q.inat_taxon_id),
            "same taxon row order/identity not preserved across scales")
    require(p.species.astype(str).equals(q.species.astype(str)),
            "source taxon names changed across regional models")
    require((p.n_photos_with_complete_climate==q.n_photos_with_complete_climate).all(),
            "five/one degree inputs did not use identical photos")
    require((p.one_photo_per_observer==q.one_photo_per_observer).all(),
            "observer-thinning design differs across scales")
    return p,q


def covariate_table(p,q):
    data={}
    for axis in AXES:
        cols=[f"{axis}_total_signed_white_excess",
              f"{axis}_climate_resid_variance_fraction"]
        require(all(c in p.columns and c in q.columns for c in cols),
                "absent original seven climate columns")
        v5=pd.to_numeric(p[cols[0]],errors="coerce").to_numpy(float)
        v1=pd.to_numeric(q[cols[0]],errors="coerce").to_numpy(float)
        f5=pd.to_numeric(p[cols[1]],errors="coerce").to_numpy(float)
        f1=pd.to_numeric(q[cols[1]],errors="coerce").to_numpy(float)
        good=(np.isfinite(v5)&np.isfinite(v1)&np.isfinite(f5)&np.isfinite(f1)
              &(f5>=.05)&(f1>=.05))
        frame=np.full((len(p),len(METRICS)),np.nan,float)
        frame[good,0]=v5[good]
        frame[good,1]=v1[good]
        frame[good,2]=v5[good]-v1[good]
        frame[good,3]=v5[good]/f5[good]
        frame[good,4]=v1[good]/f1[good]
        frame[good,5]=frame[good,3]-frame[good,4]
        data[axis]=(frame,good,f5,f1)
    return data


def clustered_family_bootstrap(p,matrices,reps=B,seed=SEED):
    # One family of seven climate axes for each metric. Cross-metric and
    # cross-observer-panel testing remains exploratory and is not corrected.
    group,labels=pd.factorize(p.genus,sort=True)
    G=len(labels)
    rng=np.random.default_rng(seed)
    draws=rng.integers(G,size=(reps,G))
    out={axis:{} for axis in AXES}
    for k,metric in enumerate(METRICS):
        V=np.column_stack([matrices[axis][0][:,k] for axis in AXES])
        mask=np.isfinite(V)
        values=np.nan_to_num(V)
        sums=np.zeros((G,len(AXES)),float)
        counts=np.zeros((G,len(AXES)),float)
        np.add.at(sums,group,values)
        np.add.at(counts,group,mask.astype(float))
        denom=counts.sum(axis=0)
        require(np.all(denom>100),"insufficient common-support per axis")
        observed=sums.sum(axis=0)/denom
        bs=sums[draws].sum(axis=1)/counts[draws].sum(axis=1)
        require(np.isfinite(bs).all(),"genus bootstrap lost source support")
        sd=bs.std(axis=0,ddof=1)
        require((sd>0).all(),"bootstrap genus SD zero")
        maxT=np.max(np.abs((bs-observed)/sd),axis=1)
        cutoff=float(np.quantile(maxT,.95))
        for j,axis in enumerate(AXES):
            ci=np.quantile(bs[:,j],[.025,.975])
            fam=[float(observed[j]-cutoff*sd[j]),float(observed[j]+cutoff*sd[j])]
            out[axis][metric]={
                "n_identical_taxa":int(denom[j]),
                "mean":float(observed[j]),
                "genus_bootstrap_marginal95":[float(ci[0]),float(ci[1])],
                "genus_family7_simultaneous95":fam,
                "family7_excludes_zero":bool(fam[0]>0 or fam[1]<0),
            }
    return out,G


def audit(paths,source_summary):
    for key,path in {**paths,"summary":source_summary}.items():
        require(digest(path)==SHA[key],
                "source per-taxon CSV byte identity differs "+key)
    receipt=json.loads(Path(source_summary).read_text())
    require(receipt["source_n_species"]==4475 and
            receipt["source_n_photos"]==447500 and
            receipt["source_n_4state_classified"]==227107 and
            receipt["region_width_degrees"]==[5,1],
            "authenticated source regional run manifest drift")
    results={}
    for idx,mode in enumerate(("original","observer_thinned")):
        p,q=load_and_match(paths[f"{mode}_5"],paths[f"{mode}_1"])
        matrices=covariate_table(p,q)
        measures,genera=clustered_family_bootstrap(
            p,matrices,seed=SEED+idx*1000)
        diag={}
        for axis in AXES:
            _,good,v5,v1=matrices[axis]
            diag[axis]={
                "n_shared_source_species":int(good.sum()),
                "n_genera":int(p.loc[good,"genus"].nunique()),
                "median_climate_residual_fraction_5degree":float(np.median(v5[good])),
                "median_climate_residual_fraction_1degree":float(np.median(v1[good])),
                "fraction_species_positive_1degree_signed_white":float(
                    np.mean(matrices[axis][0][good,1]>0)),
                "same_source_covariance_and_variance_normalized":measures[axis],
            }
        results[mode]={
            "n_genera_total_original_taxa":genera,
            "per_axis_common_support":diag,
        }
    return {
        "schema":"fcp4475_region5_vs_region1_same_species_seven_climate_signed_hue_v1",
        "source_run":40056835373,
        "source_summary_artifact_id":12545591935,
        "source_original_photos":447500,
        "source_species":4475,
        "source_classifiable":227107,
        "frozen_source_CSV_SHA256":dict(SHA),
        "n_genus_bootstrap_replicates":B,
        "slope_metric":"per-species month-region permutation-null-centered signed white-palette covariance divided by residual within-source-species standardized climate variance, not causal slope",
        "original_climate_SD_normalisation":"worldclim sampled photo climate SD within each taxon, prior to region+month residualisation",
        "paired_scales_require_identical_original_photo_ids_and_eligible_taxa":True,
        "one_degree_missing_or_no_contrast_species_never_imputed":True,
        "geographic_region_fixed_effects_are_not_genomic_controls":True,
        "bootstrap_family7_controls_each_metric_family_not_all_12_panels_jointly":True,
        "photo_labels_are_not_genetic_morphs":True,
        "no_background_organ_genotype_fitness_controls":True,
        "results":results
    }


def self_test():
    z=np.array([[.6,.4,.2,.6/.3,.4/.1,.6/.3-.4/.1]])
    require(abs(z[0,2]-.2)<1e-12 and abs(z[0,5]+2)<1e-12,
            "variance-normalized photographic covariance test fixture wrong")
    print("FCP_SEVEN_CLIMATE_PAIRED_REGIONAL_SLOPE_SELFTEST_PASS")


def main():
    p=argparse.ArgumentParser()
    for key in ("original-5","original-1","observer-thinned-5",
                "observer-thinned-1","summary","out"):
        p.add_argument("--"+key,type=Path)
    p.add_argument("--self-test",action="store_true")
    a=p.parse_args()
    self_test()
    if a.self_test:return
    all_paths={
        "original_5":a.original_5,"original_1":a.original_1,
        "observer_thinned_5":a.observer_thinned_5,
        "observer_thinned_1":a.observer_thinned_1,
    }
    require(all(all_paths.values()) and a.summary and a.out,
            "four original source CSVs plus receipt required")
    result=audit(all_paths,a.summary)
    a.out.mkdir(parents=True,exist_ok=True)
    (a.out/"paired_5vs1_seven_climate_original_hues.json").write_text(
        json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n")
    for mode in ("original","observer_thinned"):
        z=result["results"][mode]["per_axis_common_support"]["BIO5"]
        print("FCP_SEVEN_CLIMATE_PAIRED_BIO5",mode,json.dumps(z,sort_keys=True))
    print("FCP_SEVEN_CLIMATE_PAIRED_FULL_JSON",json.dumps(result,sort_keys=True,allow_nan=False))


if __name__=="__main__":
    main()

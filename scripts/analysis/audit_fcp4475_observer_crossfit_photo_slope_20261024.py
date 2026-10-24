#!/usr/bin/env python3
"""Independent source-derived audit of FCP observer-disjoint twofold colour data.

Original source of truth: prior SHA-locked GitHub Actions 447,500-photo
crossfit artifact; two 4,475-taxon photo fold directions and two heldout
photo-observer response modes. The audit verifies all crossfit PHOTO baseline
group assignments from original source counts and exact aggregated source
means/conditional intervals, and quantifies observer-fold category agreement
and reproducibility of independently estimated per-species photo-climate
slopes. It cannot independently recover observer/photo IDs from derived tables.

Any lack of photo-slope replication can reflect half-sample estimation noise.
No biological inherited morph, fitness or climatic causation is inferred.
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd

SOURCE={
 "fcp4475_observer_disjoint_crossfit_white_thermal_results.json":
     "40b8be90523f876aaf755354fd9527d7429782d8cb42414a69783dedbfa4e4c9",
 "fcp4475_source_all_taxa_twofold_observer_crossfit_thermal_white.csv":
     "079241c220a5f5993a17d361a2ea5fc276c7d13e99a168111cc14a3f51ce7253"
}
B=999
SEED=20261024
MODES=("all_heldout","one_photo_per_observer")
DIR=("forward","reverse")
CATS=("none","rare","mixed","nearly_all")


def check(ok,msg):
    if not ok:raise RuntimeError("FCP4475_INDEPENDENT_CROSSFIT_AUDIT: "+msg)


def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for ch in iter(lambda:f.read(1<<20),b""):h.update(ch)
    return h.hexdigest()


def recompute_cat(nwhite,nonwhite,cutoff):
    if nwhite==0:return "none"
    if nwhite<cutoff:return "rare"
    if nonwhite>=cutoff:return "mixed"
    return "nearly_all"


def cluster_ci_one(frame,seed):
    d=frame.loc[frame.heldout_thermal_white_slope.notna()].copy()
    x=d.heldout_thermal_white_slope.to_numpy(float)
    cod,gene=pd.factorize(d.genus,sort=True)
    G=len(gene)
    sums=np.zeros(G,float)
    n=np.bincount(cod,minlength=G)
    np.add.at(sums,cod,x)
    draw=np.random.default_rng(seed).integers(G,size=(B,G))
    boot=sums[draw].sum(axis=1)/n[draw].sum(axis=1)
    return {
        "n_source_taxa":len(x),"n_genus_clusters":G,
        "mean_species_equal_heldout_white_thermal_slope":float(x.mean()),
        "fraction_taxa_positive_heldout_slope":float(np.mean(x>0)),
        "genus_bootstrap_marginal95":[float(v) for v in np.quantile(boot,[.025,.975])],
    }


def moment_matrix(x,y,code,G):
    n=np.bincount(code,minlength=G)
    result=np.zeros((G,6),float)
    result[:,0]=n
    for i,z in enumerate((x,y,x*x,y*y,x*y),1):
        np.add.at(result[:,i],code,z)
    return result


def correlation_moments(a):
    # n, sumx, sumy, sumx2, sumy2, sumxy
    v=np.asarray(a,float)
    n=v[...,0]
    cov=v[...,5]-v[...,1]*v[...,2]/n
    xx=v[...,3]-v[...,1]**2/n
    yy=v[...,4]-v[...,2]**2/n
    den=np.sqrt(np.maximum(0,xx*yy))
    return np.divide(cov,den,out=np.full_like(cov,np.nan),where=den>1e-13)


def paired_split_diagnostics(rows,mode,seed):
    f=rows.loc[(rows.direction=="forward")&(rows.test_photo_mode==mode)].copy()
    r=rows.loc[(rows.direction=="reverse")&(rows.test_photo_mode==mode)].copy()
    ensure=["inat_taxon_id","genus","baseline_group3",
            "baseline_is_adequately_sampled","heldout_thermal_white_slope"]
    f=f[ensure].rename(columns={c:c+"_f" for c in ensure if c not in ("inat_taxon_id","genus")})
    r=r[ensure].rename(columns={c:c+"_r" for c in ensure if c not in ("inat_taxon_id","genus")})
    m=f.merge(r,on=("inat_taxon_id","genus"),validate="one_to_one")
    check(len(m)==4475,"original paired taxon set does not contain all named source taxa")
    both=(m.baseline_is_adequately_sampled_f&
          m.baseline_is_adequately_sampled_r)
    s=m.loc[both]
    agree=float(np.mean(s.baseline_group3_f==s.baseline_group3_r))
    freq1=s.baseline_group3_f.value_counts(normalize=True)
    freq2=s.baseline_group3_r.value_counts(normalize=True)
    chance=float(sum(freq1.get(key,0)*freq2.get(key,0) for key in CATS))
    kappa=(agree-chance)/(1-chance)
    score=pd.crosstab(s.baseline_group3_f,s.baseline_group3_r).reindex(
        index=CATS,columns=CATS,fill_value=0)
    present=(m.heldout_thermal_white_slope_f.notna()&
             m.heldout_thermal_white_slope_r.notna())
    d=m.loc[present].copy()
    x=d.heldout_thermal_white_slope_f.to_numpy(float)
    y=d.heldout_thermal_white_slope_r.to_numpy(float)
    g,genera=pd.factorize(d.genus,sort=True)
    G=len(genera)
    moments=moment_matrix(x,y,g,G)
    total=moments.sum(axis=0)
    pearson=float(correlation_moments(total))
    choices=np.random.default_rng(seed).integers(G,size=(B,G))
    bootstrap=correlation_moments(moments[choices].sum(axis=1))
    bootstrap=bootstrap[np.isfinite(bootstrap)]
    check(len(bootstrap)>=.9*B,"genus bootstrap correlation source variance degenerate")
    return {
       "n_taxa_both_baseline_folds_support":int(both.sum()),
       "n_baseline_white_photo_state_agree":int(np.sum(s.baseline_group3_f==s.baseline_group3_r)),
       "baseline_photo_state_agreement":agree,
       "chance_agreement_given_observed_group_marginals":chance,
       "cohens_kappa_baseline_photo_state":float(kappa),
       "baseline_white_group_confusion":{
           k:{j:int(score.loc[k,j]) for j in CATS} for k in CATS
       },
       "n_taxa_both_temperature_slope_folds_estimable":int(present.sum()),
       "n_genera_for_crossfold_photo_slope_corr":int(G),
       "species_photo_thermal_slope_crossfold_pearson":pearson,
       "slope_crossfold_pearson_genus_bootstrap95":[
           float(v) for v in np.quantile(bootstrap,[.025,.975])],
       "crossfold_thermal_slope_forward_mean_on_same_taxa":float(x.mean()),
       "crossfold_thermal_slope_reverse_mean_on_same_taxa":float(y.mean()),
       "observer_folds_are_sensitivity_splits_not_separate_genetic_populations":True,
       "lack_of_stable_slope_may_reflect_18to30_photo_sample_noise":True
    }


def source_audit(path):
    for name,known in SOURCE.items():
        check(digest(path/name)==known,"independent source image/CSV SHA mismatch "+name)
    r=json.loads((path/list(SOURCE)[0]).read_text())
    d=pd.read_csv(path/list(SOURCE)[1],low_memory=False)
    check(r["source_original_photos"]==447500 and r["source_taxa"]==4475
          and r["source_classified_nine_photo_hue"]==227107,
          "original source 447500/4475/227107 denominator drift")
    check(len(d)==17900 and not d.duplicated(
        ["inat_taxon_id","direction","test_photo_mode"]).any(),
        "source photo four-panel original 4475 taxon rows corrupted")
    check(d.n_classified_original.notna().all() and
          d.loc[(d.direction=="forward")&(d.test_photo_mode=="all_heldout"),
                "n_classified_original"].sum()==227107,
          "original source photo classifiability denominator differs from real release")
    for cutoff,key in ((3,"baseline_group3"),(5,"baseline_group5")):
        expected=[recompute_cat(int(a),int(b),cutoff)
                  for a,b in zip(d.baseline_n_white,d.baseline_n_nonwhite)]
        check((np.asarray(expected)==d[key].astype(str).to_numpy()).all(),
              "photo-white classification differs from immutable original observer fold")
    check((d.baseline_n_white+d.baseline_n_nonwhite==d.baseline_n_classified).all(),
          "original source classified photo baseline mass mismatch")
    panels={}
    for i,mode in enumerate(MODES):
        panels[mode]={"paired":paired_split_diagnostics(d,mode,SEED+i*53)}
        for j,direction in enumerate(DIR):
            key=direction+"__"+mode
            sample=d.loc[(d.direction==direction)&(d.test_photo_mode==mode)].copy()
            info=cluster_ci_one(sample,SEED+i*73+j*101)
            original=r["panels"][key]
            check(info["n_source_taxa"]==original["n_heldout_thermal_estimable"],
                  "original source per-panel photo-climate opportunity lost")
            for cat in CATS:
                c=sample.loc[sample.heldout_thermal_white_slope.notna()&
                             sample.baseline_group3.eq(cat)]
                v=original["three_image_baseline_group"]["results"][cat]
                check(len(c)==v["n_taxa"] and
                      np.isclose(c.heldout_thermal_white_slope.mean(),v["mean"],atol=1e-9),
                      "source mean mixed/none source photo class effect drift")
            panels[mode][direction]=info
    return {
        "schema":"fcp4475_observer_disjoint_photo_white_climate_source_independent_QA_v1",
        "source_photo_run_id":40080001765,
        "source_artifact_id":12554693263,
        "source_artifact_sha256":"50d5d15a2d145b0db7a176f945435701251768ba614f659d9dac496011b33177",
        "original_source_csv_and_result_sha256":SOURCE,
        "original_source_photos":447500,
        "source_taxa":4475,
        "source_classified_nine_hue_photos":227107,
        "independent_observer_photo_baseline_vs_heldout_climate":True,
        "genus_cluster_bootstrap_draws":B,
        "result":"independently source-audited original 2fold photo assignment and taxon diagnostics",
        "photo_climate_slope_photo_state_may_be_noisy_due_to_half_sample":True,
        "no_independently_verified_real_flower_petals_or_genotype_morphs":True,
        "panels":panels
    }


def main():
    a=argparse.ArgumentParser()
    a.add_argument("--source",type=Path)
    a.add_argument("--out",type=Path)
    a.add_argument("--self-test",action="store_true")
    opt=a.parse_args()
    fake=np.array([[3,3.,7.,5.,21.,10.],[4,10.,8.,30.,22.,25.]])
    check(np.all(np.isfinite(correlation_moments(fake))),"paired-photo correlation fixture broken")
    print("FCP_CROSSFIT_SOURCE_INDEPENDENT_AUDIT_SYNTHETIC_PASS")
    if opt.self_test:return
    check(opt.source and opt.out,"original SHA-authenticated 447500-photo crossfit artifact required")
    res=source_audit(opt.source)
    opt.out.mkdir(parents=True,exist_ok=True)
    (opt.out/"original_source_twofold_photo_group_stability_and_slope_reliability.json").write_text(
        json.dumps(res,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("FCP_CROSSFIT_INDEPENDENT_SOURCE_AUDIT_COMPLETE",json.dumps(res,sort_keys=True,allow_nan=False))

if __name__=="__main__":
    main()

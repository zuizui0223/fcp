#!/usr/bin/env python3
"""Source-authenticated FCP original-photo slope repeatability vs sampling support.

Are within-species temperature/white PHOTO slopes replicable across independent
observer folds when BOTH folds have more photos and informative climate range?
This is a nested original-photo opportunity sensitivity, not a new genomic
flower-colour polymorphism study and not a confirmed camera/background control.

Original inputs: SHA256-locked 447500-photo-derived source artefact from PR271.
All 4475 named taxa remain in the output, nonestimable species are not zero.
All resampling preserves genus clusters; predictor q is strictly colour-outcome
blind (min heldout classified observations * unique climate variance fraction).
"""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path

import numpy as np
import pandas as pd

EXPECTED = {
 "fcp4475_observer_disjoint_crossfit_white_thermal_results.json":
 "40b8be90523f876aaf755354fd9527d7429782d8cb42414a69783dedbfa4e4c9",
 "fcp4475_source_all_taxa_twofold_observer_crossfit_thermal_white.csv":
 "079241c220a5f5993a17d361a2ea5fc276c7d13e99a168111cc14a3f51ce7253",
}
B=1999
SEED=20261024
MODES=("all_heldout","one_photo_per_observer")
COHORTS=("discovery","validation","prospective","wave4","new_U100")
FIELDS=("n","x","y","x2","y2","xy","same_sign","pos_x","pos_y")
Q=("Q1","Q2","Q3","Q4")


def must(ok,msg):
    if not ok:raise RuntimeError("FCP4475_FOLD_RELIABILITY_SOURCE_AUDIT: "+msg)


def sha(p):
    h=hashlib.sha256()
    with Path(p).open("rb") as f:
        for z in iter(lambda:f.read(1<<20),b""):h.update(z)
    return h.hexdigest()


def load(source):
    for name,expected in EXPECTED.items():
        must(sha(source/name)==expected,"immutable original photo crossfit source SHA differs: "+name)
    root=json.loads((source/list(EXPECTED)[0]).read_text())
    d=pd.read_csv(source/list(EXPECTED)[1],low_memory=False)
    must(root["source_original_photos"]==447500
         and root["source_taxa"]==4475
         and root["source_classified_nine_photo_hue"]==227107,
         "original photo/taxon denominator not authentic")
    must(len(d)==4475*4
         and not d.duplicated(["inat_taxon_id","direction","test_photo_mode"]).any()
         and d.inat_taxon_id.nunique()==4475,
         "original 4-panel species integrity failed")
    for mode in MODES:
        for direction in ("forward","reverse"):
            a=d.loc[(d.direction==direction)&(d.test_photo_mode==mode)]
            must(len(a)==4475 and a.inat_taxon_id.is_unique,
                 "missing full original photo-taxon 4475 population")
    return root,d


def join_folds(allrows,mode):
    m=allrows.loc[allrows.test_photo_mode==mode]
    out={}
    baseline=("genus","species","source_cohort","n_classified_original")
    for name in ("forward","reverse"):
        sub=m.loc[m.direction==name].set_index("inat_taxon_id",drop=False)
        if name=="forward":
            f=sub
        else:
            r=sub
    must(f.inat_taxon_id.index.equals(r.inat_taxon_id.index),
         "observer fold-source original taxon order differs")
    for c in baseline:
        must(f[c].astype(str).equals(r[c].astype(str)),
             "original taxon identity changed between original observer folds "+c)
    cols=(
       "heldout_thermal_white_slope",
       "heldout_n_complete","heldout_n_observers","heldout_n_cells",
       "thermal_unique_var_frac","heldout_BIO5_q90_q10_C",
       "baseline_n_classified","baseline_n_observers",
       "baseline_group3","baseline_white_fraction"
    )
    d=f[["inat_taxon_id",*baseline,*cols]].join(
        r[list(cols)],rsuffix="_reverse",validate="one_to_one"
    ).reset_index(drop=True)
    # Rename the forward columns to simplify the paired reliability formulas.
    d=d.rename(columns={c:c+"_forward" for c in cols})
    d=d.rename(columns={c+"_reverse":c+"_reverse" for c in cols})
    good=np.isfinite(d.heldout_thermal_white_slope_forward.to_numpy(float))&np.isfinite(
        d.heldout_thermal_white_slope_reverse.to_numpy(float))
    d["both_folds_estimable"]=good
    d["min_climate_photo_count"]=d[
        ["heldout_n_complete_forward","heldout_n_complete_reverse"]].min(axis=1)
    d["min_climate_variance_fraction"]=d[
        ["thermal_unique_var_frac_forward","thermal_unique_var_frac_reverse"]].min(axis=1)
    d["min_BIO5_source_breadth_C"]=d[
        ["heldout_BIO5_q90_q10_C_forward","heldout_BIO5_q90_q10_C_reverse"]].min(axis=1)
    d["min_information_proxy"]=np.minimum(
        d.heldout_n_complete_forward*d.thermal_unique_var_frac_forward,
        d.heldout_n_complete_reverse*d.thermal_unique_var_frac_reverse)
    valid=d.loc[good]
    must(len(valid)>1500,"too few real photo taxa with both independently estimated slopes")
    must(np.isfinite(valid.min_information_proxy).all()
         and np.all(valid.min_information_proxy>0),
         "climate/photo-only source sampling proxy unidentifiable")
    q=pd.qcut(valid.min_information_proxy,4,labels=Q,duplicates="raise")
    d["climate_photo_information_quartile"]=pd.Series(
        q.astype(str).to_numpy(),index=valid.index).reindex(d.index)
    must(set(d.loc[good,"climate_photo_information_quartile"])==set(Q),
         "climate-only information strata missing")
    return d


def raw_stats(data):
    """n,sumx,sumy,sumx2,sumy2,sumxy,agree,positive_f,positive_r."""
    x=np.asarray(data.heldout_thermal_white_slope_forward,float)
    y=np.asarray(data.heldout_thermal_white_slope_reverse,float)
    must(np.isfinite(x).all() and np.isfinite(y).all(),"unstable source fold slope values")
    return np.stack((
        np.ones(len(x)),x,y,x*x,y*y,x*y,
        np.sign(x)==np.sign(y),x>0,y>0
    ),axis=1).astype(float)


def moment_metrics(m):
    z=np.asarray(m,float)
    n=z[...,0]
    sx=z[...,1];sy=z[...,2]
    xx=z[...,3]-sx*sx/n
    yy=z[...,4]-sy*sy/n
    cov=z[...,5]-sx*sy/n
    den=np.sqrt(np.maximum(xx*yy,0.))
    corr=np.divide(cov,den,out=np.full_like(den,np.nan),where=den>1e-14)
    posx=z[...,7]/n
    posy=z[...,8]/n
    chance=posx*posy+(1-posx)*(1-posy)
    same=z[...,6]/n
    return {
        "corr":corr,
        "sign_agreement":same,
        "sign_chance_from_marginals":chance,
        "excess_sign_agreement":same-chance,
        "mean_forward":sx/n,
        "mean_reverse":sy/n,
        "mean_difference_forward_minus_reverse":(sx-sy)/n,
        "n":n,
    }


def cluster_bootstrap(d,seed,rep=B):
    """Same genus weights across all sampling strata and between-quartile tests."""
    must(len(d)>=1500 and set(d.climate_photo_information_quartile)==set(Q),
         "cannot bootstrap unknown original climate source quartiles")
    clusters,genus=pd.factorize(d.genus,sort=True)
    G=len(genus)
    groupindex=np.array([Q.index(v) for v in d.climate_photo_information_quartile])
    tensor=np.zeros((G,5,len(FIELDS)),float)
    z=raw_stats(d)
    # Group indices: 0 pooled; 1-4 original sampled climate/photo opportunity quartiles
    np.add.at(tensor,(clusters,np.zeros(len(d),int)),z)
    np.add.at(tensor,(clusters,groupindex+1),z)
    totals=tensor.sum(axis=0)
    observed=moment_metrics(totals)
    rng=np.random.default_rng(seed)
    choice=rng.integers(G,size=(rep,G))
    weights=np.zeros((rep,G),float)
    np.add.at(weights,(np.repeat(np.arange(rep),G),choice.ravel()),1)
    totals_bs=(weights@tensor.reshape(G,-1)).reshape(rep,5,len(FIELDS))
    boot=moment_metrics(totals_bs)
    records={}
    for i,name in enumerate(("ALL",*Q)):
        ci={}
        for metric in ("corr","sign_agreement","excess_sign_agreement",
                       "mean_forward","mean_reverse","mean_difference_forward_minus_reverse"):
            bs=np.asarray(boot[metric][:,i])
            bs=bs[np.isfinite(bs)]
            must(len(bs)>.95*rep,metric+" source cluster bootstrap degenerate")
            ci[metric]={
                "point":float(observed[metric][i]),
                "genus_bootstrap95":[float(v) for v in np.quantile(bs,[.025,.975])],
            }
        records[name]={
            "n_source_taxa":int(observed["n"][i]),
            "metrics":ci,
            "chance_agreement_from_source_fold_sign_marginals":
                 float(observed["sign_chance_from_marginals"][i]),
        }
    contrasts={}
    for metric in ("corr","excess_sign_agreement"):
        v=np.asarray(boot[metric][:,4]-boot[metric][:,1])
        v=v[np.isfinite(v)]
        must(len(v)>.95*rep,"Q4 vs Q1 source opportunity comparison not estimable")
        lo,hi=np.quantile(v,[.025,.975])
        contrasts[metric+"_Q4_minus_Q1"]={
            "point":float(observed[metric][4]-observed[metric][1]),
            "genus_bootstrap95":[float(lo),float(hi)],
            "95_interval_excludes_zero":bool(lo>0 or hi<0),
        }
    return records,contrasts,G


def source_cohorts(d,seed):
    out={}
    for i,co in enumerate(COHORTS):
        x=d.loc[d.source_cohort==co]
        if len(x)<30:
            out[co]={"status":"not_estimable","n_source_taxa":len(x)}
            continue
        m=moment_metrics(raw_stats(x).sum(axis=0))
        out[co]={
            "n_source_taxa":len(x),
            "crossfold_thermal_photo_slope_Pearson":float(m["corr"]),
            "crossfold_photo_slope_sign_agreement":float(m["sign_agreement"]),
            "chance_sign_agreement":float(m["sign_chance_from_marginals"]),
            "mean_forward":float(m["mean_forward"]),
            "mean_reverse":float(m["mean_reverse"]),
            "not_independent_photo_provider_or_external_replication":True,
        }
    return out


def process(mode,whole):
    d=join_folds(whole,mode)
    v=d.loc[d.both_folds_estimable].copy()
    records,contrast,G=cluster_bootstrap(v,SEED+int(mode=="one_photo_per_observer")*710)
    original_n=1919 if mode=="all_heldout" else 1843
    must(len(v)==original_n,
         "immutable previous crossfit source slope opportunity changed")
    old=0.018082560930703307 if mode=="all_heldout" else -0.006113480036690145
    must(abs(records["ALL"]["metrics"]["corr"]["point"]-old)<1e-10,
         "original source fold correlation differs from independent PR271 audit")
    # Intermediate source-only aggregation output retains all original taxa.
    keep=d[[
        "inat_taxon_id","species","genus","source_cohort",
        "n_classified_original","both_folds_estimable",
        "heldout_thermal_white_slope_forward",
        "heldout_thermal_white_slope_reverse",
        "min_climate_photo_count",
        "min_climate_variance_fraction",
        "min_BIO5_source_breadth_C",
        "min_information_proxy","climate_photo_information_quartile",
    ]].copy()
    return {
        "n_original_source_taxa":len(d),
        "n_both_photo_slopes":len(v),
        "n_genera_in_both_source_folds":G,
        "quartile_edges_min_n_times_residual_climate_fraction":[float(x) for x in
          np.quantile(v.min_information_proxy,[0,.25,.5,.75,1])],
        "all_and_opportunity_quartiles":records,
        "higher_information_minus_lower_information":contrast,
        "source_cohort_separate_sample_direction_sensitivity":source_cohorts(v,SEED),
        "photo_information_proxy":"min across disjoint observer fold test panels of heldout N classifiable complete-climate photos times heldout unique 5-macroclimate thermal variance fraction; measured before reading hue slope outcomes",
        "quartiles_photo_colour_outcome_blind":True,
        "stratum_tests_posthoc_exploratory_four_quantiles_two_photo_panels":True,
        "no_biological_absence_inferred_from_zero_crossfold_correlation":True,
    },keep


def test():
    z=np.array([[1,1,2,1,4,2,1,1,1],
                [1,3,4,9,16,12,1,1,1]],float).sum(axis=0)
    m=moment_metrics(z)
    must(np.isclose(m["corr"],1), "sufficient-stat Pearson fixture failed")
    must(np.isclose(m["sign_agreement"],1), "fold sign agreement fixture failed")
    must(np.isclose(m["sign_chance_from_marginals"],1),
         "chance agreement 100% positive source fold fixture failed")
    values=np.array([.05,.10,.20,.35,.50,.75,.85,1.])
    labels=pd.qcut(values,4,labels=Q,duplicates="raise")
    must(list(labels)==["Q1","Q1","Q2","Q2","Q3","Q3","Q4","Q4"],
         "source photo-climate opportunity outcome-blind quartile fixture failed")
    print("FCP4475_FOLD_CLIMATE_OPPORTUNITY_RELIABILITY_SYNTHETIC_PASS")


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path)
    p.add_argument("--out",type=Path)
    p.add_argument("--self-test",action="store_true")
    a=p.parse_args()
    test()
    if a.self_test:return
    must(a.source and a.out,"SHA-locked original source crossfit artifact required")
    provenance,full=load(a.source)
    a.out.mkdir(parents=True,exist_ok=True)
    out={}
    for mode in MODES:
        panel,ledger=process(mode,full)
        ledger.to_csv(a.out/("fcp4475_"+mode+"_all4475_information_and_twofold_reliability.csv"),
                      index=False,float_format="%.12g")
        out[mode]=panel
    receipt={
        "schema":"fcp_original447500_photo_observer_fold_species_specific_thermal_slope_reliability_v1",
        "date_jst":"2026-10-24",
        "source_input_original_sha256":EXPECTED,
        "original_source_crossfit_run_id":40080001765,
        "original_crossfit_artifact_id":12554693263,
        "original_crossfit_artifact_zip_sha256":"50d5d15a2d145b0db7a176f945435701251768ba614f659d9dac496011b33177",
        "original_photos":447500,
        "original_taxa":4475,
        "original_nine_hue_classifiable_photos":227107,
        "sampling_opportunity_preset_not_response_based":True,
        "resampling":"1999 whole-genus bootstrap, equal taxon weights, SAME genus bootstrap replicates for quartile r contrasts within each original photo panel",
        "two_observer_photo_fold_directions_not_external_independent_sources":True,
        "source_hue_labels_not_genetic_morphs":True,
        "unclassified_photo_colour_not_imputed":True,
        "strict_petal_background_and_camera_bias_not_solved":True,
        "no_causal_temperature_or_pigment_cost_claim":True,
        "panels":out
    }
    (a.out/"observer_twofold_thermal_photo_slope_climate_opportunity_results.json").write_text(
        json.dumps(receipt,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("FCP4475_TWO_FOLD_THERMAL_RELIABILITY_COMPLETE",json.dumps({
        k:{"n":v["n_both_photo_slopes"],
           "edges":v["quartile_edges_min_n_times_residual_climate_fraction"],
           "all":v["all_and_opportunity_quartiles"]["ALL"],
           "high":v["all_and_opportunity_quartiles"]["Q4"],
           "low":v["all_and_opportunity_quartiles"]["Q1"],
           "contrast":v["higher_information_minus_lower_information"]}
        for k,v in out.items()},sort_keys=True,allow_nan=False))

if __name__=="__main__":main()

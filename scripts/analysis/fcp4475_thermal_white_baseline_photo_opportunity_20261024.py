#!/usr/bin/env python3
"""FCP genuine source 447500 photos: who drives shared warmth / white photo hue?

Join original frozen photo colour measurements to already separately CI-verified
per-species warmth slopes. Ask whether the across-species photographic thermal
association differs by original PHOTO white-visibility, photographic colour
class richness, classifiability and sampled WorldClim BIO5 breadth.

The raw source and source model artifacts are SHA-verified. All 4475 taxa stay
in the ledger. These are exploratory effect-modifier associations of the same
photo-derived colour and are NOT independent heritable colour morph/frequency,
pigment costs, causal selection or source-image background correction.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts"/"analysis"))
from fcp_full4475_hypothesis_census_20261021 import (
    SOURCE_SHA, BIO_ARCHIVE_SHA, TOTAL_ROWS,TOTAL_TAXA,CLASSIFIABLE_ROWS,
    load_original,sample_climate,all_coords_ok,PHOTO_LABELS
)

ARCHIVE_SHA="3f290b8028bbfcdb6fa2d1a4ae85f3155b21786d7c8d48175af4805c78c099a0"
SOURCE_MEMBER_SHA={
    "shared_thermal_moisture_nine_hue_source_result.json":
       "67a91c94a807e184178937d8e683467e5534b5c49ddfa24d83937a02eb5610ec",
    "fcp4475_original_thermal_moisture_support.csv":
       "c3b50053b7cbbd10c90958690fa6b855f4f8ff8cdc414c3140e45cd3463cbea1",
    "fcp4475_original_thermal_moisture_nine_hue_coefficients.npz":
       "77eca75f7909bcee53ec3d24766d6d6db1140d654b86b26fae81bb3bbcd15378",
    "fcp4475_observer_thinned_thermal_moisture_support.csv":
       "00ae2d01d372d56a7399c6ab12ecd755cef0dda827207d4f6209edcf19becb4d",
    "fcp4475_observer_thinned_thermal_moisture_nine_hue_coefficients.npz":
       "6a80cd28488af3fcb02ca305f5efc25ed19cb76f257fd18bd1802c7773bb86fc",
}
COHORTS=("discovery","validation","prospective","wave4","new_U100")
COARSE=("white","yellow_orange","red_pink","blue_purple")
STRATA=("no_white_label","rare_white_label_1to4","both_white_and_nonwhite_ge5",
        "nearly_all_white_label")
MIN_N_GENUS=15
BOOT=999
SEED=20261024


def require(test,reason):
    if not test:raise RuntimeError("FCP_THERMAL_WHITE_MODIFIER_SOURCE: "+reason)


def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):h.update(b)
    return h.hexdigest()


def source_measures(photo_csv,worldclim):
    full,d=load_original(photo_csv)
    require(len(full)==TOTAL_ROWS and len(d)==CLASSIFIABLE_ROWS
            and d.inat_taxon_id.nunique()==TOTAL_TAXA,
            "real original photo measurement population changed")
    groups=full.groupby("inat_taxon_id",sort=True)
    base=groups.agg(
        species=("species","first"),
        source_cohort=("fcp_source_cohort","first"),
        n_photo_original=("photo_id","size"),
    )
    require(base.n_photo_original.eq(100).all() and len(base)==TOTAL_TAXA,
            "source species x 100 photo archive changed")
    cb=d.groupby("inat_taxon_id",sort=True)
    nclass=cb.size().reindex(base.index,fill_value=0)
    colour=pd.crosstab(d.inat_taxon_id,d.morph).reindex(
        index=base.index,columns=COARSE,fill_value=0).astype(int)
    require(np.array_equal(colour.sum(axis=1).to_numpy(),nclass.to_numpy()),
            "photographic original class counts not exhaustive")
    base["n_classified_photo"]=nclass
    for key in COARSE:
        base["n_label_"+key]=colour[key]
    base["n_white_label_photo"]=colour["white"]
    base["n_nonwhite_label_photo"]=nclass-colour["white"]
    base["photo_classifiable_rate"]=nclass.to_numpy(float)/100
    p=nclass.to_numpy(float)
    nwhite=base.n_white_label_photo.to_numpy(int)
    notwhite=base.n_nonwhite_label_photo.to_numpy(int)
    labels=np.full(len(base),"",dtype="U64")
    labels[nwhite==0]="no_white_label"
    labels[(nwhite>=1)&(nwhite<5)]="rare_white_label_1to4"
    labels[(nwhite>=5)&(notwhite>=5)]="both_white_and_nonwhite_ge5"
    labels[(nwhite>=5)&(notwhite<5)]="nearly_all_white_label"
    require(set(labels)==set(STRATA) and np.all(labels!=""),
            "photographic white stratum pre-fixed categories missing")
    base["original_photo_white_presence_stratum"]=labels
    base["original_class_white_fraction"]=np.divide(
        nwhite,p,out=np.full_like(p,np.nan),where=p>0)
    base["n_original_colour_classes_with_five_each"]=(colour>=5).sum(axis=1).astype(int)
    with np.errstate(divide="ignore",invalid="ignore"):
        freq=colour.to_numpy(float)/np.maximum(1,p[:,None])
        base["original_four_photo_class_simpson"]=np.where(
            p>0,1-np.sum(freq**2,axis=1),np.nan)

    whitefraction=d.palette_norm9.map(lambda a:float(a[0]))
    q=d[["inat_taxon_id"]].copy()
    q["source_white_hue_fraction"]=whitefraction.to_numpy()
    means=q.groupby("inat_taxon_id")["source_white_hue_fraction"].agg(
        ["mean","std"])
    base["original_photo_ROI_mean_white_fraction"]=means["mean"].reindex(base.index)
    base["original_photo_ROI_sd_white_fraction"]=means["std"].reindex(base.index)

    climate=sample_climate(d,worldclim)
    ok=(all_coords_ok(climate)&climate.date.notna().to_numpy()
        &np.isfinite(climate.BIO5.to_numpy(float)))
    within=climate.loc[ok].groupby("inat_taxon_id",sort=True).BIO5
    moments=within.agg(["size","std","min","max"])
    q10=within.quantile(.10)
    q90=within.quantile(.90)
    base["n_original_valid_BIO5_photo"]=moments["size"].reindex(base.index)
    base["BIO5_sampled_original_std_C"]=moments["std"].reindex(base.index)
    base["BIO5_sampled_original_q90_minus_q10_C"]=(q90-q10).reindex(base.index)
    base["BIO5_sampled_original_range_C"]=(
        moments["max"]-moments["min"]).reindex(base.index)
    base=base.reset_index()
    require(base.n_classified_photo.sum()==CLASSIFIABLE_ROWS,
            "source photo four-colour measured denominator drift")
    return base


def source_model_archive(path):
    # GitHub's `gh run download` extracts a source ZIP rather than returning
    # the original zip bytes. Authenticate the provider digest through the
    # workflow and verify all five immutable INNER files' SHA256 here.
    if path.is_file():
        require(sha(path)==ARCHIVE_SHA,"shared thermal original archive zip SHA drift")
        with zipfile.ZipFile(path) as z:
            blobs={name:z.read(name) for name in SOURCE_MEMBER_SHA}
    else:
        require(path.is_dir(),"original source model artifact must be a zip or extracted directory")
        blobs={name:(path/name).read_bytes() for name in SOURCE_MEMBER_SHA}
    ret={}
    for name,expect in SOURCE_MEMBER_SHA.items():
        data=blobs[name]
        require(hashlib.sha256(data).hexdigest()==expect,
                "original derived source SHA mismatch "+name)
        if name.endswith(".json"):
            ret["receipt"]=json.loads(data)
        elif name.endswith(".csv"):
            panel="observer_thinned" if "observer_thinned" in name else "original"
            ret.setdefault(panel,{})["support"]=pd.read_csv(io.BytesIO(data),low_memory=False)
        elif name.endswith(".npz"):
            panel="observer_thinned" if "observer_thinned" in name else "original"
            with np.load(io.BytesIO(data),allow_pickle=False) as npz:
                ret.setdefault(panel,{})["ids"]=npz["inat_taxon_id"].copy()
                ret[panel]["source_thermal_white_slope"]=npz["conditional_slope"][:,0,0,0].copy()
                ret[panel]["source_photo_label_between_slope"]=npz["conditional_slope"][:,0,1,0].copy()
                ret[panel]["source_photo_label_within_slope"]=npz["conditional_slope"][:,0,2,0].copy()
    x=ret["receipt"]
    require(x["original_photos"]==TOTAL_ROWS and x["original_species"]==TOTAL_TAXA
            and x["source_classifiable_photos"]==CLASSIFIABLE_ROWS,
            "original photo model receipt incomplete")
    return ret


def join_source(photo,model,panel):
    meta=model[panel]
    s=meta["support"]
    require(len(s)==TOTAL_TAXA and s.inat_taxon_id.is_unique
            and np.array_equal(s.inat_taxon_id.to_numpy(),meta["ids"]),
            "old model 4475 taxon orders mismatch")
    df=photo.merge(s,how="inner",on="inat_taxon_id",
                   suffixes=("","_model"),validate="one_to_one",
                   sort=False)
    require(len(df)==TOTAL_TAXA and df.inat_taxon_id.is_unique,
            "taxonomy ID lost on original measured-photo merge")
    original_by_id=pd.Series(meta["source_thermal_white_slope"],
                            index=meta["ids"])
    between_by_id=pd.Series(meta["source_photo_label_between_slope"],
                           index=meta["ids"])
    within_by_id=pd.Series(meta["source_photo_label_within_slope"],
                          index=meta["ids"])
    for name,lookup in [
      ("thermal_photo_white_slope",original_by_id),
      ("thermal_photo_between_photo_label_slope",between_by_id),
      ("thermal_photo_within_photo_label_slope",within_by_id)]:
        df[name]=df.inat_taxon_id.map(lookup).to_numpy(float)
    require(np.array_equal(df.n_classified_photo.to_numpy(),
                           df.n_original_classifiable.to_numpy()),
            "source original photo eligibility baseline mismatched full model")
    require(int(df.n_original_classifiable.sum())==CLASSIFIABLE_ROWS,
            "original classified photo denominator drift")
    df["genus"]=df.species.str.split().str[0]
    return df


def genus_bootstrap_grouptable(frame,labels,seed=SEED,rep=BOOT):
    d=frame.loc[frame.thermal_photo_white_slope.notna()].copy()
    groups=[]
    for i,name in enumerate(labels):
        g=d.loc[d["_group"]==name]
        G=g.genus.nunique()
        y=g.thermal_photo_white_slope.to_numpy(float)
        if len(y)<25 or G<MIN_N_GENUS:
            groups.append({"label":str(name),"n_taxa":len(y),"n_genera":G,
                           "status":"too_few_taxa_or_genera"})
            continue
        codes,taxa=pd.factorize(g.genus,sort=True)
        G=len(taxa)
        sums=np.zeros(G,float)
        nums=np.bincount(codes,minlength=G)
        np.add.at(sums,codes,y)
        draw=np.random.default_rng(seed+51*i).integers(G,size=(rep,G))
        bootstrap=sums[draw].sum(axis=1)/nums[draw].sum(axis=1)
        ci=np.quantile(bootstrap,[.025,.975])
        groups.append({
           "label":str(name),"status":"photographic_source_descriptive",
           "n_taxa":len(g),"n_genera":G,
           "mean_thermal_photo_white_slope":float(y.mean()),
           "median_thermal_photo_white_slope":float(np.median(y)),
           "fraction_species_positive":float(np.mean(y>0)),
           "trimmed_mean_10pct":float(np.mean(np.sort(y)[int(len(y)*.1):max(int(len(y)*.9),1)])),
           "genus_bootstrap_marginal95":[float(v) for v in ci],
           "genus_bootstrap_reps":rep,
        })
    return groups


def standardized(array):
    y=np.asarray(array,float)
    require(np.isfinite(y).all(),"model covariates not fully measured")
    sd=float(y.std(ddof=0))
    return (y-y.mean())/max(sd,1e-10)


def adjusted_photo_group_effect(frame,seed=SEED):
    """Descriptive photo-label group contrasts controlling source opportunity.

    Genus cluster sufficient-statistic bootstrap; fixed covariate design:
    cohort, classified photo effort, spatial cell count, independently source
    sampled BIO5 climatic breadth, and residual thermal climate opportunity.

    Predictor group white/both labels is derived from the SAME photographs as
    response, so coefficients are not independent causal determinants.
    """
    d=frame.loc[frame.thermal_photo_white_slope.notna()].copy()
    for col in ("n_original_classifiable","n_original_raster_cells",
                "BIO5_sampled_original_q90_minus_q10_C",
                "thermal_level_unique_variance_fraction"):
        d=d.loc[np.isfinite(pd.to_numeric(d[col],errors="coerce"))].copy()
    require(len(d)>1200,"no original climate/photo opportunity after covariate audit")
    one=pd.get_dummies(d.original_photo_white_presence_stratum,
                       dtype=float).reindex(columns=list(STRATA)[1:],fill_value=0.)
    cohort=pd.get_dummies(d.source_cohort,
                          dtype=float).reindex(columns=list(COHORTS)[1:],fill_value=0.)
    x=np.column_stack((
        np.ones(len(d)),
        one.to_numpy(float),
        standardized(d.n_original_classifiable.to_numpy(float)),
        standardized(np.log1p(d.n_original_raster_cells.to_numpy(float))),
        standardized(np.log1p(d.BIO5_sampled_original_q90_minus_q10_C.to_numpy(float))),
        standardized(d.thermal_level_unique_variance_fraction.to_numpy(float)),
        cohort.to_numpy(float),
    ))
    colnames=["intercept","rare_white_v_no_white",
              "both_white_and_nonwhite_ge5_v_no_white",
              "near_all_white_v_no_white",
              "z_classifiable_photo_effort","z_log_raster_cells",
              "z_log_BIO5_photo_q90_q10_breadth",
              "z_thermal_macro_unique_variance_fraction",
              *["cohort_"+c for c in list(COHORTS)[1:]]]
    require(x.shape[1]==len(colnames)
          and np.linalg.matrix_rank(x)==x.shape[1],
          "photo opportunity regression terms are not full rank")
    y=d.thermal_photo_white_slope.to_numpy(float)
    bet=np.linalg.lstsq(x,y,rcond=1e-10)[0]
    label,genus=pd.factorize(d.genus,sort=True)
    G=len(genus)
    p=x.shape[1]
    XX=np.zeros((G,p,p),float)
    Xy=np.zeros((G,p),float)
    np.add.at(XX,label,np.einsum("ni,nj->nij",x,x))
    np.add.at(Xy,label,x*y[:,None])
    draw=np.random.default_rng(seed).integers(G,size=(BOOT,G))
    count=np.zeros((BOOT,G),float)
    np.add.at(count,(np.repeat(np.arange(BOOT),G),draw.ravel()),1)
    gram=np.einsum("bg,gij->bij",count,XX)
    target=count@Xy
    bootstrap=np.linalg.solve(gram,target)
    ci=np.quantile(bootstrap,[.025,.975],axis=0).T
    return {
        "n_original_species":int(len(d)),
        "n_genus_clusters":int(G),
        "base_category":"no original source photo classified as white",
        "regression_type":"descriptive OLS of per-species estimated white photo–thermal partial slope; original source photos also define group",
        "covariates":{name:{
              "point":float(bet[i]),
              "genus_bootstrap_marginal95":[float(v) for v in ci[i]],
              "marginal95_excludes_zero":bool(ci[i,0]>0 or ci[i,1]<0)}
            for i,name in enumerate(colnames)},
        "genus_bootstrap_reps":BOOT,
        "not_independent_genetic_white_morphs":True,
        "photographic_trait_and_response_derived_from_same_source_photo":True,
    }


def quartile_summary(frame,column,title,seed=SEED):
    d=frame.loc[frame.thermal_photo_white_slope.notna()].copy()
    vals=pd.to_numeric(d[column],errors="coerce")
    d=d.loc[vals.notna()].copy()
    if len(d)<200:
        return {"status":"not_enough_original_photo_taxa", "n_species":len(d)}
    groups=pd.qcut(d[column],q=4,labels=False,duplicates="drop")
    d["_group"]="Q"+(groups.astype(int)+1).astype(str)
    ordered=sorted(d._group.unique())
    out=genus_bootstrap_grouptable(d,ordered,seed)
    boundaries=[float(x) for x in np.unique(np.quantile(d[column],[0,.25,.5,.75,1]))]
    return {"title":title,"source_value":column,
            "n_original_taxa":len(d),"quantiles_edges":boundaries,
            "descriptive_by_quantile":out,
            "posthoc_effect_modifier_not_causal":True}


def one_panel(frame,mode):
    d=frame.loc[frame.thermal_photo_white_slope.notna()].copy()
    require(len(d)>2400,"unexpectedly few original thermal-model taxa")
    d["_group"]=d.original_photo_white_presence_stratum.astype(str)
    groups=genus_bootstrap_grouptable(d,STRATA,SEED+int(mode=="observer_thinned")*400)
    quart={}
    for i,(col,title) in enumerate([
        ("original_photo_ROI_mean_white_fraction","Photographed baseline continuous white hue share"),
        ("photo_classifiable_rate","Photo-classifiability rate among original 100"),
        ("BIO5_sampled_original_q90_minus_q10_C","Within-source-photo BIO5 sampled climate breadth (°C)"),
        ("n_original_raster_cells","Number of distinct source georeferenced climate raster cells"),
        ("thermal_level_unique_variance_fraction","Fraction independently identifiable shared thermal climatic variance"),
    ]):
        quart[col]=quartile_summary(d,col,title,SEED+int(mode=="observer_thinned")*100+i*23)
    adjusted=adjusted_photo_group_effect(frame,SEED+int(mode=="observer_thinned")*271)
    finite=d.thermal_photo_white_slope.to_numpy(float)
    return {
       "source_original_n_taxa":TOTAL_TAXA,
       "n_species_thermal_white_estimated":int(len(d)),
       "n_source_species_with_no_original_white_label":int((d.n_white_label_photo==0).sum()),
       "n_source_species_both_white_nonwhite_ge5":int(
           ((d.n_white_label_photo>=5)&(d.n_nonwhite_label_photo>=5)).sum()),
       "mean_thermal_photo_white_slope_all":float(finite.mean()),
       "mean_positive_fraction_species":float(np.mean(finite>0)),
       "original_photo_white_label_presence_subgroups":groups,
       "adjusted_white_label_modifier_and_photo_climate_opportunity":adjusted,
       "continuous_baseline_and_photo_opportunity_quartiles":quart,
       "source_photo_coarse_label_stratum_not_genetic_morphs":True,
    }


def synthetic_test():
    labels=np.array([0,0,2,2,5,10])
    others=np.array([10,2,12,1,6,0])
    categories=np.array([
         "no_white_label" if w==0 else
         "rare_white_label_1to4" if w<5 else
         "both_white_and_nonwhite_ge5" if o>=5 else
         "nearly_all_white_label" for w,o in zip(labels,others)])
    require(set(categories)==set(STRATA),"photographic white-category edge-case fixture failed")
    geom=pd.DataFrame({"genus":["A","A","B","B","C","C"]*5,
                       "thermal_photo_white_slope":[.2,.1,.05,.3,.1,.4]*5,
                       "_group":["one","one","two","two","one","two"]*5})
    result=genus_bootstrap_grouptable(geom,["one","two"])
    require(len(result)==2 and all(x["n_taxa"]>0 for x in result),
            "source descriptive genus bootstrap fixture failed")
    print("FCP4475_THERMAL_WHITE_SOURCE_MODIFIER_SYNTHETIC_PASS")


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-photos",type=Path)
    p.add_argument("--worldclim-bio",type=Path)
    p.add_argument("--thermal-artifact",type=Path)
    p.add_argument("--out",type=Path)
    p.add_argument("--self-test",action="store_true")
    args=p.parse_args()
    synthetic_test()
    if args.self_test:return
    require(args.source_photos and args.worldclim_bio and
            args.thermal_artifact and args.out,
            "source 447500 original photos, BIO raster and pinned model artifact missing")
    require(sha(args.source_photos)==SOURCE_SHA
            and sha(args.worldclim_bio)==BIO_ARCHIVE_SHA,
            "authentic FCP photograph or WorldClim source file changed")
    photo=source_measures(args.source_photos,args.worldclim_bio)
    models=source_model_archive(args.thermal_artifact)
    result={}
    args.out.mkdir(parents=True,exist_ok=True)
    for mode in ("original","observer_thinned"):
        joined=join_source(photo,models,mode)
        joined.to_csv(args.out/f"fcp4475_{mode}_thermal_white_original_photo_trait_opportunity.csv",
                      index=False,float_format="%.12g")
        result[mode]=one_panel(joined,mode)
    receipt={
      "schema":"fcp4475_photo_white_climate_slope_baseline_phenotype_and_opportunity_v1",
      "source_original_photo_SHA256":SOURCE_SHA,
      "source_original_BIO_climate_SHA256":BIO_ARCHIVE_SHA,
      "source_thermal_source_artifact_zip_sha256":ARCHIVE_SHA,
      "photo_model_source_inner_file_SHA256":SOURCE_MEMBER_SHA,
      "source_original_photos":TOTAL_ROWS,
      "source_species":TOTAL_TAXA,
      "source_classifiable_photos":CLASSIFIABLE_ROWS,
      "scope":"photographic 9-hue ROI and photographic white label; apparent within-species climate-associated photo colour slope, not genetic flower-colour morph",
      "same_photo_trait_baseline_and_thermal_response_not_independent":True,
      "genus_cluster_bootstrap_repetitions":BOOT,
      "classifiability_missingness_not_missing_at_random_verified":False,
      "observed_worldclim_BIO5_breadth_from_classified_valid_original_geophotos":True,
      "site_weather_not_measured_at_flower_anthesis":True,
      "no_native_status_photo_stage_or_background_fully_controlled":True,
      "original_all_4475_species_retained_in_source_data":True,
      "results":result
    }
    (args.out/"thermal_white_photo_baseline_heterogeneity.json").write_text(
       json.dumps(receipt,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("FCP4475_THERMAL_WHITE_BASELINE_PHOTO_HETEROGENEITY_COMPLETE",
          json.dumps({name:{
              "n":value["n_species_thermal_white_estimated"],
              "mean":value["mean_thermal_photo_white_slope_all"],
              "groups":value["original_photo_white_label_presence_subgroups"],
              "adjusted":{k:v for k,v in
                value["adjusted_white_label_modifier_and_photo_climate_opportunity"]["covariates"].items()
                if k in ("both_white_and_nonwhite_ge5_v_no_white",
                          "z_classifiable_photo_effort",
                          "z_log_BIO5_photo_q90_q10_breadth",
                          "z_thermal_macro_unique_variance_fraction")},
              "four_white_share_groups":value[
                "continuous_baseline_and_photo_opportunity_quartiles"][
                "original_photo_ROI_mean_white_fraction"]["descriptive_by_quantile"]
            } for name,value in result.items()},
            sort_keys=True,allow_nan=False))


if __name__=="__main__":
    main()

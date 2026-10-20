#!/usr/bin/env python3
"""Source-locked photo-only climate check within the SAME taxon, spatial cell and month.

Test BIO5 and BIO1 long-term WorldClim climatic normals vs white-photo labels,
permuting labels within observed spatial-cell x calendar-month strata.
This tests a *conditional photographic association*, not adaptation, genetic
IBD/IBE, true weather, historical connectivity or a pigment fitness effect.
No exact latitude/longitude, photo IDs or observer IDs are exported.
"""
from __future__ import annotations
import argparse
import hashlib
import io
import json
import math
from pathlib import Path
import zipfile

import numpy as np
import pandas as pd
from PIL import Image

PHOTO_SHA="8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6"
BIO_ZIP_SHA="00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf"
PHOTO_ROWS=447500
PHOTO_TAXA=4475
CLASSIFIED_ROWS=227107
SEED="FCP_NESTED_LOCAL_CLIMATE_BIO5_20261020"
REPS=199
COHORTS=("discovery","validation","prospective","wave4","new_U100")
COLOURS=("white","yellow_orange","red_pink","blue_purple")
BIO=(1,5)
FIELD=("inat_taxon_id","species","fcp_source_cohort","photo_id","observation_id",
       "observer_id","observed_on","latitude","longitude","positional_accuracy_m",
       "morph","global_classifiable")
PALETTE={
 "white":("flower_fraction_white",),
 "yellow_orange":("flower_fraction_yellow","flower_fraction_orange"),
 "red_pink":("flower_fraction_red","flower_fraction_pink","flower_fraction_magenta"),
 "blue_purple":("flower_fraction_blue","flower_fraction_purple"),
}
REQUIRED=list(dict.fromkeys([*FIELD,*(s for v in PALETTE.values() for s in v)]))
MIN_LABEL_PER_STRATUM=2
MIN_PHOTOS_PER_STRATUM=6
MIN_DISTINCT_OBSERVERS_PER_LABEL=2
MIN_BIO5_RANGE=0.5

def require(ok,msg):
    if not ok:raise RuntimeError("FCP_LOCAL_CLIMATE_SOURCE_OR_INFERENCE_INVALID: "+msg)

def sha256(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):h.update(chunk)
    return h.hexdigest()

def stable_seed(*parts):
    s="|".join(map(str,(SEED,*parts))).encode()
    return int.from_bytes(hashlib.sha256(s).digest()[:8],"little")

def photo_read(path):
    require(sha256(path)==PHOTO_SHA,"actual source CSV SHA256 mismatch")
    df=pd.read_csv(path,usecols=REQUIRED,low_memory=False)
    require(len(df)==PHOTO_ROWS and df.inat_taxon_id.nunique()==PHOTO_TAXA,
            "actual original photo/taxon census mismatch")
    require(df.photo_id.nunique()==PHOTO_ROWS,"original photo IDs not unique")
    require(df.groupby("inat_taxon_id").fcp_source_cohort.nunique().max()==1,
            "species observed in multiple sourcing cohorts")
    classification=df.global_classifiable.fillna(False).astype(str).str.strip().str.lower().isin(("true","yes","1"))
    d=df.loc[classification & df.morph.isin(COLOURS)].copy()
    require(len(d)==CLASSIFIED_ROWS,"classified original photo count changed")
    d["date"]=pd.to_datetime(d.observed_on,errors="coerce")
    d=d.loc[d.date.dt.year.between(2000,2026)].copy()
    d["month"]=d.date.dt.month.astype("int8")
    d["latitude"]=pd.to_numeric(d.latitude,errors="coerce")
    d["longitude"]=pd.to_numeric(d.longitude,errors="coerce")
    d["positional_accuracy_m"]=pd.to_numeric(d.positional_accuracy_m,errors="coerce")
    d=d.loc[d.latitude.between(-90,90)&d.longitude.between(-180,180)
      &d.positional_accuracy_m.between(0,5000)&d.observer_id.notna()
      &d.observation_id.notna()].copy()
    # Remove duplicate observation IDs using an outcome-blind photo-ID order.
    d=d.sort_values(["inat_taxon_id","observation_id","photo_id"],kind="mergesort")
    d=d.drop_duplicates(["inat_taxon_id","observation_id"]).reset_index(drop=True)
    return d

def validate_tiff(im):
    require(im.size==(2160,1080),"worldclim tif resolution changed")
    tags=im.tag_v2
    require(np.allclose(tags.get(33550,()),(1/6,1/6,0),atol=1e-14,rtol=0),"pixel scale changed")
    require(np.allclose(tags.get(33922,()),(0,0,0,-180,90,0),atol=1e-14,rtol=0),"raster tiepoint changed")
    require(42113 in tags,"nodata not defined")
    return float(str(tags[42113]).strip("\x00"))

def climate_sample(df,archive):
    require(sha256(archive)==BIO_ZIP_SHA,"WorldClim source ZIP SHA mismatch")
    lat=df.latitude.to_numpy(float);lon=df.longitude.to_numpy(float)
    rr=np.floor((90-lat)*6).astype("int64")
    cc=np.floor((lon+180)*6).astype("int64")
    valid=(rr>=0)&(rr<1080)&(cc>=0)&(cc<2160)
    out={}
    with zipfile.ZipFile(archive) as z:
        for index in BIO:
            member=f"wc2.1_10m_bio_{index}.tif"
            require(member in z.namelist(),"missing WorldClim member "+member)
            data=z.read(member)
            with Image.open(io.BytesIO(data)) as image:
                nodata=validate_tiff(image)
                # Preserve the original TIFF dtype before comparing its exact
                # nodata sentinel. Float32 -> float64 conversion otherwise
                # changes e.g. -3.4e38 and lets a gigantic false temperature
                # difference pass the scientific source gate.
                array=np.array(image)
            v=np.full(len(df),np.nan,dtype=float)
            sampled=array[rr[valid],cc[valid]]
            native_nodata=np.asarray(nodata,dtype=sampled.dtype).item()
            keep=(sampled!=native_nodata)&np.isfinite(sampled)
            v[valid]=np.where(keep,sampled.astype(float),np.nan)
            # BIO1 and BIO5 are degrees Celsius, not arbitrary climate indices.
            # Never silently treat missing or physically impossible raster
            # values as temperatures in a signed phenotype association.
            impossible=np.isfinite(v)&((v < -100.)|(v > 100.))
            require(int(np.sum(impossible))==0,
                    f"WorldClim BIO{index} has out-of-range unmasked temperature values")
            out[f"BIO{index}"]=v
    for field,values in out.items():df[field]=values
    return df.dropna(subset=["BIO1","BIO5"]).copy()

def high_purity_mask(d):
    components=np.column_stack([d[list(PALETTE[k])].sum(axis=1).to_numpy(float) for k in COLOURS])
    require(np.isfinite(components).all(),"ROI colour fraction invalid")
    indices=pd.Categorical(d.morph,categories=COLOURS).codes
    require(np.all(indices>=0),"unknown morph class")
    assigned=components[np.arange(len(d)),indices]
    components[np.arange(len(d)),indices]=-np.inf
    competing=np.max(components,axis=1)
    return (assigned>=.65)&(competing<=.20)

def make_cell(d,size):
    z=d.copy()
    xx=np.clip(np.floor((z.longitude.to_numpy(float)+180)/size).astype(int),0,round(360/size)-1)
    yy=np.clip(np.floor((z.latitude.to_numpy(float)+90)/size).astype(int),0,round(180/size)-1)
    z["spatialcell"]=yy*int(round(360/size))+xx
    if size==1:
        z=z.loc[z.positional_accuracy_m<=1000].copy()
    return z

def eligible_strata(d):
    records=[]
    diagnostics={"total_taxon_cells_monthly":0,"n_both_colours":0,"n_independent_observers":0,
                 "n_sufficient_BIO5_contrast":0,"n_qualified_strata":0}
    for (tid,cell,month),g in d.groupby(["inat_taxon_id","spatialcell","month"],sort=True):
        diagnostics["total_taxon_cells_monthly"]+=1
        labels=g.morph.eq("white").to_numpy()
        nw=int(labels.sum());nc=len(labels)-nw
        if len(g)<MIN_PHOTOS_PER_STRATUM or min(nw,nc)<MIN_LABEL_PER_STRATUM:continue
        diagnostics["n_both_colours"]+=1
        observers=g.observer_id.to_numpy()
        if (len(set(observers[labels]))<MIN_DISTINCT_OBSERVERS_PER_LABEL or
            len(set(observers[~labels]))<MIN_DISTINCT_OBSERVERS_PER_LABEL):continue
        diagnostics["n_independent_observers"]+=1
        vals=g[["BIO1","BIO5"]].to_numpy(float)
        if np.ptp(vals[:,1])<MIN_BIO5_RANGE:continue
        diagnostics["n_sufficient_BIO5_contrast"]+=1
        records.append((int(tid),int(cell),int(month),str(g.fcp_source_cohort.iloc[0]),
                        vals,labels))
    diagnostics["n_qualified_strata"]=len(records)
    return records,diagnostics

def stratum_obs_and_null(vals,labels,*,taxon,cell,month,size,tag):
    """Permutation keeps each stratum's exact white/photo counts and climate."""
    n=len(labels); k=int(labels.sum())
    sums=vals.sum(axis=0)
    multiplier=(1/k+1/(n-k))
    constant=sums/(n-k)
    s=vals[labels].sum(axis=0)
    observed=multiplier*s-constant
    rng=np.random.default_rng(stable_seed(taxon,cell,month,size,tag))
    perm=np.empty((REPS,2),float)
    for j in range(REPS):
        choice=rng.permutation(n)[:k]
        perm[j]=multiplier*vals[choice].sum(axis=0)-constant
    return observed,perm

def stats(values, nulls, coh):
    """Species-equal test. Each species may contribute multiple photo strata."""
    n=len(values)
    if n==0:return{"n_species":0,"status":"not_estimable"}
    v=np.stack(values);sims=np.stack(nulls)
    require(v.shape==(n,2) and sims.shape==(n,REPS,2),"paired species-by-null shapes invalid")
    obs=v.mean(axis=0);world=sims.mean(axis=0)
    rng=np.random.default_rng(stable_seed("bootstrap",coh,n))
    if n>=20:
        draws=rng.integers(0,n,size=(999,n))
        boot=v[draws].mean(axis=1)
        cis=np.quantile(boot,[.025,.975],axis=0)
        ci={name:[float(x) for x in cis[:,k]] for k,name in enumerate(("BIO1","BIO5"))}
    else:ci={"BIO1":None,"BIO5":None}
    p={}
    for k,name in enumerate(("BIO1","BIO5")):
        if n<20:p[name]=None
        else:p[name]=float((1+np.count_nonzero(world[:,k]>=obs[k]))/(REPS+1))
    return{"n_species":n,"status":"exploratory" if n>=20 else "underpowered_no_p_value",
           "mean_white_minus_nonwhite_degrees_C":{"BIO1":float(obs[0]),"BIO5":float(obs[1])},
           "median_species_effect_degrees_C":{"BIO1":float(np.median(v[:,0])),"BIO5":float(np.median(v[:,1]))},
           "positive_species_fraction":{"BIO1":float(np.mean(v[:,0]>0)),"BIO5":float(np.mean(v[:,1]>0))},
           "species_bootstrap95_conditional":ci,
           "null_upper_p_199_one_sided_white_warmer":p,
           "n_null_worlds":REPS}

def analyse(d,size,purity,tag,restrict_keys=None):
    geo=make_cell(d,size)
    if purity:geo=geo.loc[high_purity_mask(geo)].copy()
    retained_rows=int(len(geo))
    groups,diagnostics=eligible_strata(geo)
    unrestricted_n=len(groups)
    if restrict_keys is not None:
        groups=[r for r in groups if r[:3] in restrict_keys]
        require(len(groups)==len(restrict_keys),
                "selected strict-quality taxon/cell/month support changed in broad photos")
        diagnostics["n_qualified_strata_before_matched_quality"]=unrestricted_n
        diagnostics["n_qualified_strata"]=len(groups)
    by_taxon={}
    for tid,cell,month,cohort,vals,labels in groups:
        observed,sim=stratum_obs_and_null(vals,labels,
                    taxon=tid,cell=cell,month=month,size=size,tag=tag)
        x=by_taxon.setdefault(tid,{"cohort":cohort,"observed":[],"null":[]})
        require(x["cohort"]==cohort,"mismatched taxon source cohort")
        x["observed"].append(observed);x["null"].append(sim)
    ids=sorted(by_taxon)
    observed=[np.stack(by_taxon[k]["observed"]).mean(axis=0) for k in ids]
    nulls=[np.stack(by_taxon[k]["null"]).mean(axis=0) for k in ids]
    cohort=np.array([by_taxon[k]["cohort"] for k in ids])
    full=stats(observed,nulls,"all:"+tag)
    by_cohort={}
    for c in COHORTS:
        idx=np.flatnonzero(cohort==c)
        by_cohort[c]=stats([observed[k] for k in idx],[nulls[k] for k in idx],c+":"+tag)
    return{
        "cell_size_degrees":size,"condition_on_calendar_month_of_year":True,
        "matched_exact_taxon_cell_month_between_photo_quality_conditions":restrict_keys is not None,
        "photo_quality":"strict_ROI_65pct_other_le20pct" if purity else "all_four_state_classified",
        "max_source_coordinate_accuracy_m":1000 if size==1 else 5000,
        "minimum_each_white_nonwhite":MIN_LABEL_PER_STRATUM,
        "minimum_observers_each_colour":MIN_DISTINCT_OBSERVERS_PER_LABEL,
        "minimum_n_photos_each_stratum":MIN_PHOTOS_PER_STRATUM,
        "minimum_BIO5_contrast_range_C":MIN_BIO5_RANGE,
        "source_photo_rows_retained_before_stratum_gates":retained_rows,
        "opportunity_gates":diagnostics,
        "species_equal_pooled":full,
        "species_disjoint_cohorts":by_cohort,
        "n_species_with_two_plus_eligible_strata":int(sum(len(by_taxon[t]["observed"])>=2 for t in ids)),
        "no_genetic_IBE_or_population_polymorphism_inference":True,
    }

def test_float32_nodata_cast():
    raw=np.asarray([-3.4e38,20.0],dtype=np.float32)
    tag=-3.4e38
    expected=np.asarray(tag,dtype=raw.dtype).item()
    require(raw[0]==expected and np.isfinite(raw[1]),
            "original-dtype nodata sentinel mismatch")
    require(raw.astype(float)[0]!=tag,
            "test must detect float32-to-float64 nodata representation drift")


def synthetic_test():
    test_float32_nodata_cast()
    vals=np.array([[20.,29.],[20.,29.],[21.,30.],[21.,30.],
                   [18.,25.],[18.,25.]],float)
    labels=np.array([True,True,False,False,False,False])
    obs,null=stratum_obs_and_null(vals,labels,taxon=1,cell=2,month=3,size=1,tag="test")
    require(obs.shape==(2,) and null.shape==(REPS,2),"synthetic shape error")
    require(np.allclose(obs,[.5,1.5]),"true signed temperature contrast error")
    require(np.allclose(np.mean(null,axis=0),[0,0],atol=.3),"null unexpected")
    # A cell with no actual climate difference is not evidence for local sorting.
    demo=pd.DataFrame({
      "inat_taxon_id":[1]*6,"spatialcell":[4]*6,"month":[6]*6,
      "morph":["white"]*2+["blue_purple"]*4,
      "observer_id":list(range(6)),
      "BIO1":[10.]*6,"BIO5":[20.]*6,
      "fcp_source_cohort":["discovery"]*6})
    groups,d=eligible_strata(demo)
    require(len(groups)==0 and d["n_sufficient_BIO5_contrast"]==0,
            "climate contrast gate failure")
    print("FCP_NESTED_LOCAL_CLIMATE_SYNTHETIC_TESTS_PASS")

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--photos",type=Path)
    a.add_argument("--worldclim-zip",type=Path)
    a.add_argument("--out",type=Path)
    a.add_argument("--self-test",action="store_true")
    z=a.parse_args()
    if z.self_test:synthetic_test();return
    require(z.photos and z.worldclim_zip and z.out,"source and output paths required")
    synthetic_test()
    d=photo_read(z.photos)
    start_n=len(d)
    d=climate_sample(d,z.worldclim_zip)
    outputs={}
    for cell_size in (1,5):
        strict_geo=make_cell(d,cell_size)
        strict_geo=strict_geo.loc[high_purity_mask(strict_geo)].copy()
        strict_groups,_=eligible_strata(strict_geo)
        exact_strict_opportunity={(r[0],r[1],r[2]) for r in strict_groups}
        for strict in (False,True):
            tag=f"cell{cell_size}_{'strict' if strict else 'all'}"
            outputs[tag]=analyse(d,cell_size,strict,tag)
        tag=f"cell{cell_size}_all_matched_to_strict_opportunity"
        if exact_strict_opportunity:
            outputs[tag]=analyse(d,cell_size,False,tag,restrict_keys=exact_strict_opportunity)
            require(outputs[tag]["species_equal_pooled"]["n_species"]==
                    outputs[f"cell{cell_size}_strict"]["species_equal_pooled"]["n_species"],
                    "same species/cell/month photo-quality comparison failed")
        for tag in (f"cell{cell_size}_all",f"cell{cell_size}_strict",
                    f"cell{cell_size}_all_matched_to_strict_opportunity"):
            if tag not in outputs:continue
            summary=outputs[tag]["species_equal_pooled"]
            print("FCP_NESTED_LOCAL_CLIMATE_PANEL",tag,json.dumps({
              "strata":outputs[tag]["opportunity_gates"]["n_qualified_strata"],
              "n_taxa":summary["n_species"],
              "mean_bio5_degC":summary.get("mean_white_minus_nonwhite_degrees_C",{}).get("BIO5"),
              "p":summary.get("null_upper_p_199_one_sided_white_warmer",{}).get("BIO5")},sort_keys=True))
    result={
      "schema":"fcp_source_locked_within_cell_month_climate_white_photo_v1",
      "source_photo_sha256":PHOTO_SHA,
      "worldclim_zip_sha256":BIO_ZIP_SHA,
      "input_source_measured_photos":PHOTO_ROWS,
      "input_source_measured_species":PHOTO_TAXA,
      "eligible_dated_independent_obs_classifiable_before_climate":start_n,
      "eligible_dated_independent_obs_with_climate":int(len(d)),
      "date_month_conditioning":"calendar month-of-year, photographs may be from different years",
      "climate_years":"WorldClim 2.1 long-term 1970–2000 normals (not observation-year weather)",
      "primary_environment":"BIO5 maximum temperature of warmest month",
      "secondary_environment":"BIO1 annual mean temperature",
      "analysis":"same taxon, photo-lat/lon cell and calendar month, conditional white/photo label permutation",
      "no_source_species_biological_inheritance_validation":True,
      "no_native_or_cultivation_filter":True,
      "no_genetic_IBE_history_or_local_adaptation_identification":True,
      "exploratory_posthoc_multiple_comparisons":True,
      "panels":outputs,
    }
    z.out.mkdir(parents=True,exist_ok=True)
    (z.out/"local_climate_stratified_original_photo_result.json").write_text(
      json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("FCP_NESTED_LOCAL_CLIMATE_FINAL",json.dumps({
      "n_classifiable_date_eligible":start_n,
      "with_climate":len(d),
      "panels":{key:{
          "n_species":p["species_equal_pooled"]["n_species"],
          "n_strata":p["opportunity_gates"]["n_qualified_strata"],
          "bio5_C":p["species_equal_pooled"].get("mean_white_minus_nonwhite_degrees_C",{}).get("BIO5"),
          "bio5_p":p["species_equal_pooled"].get("null_upper_p_199_one_sided_white_warmer",{}).get("BIO5")
       } for key,p in outputs.items()}},sort_keys=True,allow_nan=False))

if __name__=="__main__":main()

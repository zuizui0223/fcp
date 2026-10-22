#!/usr/bin/env python3
"""Complete 4475-original-photo WorldClim 7-axis phenotype contrast.

All original photo species remain present. Within each species compare the
mean seven historical climate normals at white-labelled vs nonwhite-labelled
photo locations (also scaled by the species' observed climate SD).
Taxa with only one observed colour remain in the same 4475-row table, with
null, NEVER zero effect. This is a post hoc photographic descriptive result,
not a new physically controlled climate-selection experiment.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
import zipfile

import numpy as np
import pandas as pd
from PIL import Image

from fcp_full4475_hypothesis_census_20261021 import (
    SOURCE_SHA,BIO_ARCHIVE_SHA,TOTAL_TAXA,TOTAL_ROWS,CLASSIFIABLE_ROWS,
    load_original,sample_climate,check,
)

SRAD_ZIP_SHA="c72ee7f4f9a0eb4b5f6cd7a003eddc05bda22a2e5666d968ed4e817fd36b9026"
AXES=("BIO1","BIO4","BIO5","BIO12","BIO14","BIO15","SRADmean")
B=1999
SEED=20261021


def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for part in iter(lambda:f.read(1<<20),b""):h.update(part)
    return h.hexdigest()


def add_srad(d,archive):
    check(digest(archive)==SRAD_ZIP_SHA,"solar radiation source archive SHA mismatch")
    la=d.latitude.to_numpy(float);lo=d.longitude.to_numpy(float)
    valid=np.isfinite(la)&np.isfinite(lo)&(la>=-90)&(la<90)&(lo>=-180)&(lo<180)
    rr=np.zeros(len(la),dtype=np.int64); cc=np.zeros(len(la),dtype=np.int64)
    rr[valid]=np.floor((90-la[valid])*6).astype(np.int64)
    cc[valid]=np.floor((lo[valid]+180)*6).astype(np.int64)
    valid&=(rr>=0)&(rr<1080)&(cc>=0)&(cc<2160)
    monthly=np.empty((len(d),12),float); monthly[:]=np.nan
    with zipfile.ZipFile(archive) as z:
        for month in range(1,13):
            name=f"wc2.1_10m_srad_{month:02d}.tif"
            check(name in z.namelist(),"missing solar radiation TIFF "+name)
            with Image.open(io.BytesIO(z.read(name))) as im:
                check(im.size==(2160,1080),"srad grid resolution changed")
                tag=im.tag_v2
                check(np.allclose(tag.get(33550,()),(1/6,1/6,0),atol=1e-14,rtol=0),
                      "wrong solar radiation cell scale")
                check(np.allclose(tag.get(33922,()),(0,0,0,-180,90,0),atol=1e-14,rtol=0),
                      "wrong solar radiation grid origin")
                check(42113 in tag,"missing solar radiation nodata")
                array=np.array(im)
                nodata=np.asarray(float(str(tag[42113]).strip("\x00")),dtype=array.dtype).item()
            sample=np.full(len(d),np.nan,float)
            chunk=array[rr[valid],cc[valid]]
            sample[valid]=np.where((chunk!=nodata)&np.isfinite(chunk),chunk.astype(float),np.nan)
            check(not np.any(np.isfinite(sample)&((sample<0)|(sample>150000))),
                  "unmasked invalid solar radiation physical units")
            monthly[:,month-1]=sample
    ok=np.isfinite(monthly).all(axis=1)
    result=np.full(len(d),np.nan)
    result[ok]=np.mean(monthly[ok],axis=1)
    d["SRADmean"]=result
    return d


def per_species(g):
    ys=g.morph.eq("white").to_numpy()
    result={
      "inat_taxon_id":int(g.inat_taxon_id.iloc[0]),
      "species":str(g.species.iloc[0]),
      "cohort":str(g.fcp_source_cohort.iloc[0]),
      "n_classifiable_photo":len(g),
      "n_white_photo":int(np.sum(ys)),
      "n_nonwhite_photo":int(np.sum(~ys)),
    }
    for k in AXES:
        v=g[k].to_numpy(float)
        valid=np.isfinite(v)
        w=ys&valid
        c=(~ys)&valid
        nW=int(w.sum());nC=int(c.sum())
        result[f"{k}_n_white"]=nW
        result[f"{k}_n_nonwhite"]=nC
        result[f"{k}_signed_white_minus_nonwhite_native"]=None
        result[f"{k}_signed_white_minus_nonwhite_SD"]=None
        result[f"{k}_status"]="unidentified_no_both_colour_classes_with_normal"
        if nW<1 or nC<1:continue
        raw=float(v[w].mean()-v[c].mean())
        sd=float(np.std(v[valid]))
        result[f"{k}_signed_white_minus_nonwhite_native"]=raw
        if sd>max(1e-10,32*np.finfo(float).eps*max(1,float(np.max(np.abs(v[valid]))))):
            result[f"{k}_signed_white_minus_nonwhite_SD"]=raw/sd
            result[f"{k}_status"]="estimable"
        else:
            result[f"{k}_status"]="zero_or_tiny_within_species_normal_variance"
    return result


def taxonomy_genus_cluster_bootstrap(records):
    data=pd.DataFrame(records)
    genera=data.species.astype(str).str.split().str[0]
    unique=sorted(set(genera))
    idx={k:i for i,k in enumerate(unique)}
    # Coefficient for each species can be absent on one or many axes.
    agg_sum=np.zeros((len(unique),len(AXES)),float)
    agg_count=np.zeros((len(unique),len(AXES)),float)
    original=np.empty((len(data),len(AXES)),float)
    for j,k in enumerate(AXES):
        z=pd.to_numeric(data[f"{k}_signed_white_minus_nonwhite_SD"],errors="coerce").to_numpy(float)
        original[:,j]=z
        for genus, val in zip(genera,z):
            if np.isfinite(val):
                i=idx[genus];agg_sum[i,j]+=val;agg_count[i,j]+=1
    denom=agg_count.sum(axis=0)
    check(np.all(denom>10),"too few taxon-axis combinations for global 7-axis resampling")
    obs=agg_sum.sum(axis=0)/denom
    rng=np.random.default_rng(SEED)
    samples=np.empty((B,len(AXES)),float)
    for draw in range(B):
        choose=rng.integers(0,len(unique),size=len(unique))
        a=agg_sum[choose].sum(axis=0)
        b=agg_count[choose].sum(axis=0)
        samples[draw]=a/b
    se=np.std(samples,axis=0,ddof=1)
    check(np.all(se>0),"bootstrap SE unexpectedly degenerate")
    center=(samples-obs)/se
    maxAbs=np.max(np.abs(center),axis=1)
    threshold=float(np.quantile(maxAbs,.95))
    output=[]
    for j,k in enumerate(AXES):
        values=original[:,j];use=values[np.isfinite(values)]
        check(len(use)==int(denom[j]),"per-species estimable axis denominator differs")
        interval=[float(x) for x in np.quantile(samples[:,j],[.025,.975])]
        family=[float(obs[j]-threshold*se[j]),float(obs[j]+threshold*se[j])]
        output.append({
            "axis":k,
            "n_species_estimable":int(len(use)),
            "mean_species_equal_white_minus_nonwhite_SD":float(obs[j]),
            "median_species_white_minus_nonwhite_SD":float(np.median(use)),
            "fraction_positive_species_given_estimable":float(np.mean(use>0)),
            "genus_bootstrap_marginal95":interval,
            "genus_bootstrap_family7_maxT_95":family,
            "family7_excludes_zero":bool(family[0]>0 or family[1]<0),
        })
    return output,threshold,len(unique)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--photos",type=Path)
    ap.add_argument("--worldclim-bio",type=Path)
    ap.add_argument("--worldclim-srad",type=Path)
    ap.add_argument("--out",type=Path)
    ap.add_argument("--self-test",action="store_true")
    o=ap.parse_args()
    if o.self_test:
        g=pd.DataFrame({"inat_taxon_id":[1]*4,"species":["X y"]*4,"fcp_source_cohort":["discovery"]*4,
                        "morph":["white"]*2+["red_pink"]*2,
                        **{k:[1.,2.,3.,4.] for k in AXES}})
        z=per_species(g)
        check(all(z[f"{k}_status"]=="estimable" for k in AXES),"7-axis synthetic estimability")
        check(abs(z["BIO5_signed_white_minus_nonwhite_native"]+2.)<1e-10,
              "H3 signed phenotype climate separation wrong")
        print("FCP_ALL4475_SEVEN_WORLDCLIM_AXES_SYNTHETIC_PASS");return
    check(o.photos and o.worldclim_bio and o.worldclim_srad and o.out,
          "missing source path")
    full,photos=load_original(o.photos)
    photos=sample_climate(photos,o.worldclim_bio)
    photos=add_srad(photos,o.worldclim_srad)
    rec=[per_species(g) for _,g in photos.groupby("inat_taxon_id",sort=True)]
    check(len(rec)==TOTAL_TAXA and sum(x["n_classifiable_photo"] for x in rec)==CLASSIFIABLE_ROWS,
          "the full 4475 species/source photo table lost rows")
    axes,mult,gen=taxonomy_genus_cluster_bootstrap(rec)
    cohorts={}
    for c,g in pd.DataFrame(rec).groupby("cohort",sort=True):
        result={}
        for a in AXES:
            x=pd.to_numeric(g[f"{a}_signed_white_minus_nonwhite_SD"],errors="coerce").dropna()
            result[a]={"n_species_estimable":int(len(x)),"mean_species_equal":float(x.mean()) if len(x) else None}
        cohorts[str(c)]={"n_total_species":int(len(g)),"axes":result}
    result={
      "schema":"fcp_all4475_source_original_seven_bioclim_srad_white_label_photo_contrasts_v1",
      "full_source_447500_photo_measurements_sha256":SOURCE_SHA,
      "complete_source_4475_species":True,
      "classified_source_photos":CLASSIFIABLE_ROWS,
      "worldclim_BIO_zip_sha256":BIO_ARCHIVE_SHA,
      "worldclim_SRAD_zip_sha256":SRAD_ZIP_SHA,
      "seven_axis_multiplicity":"exploratory genus-bootstrap 7-axis simultaneous max abs studentized statistic",
      "simultaneous_family7_cutoff":mult,
      "n_genera_in_bootstrap":gen,
      "n_genera_bootstrap_draws":B,
      "axes":axes,
      "source_cohorts":cohorts,
      "no_missing_phenotype_imputation":True,
      "no_genetic_IBE_or_local_adaptation_claim":True,
      "no_annual_climate_photo_date_mismatch_claim":True,
      "no_fitness_and_biosynthetic_cost_measured":True,
    }
    o.out.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(rec).to_csv(o.out/"all4475_species_seven_environment_axes_NO_COORDINATES.csv",
                              index=False,float_format="%.12g")
    (o.out/"all4475_seven_climate_axes_result.json").write_text(
        json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("FCP_REAL_ALL4475_SEVEN_ENVIRONMENT_AXES_COMPLETE",json.dumps(result,sort_keys=True,allow_nan=False))


if __name__=="__main__":
    main()

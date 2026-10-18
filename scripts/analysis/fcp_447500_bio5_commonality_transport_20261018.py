#!/usr/bin/env python3
"""Exploratory BIO5 transport on the complete authenticated 447500 FCP photos.

Distinct from the original 149900-photo prospective white-environment filter:
no 2026-era image high-clip correction is available for the new fifth-wave
photos, and BIO5 is historical 1970-2000 *climate*, not observation weather.
This is an outcome-opened, SAME-provider photo-phenotype stress test, not
causal local adaptation or genetically confirmed population polymorphism.

Use species as equal inferential units, not photograph rows as replicates.
Always report original primary, observer-paired, and cross-cohort negative
results alongside any new U100 BIO5 association.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from scipy.stats import binomtest, spearmanr, wilcoxon

ORIGINAL_PHOTO_SHA256="8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6"
BIO5_RASTER_SHA256="32e02d85868734c32547da70cf9c620f4739a3cf5dbc4f2b751e16ebb2338ae1"
BIO5_ARCHIVE_SHA256="00513224583665ec0f2f955a4ec252730c4deb2004cce9e793492a3f26df4dcf"
SPECIES_PHOTOS=447500
CLASSIFIABLE_PHOTOS=227107
SPECIES=4475
MORPHS=("white","yellow_orange","red_pink","blue_purple")
COHORTS=("discovery","validation","prospective","wave4","new_U100")
SOURCE_COLUMNS=("inat_taxon_id","species","photo_id","fcp_source_cohort","observed_on",
                "observer_id","latitude","longitude","morph","global_classifiable",
                "background_palette_count_white","background_palette_count_yellow",
                "background_palette_count_orange","background_palette_count_red",
                "background_palette_count_pink","background_palette_count_magenta",
                "background_palette_count_purple","background_palette_count_blue",
                "background_palette_count_bronze","background_palette_count_green",
                "background_palette_count_brown","background_palette_count_black")


def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as fp:
        for block in iter(lambda:fp.read(1<<20),b""):
            h.update(block)
    return h.hexdigest()


def check(flag,why):
    if not flag:
        raise ValueError("FCP_BIO5_NOT_AUTHENTICATED: "+why)


def sample_climate(path, longitude, latitude):
    result=np.full(len(longitude),np.nan,dtype=float)
    lon=np.asarray(longitude,float)
    lat=np.asarray(latitude,float)
    ok=np.isfinite(lon)&np.isfinite(lat)&(np.abs(lat)<=90)&(np.abs(lon)<=180)
    ids=np.flatnonzero(ok)
    with rasterio.open(path) as src:
        check(src.count==1 and src.crs is not None,"no one-band georeferenced WorldClim BIO5")
        # Piecewise sampling from the exact original WorldClim raster;
        # reject nonfinite, nodata and absurd climate values.
        for indices in np.array_split(ids,max(1,int(np.ceil(len(ids)/15000)))):
            if not len(indices):
                continue
            vals=np.fromiter((float(x[0]) for x in src.sample(
                list(zip(lon[indices],lat[indices])))),float,count=len(indices))
            if src.nodata is not None:
                vals[np.isclose(vals,src.nodata)]=np.nan
            vals[(vals < -60) | (vals > 65)] = np.nan
            result[indices]=vals
    return result


def summarise_deltas(values):
    x=np.asarray(values,float)
    x=x[np.isfinite(x)]
    if len(x)==0:
        return {"estimable":False,"n_species":0,"p_wilcoxon":None,
                "p_sign_two_sided":None}
    nonzero=x[x!=0]
    wp=float(wilcoxon(x,zero_method="wilcox").pvalue) if len(nonzero) else 1.0
    signs=int(np.count_nonzero(x>0))
    signp=float(binomtest(signs,len(nonzero),.5).pvalue) if len(nonzero) else 1.0
    return {
        "estimable":True,"n_species":len(x),
        "median_white_minus_chromatic_BIO5_within_species_SD":float(np.median(x)),
        "mean_white_minus_chromatic_BIO5_within_species_SD":float(np.mean(x)),
        "fraction_positive_among_all_evaluable_species":float(np.mean(x>0)),
        "n_exact_zero":int(np.count_nonzero(x==0)),
        "p_wilcoxon":wp,"p_sign_two_sided":signp,
    }


def paired_species_strata(d,columns):
    chosen=d.dropna(subset=list(columns)+["bio5_z"]).copy()
    p=chosen.groupby([*columns,"white"]).bio5_z.mean().unstack()
    if 0 not in p.columns or 1 not in p.columns:
        return pd.Series(dtype=float),0
    p=p.dropna(subset=[0,1])
    if len(p)==0:
        return pd.Series(dtype=float),0
    value=p[1]-p[0]
    # All strata from the same taxon collapse to ONE species weight.
    return value.groupby(level=0).mean(),len(value)


def equal_genus(deltas,species_to_genus):
    df=deltas.to_frame("delta").join(species_to_genus.rename("genus"),how="inner")
    return summarise_deltas(df.groupby("genus").delta.mean())


def background_white_fraction(d):
    cols=[col for col in d.columns if col.startswith("background_palette_count_")]
    denominator=d[cols].sum(axis=1)
    d["background_white_fraction"]=(
        d["background_palette_count_white"].div(denominator.replace(0,np.nan))
    )
    return d


def analyze(photo,climate):
    check(len(photo)==SPECIES_PHOTOS and photo.inat_taxon_id.nunique()==SPECIES
          and photo.photo_id.nunique()==SPECIES_PHOTOS,
          "full original measured terminal/photo ID denominator drift")
    good=photo.global_classifiable.fillna(False).astype(str).str.lower().isin(("true","1","yes"))
    d=photo.loc[good & photo.morph.isin(MORPHS)].copy()
    check(len(d)==CLASSIFIABLE_PHOTOS,
          "four-state actual classified measured rows drift")
    d["bio5"]=sample_climate(climate,d.longitude,d.latitude)
    n_before=len(d)
    d=d.dropna(subset=["bio5"]).copy()
    d["white"]=d.morph.eq("white").astype(int)
    d["calendar_yearmonth"]=pd.to_datetime(d.observed_on,errors="coerce").dt.to_period("M").astype("string")
    d.loc[d.calendar_yearmonth.eq("NaT"),"calendar_yearmonth"]=pd.NA
    sigma=d.groupby("inat_taxon_id").bio5.transform(lambda x: x.std(ddof=0))
    mean=d.groupby("inat_taxon_id").bio5.transform("mean")
    d["bio5_z"]=(d.bio5-mean)/sigma.replace(0,np.nan)
    d=background_white_fraction(d)

    summary={
        "schema":"fcp_447500_unfrozen_BIO5_cross_cohort_transport_v1",
        "original_genuine_photo_rows":SPECIES_PHOTOS,
        "true_classified_photo_rows":CLASSIFIABLE_PHOTOS,
        "classified_photos_with_valid_BIO5":len(d),
        "classifiable_photos_without_BIO5":n_before-len(d),
        "species_with_nonconstant_sourced_BIO5":int(d.dropna(subset=["bio5_z"]).inat_taxon_id.nunique()),
        "WorldClim_2_1_10m_BIO5_1970_2000_not_photo_date_temperature":True,
        "not_independent_provider":True,
        "not_high_clip_corrected_on_new_original_photo_waves":True,
        "exploratory_post_outcome":True,
        "primary_and_historical_frozen_inference_not_overwritten":True,
    }
    d=d.dropna(subset=["bio5_z"]).copy()
    group=d.groupby("inat_taxon_id").white.agg(["sum","count"])
    eligible=group.index[(group["sum"]>=5)&((group["count"]-group["sum"])>=5)]
    e=d.loc[d.inat_taxon_id.isin(eligible)].copy()
    check(len(eligible)>5,"insufficient five-cohort within-species photographic white/chromatic support")
    contrasts=e.groupby(["inat_taxon_id","white"]).bio5_z.mean().unstack()
    check(0 in contrasts.columns and 1 in contrasts.columns,"white/nonwhite missing")
    values=(contrasts[1]-contrasts[0]).rename("white_minus_nonwhite_bio5_z")
    one=e.groupby("inat_taxon_id").fcp_source_cohort.first()
    genus=e.groupby("inat_taxon_id").species.first().map(lambda s:s.split()[0])
    report={"pooled":summarise_deltas(values),"by_cohort":{},
            "original_3_cohort_only":summarise_deltas(values.loc[one.isin(COHORTS[:3])]),
            "later_wave4_plus_new_U100_only":summarise_deltas(values.loc[one.isin(COHORTS[3:])]),
            "genus_equal_weight":equal_genus(values,genus)}
    for cohort in COHORTS:
        report["by_cohort"][cohort]=summarise_deltas(values.loc[one.eq(cohort)])
    obspaired,nobserver=paired_species_strata(e,("inat_taxon_id","observer_id"))
    monthpaired,nmonth=paired_species_strata(
        e,("inat_taxon_id","calendar_yearmonth"))
    report["observer_paired"]={
        **summarise_deltas(obspaired),
        "n_species_observer_strata_with_both_colour_states":nobserver
    }
    report["yearmonth_paired"]={
        **summarise_deltas(monthpaired),
        "n_species_yearmonth_strata_with_both_colour_states":nmonth
    }
    # A negative-control correlation from classified PHOTO BACKGROUNDS; this is
    # not a substitute for the original response-blind high-clip/ROI validation.
    bgcolor=e.groupby(["inat_taxon_id","white"]).background_white_fraction.mean().unstack()
    if 0 in bgcolor.columns and 1 in bgcolor.columns:
        bgdelta=(bgcolor[1]-bgcolor[0]).rename("delta_background_white")
        m=pd.concat([values,bgdelta],axis=1).dropna()
        if len(m)>5:
            rho,p=spearmanr(m.iloc[:,0],m.iloc[:,1])
            report["background_white_control"]={
                "n_species_with_background_pixels":len(m),
                "rho_BIO5_delta_vs_background_whiteness_delta":float(rho) if np.isfinite(rho) else None,
                "p_spearman_exploratory":float(p) if np.isfinite(p) else None,
            }
    report["eligibility"]={
        "n_white_chromatic_species_ge5_each_and_BIO5_varies":len(eligible),
        "n_classified_photos_in_evaluable_species":len(e),
        "n_observer_paired_species":len(obspaired),
        "n_yearmonth_paired_species":len(monthpaired),
        "n_genera_in_main_analysis":genus.nunique(),
    }
    summary["outcomes"]=report
    summary["causal_scope"]={
        "BIO5_alone_does_not_identify_heat_selection":True,
        "BIO5_cannot_identify_genetically_inherited_FCP_morphs":True,
        "observer_and_yearmonth_paired_require_separate_replication":True,
        "same_source_extended_results_cannot_be_called_original_prospective_H2":True,
        "no_molecular_genotype_or_fitness_join_in_this_dataset":True,
    }
    # No private photo-level latitude/longitude, identifiable observer IDs or
    # observation IDs are released. Taxon-level profiles are data products.
    profiles=pd.DataFrame({
        "taxon_id":values.index.astype(int),
        "species":e.groupby("inat_taxon_id").species.first().reindex(values.index).to_numpy(),
        "cohort":one.reindex(values.index).to_numpy(),
        "genus":genus.reindex(values.index).to_numpy(),
        "delta_bio5_z":values.to_numpy(),
        "n_classified_photos":group.loc[values.index,"count"].to_numpy(),
        "n_white_photos":group.loc[values.index,"sum"].to_numpy(),
    })
    check(not {"latitude","longitude","observer_id","photo_id","observation_id"}.intersection(profiles.columns),
          "source observation identity/location cannot leak from comparative result")
    return summary,profiles


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--photos",required=True,type=Path)
    ap.add_argument("--bio5",required=True,type=Path)
    ap.add_argument("--out",required=True,type=Path)
    a=ap.parse_args()
    check(a.photos.is_file() and digest(a.photos)==ORIGINAL_PHOTO_SHA256,
          "photograph source SHA differs from actual completed original study")
    check(a.bio5.is_file() and digest(a.bio5)==BIO5_RASTER_SHA256,
          "climate raster SHA differs from frozen WorldClim 2.1 10m source")
    photo=pd.read_csv(a.photos,usecols=list(SOURCE_COLUMNS),low_memory=False)
    result,taxa=analyze(photo,a.bio5)
    result["source_photo_sha256"]=ORIGINAL_PHOTO_SHA256
    result["source_BIO5_tif_sha256"]=BIO5_RASTER_SHA256
    a.out.mkdir(parents=True,exist_ok=True)
    (a.out/"fcp_real_447500_BIO5_unfrozen_transport.json").write_text(
        json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n",encoding="utf-8")
    taxa.to_csv(a.out/"fcp_real_447500_BIO5_species_deltas_no_coordinates.csv",index=False)
    print("FCP_REAL_FULL447500_BIO5_TRANSPORT_TERMINAL",
          json.dumps(result["outcomes"],sort_keys=True))


if __name__=="__main__":
    main()

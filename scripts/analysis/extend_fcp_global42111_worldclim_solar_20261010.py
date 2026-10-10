#!/usr/bin/env python3
"""Source-exact WorldClim 2.1 solar radiation and expanded thermal/rain climate.

Join ONLY original 42111 measured-photo species and real source photo coords.
Read source WC2.1 historical 1970-2000 10-arc-minute BIO and srad rasters.
Srad month 1-12 means seasonal long-term solar exposure, NOT local canopy light.
Preserve every original species/photo colour state including 23654 unclassified.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from attach_fcp_global42111_true_site_climate_soil_20261009 import sample_grid

ORIGINAL_N=42111
CLASSIFIED_N=18457
SOURCE_GEO="VALID_PUBLIC_ORIGINAL_PHOTO_POINT"
ADDITIONAL_BIO=tuple(i for i in range(1,20) if i not in (1,5,12,15))
S_RAD_MONTHS=tuple(range(1,13))
ADDITIONAL_MONTHLY=("srad","wind","vapr")
SOIL_NEW_PROPERTIES=("cec","sand","silt","bdod","cfvo")
SOIL_DEPTHS=(("0-5cm",5),("5-15cm",10),("15-30cm",15))
NEW_FEATURES=tuple("wc_bio"+str(i) for i in ADDITIONAL_BIO)+tuple("wc_"+p+"_"+q for p in ADDITIONAL_MONTHLY for q in ("annual_mean","monthly_cv"))+tuple("soil_"+p+"_0_30cm_source_raw" for p in SOIL_NEW_PROPERTIES)
SCHEMA="fcp_42111_original_photo_expanded_worldclim_solar_v1"


def sha256(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):h.update(b)
    return h.hexdigest()


def extract(source:pd.DataFrame,biodir:Path,solardir:Path,soil_dir:Path,*,winddir:Path|None=None,vaprdir:Path|None=None,role:str="breadth",sampler=sample_grid)->tuple[pd.DataFrame,dict]:
    required={"inat_taxon_id","observation_id","photo_id","latitude","longitude","site_geo_status",
              "wc_bio1","wc_bio5","wc_bio12","wc_bio15","wc_elevation_m",
              "soil_pH","soil_SOC","soil_N","soil_clay","soil_available_water_proxy"}
    if not required.issubset(source):
        raise ValueError("Original unchanged species/photo/abiotic column schema missing")
    if source.photo_id.duplicated().any() or source.inat_taxon_id.nunique()!=ORIGINAL_N:
        raise ValueError("Original 42111 taxon photo-ID identity changed")
    if role=="breadth":
        if len(source)!=ORIGINAL_N or not {"morph","measurement_status"}.issubset(source):
            raise ValueError("Original one-per-species photograph breadth source missing")
        classified=(source.measurement_status.eq("classified_four_state_morph") &
                    source.morph.isin(("white","yellow_orange","red_pink","blue_purple")))
        if int(classified.sum())!=CLASSIFIED_N:
            raise ValueError("Original 18457 photo-colour classified source changed")
    elif role=="all_photo_sites":
        if len(source)!=100543 or not {"present_in_breadth","present_in_taxon_cell"}.issubset(source):
            raise ValueError("Original 100543 photo-site source missing")
        if int(source.present_in_breadth.sum())!=42111 or int(source.present_in_taxon_cell.sum())!=85337:
            raise ValueError("Original 42111/85337 source photo roles drifted")
        classified=pd.Series(False,index=source.index)
    else:
        raise ValueError("Unknown original source phenotype grain")
    good=source.site_geo_status.eq(SOURCE_GEO)
    lat=pd.to_numeric(source.latitude,errors="coerce").to_numpy(float)
    lon=pd.to_numeric(source.longitude,errors="coerce").to_numpy(float)
    if not (np.isfinite(lat[good]).all() and np.isfinite(lon[good]).all() and
            (np.abs(lat[good])<=90).all() and (np.abs(lon[good])<=180).all()):
        raise ValueError("The original valid photo coordinates are not genuine")
    if np.isfinite(lat[~good]).any() or np.isfinite(lon[~good]).any():
        raise ValueError("Unlocated original photograph carried a fabricated source site")
    d=source.copy()
    for i in ADDITIONAL_BIO:
        key="wc_bio"+str(i)
        if key in d: raise ValueError("Original measured environmental features cannot be overwritten")
        vals=sampler(biodir/f"wc2.1_10m_bio_{i}.tif",lon,lat)
        if len(vals)!=len(d):raise ValueError("Source historical BIO raster lost photo records")
        d[key]=vals
    monthly_completeness={}
    for prop in ADDITIONAL_MONTHLY:
        monthly=[]
        directory={"srad":solardir,"wind":winddir or solardir.parent/"wind","vapr":vaprdir or solardir.parent/"vapr"}[prop]
        for month in S_RAD_MONTHS:
            vals=np.asarray(sampler(directory/f"wc2.1_10m_{prop}_{month}.tif",lon,lat),float)
            if len(vals)!=len(d):
                raise ValueError("Original source monthly "+prop+" raster mismatch")
            if np.isfinite(vals).any() and (vals[np.isfinite(vals)]<0).any():
                raise ValueError("Negative physical monthly "+prop+" value")
            monthly.append(vals)
        matrix=np.stack(monthly,axis=1)
        mean_of_months=matrix.mean(axis=1)
        valid=np.isfinite(matrix).all(axis=1)&(mean_of_months>0)
        annual=np.full(len(d),np.nan,float)
        cv=np.full(len(d),np.nan,float)
        annual[valid]=mean_of_months[valid]
        cv[valid]=np.std(matrix[valid],axis=1)/annual[valid]
        d["wc_"+prop+"_annual_mean"]=annual
        d["wc_"+prop+"_monthly_cv"]=cv
        monthly_completeness[prop]=int(valid.sum())
    # CEC, sand/silt, bulk density and coarse fragments are raw provider
    # integer-scaled modelled properties, NOT direct rhizosphere measurements.
    # Keep raw original scale rather than guess unverified physical unit factors.
    for prop in SOIL_NEW_PROPERTIES:
        sublayers=[]
        for depth,weight in SOIL_DEPTHS:
            source_path=soil_dir/prop/f"{prop}_{depth}_mean_5000.tif"
            val=np.asarray(sampler(source_path,lon,lat),float)
            if len(val)!=len(d):raise ValueError("SoilGrids source layer length differs from photographed source")
            sublayers.append(val)
        soil=np.column_stack(sublayers)
        complete=np.isfinite(soil).all(axis=1)
        averaged=np.full(len(d),np.nan,float)
        weights=np.array([w for _,w in SOIL_DEPTHS],float)
        averaged[complete]=soil[complete]@weights/weights.sum()
        if np.isfinite(averaged).any() and (averaged[np.isfinite(averaged)]<0).any():
            raise ValueError("Negative source-raw soil property")
        d["soil_"+prop+"_0_30cm_source_raw"]=averaged
    if not d.loc[~good,list(NEW_FEATURES)].isna().all().all():
        raise ValueError("Missing original source coordinate was assigned artificial solar/bioclim")
    n_sun=monthly_completeness["srad"]
    all_soil_new=d[["soil_"+x+"_0_30cm_source_raw" for x in SOIL_NEW_PROPERTIES]].notna().all(axis=1)
    n_sun_classified=int((classified&d["wc_srad_annual_mean"].notna()).sum())
    n_new_complete=int((classified & d[list(NEW_FEATURES)].notna().all(axis=1)).sum())
    report={
        "schema":SCHEMA,"date_jst":"2026-10-10",
        "original_source_species":int(d.inat_taxon_id.nunique()),
        "original_source_photo_rows":len(d),
        "original_source_grain":role,
        "original_classified_photo_colours":int(classified.sum()),
        "original_photo_unclassified":int((~classified).sum()) if role=="breadth" else None,
        "new_bioclim_variable_numbers":list(ADDITIONAL_BIO),
        "new_solar_monthly_rasters":list(S_RAD_MONTHS),
        "monthly_variables":list(ADDITIONAL_MONTHLY),
        "monthly_variable_complete_source_photo_counts":monthly_completeness,
        "annual_solar_mean_unit":"kJ m^-2 day^-1 (WorldClim v2.1 source)",
        "solar_cv":"sample-location mean annual cycle, monthly population SD / annual monthly mean",
        "monthly_wind_unit":"m/s WorldClim v2.1 source",
        "monthly_vapor_pressure_unit":"kPa WorldClim v2.1 source",
        "n_original_public_photo_sites":int(good.sum()),
        "n_solar_complete_source_taxa":n_sun,
        "n_solar_complete_classified_original_taxa":n_sun_classified,
        "n_additional_sun_bio_soil_complete_classified":n_new_complete,
        "n_additional_new_soil_properties_available":int(all_soil_new.sum()),
        "source_extended_soil_properties":list(SOIL_NEW_PROPERTIES),
        "source_extended_soil_depth_weights":[5,10,15],
        "extended_soil_units":"original SoilGrids mean TIFF source stored scales (no guessed rescaling)",
        "solar_annual_mean_sample_median":float(d["wc_srad_annual_mean"].median()) if n_sun else None,
        "solar_cv_sample_median":float(d["wc_srad_monthly_cv"].median()) if n_sun else None,
        "all_original_source_photo_rows_retained":True,
        "original_photo_colours_untouched":True,
        "source_cell_centroid_never_used_as_photo_location":True,
        "no_claim_sunshine_at_photo_date_or_under_canopy":True,
    }
    return d,report


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-breadth-abiotic",required=True,type=Path)
    p.add_argument("--original-photo-sites",required=True,type=Path)
    p.add_argument("--bio-dir",required=True,type=Path)
    p.add_argument("--solar-dir",required=True,type=Path)
    p.add_argument("--soil-dir",required=True,type=Path)
    p.add_argument("--wind-dir",required=True,type=Path)
    p.add_argument("--vapr-dir",required=True,type=Path)
    p.add_argument("--bio-zip",required=True,type=Path)
    p.add_argument("--solar-zip",required=True,type=Path)
    p.add_argument("--wind-zip",required=True,type=Path)
    p.add_argument("--vapr-zip",required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    a=p.parse_args()
    orig=pd.read_csv(a.original_breadth_abiotic,low_memory=False)
    all_sites=pd.read_csv(a.original_photo_sites,low_memory=False)
    site_features,report=extract(all_sites,a.bio_dir,a.solar_dir,a.soil_dir,
                                 winddir=a.wind_dir,vaprdir=a.vapr_dir,role="all_photo_sites")
    if len(orig)!=ORIGINAL_N or orig.inat_taxon_id.nunique()!=ORIGINAL_N:
        raise RuntimeError("Original 42111 breadth photo source denominator drifted")
    classified=orig.measurement_status.eq("classified_four_state_morph") & orig.morph.isin(
        ("white","yellow_orange","red_pink","blue_purple"))
    if int(classified.sum())!=CLASSIFIED_N:
        raise RuntimeError("Original 18457 classifiable photo source drifted")
    key=["inat_taxon_id","observation_id","photo_id"]
    result=orig.merge(site_features[key+list(NEW_FEATURES)],on=key,how="left",
                      validate="one_to_one",indicator=True)
    if len(result)!=ORIGINAL_N or not result["_merge"].eq("both").all():
        raise RuntimeError("Original 42111 breadth photo original-site extension lost original images")
    result=result.drop(columns=["_merge"])
    report["original_classified_photo_colours"]=CLASSIFIED_N
    report["original_photo_unclassified"]=ORIGINAL_N-CLASSIFIED_N
    report["n_solar_complete_classified_original_taxa"]=int((classified & result["wc_srad_annual_mean"].notna()).sum())
    report["n_additional_sun_bio_soil_complete_classified"]=int(
        (classified & result[list(NEW_FEATURES)].notna().all(axis=1)).sum())
    report["original_source_42111_breadth_and_85337_taxon_cell_photos_retained"]=True
    report["worldclim_bioclim_zip_sha256"]=sha256(a.bio_zip)
    report["worldclim_srad_zip_sha256"]=sha256(a.solar_zip)
    report["worldclim_wind_zip_sha256"]=sha256(a.wind_zip)
    report["worldclim_vapr_zip_sha256"]=sha256(a.vapr_zip)
    report["source_urls"]={
        "bio":"https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_10m_bio.zip",
        "srad":"https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_10m_srad.zip",
        "wind":"https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_10m_wind.zip",
        "vapr":"https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_10m_vapr.zip"}
    a.outdir.mkdir(parents=True,exist_ok=True)
    site_features.to_csv(a.outdir/"all_original_100543_photo_sites_expanded_environment.csv.gz",index=False,
                         compression={"method":"gzip","compresslevel":9,"mtime":0})
    result.to_csv(a.outdir/"source_original_42111_species_expanded_solar_bio_soil.csv.gz",index=False,
                  compression={"method":"gzip","compresslevel":9,"mtime":0})
    (a.outdir/"result.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print(json.dumps(report,sort_keys=True),flush=True)


if __name__=="__main__":
    main()

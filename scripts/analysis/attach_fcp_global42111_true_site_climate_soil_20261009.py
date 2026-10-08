#!/usr/bin/env python3
"""Attach immutable WorldClim 2.1 and SoilGrids 5km to real FCP photo points.

Reads old flower-photo results; never infers any flower-colour class from
environment, never substitutes a 162-cell centroid for a photo position.
Each old photo contributes at most one environmental record. Preserve all
source 42111 one-photo species and 85337 observed taxon-cell cases, including
unclassified, soil-missing and geolocationally missing source rows.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

CLIMATE_VARS=("bio1","bio5","bio12","bio15")
SOIL_PROPERTIES=("phh2o","soc","nitrogen","clay","wv0033","wv1500")
SOIL_DEPTHS=(("0-5cm",5),("5-15cm",10),("15-30cm",15))
SOIL_SCALE={"phh2o":10.0,"soc":10.0,"nitrogen":100.0,
            "clay":10.0,"wv0033":10.0,"wv1500":10.0}
PHOTO_STATUS="VALID_PUBLIC_ORIGINAL_PHOTO_POINT"
BREADTH_ORIGINAL=42111
CELL_ORIGINAL=85337
MEASURED_SOURCE_SHA={
    "breadth":"38aa42123b4e9b05753020ff1a3b050f4d14dbd3de3194557ead90b75c0cc605",
    "taxon_cell":"7fddcf3449fbcdbaba63d857e4028636e30d2b330a6b39c01c67d42b0470140e",
}
KEY=("inat_taxon_id","observation_id","photo_id")
ALL_FEATURES=tuple("wc_"+x for x in CLIMATE_VARS)+(
    "wc_elevation_m", "soil_pH", "soil_SOC", "soil_N", "soil_clay", "soil_available_water_proxy",
)
NO_MORPH_RECLASSIFICATION=True


def sha(path:Path)->str:
    d=hashlib.sha256()
    with path.open("rb") as f:
        for part in iter(lambda:f.read(1<<20),b""):d.update(part)
    return d.hexdigest()


def check_original(path:Path,kind:str)->pd.DataFrame:
    if sha(path)!=MEASURED_SOURCE_SHA[kind]:
        raise ValueError(f"{kind}: frozen original measured 2026-09 source changed")
    return pd.read_csv(path,low_memory=False)


def check_point_identity(points:pd.DataFrame)->pd.DataFrame:
    required=set((*KEY,"site_geo_status","latitude","longitude",
                  "present_in_breadth","present_in_taxon_cell"))
    if not required.issubset(points):
        raise ValueError(f"Missing sealed original photo coordinate schema: {required-set(points)}")
    if points.photo_id.duplicated().any():
        raise ValueError("Original source photo IDs cannot be duplicated")
    if points.inat_taxon_id.nunique()!=BREADTH_ORIGINAL:
        raise ValueError("The old 42111 global taxa are incomplete")
    if int(points.present_in_breadth.sum())!=BREADTH_ORIGINAL or int(points.present_in_taxon_cell.sum())!=CELL_ORIGINAL:
        raise ValueError("The original 42111/85337 source usages changed")
    x=points.copy()
    good=x.site_geo_status.eq(PHOTO_STATUS)
    lat=pd.to_numeric(x.latitude,errors="coerce")
    lon=pd.to_numeric(x.longitude,errors="coerce")
    if not (lat[good].between(-90,90)&lon[good].between(-180,180)).all():
        raise RuntimeError("A valid source public point lacks a proper coordinate")
    if (lat[~good].notna()|lon[~good].notna()).any():
        raise RuntimeError("An invalid source location received fabricated coordinates")
    x["latitude"]=lat
    x["longitude"]=lon
    return x


def sample_grid(path:Path,lon:np.ndarray,lat:np.ndarray)->np.ndarray:
    import rasterio
    from rasterio.warp import transform as rio_transform
    out=np.full(len(lat),np.nan,dtype=float)
    finite=np.flatnonzero(np.isfinite(lon)&np.isfinite(lat))
    if len(finite)==0:return out
    with rasterio.open(path) as r:
        xs=lon[finite].astype(float).tolist()
        ys=lat[finite].astype(float).tolist()
        if r.crs is not None and r.crs.to_string().upper() not in ("EPSG:4326","OGC:CRS84"):
            xs,ys=rio_transform("EPSG:4326",r.crs,xs,ys)
        for i,value in zip(finite,r.sample(zip(xs,ys))):
            v=float(value[0])
            if r.nodata is not None and np.isclose(v,r.nodata):
                continue
            if np.isfinite(v):
                out[i]=v
    return out


def attach(points:pd.DataFrame,bio_dir:Path,elev:Path,
           soil_dir:Path,*,sampler=sample_grid)->pd.DataFrame:
    x=check_point_identity(points)
    lon=x.longitude.to_numpy(float)
    lat=x.latitude.to_numpy(float)
    for var in CLIMATE_VARS:
        x["wc_"+var]=sampler(bio_dir/f"wc2.1_10m_{var.replace('bio','bio_')}.tif",lon,lat)
    x["wc_elevation_m"]=sampler(elev,lon,lat)
    raw={}
    for prop in SOIL_PROPERTIES:
        vals=[]
        for depth,thickness in SOIL_DEPTHS:
            a=sampler(soil_dir/prop/f"{prop}_{depth}_mean_5000.tif",lon,lat)
            vals.append(a)
        z=np.column_stack(vals)
        good=np.isfinite(z).all(axis=1)
        final=np.full(len(x),np.nan)
        weights=np.array([v for _,v in SOIL_DEPTHS],float)
        final[good]=(z[good]@weights)/weights.sum()/SOIL_SCALE[prop]
        raw[prop]=final
    x["soil_pH"]=raw["phh2o"]
    x["soil_SOC"]=raw["soc"]
    x["soil_N"]=raw["nitrogen"]
    x["soil_clay"]=raw["clay"]
    diff=raw["wv0033"]-raw["wv1500"]
    x["soil_available_water_proxy"]=np.where(np.isfinite(diff)&(diff>=0),diff,np.nan)
    # Sanity only; retain missing values rather than choosing nearby pixels.
    for c in ALL_FEATURES:
        if c not in x:raise RuntimeError(f"Missing original fixed environmental covariate {c}")
        if not x.loc[x.site_geo_status.ne(PHOTO_STATUS),c].isna().all():
            raise RuntimeError("Non-geolocated source photo was assigned fake environmental data")
    x["environment_climate_complete"]=x[list("wc_"+v for v in CLIMATE_VARS)].notna().all(axis=1)
    x["environment_soil_complete"]=x[["soil_pH","soil_SOC","soil_N","soil_clay","soil_available_water_proxy"]].notna().all(axis=1)
    x["environment_all_complete"]=x.environment_climate_complete & x.environment_soil_complete & x.wc_elevation_m.notna()
    return x


def left_join_measurement(measured:pd.DataFrame,photo_env:pd.DataFrame,role:str)->pd.DataFrame:
    if len(measured)!=(BREADTH_ORIGINAL if role=="breadth" else CELL_ORIGINAL):
        raise ValueError("FCP source measured colour table row count changed")
    if not set(KEY).issubset(measured):
        raise ValueError("Measured photo key missing")
    f=photo_env[list(KEY)+["site_geo_status","latitude","longitude",*ALL_FEATURES,
                          "environment_climate_complete","environment_soil_complete",
                          "environment_all_complete"]].copy()
    out=measured.merge(f,on=list(KEY),how="left",validate="one_to_one",indicator=True)
    if len(out)!=len(measured) or not out["_merge"].eq("both").all():
        raise RuntimeError("Full source photo-ID environmental join dropped historic colour observations")
    return out.drop(columns=["_merge"])


def aggregate(breadth:pd.DataFrame,cell:pd.DataFrame,photo_env:pd.DataFrame)->dict:
    return {
        "schema":"fcp_global42111_original_photo_climate_soil_attachment_v1",
        "status":"ENVIRONMENTAL_COVARIATES_ATTACHED_NO_CAUSAL_COLOUR_MODEL",
        "original_species_count":BREADTH_ORIGINAL,
        "original_taxon_cell_count":CELL_ORIGINAL,
        "n_unique_original_photo_IDs":len(photo_env),
        "n_source_photo_positions_valid":int(photo_env.site_geo_status.eq(PHOTO_STATUS).sum()),
        "source_photo_environment":{k:int(photo_env[k].sum()) for k in
            ("environment_climate_complete","environment_soil_complete","environment_all_complete")},
        "breadth_classifiable_by_covariate": {
            name:int((breadth.morph.isin(("white","yellow_orange","red_pink","blue_purple")) &
                       breadth.measurement_status.eq("classified_four_state_morph") &
                       breadth[name]).sum())
            for name in ("environment_climate_complete","environment_soil_complete","environment_all_complete")
        },
        "taxon_cell_classifiable_by_covariate": {
            name:int((cell.morph.isin(("white","yellow_orange","red_pink","blue_purple")) &
                       cell.measurement_status.eq("classified_four_state_morph") &
                       cell[name]).sum())
            for name in ("environment_climate_complete","environment_soil_complete","environment_all_complete")
        },
        "source_worldclim_bio_features":list(CLIMATE_VARS),
        "soil_0_to_30cm_thickness_weights":[5,10,15],
        "soil_property_5km_source":list(SOIL_PROPERTIES),
        "soil_pH_and_nitrogen_sublayer_unit_conversions_preserved":True,
        "geographic_cell_centroid_imputation":False,
        "all_original_unclassified_photo_rows_preserved":True,
        "new_colour_label_or_photo_pixel_opened":False,
        "original_manuscript_or_prospective_2000_730_changed":False,
        "inferential_boundary":"Associative, sampled-visible-phenotype and source-missingness limited; no causal flower-colour selection inferred",
    }


def main():
    p=argparse.ArgumentParser()
    for v in ("original-coordinates","breadth","taxon-cell","bio-dir","elevation","soil-dir","outdir"):
        p.add_argument("--"+v,required=True,type=Path)
    a=p.parse_args()
    locations=pd.read_csv(a.original_coordinates,low_memory=False)
    env=attach(locations,a.bio_dir,a.elevation,a.soil_dir)
    b=left_join_measurement(check_original(a.breadth,"breadth"),env,"breadth")
    c=left_join_measurement(check_original(a.taxon_cell,"taxon_cell"),env,"taxon_cell")
    report=aggregate(b,c,env)
    a.outdir.mkdir(parents=True,exist_ok=True)
    for kind,d in (("breadth",b),("taxon_cell",c)):
        d.to_csv(a.outdir/f"{kind}_flower_colour_climate_soil.csv.gz",index=False,
                 compression={"method":"gzip","compresslevel":9,"mtime":0})
    env.to_csv(a.outdir/"all_original_photo_site_abiotic_coverage.csv.gz",index=False,
               compression={"method":"gzip","compresslevel":9,"mtime":0})
    (a.outdir/"result.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print(json.dumps(report,indent=2,sort_keys=True),flush=True)
if __name__=="__main__":
    main()

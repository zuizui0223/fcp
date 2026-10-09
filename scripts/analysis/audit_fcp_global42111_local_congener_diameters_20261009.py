#!/usr/bin/env python3
"""Distance-bounded sensitivity for 42,111-source local congeneric flower photos.

The original 162-cell support can span thousands of km. On fixed source
same-genus same-cell source groups with >=3 different species, compute
TRUE photographed-site maximum group diameter (spherical great-circle).
Audit 100, 250, 500-km nested subsets regardless of their observed flower
colour outcome. Do not choose thresholds after seeing prediction scores.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd

from audit_fcp_global42111_local_congeners_20261009 import (
    photo_populations,source_group_coverage,test_local_congeners,MIN_GROUP_TAXA
)

THRESHOLDS_KM=(100.0,250.0,500.0)
EARTH_KM=6371.0088


def maximum_great_circle_km(lat:np.ndarray,lon:np.ndarray)->float:
    lat=np.asarray(lat,float);lon=np.asarray(lon,float)
    if not (len(lat)==len(lon) and len(lat)>=2):
        raise ValueError("Group needs two or more actual source positions")
    if not (np.isfinite(lat).all() and np.isfinite(lon).all() and
            np.all(np.abs(lat)<=90) and np.all(np.abs(lon)<=180)):
        raise ValueError("Unverified original photograph source position")
    a=np.deg2rad(lat)
    b=np.deg2rad(lon)
    xyz=np.column_stack([np.cos(a)*np.cos(b),np.cos(a)*np.sin(b),np.sin(a)])
    cos=np.clip(xyz@xyz.T,-1,1)
    return float(EARTH_KM*np.arccos(np.min(cos)))


def group_geometry(pop:pd.DataFrame)->pd.DataFrame:
    if pop.inat_taxon_id.duplicated().any():
        raise ValueError("Repeated source nominal species")
    required={"genus_cell_id","genus","photo_cell_162","latitude","longitude","inat_taxon_id"}
    if not required.issubset(pop):
        raise ValueError("Missing original genus, photo ID and genuine location")
    records=[]
    for group,g in pop.groupby("genus_cell_id",sort=True):
        if len(g)<MIN_GROUP_TAXA:continue
        if len(set(g.inat_taxon_id))!=len(g) or g.photo_cell_162.nunique()!=1 or g.genus.nunique()!=1:
            raise RuntimeError("Source group merged different species/photographed regions or genera")
        diam=maximum_great_circle_km(g.latitude.to_numpy(float),g.longitude.to_numpy(float))
        records.append({"genus_cell_id":group,"source_genera":str(g.genus.iloc[0]),
                        "photo_cell_162":int(g.photo_cell_162.iloc[0]),
                        "n_distinct_original_species":int(len(g)),
                        "max_source_photo_distance_km":diam})
    return pd.DataFrame(records,columns=[
        "genus_cell_id","source_genera","photo_cell_162",
        "n_distinct_original_species","max_source_photo_distance_km"])


def assess(pop:pd.DataFrame,*,include_soil:bool)->tuple[dict,pd.DataFrame]:
    g=group_geometry(pop)
    info=source_group_coverage(pop)
    summaries={}
    for km in THRESHOLDS_KM:
        matching=g.loc[g.max_source_photo_distance_km.le(km+1e-8)]
        select=pop.genus_cell_id.isin(matching.genus_cell_id)
        chosen=pop.loc[select].copy()
        model=test_local_congeners(chosen,include_soil=include_soil)
        if len(chosen)!=int(matching.n_distinct_original_species.sum()):
            raise RuntimeError("Photo source genus-cell diameter selection changed species denominator")
        if model["original_photo_support"]["n_source_species_lacking_three_local_congeners"]!=0:
            raise RuntimeError("Original photo source diameter bound selected unsupported group")
        summaries[str(int(km))]={
            "max_allowed_original_photo_group_diameter_km":km,
            "n_original_source_genus_cell_groups":int(len(matching)),
            "n_original_source_species":int(len(chosen)),
            "n_distinct_original_genera":int(matching.source_genera.nunique()),
            "n_original_source_photo_cells":int(matching.photo_cell_162.nunique()),
            "model":model,
            "sampling_rule":"Select entire original same-genus/same-162-cell photographed group iff >=3 taxa and EVERY ORIGINAL PHOTO PAIR within diameter; no outcome- or climate-based selection",
        }
    total=int(g.n_distinct_original_species.sum()) if len(g) else 0
    return {
        "schema":"fcp_42111_true_local_congener_distance_sensitivity_v1",
        "n_original_complete_case_photo_species":len(pop),
        "n_original_genus_cell_groups_min3":len(g),
        "n_photo_species_from_original_genus_cell_groups_min3":total,
        "distance_cutoff_km_prespecified_all_reported":list(THRESHOLDS_KM),
        "original_group_all_member_pairwise_diameter":True,
        "source_colour_outcome_not_used_for_threshold_group_selection":True,
        "distance_bounded_studies":summaries,
    },g


def run(source:pd.DataFrame,*,strict=True)->tuple[dict,pd.DataFrame]:
    pops,coverage=photo_populations(source,strict=strict)
    climate,geo=assess(pops["CLIMATE_ALL"],include_soil=False)
    soil,soil_geo=assess(pops["CLIMATE_SOIL_COMPLETE"],include_soil=True)
    return {
        "schema":"fcp_global42111_same_genus_cell_distance_bounded_colour_climate_v1",
        "date_jst":"2026-10-09",
        "status":"SOURCE_PHOTO_LOCAL_GEODESIC_FEASIBILITY_AND_CONDITIONAL_PREDICTION",
        "original_species_source_denominator":len(source),
        "original_source_photo_coverage":coverage,
        "climate_photo_population":climate,
        "soil_photo_population":soil,
        "no_true_photo_colour_reclassification_or_new_raster_download":True,
        "old_1499_and_independent_2000_730_untouched":True,
        "interpretation_limits":[
            "A 100-500 km genus-cell historical source diameter is not proof of common population or adaptation",
            "Requiring all source group members inside fixed diameter excludes wide-ranging genera and is observation effort selected",
            "The heldout design knows other source congeneric species from the same area by construction; not transfer to an unobserved region",
            "All three 100/250/500 km results retained even if insufficient precision or missing group support",
            "Retrospective photographed colour comparisons do not establish inheritance or direct pigment selection",
        ],
    },geo


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-breadth-abiotic",required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    a=p.parse_args()
    d=pd.read_csv(a.original_breadth_abiotic,low_memory=False)
    result,geometry=run(d)
    a.outdir.mkdir(parents=True,exist_ok=True)
    (a.outdir/"result.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    geometry.to_csv(a.outdir/"original_climate_source_genus_cell_photo_group_diameters.csv",index=False)
    print(json.dumps(result,indent=2,sort_keys=True),flush=True)

if __name__=="__main__":
    main()

#!/usr/bin/env python3
"""Preflight for 42,111-species flower-colour x climate/soil geography join.

Reads ONLY previously measured, hash-frozen 2026-09 source tables. No new
image acquisition, no environmental rasters, no result-conditioned filtering.
Missing or ocean-centred locations remain missing, NEVER invented locations.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

INPUT_SHA256 = {
    "breadth": "38aa42123b4e9b05753020ff1a3b0504d14dbd3de3194557ead90b75c0cc605",
    "taxon_cell": "7fddcf3449fbcdbaba63d857e4028636e30d2b330a6b39c01c67d42b0470140e",
    "pairs": "26755316b2cbd07e1e2241eb4ddc1cddc85219fd7e96424df8c14df32575a4d2",
}
# The literal above is checked against a terminal result -- FAIL if edited.
EXPECTED = {"breadth": 42111, "taxon_cell": 85337, "pairs": 13416}
COLOURS = ("white", "yellow_orange", "red_pink", "blue_purple")
COORD_COLUMNS = (
    ("latitude","longitude"),("lat","lon"),("observed_latitude","observed_longitude"),
    ("latitude_observed","longitude_observed"),("decimalLatitude","decimalLongitude")
)


def sha256(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
    return h.hexdigest()


def check_input(path:Path,role:str)->pd.DataFrame:
    if sha256(path)!=INPUT_SHA256[role]:
        raise ValueError(f"{role}: original 2026-09 measured-source SHA256 mismatch")
    data=pd.read_csv(path,low_memory=False)
    if len(data)!=EXPECTED[role]:
        raise ValueError(f"{role}: original row denominator altered")
    return data


def coordinate_contract(data:pd.DataFrame)->dict:
    present=[(a,b) for a,b in COORD_COLUMNS if a in data and b in data]
    if not present:
        return {
            "status":"NO_VERIFIED_OBSERVATION_COORDINATE_COLUMNS",
            "column_names":None, "n_valid":0, "n_missing":len(data),
            "n_unique_coordinate_pairs":0,"n_conflicting_coord_alias_rows":None,
            "site_level_soil_or_climate_sampling_permitted":False,
        }
    a,b=present[0]
    lat=pd.to_numeric(data[a],errors="coerce")
    lon=pd.to_numeric(data[b],errors="coerce")
    good=(lat.between(-90,90)&lon.between(-180,180)&lat.notna()&lon.notna())
    conflict=0
    for c,d in present[1:]:
        altlat=pd.to_numeric(data[c],errors="coerce")
        altlon=pd.to_numeric(data[d],errors="coerce")
        overlapping=good & altlat.notna() & altlon.notna()
        conflict+=int(((lat[overlapping]-altlat[overlapping]).abs()>1e-5).sum()+
                      ((lon[overlapping]-altlon[overlapping]).abs()>1e-5).sum())
    valid=pd.DataFrame({"lat":lat[good].round(5),"lon":lon[good].round(5)})
    return {
        "status":"SOURCE_OBSERVATION_COORDINATES_AVAILABLE" if good.any() else "ALL_OBSERVATION_COORDINATES_MISSING",
        "column_names":[a,b],"n_valid":int(good.sum()),
        "n_missing":int((~good).sum()),
        "fraction_valid":float(good.mean()) if len(data) else None,
        "n_unique_coordinate_pairs":int(valid.drop_duplicates().shape[0]),
        "n_conflicting_coord_alias_rows":conflict,
        "site_level_soil_or_climate_sampling_permitted":bool(good.any() and conflict==0),
    }


def audit(breadth:pd.DataFrame,geo:pd.DataFrame,pairs:pd.DataFrame)->dict:
    if (len(breadth),len(geo),len(pairs)) != (42111,85337,13416):
        raise ValueError("Original denominators drifted")
    for name,d in (("breadth",breadth),("taxon_cell",geo),("pairs",pairs)):
        if "inat_taxon_id" not in d:
            raise ValueError(f"{name}: taxon ID unavailable")
    if breadth.inat_taxon_id.duplicated().any():
        raise ValueError("Breadth must contain one unique original taxon")
    if not {"cell_id","measurement_status","morph"}.issubset(geo):
        raise ValueError("Regional measured colour schema not found")
    if geo[["inat_taxon_id","cell_id"]].duplicated().any():
        raise ValueError("Expected one original flower-photo outcome per taxon-cell")
    if not set(geo.inat_taxon_id).issubset(set(breadth.inat_taxon_id)):
        raise ValueError("Taxon-cell rows not found within original global frame")
    if not (geo.cell_id.astype(int).between(0,161)).all():
        raise ValueError("The source 162-cell equal-area grid has changed")
    if len(set(pairs.inat_taxon_id))!=13416:
        raise ValueError("One cross-cell pair per original source species required")
    if not {"cell_id_1","cell_id_2","pair_state"}.issubset(pairs):
        raise ValueError("Cross-cell pair provenance missing")
    classified=(geo.morph.isin(COLOURS)&
                geo.measurement_status.astype(str).eq("classified_four_state_morph"))
    breadth_class=(breadth.morph.isin(COLOURS)&
                   breadth.measurement_status.astype(str).eq("classified_four_state_morph"))
    if int(classified.sum())!=39075 or int(breadth_class.sum())!=18457:
        raise ValueError("Historical four-state colour classification total drift")
    if pairs.cell_id_1.eq(pairs.cell_id_2).any():
        raise ValueError("Cross-cell contrast contains within-cell pair")
    rows=[]
    for role,d,n in [("breadth",breadth,42111),("taxon_cell",geo,85337),("pairs",pairs,13416)]:
        coord=coordinate_contract(d)
        rows.append({
            "source_role":role,"n_rows":n,
            "n_distinct_taxa":int(d.inat_taxon_id.nunique()),
            "available_field_names":sorted(map(str,d.columns)),
            "coordinate_audit":coord,
            "geographic_cell_id_available":"cell_id" in d,
            "photo_id_available":"photo_id" in d,
            "observation_id_available":"observation_id" in d,
            "year_or_date_columns":[c for c in d.columns if any(k in c.lower() for k in ("date","observed_on","year"))],
        })
    return {
        "schema":"fcp_global42111_climate_soil_join_preflight_v1",
        "status":"SOURCE_SCHEMA_AND_SPATIAL_OPPORTUNITY_AUDIT_ONLY",
        "original_source_commit":"2b5390ea74f8d196012499d2d35dd477f9795938",
        "source_sha256":INPUT_SHA256,
        "n_original_taxa":42111,"n_original_taxon_cell_records":85337,
        "n_original_observer_disjoint_pairs":13416,
        "n_breadth_colour_classified":18457,
        "n_taxon_cell_colour_classified":39075,
        "n_breadth_unclassifiable":23654,
        "n_taxon_cell_unclassifiable":46262,
        "n_world_cells":162,"n_occupied_cells":int(geo.cell_id.nunique()),
        "layers":rows,
        "global_climate_soil_covariates_computed":False,
        "historical_fcp_manuscript_modified":False,
        "prospective_new_2000_730_outcome_opened":False,
        "interpretation_rules":[
            "exact 42111 species remains denominator; photo-unclassifiable is missing, never a morph",
            "species-wide single photograph colour is NOT species natural modal phenotype",
            "taxon-cell regional observations are repeated taxa, not 85337 independent species",
            "climate/soil extraction from a broad equal-area cell centroid alone is NOT site-local ecological covariate",
            "environmental join requires actual specimen/photo latitude/longitude, or else an explicit cell-level aggregate with land-area mask and no site-scale claims",
            "soil structural missingness must be audited by latitude, habitat and original colour-classifiability before regressions",
            "ecological covariation does not establish pigment selection, plasticity, heritable morphology or causal fitness",
        ],
    }


def main():
    p=argparse.ArgumentParser()
    for role in INPUT_SHA256:
        p.add_argument("--"+role.replace("_","-"),required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    a=p.parse_args()
    dfs={name:check_input(getattr(a,name),name) for name in INPUT_SHA256}
    report=audit(dfs["breadth"],dfs["taxon_cell"],dfs["pairs"])
    a.outdir.mkdir(parents=True,exist_ok=True)
    (a.outdir/"result.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in report.items() if k!="layers"},sort_keys=True))
    for layer in report["layers"]:
        print(json.dumps(layer,sort_keys=True))
if __name__=="__main__":
    main()

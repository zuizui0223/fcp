#!/usr/bin/env python3
"""Global 162-cell flower-colour/environment atlas from 85,337 original photo sites.

Cell geometry is for geographic display only; climate and SoilGrids covariates
were previously sampled at true public original PHOTO coordinates. Include
all 162 geographical cells and all 85,337 old taxon-cell measurements, even
source-empty cells, unclassified flowers and soil-masked sites.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd

N_SOURCE=85337
N_CLASSIFIED=39075
N_TOTAL_CELLS=162
N_LONG=18
N_SINLAT=9
COLOURS=("white","yellow_orange","red_pink","blue_purple")
FEATURES=("wc_bio1","wc_bio5","wc_bio12","wc_bio15",
          "wc_elevation_m","soil_pH","soil_SOC","soil_N","soil_clay","soil_available_water_proxy")
VALID_SITE="VALID_PUBLIC_ORIGINAL_PHOTO_POINT"


def centre_and_bounds(cell_id:int)->dict:
    if not 0<=cell_id<N_TOTAL_CELLS:raise ValueError("Invalid original 162-cell range")
    r,c=divmod(cell_id,N_LONG)
    left=-180+20*c
    right=left+20
    latlo=float(np.rad2deg(np.arcsin(-1+2*r/N_SINLAT)))
    lathi=float(np.rad2deg(np.arcsin(-1+2*(r+1)/N_SINLAT)))
    midlat=float(np.rad2deg(np.arcsin(-1+2*(r+.5)/N_SINLAT)))
    return {"cell_id":cell_id,"display_only_cell_centroid_lat":midlat,
            "display_only_cell_centroid_lon":left+10,
            "display_only_cell_latitude_low":latlo,
            "display_only_cell_latitude_high":lathi,
            "display_only_cell_lon_low":left,
            "display_only_cell_lon_high":right}


def color_and_environment_atlas(source:pd.DataFrame)->tuple[dict,pd.DataFrame,pd.DataFrame,dict]:
    must={"inat_taxon_id","cell_id","morph","measurement_status","site_geo_status",
          "environment_climate_complete","environment_soil_complete","environment_all_complete",*FEATURES}
    if not must.issubset(source):raise ValueError(f"Missing real-site source fields {sorted(must-set(source))}")
    if len(source)!=N_SOURCE or source[["inat_taxon_id","cell_id"]].duplicated().any():
        raise ValueError("Original 85337 source taxon-cell denominator or unique identity failed")
    if not pd.to_numeric(source.cell_id,errors="raise").between(0,161).all():
        raise ValueError("Original geographically fixed cell was altered")
    classified=(source.morph.isin(COLOURS)&
                source.measurement_status.eq("classified_four_state_morph"))
    if int(classified.sum())!=N_CLASSIFIED:
        raise ValueError("Original 39075 source-image classified states drifted")
    geolocated=source.site_geo_status.eq(VALID_SITE)
    if source.loc[~geolocated,list(FEATURES)].notna().any().any():
        raise ValueError("Non-geolocated photo was assigned artificial climate or soil")
    climate=source.environment_climate_complete.astype(bool)
    soil=source.environment_soil_complete.astype(bool)
    allcomplete=source.environment_all_complete.astype(bool)
    if (allcomplete&(~climate|~soil|~geolocated)).any():
        raise ValueError("Impossible claimed complete photo abiotic site")
    source=source.copy()
    source["four_state_photo_classifiable"]=classified
    source["real_public_geolocation"]=geolocated
    source["classified_all_abiotic"]=classified&allcomplete
    rows=[]
    for cell in range(N_TOTAL_CELLS):
        region=source.loc[source.cell_id==cell]
        valid=region.loc[region.four_state_photo_classifiable]
        v=len(valid)
        geom=centre_and_bounds(cell)
        counts=valid.morph.value_counts().reindex(COLOURS,fill_value=0)
        row={**geom,
             "source_taxon_cell_rows":int(len(region)),
             "n_unique_taxa":int(region.inat_taxon_id.nunique()),
             "n_four_state_classified":int(v),
             "n_unclassified":int(len(region)-v),
             "fraction_classifiable":float(v/len(region)) if len(region) else None,
             "n_public_site_coordinates":int(region.real_public_geolocation.sum()),
             "n_source_climate_complete":int(climate.loc[region.index].sum()),
             "n_source_soil_complete":int(soil.loc[region.index].sum()),
             "n_source_all_abiotic_complete":int(allcomplete.loc[region.index].sum()),
             "n_colour_classified_and_full_abiotic":int(region.classified_all_abiotic.sum()),
             }
        for color in COLOURS:
            row["count_"+color]=int(counts[color])
            row["fraction_"+color+"_given_classified"]=float(counts[color]/v) if v else None
        for x in FEATURES:
            n=region[x].dropna()
            # Distributions among ALL original source public-site photo metadata:
            # avoid selecting climate sites based on colour-classifiability.
            row["median_real_photo_"+x]=float(n.median()) if len(n)>0 else None
            row["n_real_photo_"+x]=int(len(n))
        rows.append(row)
    grid=pd.DataFrame(rows).sort_values("cell_id").reset_index(drop=True)
    if len(grid)!=162 or int(grid.source_taxon_cell_rows.sum())!=85337:
        raise RuntimeError("Global cell atlas has missing source records")
    if int(grid.n_four_state_classified.sum())!=39075:
        raise RuntimeError("Classified photograph denominator does not reproduce original global result")
    if int(grid.n_source_all_abiotic_complete.sum())!=int(allcomplete.sum()):
        raise RuntimeError("Real-site environmental coverage changed while mapping")
    region=grid.copy()
    region["absolute_latitude_region"]=pd.cut(region.display_only_cell_centroid_lat.abs(),
           [0,30,60,90],labels=["0_30","30_60","60_90"],
           include_lowest=True,right=False).astype(str)
    summaries=[]
    for band,g in region.groupby("absolute_latitude_region",sort=True):
        total=int(g.source_taxon_cell_rows.sum())
        n=int(g.n_four_state_classified.sum())
        z={
            "absolute_latitude_region":band,"n_original_global_cells":int(len(g)),
            "n_occupied_cells":int(g.source_taxon_cell_rows.gt(0).sum()),
            "source_taxon_cell_photo_records":total,
            "four_state_classifiable":n,
            "classifiability_fraction":n/total if total else None,
            "soil_complete_photo_records":int(g.n_source_soil_complete.sum()),
            "full_environment_photo_records":int(g.n_source_all_abiotic_complete.sum()),
        }
        for col in COLOURS:
            c=int(g["count_"+col].sum())
            z["count_"+col]=c
            z["fraction_"+col+"_given_classified"]=c/n if n else None
        summaries.append(z)
    broad=pd.DataFrame(summaries)
    polygons=[]
    for r in rows:
        a,b=r["display_only_cell_lon_low"],r["display_only_cell_lon_high"]
        c,d=r["display_only_cell_latitude_low"],r["display_only_cell_latitude_high"]
        props={k:(None if pd.isna(v) else v) for k,v in r.items()
               if not k.startswith("median_real_photo_") or v is not None}
        polygons.append({"type":"Feature","id":r["cell_id"],"properties":props,
            "geometry":{"type":"Polygon","coordinates":[[[a,c],[b,c],[b,d],[a,d],[a,c]]]}})
    geojson={"type":"FeatureCollection","features":polygons}
    report={
        "schema":"fcp_global42111_original_taxon_cell_colour_and_abiotic_atlas_v1",
        "status":"GEOGRAPHIC_ORIGINAL_PHOTO_COLOUR_COMPOSITION_AND_TRUE_SITE_ABIOTIC_DESCRIPTIVE",
        "original_42111_species_breadth_denominator_kept_separate":True,
        "source_taxon_cell_denominator":85337,
        "source_four_state_colour_classifiable":39075,
        "n_world_equal_area_cells":162,
        "n_original_source_occupied_cells":int(grid.source_taxon_cell_rows.gt(0).sum()),
        "original_taxon_cell_photo_colour_counts":{c:int(grid["count_"+c].sum()) for c in COLOURS},
        "n_source_real_geo_site":int(grid.n_public_site_coordinates.sum()),
        "n_source_climate_complete":int(grid.n_source_climate_complete.sum()),
        "n_source_soil_complete":int(grid.n_source_soil_complete.sum()),
        "n_source_all_abiotic_complete":int(grid.n_source_all_abiotic_complete.sum()),
        "n_classified_and_full_abiotic":int(grid.n_colour_classified_and_full_abiotic.sum()),
        "site_rasters_sampled_at_true_photo_coordinates_only":True,
        "original_equal_area_centroids_display_only":True,
        "all_unclassifiable_photo_records_preserved":True,
        "causal_adaptation_claim":False,
        "new_flower_pixels_or_colors_reclassified":False,
        "hard_nonclaims":[
            "Cell composition is one photo per original taxon-cell, not a natural within-species colour-frequency map",
            "Region totals include repeated source species; they are NOT equal-species global flower-colour prevalence",
            "Environment medians are from true photographed locations in each cell, not geometric centroids",
            "No photo classification or taxonomy/observer bias is eliminated by making a map",
            "Cell-colour associations can reflect species turnover, missing photos and unequal recording effort",
        ]
    }
    return report,grid,broad,geojson


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source-taxon-cell-abiotic",required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    a=p.parse_args()
    r,grid,bands,geojson=color_and_environment_atlas(
        pd.read_csv(a.source_taxon_cell_abiotic,low_memory=False))
    a.outdir.mkdir(parents=True,exist_ok=True)
    (a.outdir/"result.json").write_text(json.dumps(r,indent=2,sort_keys=True)+"\n")
    grid.to_csv(a.outdir/"all_162_cells_photo_colour_climate_soil.csv",index=False)
    bands.to_csv(a.outdir/"three_abs_latitude_regions_colour_and_abiotic_coverage.csv",index=False)
    (a.outdir/"original_162_equal_area_photo_colour_map.geojson").write_text(
        json.dumps(geojson,ensure_ascii=False)+"\n")
    print(json.dumps(r,sort_keys=True),flush=True)
if __name__=="__main__":
    main()

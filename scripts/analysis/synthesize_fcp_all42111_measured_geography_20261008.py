#!/usr/bin/env python3
"""Combine already-MEASURED FCP 42,111-species breadth and geographic atlas.

2026-09 completed source data, not a new photo/colour processing run.
One per species breadth phenotype and one per taxon-cell geographic anchor
are DISTINCT denominators. Repeated taxon-cell photos cannot be pooled
to estimate worldwide species-normalised flower-colour composition.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

FROZEN_SOURCE_COMMIT="2b5390ea74f8d196012499d2d35dd477f9795938"
FILES_SHA={
    "breadth":"38aa42123b4e9b05753020ff1a3b050f4d14dbd3de3194557ead90b75c0cc605",
    "taxon_cell":"7fddcf3449fbcdbaba63d857e4028636e30d2b330a6b39c01c67d42b0470140e",
    "crosscell_pairs":"a5d8e677a87fc3b0771413bbb5b2aa40a1b2242acb4156813ae986518ffa7e4f",
}
MORPHS=("white","yellow_orange","red_pink","blue_purple")
N_WORLD_SPECIES=42111
N_WORLD_CELL_ROWS=85337
N_CELLS=162
N_OBSERVED_CELLS=128
N_CROSSCELL_PAIRS=13416
N_LON=18
N_SINLAT=9
REGIONS={
    "low_0_30":(0.0,30.0),
    "middle_30_60":(30.0,60.0),
    "high_60_90":(60.0,90.0000001),
}
BOOTSTRAPS=1999
SEED=20261008136


def file_sha(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for v in iter(lambda:f.read(1<<20),b""):h.update(v)
    return h.hexdigest()


def row_cell_latitude(cell_ids:np.ndarray)->np.ndarray:
    ci=np.asarray(cell_ids,dtype=int)
    if np.any(ci<0)|np.any(ci>=N_CELLS):
        raise ValueError("Cell IDs beyond original 162 equal-area globe")
    r=ci//N_LON
    frac=-1.0+(r+0.5)*(2.0/N_SINLAT)
    return np.rad2deg(np.arcsin(np.clip(frac,-1,1)))


def row_cell_lon(cell_ids:np.ndarray)->np.ndarray:
    ci=np.asarray(cell_ids,int)
    if np.any((ci<0)|(ci>=N_CELLS)):raise ValueError("Invalid equal area cell")
    return -180.0+(ci%N_LON+0.5)*(360.0/N_LON)


def classify_region(lat:np.ndarray)->np.ndarray:
    x=np.abs(np.asarray(lat,float))
    result=np.full(len(x),"UNMAPPED",dtype=object)
    for name,(low,high) in REGIONS.items():
        result[(x>=low)&(x<high)]=name
    if "UNMAPPED" in result:
        raise ValueError("Original equal-area cell centre latitude impossible")
    return result


def pair_mid_latitude(c1:np.ndarray,c2:np.ndarray)->np.ndarray:
    l1=row_cell_latitude(c1)
    l2=row_cell_latitude(c2)
    o1=row_cell_lon(c1)
    o2=row_cell_lon(c2)
    def xyz(lat,lon):
        a=np.deg2rad(lat);b=np.deg2rad(lon)
        return np.column_stack((np.cos(a)*np.cos(b),np.cos(a)*np.sin(b),np.sin(a)))
    v=xyz(l1,o1)+xyz(l2,o2)
    mag=np.linalg.norm(v,axis=1)
    # Antipodal cell centre pairs do not have a unique geographical midpoint;
    # label missing rather than assign an invented latitude.
    return np.where(mag>1e-12,np.rad2deg(np.arcsin(
        np.clip(v[:,2]/np.maximum(mag,1e-12),-1,1))),np.nan)


def region_cell_stats(frame:pd.DataFrame)->tuple[list[dict],pd.DataFrame]:
    must={"inat_taxon_id","cell_id","morph","measurement_status"}
    if not must.issubset(frame):raise ValueError(f"Missing actual taxon-cell data {must-set(frame)}")
    d=frame.copy()
    if d[["inat_taxon_id","cell_id"]].duplicated().any():
        raise RuntimeError("A taxon-cell observation appears more than once")
    d["cell_latitude"]=row_cell_latitude(d.cell_id.to_numpy(dtype=int))
    d["absolute_latitude_region"]=classify_region(d.cell_latitude.to_numpy(float))
    d["classified"]=(d.morph.astype(str).isin(MORPHS)
        &d.measurement_status.astype(str).eq("classified_four_state_morph"))
    summary=[]
    for region in REGIONS:
        g=d.loc[d.absolute_latitude_region==region]
        valid=g.loc[g.classified]
        n=len(g); v=len(valid)
        counts=valid.morph.value_counts().reindex(MORPHS,fill_value=0)
        summary.append({
            "region":region,
            "taxon_cell_opportunities":int(n),
            "distinct_observed_species":int(g.inat_taxon_id.nunique()),
            "classified_taxon_cell":int(v),
            "unclassified_taxon_cell":int(n-v),
            "classified_fraction":float(v/n) if n else None,
            "occupied_equal_area_cells":int(g.cell_id.nunique()),
            "colours_among_classified":{k:int(counts[k]) for k in MORPHS},
            "species_equal_within_region_photo_colour_fraction_conditional":{
                k:float(counts[k]/v) if v else None for k in MORPHS
            },
            "population_claim_boundary":"taxon-cell anchors are NOT equal-weight species globally; some species appear in multiple cells",
        })
    cell=d.groupby(["cell_id","absolute_latitude_region"],sort=True).agg(
        n_taxa=("inat_taxon_id","size"),n_classified=("classified","sum")).reset_index()
    return summary,cell


def pair_metrics(pairs:pd.DataFrame)->tuple[list[dict],dict]:
    req={"inat_taxon_id","cell_id_1","cell_id_2","both_endpoints_classifiable",
         "pair_state","observer_id_1","observer_id_2"}
    if not req.issubset(pairs):raise ValueError(f"Missing real observer-disjoint pairs {req-set(pairs)}")
    d=pairs.copy()
    if d.inat_taxon_id.duplicated().any():
        raise RuntimeError("Fixed observer-disjoint sample must be exactly one pair per species")
    if (d.cell_id_1==d.cell_id_2).any():
        raise RuntimeError("Same-cell pairs found in cross-cell data")
    if (d.observer_id_1.astype(str)==d.observer_id_2.astype(str)).any():
        raise RuntimeError("Same-observer contamination")
    valid=d.both_endpoints_classifiable.astype(str).str.strip().str.casefold().isin(("true","1","yes"))
    if not (d.pair_state.astype(str).eq("unclassifiable")==~valid).all():
        raise RuntimeError("Fixed endpoint classification and pair-state status conflict")
    d["mid_lat"]=pair_mid_latitude(d.cell_id_1.to_numpy(int),
                                   d.cell_id_2.to_numpy(int))
    d["region"]=classify_region(np.where(d.mid_lat.isna(),0,d.mid_lat.to_numpy(float)))
    d.loc[d.mid_lat.isna(),"region"]="UNMAPPED_ANTIPODAL"
    rng=np.random.default_rng(SEED)
    summaries=[]
    for region in list(REGIONS)+["ALL"]:
        g=d if region=="ALL" else d.loc[d.region==region]
        n=len(g)
        same=int(g.pair_state.eq("same").sum())
        diff=int(g.pair_state.eq("discordant").sum())
        missing=int(g.pair_state.eq("unclassifiable").sum())
        if n!=same+diff+missing:raise RuntimeError("Pair state accounting broken")
        identified=same+diff
        y=np.r_[np.ones(diff),np.zeros(same)]
        if identified>1:
            boot=y[rng.integers(0,identified,size=(BOOTSTRAPS,identified))].mean(axis=1)
            bounds=[float(z) for z in np.quantile(boot,[.025,.975])]
        else:
            bounds=None
        summaries.append({
            "region":region,"fixed_observer_disjoint_species_pairs":int(n),
            "both_photos_classified":identified,
            "same_photo_colour":same,
            "discordant_photo_colour":diff,
            "one_or_both_unclassifiable":missing,
            "discordance_fraction_among_both_classified":float(diff/identified) if identified else None,
            "species_bootstrap_95CI_conditional":bounds,
            "no_assumption_full_fixed_pair_discordance_bounds":[float(diff/n),float((diff+missing)/n)] if n else None,
            "n_observed_species_in_region":int(g.inat_taxon_id.nunique()),
            "eligible_for_representative_regional_comparison":bool(identified>=100),
            "geographic_region_assignment":"great circle midpoint of TWO different frozen equal-area cell centroids",
        })
    return summaries,{
        "n_not_assignable_pair_midpoints":int(d.mid_lat.isna().sum()),
        "n_fixed_pairs":len(d),
        "n_identified_pairs":int((valid).sum()),
        "n_discordant_pairs":int(d.pair_state.eq("discordant").sum()),
    }


def analyze(breadth:pd.DataFrame,geo:pd.DataFrame,pairs:pd.DataFrame)->tuple[dict,pd.DataFrame]:
    if len(breadth)!=N_WORLD_SPECIES or breadth.inat_taxon_id.nunique()!=N_WORLD_SPECIES:
        raise RuntimeError("Original one-photo species breadth cannot be altered")
    if len(geo)!=N_WORLD_CELL_ROWS or geo[["inat_taxon_id","cell_id"]].duplicated().any():
        raise RuntimeError("Original 85337 taxon-cell outcome denominator drift")
    if len(pairs)!=N_CROSSCELL_PAIRS or pairs.inat_taxon_id.nunique()!=N_CROSSCELL_PAIRS:
        raise RuntimeError("Original cross-cell observer-disjoint species count drift")
    if not set(geo.inat_taxon_id.astype(int)).issubset(set(breadth.inat_taxon_id.astype(int))):
        raise RuntimeError("Some taxon-cell species missing from original world 42111 set")
    validbreadth=(breadth.morph.isin(MORPHS) &
                  breadth.measurement_status.eq("classified_four_state_morph"))
    counts=breadth.loc[validbreadth,"morph"].value_counts().reindex(MORPHS,fill_value=0)
    if int(validbreadth.sum())!=18457 or counts.to_dict()!={
        "white":8597,"yellow_orange":6323,"red_pink":2001,"blue_purple":1536}:
        raise RuntimeError("Previously published source breadth phenotype counts drift")
    rs,cells=region_cell_stats(geo)
    pair_regions,pair_audit=pair_metrics(pairs)
    if sum(x["classified_taxon_cell"] for x in rs)!=39075:
        raise RuntimeError("Cell group sums do not reproduce 39075 classified anchor observations")
    totals=next(x for x in pair_regions if x["region"]=="ALL")
    if (totals["both_photos_classified"],totals["discordant_photo_colour"],
        totals["one_or_both_unclassifiable"])!=(3196,797,10220):
        raise RuntimeError("Previous source pair totals drift")
    return {
        "schema":"fcp_all42111_existing_measured_phenotype_geography_synthesis_v1",
        "date_jst":"2026-10-08",
        "status":"EXISTING_REAL_IMAGE_MEASUREMENT_GLOBAL_AND_REGIONAL_AUDITED",
        "role":"integrated_source-frozen_discovery_not_independent_new_image_confirmation",
        "origin_commit":FROZEN_SOURCE_COMMIT,
        "input_SHA256":FILES_SHA,
        "original_species_denominator":N_WORLD_SPECIES,
        "breadth_single_original_photo_per_species":{
            "species_source_rows":N_WORLD_SPECIES,
            "species_four_colour_classifiable":int(validbreadth.sum()),
            "unclassifiable":N_WORLD_SPECIES-int(validbreadth.sum()),
            "colour_counts":{k:int(counts[k]) for k in MORPHS},
            "conditional_colour_fractions":{k:float(counts[k]/counts.sum()) for k in MORPHS}
        },
        "regional_taxon_cell_measurement":{
            "all_taxon_cell_rows":len(geo),
            "classified_taxon_cell_rows":int(sum(x["classified_taxon_cell"] for x in rs)),
            "occupied_cells":int(cells.cell_id.nunique()),
            "original_equal_area_grid":{"longitude_sections":N_LON,"sin_latitude_sections":N_SINLAT,"total_world_cells":N_CELLS},
            "abs_latitude_regions":rs,
            "note":"colour fractions across repeat-cell photos are not a worldwide species-equal estimator"
        },
        "observer_disjoint_crosscell_morph_comparison":{
            **pair_audit,
            "regional_and_global_descriptive":pair_regions,
            "prior_selected_pair_frame":True,
            "all_intervals_conditional_on_both_photos_classified":True
        },
        "hard_nonclaims":[
            "only coarse photographed floral colours; not genetic or biochemical polymorphism prevalence",
            "one original photo/species identifies a visible state, not dominant species morph",
            "the repeated cell table has different species weight per cell and cannot replace single-photo species-equal global composition",
            "only 3196 of 13416 observer-disjoint pairs were both classifiable; absence of a colour difference among the other 10220 is unknown",
            "regional differences can reflect species turnover, observer effort, taxonomic composition, image exposure, and biased classification",
            "the source data overlap prior 42111 and 1499 discovery work; these are NOT a prospectively independent confirmation set",
            "geographic distance and photo differences cannot identify natural selection, fitness balancing or white pigment synthesis costs",
            "the 42111 species are a citizen-science discovery frame, not every flowering species globally"
        ],
        "confirmatory_decisions_changed":False,
    },cells


def main()->None:
    p=argparse.ArgumentParser()
    for key in FILES_SHA:
        p.add_argument("--"+key.replace("_","-"),required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    args=p.parse_args()
    paths={k:getattr(args,k) for k in FILES_SHA}
    for name,filename in paths.items():
        if file_sha(filename)!=FILES_SHA[name]:
            raise RuntimeError(f"Previously completed 2026-09 FCP biological measurement source changed: {name}")
    breadth=pd.read_csv(paths["breadth"],low_memory=False)
    geo=pd.read_csv(paths["taxon_cell"],low_memory=False)
    pairs=pd.read_csv(paths["crosscell_pairs"],low_memory=False)
    result,cells=analyze(breadth,geo,pairs)
    args.outdir.mkdir(parents=True,exist_ok=True)
    cell_path=args.outdir/"cell_region_support.csv"
    cells.to_csv(cell_path,index=False,lineterminator="\n")
    result["cell_support_csv_SHA256"]=file_sha(cell_path)
    target=args.outdir/"result.json"
    target.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False))


if __name__=="__main__":
    main()

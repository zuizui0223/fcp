#!/usr/bin/env python3
"""Recover original photograph coordinates by exact frozen photo/observation IDs.

Only original 2026-09 measured photo IDs and 2026-09 discovered photo-location
indexes. Never impute cell centroid as plant environment. Every ambiguous
metadata match is MISSING, not borrowed from another photo or species.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

MEASURED_SHA = {
    "breadth": "38aa42123b4e9b05753020ff1a3b050f4d14dbd3de3194557ead90b75c0cc605",
    "taxon_cell": "7fddcf3449fbcdbaba63d857e4028636e30d2b330a6b39c01c67d42b0470140e",
}
INDEX_SHA = {
    "index_v1": "b60a8b1b98fcf6745344d6acc394bb6ff905b9d5bded134031402892ff2efa2e",
    "index_v2": "d7f08601ccc4d7dc9ff1943ea3e58ed903b9759277f7b536bc7423bc08189b3a",
}
KEY = ("inat_taxon_id","observation_id","photo_id")
COORD_ALIASES = (
    ("latitude","longitude"),("lat","lon"),
    ("observed_latitude","observed_longitude"),
    ("decimalLatitude","decimalLongitude"),
)


def sha(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for x in iter(lambda:f.read(1<<20),b""):h.update(x)
    return h.hexdigest()


def original(path:Path,expected:str)->pd.DataFrame:
    if sha(path)!=expected:raise ValueError(f"historical original byte SHA mismatch: {path.name}")
    return pd.read_csv(path,low_memory=False)


def actual_coordinate_pair(df:pd.DataFrame)->tuple[str,str]|None:
    for a,b in COORD_ALIASES:
        if {a,b}.issubset(df.columns):
            return a,b
    return None


def source_index(v1:pd.DataFrame,v2:pd.DataFrame)->tuple[pd.DataFrame,dict]:
    outputs=[]
    audit={}
    for tag,z in (("v1",v1),("v2",v2)):
        if not set((*KEY,"cell_id")).issubset(z.columns):
            raise ValueError(f"{tag}: taxon, observation, photograph or source cell identity missing")
        latlon=actual_coordinate_pair(z)
        out=z[list(KEY)+["cell_id"]].copy()
        if latlon:
            out["latitude"]=pd.to_numeric(z[latlon[0]],errors="coerce")
            out["longitude"]=pd.to_numeric(z[latlon[1]],errors="coerce")
        else:
            out["latitude"]=np.nan
            out["longitude"]=np.nan
        out["source_index_tag"]=tag
        audit[tag]={
            "n_index_photo_rows":len(z),
            "columns":list(map(str,z.columns)),
            "coordinate_source_fields":list(latlon) if latlon else None,
            "n_valid_coordinate_values":int(
                (out.latitude.between(-90,90)&out.longitude.between(-180,180)).sum()),
        }
        outputs.append(out)
    x=pd.concat(outputs,ignore_index=True)
    for c in KEY:
        x[c]=pd.to_numeric(x[c],errors="raise").astype("int64")
    x["cell_id"]=pd.to_numeric(x.cell_id,errors="raise").astype("int64")
    if not x.cell_id.between(0,161).all():
        raise ValueError("Source cell outside fixed global equal-area design")
    x["valid_coord"]=x.latitude.between(-90,90)&x.longitude.between(-180,180)
    dups=x.duplicated(list(KEY),keep=False)
    unstable=[]
    for _,g in x.loc[dups].groupby(list(KEY),sort=False):
        if g.cell_id.nunique()!=1:
            unstable.extend(g.index)
        elif g.valid_coord.any():
            coords=g.loc[g.valid_coord,["latitude","longitude"]].round(5).drop_duplicates()
            if len(coords)>1:
                unstable.extend(g.index)
    # ambiguous same original photo metadata is unusable for point sampling
    x=x.drop(index=unstable)
    x=x.sort_values("source_index_tag",kind="stable").drop_duplicates(list(KEY),keep="first")
    audit["combined"]={
        "n_source_index_rows":len(v1)+len(v2),
        "n_unique_stable_photo_identity_triples":len(x),
        "n_ambiguous_rows_unusable":len(unstable),
        "n_rows_with_valid_coordinate":int(x.valid_coord.sum())
    }
    return x,audit


def grid_cell_for(lat:np.ndarray,lon:np.ndarray)->np.ndarray:
    a=np.asarray(lat,float); b=np.asarray(lon,float)
    i=np.minimum(8,np.maximum(0,np.floor((np.sin(np.deg2rad(a))+1.0)*4.5).astype(int)))
    j=np.minimum(17,np.maximum(0,np.floor((b+180.0)/20.0).astype(int)))
    return 18*i+j


def join_immutable(measured:pd.DataFrame,idx:pd.DataFrame,role:str)->tuple[pd.DataFrame,dict]:
    if not set(KEY).issubset(measured.columns):
        raise ValueError(f"{role}: measured original photo and observation identities missing")
    m=measured.copy()
    for c in KEY:
        m[c]=pd.to_numeric(m[c],errors="raise").astype("int64")
    if m.duplicated(list(KEY)).any():
        raise ValueError(f"{role}: exact photo identity is duplicated")
    ref=idx[list(KEY)+["cell_id","latitude","longitude","valid_coord","source_index_tag"]].rename(
        columns={"cell_id":"index_cell_id"})
    j=m.merge(ref,on=list(KEY),how="left",validate="one_to_one",indicator=True)
    if len(j)!=len(m):
        raise RuntimeError("Non-injective original photo metadata merge")
    hits=j["_merge"].eq("both")
    if "cell_id" in j:
        matched_cell=(hits & j.cell_id.eq(j.index_cell_id))
    else:
        matched_cell=hits
    good=matched_cell & j.valid_coord.fillna(False).astype(bool)
    valid=np.flatnonzero(good.to_numpy())
    geographic_mismatch=0
    if len(valid):
        target=j.iloc[valid]
        predicted=grid_cell_for(target.latitude.to_numpy(float),target.longitude.to_numpy(float))
        geographic_mismatch=int(np.sum(predicted!=target.index_cell_id.to_numpy(int)))
        good.iloc[valid]=predicted==target.index_cell_id.to_numpy(int)
    j["original_photo_site_locatable"]=good
    j.loc[~good,["latitude","longitude"]]=np.nan
    needed=("inat_taxon_id","species","observation_id","photo_id",
            "latitude","longitude","original_photo_site_locatable")
    output=j[list(needed)+(["cell_id"] if "cell_id" in j else [])+
             (["morph","measurement_status"] if {"morph","measurement_status"}.issubset(j) else [])].copy()
    return output, {
        "n_measured_photo_records":len(m),
        "n_matched_index_id_triples":int(hits.sum()),
        "n_matched_correct_taxon_cell":int(matched_cell.sum()),
        "n_recorded_coordinate_valid_and_original_cell_consistent":int(good.sum()),
        "n_not_point_sampleable":int((~good).sum()),
        "n_source_cell_coordinate_mismatch":geographic_mismatch,
        "n_measured_four_state_and_point_sampleable":int((
            good & j.morph.isin(("white","yellow_orange","red_pink","blue_purple")) &
            j.measurement_status.eq("classified_four_state_morph")
        ).sum()),
    }


def audit(breadth:pd.DataFrame,cell:pd.DataFrame,v1:pd.DataFrame,v2:pd.DataFrame)->tuple[dict,dict]:
    if len(breadth)!=42111 or len(cell)!=85337:
        raise ValueError("Global observed-measured denominators changed")
    idx,index_report=source_index(v1,v2)
    b,bstats=join_immutable(breadth,idx,"breadth")
    c,cstats=join_immutable(cell,idx,"taxon_cell")
    result={
        "schema":"fcp_global42111_original_photo_coordinate_recovery_v1",
        "role":"source_identity_and_coordinate_coverage_only_no_climate_soil_effect",
        "source_historical_photo_index_sha":INDEX_SHA,
        "source_historical_measured_sha":MEASURED_SHA,
        "status":"ORIGINAL_SOURCE_POINT_LOCATION_RECOVERY_AUDITED",
        "source_index_metadata":index_report,
        "breadth":bstats,
        "taxon_cell":cstats,
        "n_breadth_species_denominator":42111,
        "n_taxon_cell_denominator":85337,
        "confirmatory_decisions_changed":False,
        "no_environment_raster_sampled":True,
        "no_new_pixels_or_photo_labels":True,
        "hard_nonclaims":[
            "An equal-area cell centroid must never be attached as the observed plant soil location",
            "If archived metadata lacks coordinates, no site-scale climate or soil estimates can be reported",
            "Photo coordinates when valid are citizen-science sites, not measured root-zone edaphic conditions",
            "No inferred plant genotype, population allele frequency or environmental adaptation",
            "Original 42111 one-photo per species and 85337 repeated taxon-cell denominators remain separate",
        ],
    }
    return result,{"breadth":b,"taxon_cell":c}


def main():
    p=argparse.ArgumentParser()
    for k in MEASURED_SHA:
        p.add_argument("--"+k.replace("_","-"),required=True,type=Path)
    for k in INDEX_SHA:
        p.add_argument("--"+k.replace("_","-"),required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    a=p.parse_args()
    dfs={r:original(getattr(a,r),sha) for r,sha in {**MEASURED_SHA,**INDEX_SHA}.items()}
    result,position=audit(dfs["breadth"],dfs["taxon_cell"],dfs["index_v1"],dfs["index_v2"])
    a.outdir.mkdir(parents=True,exist_ok=True)
    for k,d in position.items():
        d.to_csv(a.outdir/f"{k}_source_matched_photo_positions.csv.gz",index=False,
                 compression={"method":"gzip","mtime":0})
    (a.outdir/"result.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,indent=2,sort_keys=True))
if __name__=="__main__":
    main()

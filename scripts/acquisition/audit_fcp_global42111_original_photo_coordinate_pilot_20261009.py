#!/usr/bin/env python3
"""Bounded original-ID iNaturalist coordinate pilot for the 42111 colour atlas.

Only source-frozen public observation IDs. No images, colour classifications,
genetic/fitness data, or future 2000+730 source observations are queried.
The sample is keyed by geographic opportunity and hash, never flower colour.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import time
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request,urlopen

import numpy as np
import pandas as pd

PHOTO_ORIGIN_RUN=37855050150
MAX_RECORDS=1000
MAX_PER_GEO_CELL=6
MAX_BATCH=100
RATE_INTERVAL_SEC=1.2
TIMEOUT=35.0
ENDPOINT="https://api.inaturalist.org/v1/observations"
HEADERS={"User-Agent":"fcp-42111-original-observation-coordinate-pilot/1.0","Accept":"application/json"}


def seedrank(text:str)->str:
    return hashlib.sha256(("FCP_20261009_OUTCOME_BLIND_GEO_PILOT_V1|"+text).encode()).hexdigest()


def build_source_sample(breadth:pd.DataFrame,cells:pd.DataFrame)->pd.DataFrame:
    cols={"inat_taxon_id","observation_id","photo_id"}
    if not cols.issubset(breadth) or not cols|{"cell_id"} <= set(cells):
        raise ValueError("Unlinked original measured IDs or species-cell source")
    if len(breadth)!=42111 or len(cells)!=85337:
        raise ValueError("Original measured photo denominators changed")
    first=cells[list(cols)+["cell_id"]].copy().drop_duplicates(
        ["inat_taxon_id","observation_id","photo_id","cell_id"])
    first["rank"]=[seedrank(f"{int(k)}|{int(o)}|{int(p)}")
                   for k,o,p in zip(first.inat_taxon_id,first.observation_id,first.photo_id)]
    chosen=[]
    for _,g in first.groupby("cell_id",sort=True):
        chosen.append(g.sort_values("rank",kind="stable").head(MAX_PER_GEO_CELL))
    regional=pd.concat(chosen,ignore_index=True)
    regional["origin"]="taxon_cell"
    secondary=breadth[list(cols)].copy()
    secondary["cell_id"]=-1
    secondary["rank"]=[seedrank(f"{int(k)}|{int(o)}|{int(p)}")
                       for k,o,p in zip(secondary.inat_taxon_id,secondary.observation_id,secondary.photo_id)]
    secondary["origin"]="breadth"
    combined=pd.concat([regional,secondary.sort_values("rank",kind="stable")],ignore_index=True)
    combined=combined.drop_duplicates(["observation_id"],keep="first").head(MAX_RECORDS).copy()
    if combined.observation_id.duplicated().any() or len(combined)>MAX_RECORDS:
        raise RuntimeError("Repeated original observation or exceeded API cap")
    if not combined[["inat_taxon_id","observation_id","photo_id"]].notna().all().all():
        raise ValueError("No blank original IDs may be queried")
    return combined.sort_values(["origin","cell_id","rank"],kind="stable").reset_index(drop=True)


def public_position(record:dict)->tuple[float,float]|None:
    if record.get("coordinates_obscured") is True:
        return None
    if str(record.get("geoprivacy") or "").lower() in {"private","obscured"}:
        return None
    pos=record.get("geojson")
    if isinstance(pos,dict) and isinstance(pos.get("coordinates"),list) and len(pos["coordinates"])==2:
        lon,lat=pos["coordinates"]
    elif isinstance(record.get("location"),str) and "," in record["location"]:
        lat,lon=record["location"].split(",",1)
    else:
        return None
    try: lat=float(lat);lon=float(lon)
    except (TypeError,ValueError):return None
    if not np.isfinite(lat) or not np.isfinite(lon) or not -90<=lat<=90 or not -180<=lon<=180:
        return None
    return lat,lon


def cell_id_for(lat:float,lon:float)->int:
    row=min(8,max(0,int(np.floor((np.sin(np.deg2rad(lat))+1)*4.5))))
    col=min(17,max(0,int(np.floor((lon+180)/20.0))))
    return 18*row+col


def fetch_batch(obs_ids:list[int])->dict:
    if not 1<=len(obs_ids)<=MAX_BATCH or len(obs_ids)!=len(set(obs_ids)):
        raise ValueError("Duplicate or excessive public observation-ID request")
    url=ENDPOINT+"?"+urlencode({"per_page":200,"id":",".join(str(int(i)) for i in obs_ids)})
    with urlopen(Request(url,headers=HEADERS),timeout=TIMEOUT) as response:
        result=json.load(response)
    if not isinstance(result,dict) or not isinstance(result.get("results"),list):
        raise ValueError("Invalid source API response")
    return result


def summarize_one(original:dict, obs:dict|None)->dict:
    out={
        "inat_taxon_id":int(original["inat_taxon_id"]),
        "observation_id":int(original["observation_id"]),
        "photo_id":int(original["photo_id"]),
        "source_cell_id":int(original["cell_id"]),
        "source_sample_role":original["origin"],
        "latitude":None,"longitude":None,
        "positional_accuracy_m":None,
        "status":"UNKNOWN"
    }
    if obs is None:
        out["status"]="OBSERVATION_ID_NOT_RETURNED";return out
    if int(obs.get("id") or 0)!=out["observation_id"]:
        raise ValueError("Source observation response ID mismatch")
    taxon=obs.get("taxon")
    if not isinstance(taxon,dict) or int(taxon.get("id") or 0)!=out["inat_taxon_id"]:
        out["status"]="TAXON_IDENTITY_CHANGED";return out
    photos=obs.get("photos") or []
    if not any(isinstance(p,dict) and int(p.get("id") or 0)==out["photo_id"] for p in photos):
        out["status"]="ORIGINAL_PHOTO_ID_MISSING";return out
    if bool(obs.get("captive")):
        out["status"]="CAPTIVE_OR_CULTIVATED";return out
    loc=public_position(obs)
    if loc is None:
        out["status"]="NO_UNOBSCURED_PUBLIC_COORDINATES";return out
    lat,lon=loc
    if out["source_cell_id"]>=0 and cell_id_for(lat,lon)!=out["source_cell_id"]:
        out["status"]="LOCATION_CHANGED_OR_SOURCE_CELL_MISMATCH";return out
    acc=obs.get("positional_accuracy")
    try: acc=float(acc)
    except (TypeError,ValueError):acc=None
    out["latitude"]=lat
    out["longitude"]=lon
    out["positional_accuracy_m"]=acc if acc is not None and np.isfinite(acc) else None
    out["status"]="EXACT_SOURCE_PHOTO_PUBLIC_COORDINATE_AVAILABLE"
    return out


def pilot(sample:pd.DataFrame, *, client=fetch_batch, pause=time.sleep)->tuple[dict,pd.DataFrame]:
    results=[]
    errors=0
    for start in range(0,len(sample),MAX_BATCH):
        g=sample.iloc[start:start+MAX_BATCH].copy()
        ids=g.observation_id.astype(int).tolist()
        try:
            payload=client(ids)
            received={int(o["id"]):o for o in payload.get("results",[]) if isinstance(o,dict) and o.get("id")}
            if len(received)>MAX_BATCH:
                raise ValueError("API response exceeded frozen ID request")
            part=[summarize_one(z,received.get(int(z["observation_id"]))) for z in g.to_dict("records")]
        except Exception as e:
            errors+=1
            part=[{**{
                "inat_taxon_id":int(z["inat_taxon_id"]),"observation_id":int(z["observation_id"]),
                "photo_id":int(z["photo_id"]),"source_cell_id":int(z["cell_id"]),
                "source_sample_role":z["origin"],"latitude":None,"longitude":None,
                "positional_accuracy_m":None},"status":"API_ERROR_UNRESOLVED",
                "error_class":type(e).__name__} for z in g.to_dict("records")]
        results.extend(part)
        if start+MAX_BATCH<len(sample):
            pause(RATE_INTERVAL_SEC)
    out=pd.DataFrame(results)
    if len(out)!=len(sample) or out.observation_id.duplicated().any():
        raise RuntimeError("Pilot changed original observation denominator")
    status=out.status.value_counts().to_dict()
    return {
        "schema":"fcp_global42111_photo_ID_public_coordinates_pilot_v1",
        "status":"PUBLIC_METADATA_LOCATION_RECOVERY_NOT_ENVIRONMENT_EFFECT",
        "n_original_measured_source_species":42111,
        "n_original_measured_taxon_cell":85337,
        "n_photo_metadata_targets":len(sample),
        "n_requests":int(np.ceil(len(sample)/MAX_BATCH)),
        "n_failed_request_batches":errors,
        "n_exact_source_photo_coordinates":int((out.status=="EXACT_SOURCE_PHOTO_PUBLIC_COORDINATE_AVAILABLE").sum()),
        "by_photo_metadata_status":{str(k):int(v) for k,v in status.items()},
        "n_pilot_original_cells":int(sample.loc[sample.cell_id>=0,"cell_id"].nunique()),
        "original_162_equal_area_cell_ids_not_centroind_imputed":True,
        "no_photo_pixels_opened":True,
        "no_new_species_or_photo_identity_selection":True,
        "no_climate_or_soil_effect_estimated":True,
        "outcome_labels_not_used_for_selection":True,
        "confirmatory_decisions_changed":False,
    },out


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--breadth",required=True,type=Path)
    p.add_argument("--taxon-cell",required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    a=p.parse_args()
    breadth=pd.read_csv(a.breadth,low_memory=False)
    cells=pd.read_csv(a.taxon_cell,low_memory=False)
    sample=build_source_sample(breadth,cells)
    report,rows=pilot(sample)
    a.outdir.mkdir(parents=True,exist_ok=True)
    sample.to_csv(a.outdir/"original_photo_ID_pilot_prequery_manifest.csv",index=False)
    rows.to_csv(a.outdir/"photo_coordinate_pilot_public_metadata.csv",index=False)
    (a.outdir/"result.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print(json.dumps(report,indent=2,sort_keys=True),flush=True)


if __name__=="__main__":
    main()

#!/usr/bin/env python3
"""Bounded, metadata-only iNaturalist audit of 22 FCP Validation gap-one targets.

Never fetches photo URLs or pixels. NEVER upgrades original 10-km HOLD.
All source IDs/date/observers and 42,111-photo breadth exclusions are fixed.
Results from this API pilot may be changed by current taxonomy, licensing,
obscured locations and pagination; they are not biological outcomes.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
import sys

# Standalone CLI must locate sibling analysis modules without PYTHONPATH setup.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "analysis"))

import numpy as np
import pandas as pd

from audit_fcp_siteyear_same_month_opportunity_20261008 import (
    SHA256, gc_matrix_km, load_photo_opportunity
)

QUEUE_SHA256 = "a69df5ad49cc8cabde2217fdd7e792b57e9ce972abebefd37b044b3fd5e60bb0"
BREADTH_SHA256 = "38aa42123b4e9b05753020ff1a3b050f4d14dbd3de3194557ead90b75c0cc605"
API_BASE = "https://api.inaturalist.org/v1/observations"
ALLOWED_LICENSES = {"cc0","cc-by","cc-by-sa","cc-by-nc","cc-by-nc-sa"}
COHORT = "validation"
PRIORITY = 1
MAX_REQUESTS = 22
PER_PAGE = 200
MAX_SITE_RADIUS_KM = 5.0
MAX_POSITIONAL_ACCURACY_M = 5000.0
MIN_INTERVAL_SECONDS = 1.15
TIMEOUT_SECONDS = 25
USER_AGENT = "fcp-flower-colour-siteyear-meta-pilot/1.0 (metadata-only)"

def sha(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as stream:
        for b in iter(lambda: stream.read(1<<20),b""):
            h.update(b)
    return h.hexdigest()

def source_exclusions(paths: dict[str,Path], breadth:Path) -> tuple[set[str],set[str]]:
    observations=set()
    photos=set()
    for c,f in paths.items():
        if sha(f)!=SHA256[c]:
            raise RuntimeError(f"Source photo SHA mismatch: {c}")
        d=pd.read_csv(f, low_memory=False, usecols=["observation_id","photo_id"])
        observations.update(d.observation_id.dropna().astype("int64").astype(str))
        photos.update(d.photo_id.dropna().astype("int64").astype(str))
    if sha(breadth)!=BREADTH_SHA256:
        raise RuntimeError("Frozen 42111 breadth photo-ID exclusion SHA mismatch")
    b=pd.read_csv(breadth,compression="gzip",usecols=["observation_id","photo_id"],low_memory=False)
    observations.update(b.observation_id.dropna().astype("int64").astype(str))
    photos.update(b.photo_id.dropna().astype("int64").astype(str))
    return observations, photos

def anchor_context(d:pd.DataFrame, row:pd.Series)->tuple[float,float,set[str]]:
    focal=d.loc[d.inat_taxon_id.astype(str)==str(row.inat_taxon_id)].copy()
    anchor=focal.loc[focal.photo_id.astype(str)==str(row.source_anchor_photo_id)]
    if len(anchor)!=1:
        raise RuntimeError(f"Anchor photo missing or duplicated for taxon {row.inat_taxon_id}")
    lat=float(anchor.latitude.iloc[0])
    lon=float(anchor.longitude.iloc[0])
    month=int(row.calendar_month); year=int(row.target_year)
    match=focal.loc[(focal.month==month)&(focal.year==year)&(focal.observer!="")].copy()
    if len(match):
        dist=gc_matrix_km(np.r_[lat,match.latitude.to_numpy(float)],np.r_[lon,match.longitude.to_numpy(float)])[0,1:]
        match=match.loc[dist<=MAX_SITE_RADIUS_KM+1e-8]
    existing=set(match.observer.astype(str))
    if len(existing)!=int(row.existing_distinct_observers_in_target_year):
        raise RuntimeError(f"Original observer counts drifted for {row.inat_taxon_id}/{year}/{month}")
    if int(row.required_additional_distinct_observer_photos)!=1 or len(existing)!=1:
        raise RuntimeError("Pilot can inspect only gap-one source-year cells with exactly one existing observer")
    return lat,lon,existing

def public_coordinates(item:dict)->tuple[float,float]|None:
    if item.get("coordinates_obscured") is True or str(item.get("geoprivacy") or "").lower() in {"private","obscured"}:
        return None
    coords=item.get("geojson")
    if isinstance(coords,dict) and isinstance(coords.get("coordinates"),list) and len(coords["coordinates"])==2:
        lon,lat=coords["coordinates"]
    else:
        s=item.get("location")
        if not s or not isinstance(s,str) or "," not in s:
            return None
        lat,lon=s.split(",",1)
    try:
        lat=float(lat);lon=float(lon)
    except (TypeError,ValueError):
        return None
    return (lat,lon) if np.isfinite(lat) and np.isfinite(lon) and -90<=lat<=90 and -180<=lon<=180 else None

def qualify(item:dict,*,taxon_id:int,year:int,month:int,
            anchor_lat:float,anchor_lon:float,existing_observers:set[str],
            banned_obs:set[str],banned_photos:set[str])->list[dict]:
    if not isinstance(item,dict) or not isinstance(item.get("taxon"),dict) or str(item["taxon"].get("id"))!=str(taxon_id):
        return []
    if item.get("captive") is True:
        return []
    date=str(item.get("observed_on") or "")
    if not date.startswith(f"{year:04}-{month:02}-"):
        return []
    user=item.get("user") or {}
    observer=str(user.get("id") or "")
    if not observer or observer in existing_observers:
        return []
    obsid=str(item.get("id") or "")
    if not obsid or obsid in banned_obs:
        return []
    position=public_coordinates(item)
    if position is None:
        return []
    try:
        accuracy=float(item.get("positional_accuracy"))
    except (TypeError,ValueError):
        return []
    if not np.isfinite(accuracy) or accuracy<0 or accuracy>MAX_POSITIONAL_ACCURACY_M:
        return []
    lat,lon=position
    d=float(gc_matrix_km(np.array([anchor_lat,lat]),np.array([anchor_lon,lon]))[0,1])
    if d>MAX_SITE_RADIUS_KM+1e-8:
        return []
    accepted=[]
    for p in item.get("photos") or []:
        pid=str(p.get("id") or "")
        license_code=str(p.get("license_code") or p.get("license") or "").lower()
        if pid and pid not in banned_photos and license_code in ALLOWED_LICENSES:
            accepted.append({"observation_id":obsid,"photo_id":pid,
                             "observer_id":observer,"public_distance_km":d,
                             "position_accuracy_m":accuracy,"photo_license":license_code})
    return accepted

def fetch_first_page(params:dict)->dict:
    url=API_BASE+"?"+urllib.parse.urlencode(params)
    req=urllib.request.Request(url,headers={"User-Agent":USER_AGENT,"Accept":"application/json"})
    with urllib.request.urlopen(req,timeout=TIMEOUT_SECONDS) as response:
        if getattr(response,"status",200)!=200:
            raise RuntimeError("Non-200 API response")
        data=json.load(response)
    if not isinstance(data,dict) or not isinstance(data.get("results"),list) or not isinstance(data.get("total_results"),int):
        raise RuntimeError("Unexpected iNaturalist JSON envelope")
    return data

def run(queue:Path,sources:dict[str,Path],breadth:Path,output:Path,fetcher=fetch_first_page,sleep=time.sleep)->dict:
    if sha(queue)!=QUEUE_SHA256:
        raise RuntimeError("Frozen label-blind photo gap queue SHA256 drift")
    t=pd.read_csv(queue,low_memory=False)
    t=t.loc[(t.cohort==COHORT)&(t.priority_tier==PRIORITY)].copy()
    if len(t)!=MAX_REQUESTS or t.inat_taxon_id.nunique()!=MAX_REQUESTS:
        raise RuntimeError("The frozen 22 Validation gap-one targets changed")
    if not t.live_metadata_status.eq("NOT_CHECKED").all():
        raise RuntimeError("Queue already claims live metadata outcome")
    banned_obs,banned_photos=source_exclusions(sources,breadth)
    d,_=load_photo_opportunity(sources[COHORT],COHORT)
    output.mkdir(parents=True,exist_ok=True)
    results=[]
    for i,row in enumerate(t.sort_values(["inat_taxon_id","target_year"],kind="stable").itertuples(index=False)):
        lat,lon,existing=anchor_context(d,row)
        params={"taxon_id":int(row.inat_taxon_id),"year":int(row.target_year),
                "month":int(row.calendar_month),"lat":round(lat,6),"lng":round(lon,6),
                "radius":MAX_SITE_RADIUS_KM,"photos":"true","geo":"true",
                "per_page":PER_PAGE,"page":1}
        result={"cohort":COHORT,"inat_taxon_id":str(row.inat_taxon_id),
                "target_year":int(row.target_year),"month":int(row.calendar_month),
                "original_anchor_photo_id":str(row.source_anchor_photo_id),
                "source_existing_observer_count":len(existing),
                "source_complete_10km_hold_unchanged":True,
                "query_status":"REQUEST_NOT_STARTED","candidate_metadata_count":0,
                "API_reported_total":None,"API_page_limit":PER_PAGE,
                "eligible_photo_observation_ids":[]}
        try:
            response=fetcher(params)
            total=int(response["total_results"])
            result["API_reported_total"]=total
            found=[]
            for item in response["results"]:
                found.extend(qualify(item,taxon_id=int(row.inat_taxon_id),
                    year=int(row.target_year),month=int(row.calendar_month),
                    anchor_lat=lat,anchor_lon=lon,existing_observers=existing,
                    banned_obs=banned_obs,banned_photos=banned_photos))
            # Distinct-observer and observation identity; deterministic ID ordering.
            unique={}
            for item in sorted(found,key=lambda v:(v["observer_id"],v["observation_id"],v["photo_id"])):
                unique.setdefault(item["observer_id"],item)
            result["eligible_photo_observation_ids"]=list(unique.values())
            result["candidate_metadata_count"]=len(unique)
            result["query_status"]=(
                "ELIGIBLE_PUBLIC_METADATA_CANDIDATE_EXISTS"
                if unique else
                "TRUNCATED_FIRST_PAGE_UNKNOWN" if total>PER_PAGE else
                "NO_QUALIFYING_METADATA_ON_COMPLETE_FIRST_PAGE"
            )
            if total>PER_PAGE:
                result["pagination_incomplete"]=True
        except (urllib.error.HTTPError,urllib.error.URLError,TimeoutError,RuntimeError,ValueError) as e:
            result["query_status"]="API_ERROR_UNRESOLVED"
            result["error_class"]=type(e).__name__
        results.append(result)
        if i+1<len(t):
            sleep(MIN_INTERVAL_SECONDS)
    n=sum(x["query_status"]=="ELIGIBLE_PUBLIC_METADATA_CANDIDATE_EXISTS" for x in results)
    errors=sum(x["query_status"]=="API_ERROR_UNRESOLVED" for x in results)
    pending=sum(x["query_status"]=="TRUNCATED_FIRST_PAGE_UNKNOWN" for x in results)
    report={"schema":"fcp_validation_gap1_live_observation_metadata_pilot_v1",
            "date_jst":"2026-10-08","status":"POSTHOC_LIVE_METADATA_ONLY_NOT_IMAGE_CLASSIFICATION",
            "original_10km_gate":"HOLD_UNCHANGED",
            "n_query_target_species":MAX_REQUESTS,
            "n_api_requests_attempted":len(results),"n_metadata_candidate_species":n,
            "n_api_error_species":errors,"n_truncated_unknown_species":pending,
            "n_first_page_no_candidate_species":MAX_REQUESTS-n-errors-pending,
            "eligible_metadata_is_NOT_classifiable_flower_photo":True,
            "queue_sha256":QUEUE_SHA256,"breadth_exclusion_sha256":BREADTH_SHA256,
            "source_cohort_sha256":SHA256,"max_requests":MAX_REQUESTS,
            "strict_query":{"same_species_taxon_ID":True,"same_year_month":True,
                            "max_anchor_distance_km":MAX_SITE_RADIUS_KM,
                            "max_position_accuracy_m":MAX_POSITIONAL_ACCURACY_M,
                            "observer_must_differ_from_original_year":True,
                            "exclude_all_original_three_cohorts_and_42111_breadth_photo_ids":True,
                            "live_photo_pixels_opened":False},
            "hard_nonclaims":[
                "live metadata is not a usable flowering image nor a genetically segregating morph",
                "the one-page capped search cannot rule out unused photos on later pages",
                "photo licenses and current taxon assignment may have changed since historical freezing",
                "no photographs or colour labels are opened, and no 2000+730 prospective taxa are queried",
                "photographic site/year does not identify the same plant nor adaptation",
                "original precommitted 10km sample coverage HOLD remains in force",
            ]}
    (output/"result.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    pd.DataFrame(results).drop(columns=["eligible_photo_observation_ids"]).to_csv(output/"per_species_technical_status.csv",index=False)
    # Keep candidate IDs as JSON for subsequent independent quality checks.
    (output/"public_metadata_candidates.json").write_text(json.dumps(results,indent=2,sort_keys=True)+"\n")
    print(json.dumps(report,indent=2,sort_keys=True))
    return report

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--queue",type=Path,required=True)
    for c in ("discovery","validation","third"):
        p.add_argument("--"+c,type=Path,required=True)
    p.add_argument("--breadth",type=Path,required=True)
    p.add_argument("--outdir",type=Path,required=True)
    a=p.parse_args()
    run(a.queue,{c:getattr(a,c) for c in ("discovery","validation","third")},a.breadth,a.outdir)

if __name__=="__main__":
    main()

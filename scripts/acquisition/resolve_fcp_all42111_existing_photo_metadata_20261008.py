#!/usr/bin/env python3
"""Resolve URL/licence for exactly one ORIGINAL photo per all 42,111 FCP taxa.

No photo pixels downloaded, no new source photos substituted. iNaturalist API
allows comma-separated observation IDs in one request; use batches of 100 and
<=1 request/sec per API recommended practices. All taxon rows, including
missing, revoked-licence and transport failures, remain in the denominator.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import pandas as pd

FROZEN_IDENTITIES_SHA256="aa179ac271be2731cf1fdd00f305999b765cd7680b358554b75fd9d9372091ba"
EXPECTED_SPECIES=42111
IDS_PER_REQUEST=100
REQUEST_START_MIN_INTERVAL_SECONDS=1.2
TIMEOUT_SECONDS=45
RETRIES=1
ALLOWED_LICENSES=frozenset(["cc0","cc-by","cc-by-sa","cc-by-nc","cc-by-nc-sa"])
UA="zuizui0223-fcp-whole42111-original-photo-metadata/1.0 (github.com/zuizui0223/fcp)"
URL_PREFIX="https://api.inaturalist.org/v1/observations"
TRANSPORT_RECOVERY_OF_RUN=37762416621
# The /observations/{comma-separated-ids} route may reject >30 IDs (422).
# The documented ID search query accepts longer lists with per_page=200.
ALLOWED_STATUSES=frozenset([
    "VALID_SOURCE_PHOTO_URL_AND_LICENSE",
    "OBSERVATION_ID_NOT_RETURNED",
    "ORIGINAL_PHOTO_ID_MISSING",
    "PHOTO_LICENSE_NOT_ALLOWED",
    "PHOTO_URL_MISSING",
    "SOURCE_API_QUERY_FAILURE",
    "TAXON_IDENTITY_MISMATCH"
])


def sha(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for buf in iter(lambda:f.read(1<<20),b""):h.update(buf)
    return h.hexdigest()


def original_photo_manifest(path:Path)->pd.DataFrame:
    if sha(path)!=FROZEN_IDENTITIES_SHA256:
        raise RuntimeError("All-42111 historically selected original photo identities have drifted")
    d=pd.read_csv(path,low_memory=False)
    if len(d)!=EXPECTED_SPECIES:
        raise RuntimeError("Not exactly 42111 original source photo IDs")
    for field in ("inat_taxon_id","observation_id","photo_id","species"):
        if field not in d.columns:
            raise RuntimeError(f"Original photo identity missing: {field}")
    if any(d[x].duplicated().any() for x in ("inat_taxon_id","observation_id","photo_id")):
        raise RuntimeError("Source one-photo-per-species IDs were reused")
    for field in ("inat_taxon_id","observation_id","photo_id"):
        d[field]=pd.to_numeric(d[field],errors="raise").astype("int64")
    return d.sort_values("inat_taxon_id",kind="stable").reset_index(drop=True)


def large_url(url:str)->str:
    """Change only the photo-size token and retain its original extension."""
    text=str(url).strip()
    for size in ("square","small","medium","thumb"):
        token=f"/{size}."
        if token in text:
            return text.replace(token,"/large.")
    return text


def decode_batch(source:pd.DataFrame,payload:dict)->list[dict]:
    if not isinstance(payload,dict):
        raise ValueError("Observation API returned a non-JSON-object")
    results=payload.get("results")
    if not isinstance(results,list):
        raise ValueError("Observation API returned no results array")
    obsbyid={}
    for row in results:
        if not isinstance(row,dict):
            continue
        try: k=int(row.get("id"))
        except (TypeError,ValueError):continue
        if k in obsbyid:
            raise ValueError(f"Duplicate observation ID in one API response: {k}")
        obsbyid[k]=row
    out=[]
    for row in source.itertuples(index=False):
        obs=obsbyid.get(int(row.observation_id))
        item={
            "inat_taxon_id":int(row.inat_taxon_id),
            "species":str(row.species),
            "observation_id":int(row.observation_id),
            "photo_id":int(row.photo_id),
            "historical_source_photo_capacity":int(row.after_observer_cap),
            "photo_url_large":"",
            "source_photo_license":"",
            "returned_taxon_id":None,
            "source_taxon_identity_match":False,
            "photo_url_status":"",
        }
        if obs is None:
            item["photo_url_status"]="OBSERVATION_ID_NOT_RETURNED"
            out.append(item);continue
        taxon=obs.get("taxon")
        if isinstance(taxon,dict):
            try:
                item["returned_taxon_id"]=int(taxon.get("id"))
            except (TypeError,ValueError):
                pass
        item["source_taxon_identity_match"]=item["returned_taxon_id"]==int(row.inat_taxon_id)
        # A different current source taxon must NOT be counted as
        # verified colour evidence for the immutable historical species.
        if not item["source_taxon_identity_match"]:
            item["photo_url_status"]="TAXON_IDENTITY_MISMATCH"
            out.append(item);continue
        photos=obs.get("photos") or []
        matching=[]
        for photo in photos:
            if not isinstance(photo,dict):continue
            try:
                if int(photo.get("id"))==int(row.photo_id):
                    matching.append(photo)
            except (TypeError,ValueError):
                continue
        if len(matching)!=1:
            item["photo_url_status"]="ORIGINAL_PHOTO_ID_MISSING"
        else:
            photo=matching[0]
            item["source_photo_license"]=str(photo.get("license_code") or "").strip().casefold()
            url=photo.get("url") or ""
            if item["source_photo_license"] not in ALLOWED_LICENSES:
                item["photo_url_status"]="PHOTO_LICENSE_NOT_ALLOWED"
            elif not isinstance(url,str) or not url.startswith("https://"):
                item["photo_url_status"]="PHOTO_URL_MISSING"
            else:
                item["photo_url_large"]=large_url(url)
                item["photo_url_status"]="VALID_SOURCE_PHOTO_URL_AND_LICENSE"
        out.append(item)
    return out


def fetch_api_batch(observation_ids:list[int], *,rate_sleep=time.sleep)->dict:
    if not 1<=len(observation_ids)<=IDS_PER_REQUEST:
        raise ValueError("Unexpected iNaturalist batch size")
    if len(set(observation_ids))!=len(observation_ids):
        raise ValueError("No duplicate observation IDs may enter API request")
    url=URL_PREFIX + "?" + urlencode({"per_page":200, "id":",".join(str(int(x)) for x in observation_ids)})
    err=None
    for attempt in range(RETRIES+1):
        if attempt:
            rate_sleep(5*attempt)
        req=Request(url,headers={"Accept":"application/json","User-Agent":UA})
        try:
            with urlopen(req,timeout=TIMEOUT_SECONDS) as reply:
                return json.loads(reply.read().decode("utf-8"))
        except (HTTPError,URLError,TimeoutError,ValueError) as exc:
            err=exc
            if isinstance(exc,HTTPError) and exc.code in (401,403,404,422):
                break
            if isinstance(exc,HTTPError) and exc.code==429:
                rate_sleep(60)
    raise RuntimeError("Batch source fetch failed after frozen bounded retries: "+type(err).__name__)


def resolve_all(identities:pd.DataFrame,client=fetch_api_batch,*,pause=time.sleep)->tuple[pd.DataFrame,dict]:
    if len(identities)!=EXPECTED_SPECIES:
        raise RuntimeError("Source photo selection must cover all 42111")
    parts=[];failures=0;calls=0
    last_request_start=None
    for offset in range(0,len(identities),IDS_PER_REQUEST):
        chunk=identities.iloc[offset:offset+IDS_PER_REQUEST]
        if last_request_start is not None:
            elapsed=time.monotonic()-last_request_start
            if elapsed<REQUEST_START_MIN_INTERVAL_SECONDS:
                pause(REQUEST_START_MIN_INTERVAL_SECONDS-elapsed)
        last_request_start=time.monotonic()
        calls+=1
        try:
            payload=client(chunk.observation_id.astype(int).tolist())
            rows=decode_batch(chunk,payload)
        except Exception as ex:
            failures+=1
            rows=[{
                "inat_taxon_id":int(x.inat_taxon_id),
                "species":str(x.species),
                "observation_id":int(x.observation_id),
                "photo_id":int(x.photo_id),
                "historical_source_photo_capacity":int(x.after_observer_cap),
                "photo_url_large":"","source_photo_license":"",
                "returned_taxon_id":None,
                "source_taxon_identity_match":False,
                "photo_url_status":"SOURCE_API_QUERY_FAILURE",
            } for x in chunk.itertuples(index=False)]
        parts.extend(rows)
    out=pd.DataFrame(parts)
    if len(out)!=EXPECTED_SPECIES or not out.inat_taxon_id.is_unique:
        raise RuntimeError("At least one originally archived species was omitted during URL resolution")
    if out.observation_id.duplicated().any() or out.photo_id.duplicated().any():
        raise RuntimeError("Reused a source photo across two plant taxa")
    if not out.photo_url_status.isin(ALLOWED_STATUSES).all():
        raise RuntimeError("Unhandled source URL / licensing state")
    counts=out.photo_url_status.value_counts().to_dict()
    return out,{
        "n_original_source_species":len(out),
        "n_observation_batches":calls,
        "expected_n_batches":(EXPECTED_SPECIES+IDS_PER_REQUEST-1)//IDS_PER_REQUEST,
        "n_failed_batches":failures,
        "species_by_current_original_photo_URL_status":{k:int(counts.get(k,0)) for k in sorted(ALLOWED_STATUSES)},
        "n_url_and_license_verified":int(out.photo_url_status.eq("VALID_SOURCE_PHOTO_URL_AND_LICENSE").sum()),
        "n_url_and_license_verified_original_capacity_zero":int(
            (out.photo_url_status.eq("VALID_SOURCE_PHOTO_URL_AND_LICENSE") &
             out.historical_source_photo_capacity.eq(0)).sum()),
        "n_name_id_mismatches_in_source_taxonomy":int((~out.source_taxon_identity_match).sum())
    }


def main()->None:
    p=argparse.ArgumentParser()
    p.add_argument("--frozen-photo-identities",type=Path,required=True)
    p.add_argument("--outdir",type=Path,required=True)
    args=p.parse_args()
    if os.environ.get("GITHUB_ACTIONS")=="true" and os.environ.get("GITHUB_RUN_ATTEMPT","1")!="1":
        raise RuntimeError("No rerun of URL/permission snapshot after viewing response outcomes")
    manifest=original_photo_manifest(args.frozen_photo_identities)
    rows,stats=resolve_all(manifest)
    args.outdir.mkdir(parents=True,exist_ok=True)
    target=args.outdir/"whole_42111_original_photo_URL_license_audit.csv.gz"
    rows.to_csv(target,index=False,lineterminator="\n",
                compression={"method":"gzip","compresslevel":9,"mtime":0})
    receipt={
        "schema":"fcp_all42111_first_historical_photo_source_URL_resolution_v1",
        "date_jst":"2026-10-08",
        "status":"WHOLE_SOURCE_PHOTO_URL_LICENCE_RESOLVED_WITH_MISSINGNESS",
        "original_single_photo_species":EXPECTED_SPECIES,
        "original_archive_single_photo_identity_sha256":FROZEN_IDENTITIES_SHA256,
        "api_endpoint":"GET /v1/observations?per_page=200&id=id1,id2,... in batches of 100",
        "technical_only_transport_recovery_of_failed_run":TRANSPORT_RECOVERY_OF_RUN,
        "technical_failure_of_previous_endpoint":"421 of 422 source batches failed on the comma-separated IDs path; route replaced without changing selected photo identities",
        "api_requests_rate_policy":"<=1 batch per 1.2 seconds; bounded 1 retry; follows bulk observation-ID recommendations",
        "original_photo_identity_unchanged":True,
        "downloaded_image_pixels":False,
        "classified_flower_colour":False,
        "new_external_photo_substitution":False,
        "source_photo_URL_workload":stats,
        "URL_licence_ledger_SHA256":sha(target),
        "external_photo_licence_consent_required_per_original_asset":True,
        "historical_photo_not_independent_species_confirmation":True,
        "run_completion_does_not_authorize_image_pixel_download":True,
        "confirmatory_decisions_changed":False
    }
    (args.outdir/"result.json").write_text(json.dumps(receipt,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print(json.dumps(receipt,indent=2,ensure_ascii=False,allow_nan=False))


if __name__=="__main__":
    main()

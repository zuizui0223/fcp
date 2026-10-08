#!/usr/bin/env python3
"""Bounded posthoc iNaturalist metadata feasibility for 60 Validation gap-two taxa.

Reuses original 10-km conservative anchored photo opportunity and the frozen
outcome-blind queue. A species requires the full gap of two DISTINCT observer
photo slots, in the fixed photographed year(s) and same calendar month.
This makes no photo pixel request and does not establish a colour phenotype.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import time

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "analysis"))
from audit_fcp_validation_gap1_live_metadata_pilot_20261008 import (
    SHA256, QUEUE_SHA256, BREADTH_SHA256, PER_PAGE,
    COHORT, MAX_SITE_RADIUS_KM, MAX_POSITIONAL_ACCURACY_M,
    MIN_INTERVAL_SECONDS, sha, source_exclusions, anchor_context,
    fetch_first_page, qualify,
)
from audit_fcp_siteyear_same_month_opportunity_20261008 import load_photo_opportunity

TARGET_GAP = 2
TARGET_SPECIES = 60
TOTAL_MISSING_SLOTS = 120
MAX_REQUESTS = 120


def decide_species(requirements: list[dict]) -> str:
    """Require complete metadata eligibility for EVERY deficient source year.

    Possible later-page candidates and failed API requests are UNKNOWN, never
    counted as absence or success. Existence at one year cannot rescue the
    opposite year's missing observer quota.
    """
    if not requirements:
        raise ValueError("Species missing required year/month slots")
    if sum(int(x["required"]) for x in requirements) != TARGET_GAP:
        raise ValueError("Not a frozen two-slot species")
    if len(set(int(x["target_year"]) for x in requirements)) != len(requirements):
        raise ValueError("Repeated target year")
    if any(x.get("api_error") for x in requirements):
        return "API_ERROR_UNRESOLVED"
    if all(int(x["found"]) >= int(x["required"]) for x in requirements):
        return "POSSIBLE_COMPLETE_METADATA_ONLY"
    if any(x.get("pagination_incomplete") for x in requirements):
        return "PAGINATION_INCOMPLETE_UNKNOWN"
    return "NO_COMPLETE_METADATA_CANDIDATE_FIRST_PAGES"


def check_queue(queue: pd.DataFrame) -> pd.DataFrame:
    expected={
        "cohort","priority_tier","inat_taxon_id","calendar_month",
        "target_year","source_anchor_photo_id","existing_distinct_observers_in_target_year",
        "required_additional_distinct_observer_photos","live_metadata_status"
    }
    if not expected.issubset(queue):
        raise ValueError("Missing photo-opportunity query columns")
    q=queue.loc[
        (queue.cohort==COHORT) & (queue.priority_tier==TARGET_GAP)
    ].copy()
    if q.inat_taxon_id.nunique()!=TARGET_SPECIES or not TARGET_SPECIES<=len(q)<=MAX_REQUESTS:
        raise ValueError("Frozen Validation gap-two taxa or request bound changed")
    if int(q.required_additional_distinct_observer_photos.sum())!=TOTAL_MISSING_SLOTS:
        raise ValueError("Source two-photo requirement drift")
    if q.duplicated(["inat_taxon_id","target_year"]).any():
        raise ValueError("Repeat query for same species-year")
    if not q.live_metadata_status.eq("NOT_CHECKED").all():
        raise ValueError("Future outcomes leaked into frozen queue")
    for taxon, group in q.groupby("inat_taxon_id"):
        if int(group.required_additional_distinct_observer_photos.sum())!=TARGET_GAP:
            raise ValueError(f"{taxon}: requirement is not exactly two")
        if group.calendar_month.nunique()!=1 or group.source_anchor_photo_id.nunique()!=1:
            raise ValueError("Changing source anchored month or place across years")
        if (group.required_additional_distinct_observer_photos.astype(int)<=0).any():
            raise ValueError("Only genuinely missing year cells may be queried")
    return q.sort_values(["inat_taxon_id","target_year"],kind="stable")


def run(queue_file:Path,sources:dict[str,Path],breadth:Path,outdir:Path,
        *,fetcher=fetch_first_page,sleep=time.sleep) -> dict:
    if sha(queue_file)!=QUEUE_SHA256:
        raise RuntimeError("Validated earlier observer-photo queue SHA mismatch")
    q=check_queue(pd.read_csv(queue_file,low_memory=False))
    old_obs,old_photos=source_exclusions(sources,breadth)
    source,_=load_photo_opportunity(sources[COHORT],COHORT)
    # Persist every metadata response before aggregation, so late analysis errors
    # never discard API work or turn queried photo IDs into unknown negatives.
    outdir.mkdir(parents=True,exist_ok=True)
    progress=outdir/"technical_query_progress.jsonl"
    if progress.exists():
        raise RuntimeError("Pre-existing partial metadata output: do not silently repeat requests")
    technical=[]
    for i,row in enumerate(q.itertuples(index=False)):
        lat,lon,existing=anchor_context_gap2(source,row)
        requirement=int(row.required_additional_distinct_observer_photos)
        params={
            "taxon_id":int(row.inat_taxon_id),"year":int(row.target_year),
            "month":int(row.calendar_month),"lat":round(lat,6),
            "lng":round(lon,6),"radius":MAX_SITE_RADIUS_KM,
            "photos":"true","geo":"true","per_page":PER_PAGE,"page":1
        }
        z={
            "cohort":COHORT, "inat_taxon_id":str(row.inat_taxon_id),
            "source_anchor_photo_id":str(row.source_anchor_photo_id),
            "target_year":int(row.target_year),"month":int(row.calendar_month),
            "existing_observers_in_year":len(existing),
            "required":requirement, "found":0,
            "api_error":False,"pagination_incomplete":False,
            "api_total_results":None,"query_status":"REQUEST_NOT_STARTED",
            "metadata_candidates":[],
        }
        try:
            response=fetcher(params)
            z["api_total_results"]=int(response["total_results"])
            candidates=[]
            for item in response["results"]:
                candidates.extend(qualify(
                    item,taxon_id=int(row.inat_taxon_id),
                    year=int(row.target_year),month=int(row.calendar_month),
                    anchor_lat=lat,anchor_lon=lon,
                    existing_observers=existing,banned_obs=old_obs,banned_photos=old_photos))
            per_observer={}
            for c in sorted(candidates,key=lambda x:(x["observer_id"],x["observation_id"],x["photo_id"])):
                per_observer.setdefault(c["observer_id"],c)
            z["metadata_candidates"]=list(per_observer.values())
            z["found"]=len(per_observer)
            z["pagination_incomplete"]=z["api_total_results"]>PER_PAGE
            z["query_status"]=(
                "SUFFICIENT_YEAR_METADATA_CANDIDATES" if z["found"]>=requirement
                else "INCOMPLETE_API_FIRST_PAGE" if z["pagination_incomplete"]
                else "INSUFFICIENT_ON_COMPLETE_FIRST_PAGE"
            )
        except Exception as e:
            # Do not treat transient errors or schema mismatch as no candidates.
            z["api_error"]=True
            z["query_status"]="API_ERROR_UNRESOLVED"
            z["error_class"]=type(e).__name__
        technical.append(z)
        with progress.open("a",encoding="utf-8") as fp:
            fp.write(json.dumps(z,sort_keys=True)+"\\n")
        if i+1<len(q):
            sleep(MIN_INTERVAL_SECONDS)
    if len(technical)!=len(q):
        raise RuntimeError("Incomplete metadata collection cannot be aggregated")
    (outdir/"technical_all_queries_completed.json").write_text(
        json.dumps(technical,indent=2,sort_keys=True)+"\\n",encoding="utf-8")
    print("FCP_GAP2_ALL_METADATA_QUERY_RESPONSES_DURABLY_WRITTEN",flush=True)
    by_species={}
    for tid in sorted(q.inat_taxon_id.astype(str).unique(),key=int):
        slots=[x for x in technical if x["inat_taxon_id"]==tid]
        decision=decide_species(slots)
        by_species[tid]={
            "inat_taxon_id":tid,"source_year_cells_checked":len(slots),
            "required_total":sum(z["required"] for z in slots),
            "found_total":sum(min(z["required"],z["found"]) for z in slots),
            "status":decision,
        }
    successes=sum(x["status"]=="POSSIBLE_COMPLETE_METADATA_ONLY" for x in by_species.values())
    errors=sum(x["status"]=="API_ERROR_UNRESOLVED" for x in by_species.values())
    partial=sum(x["status"]=="PAGINATION_INCOMPLETE_UNKNOWN" for x in by_species.values())
    report={
        "schema":"fcp_validation_gap2_live_observation_metadata_pilot_v1",
        "date_jst":"2026-10-09",
        "status":"POSTHOC_METADATA_ONLY_NOT_BIOLOGICAL_OUTCOME",
        "source_validation_species":363,"previously_qualified_anchor10_species":10,
        "frozen_validation_gap1_metadata_candidate_species":15,
        "original_10km_gate":"HOLD_UNCHANGED",
        "n_species_target":TARGET_SPECIES,
        "n_api_requests_attempted":len(technical),
        "max_api_requests":MAX_REQUESTS,
        "required_missing_observer_photo_slots":TOTAL_MISSING_SLOTS,
        "n_species_both_year_requirements_metadata_possible":successes,
        "n_species_api_error_unresolved":errors,
        "n_species_pagination_unknown":partial,
        "n_species_first_pages_insufficient":TARGET_SPECIES-successes-errors-partial,
        "max_metadata_only_qualified_validation_species_if_all_images_valid":10+15+successes,
        "frozen_validation_species_threshold":30,
        "could_cross_30_using_only_found_metadata":10+15+successes>=30,
        "queue_sha256":QUEUE_SHA256,"breadth_exclusion_sha256":BREADTH_SHA256,
        "source_sha256":SHA256,
        "no_photo_pixels_fetched":True,
        "confirmatory_decisions_changed":False,
        "nonclaims":[
            "Current metadata candidates are not independently validated flowering images or genetic colour morphs",
            "The original 10km source-only sample qualification remains HOLD",
            "A source-photo-year/month anchor is not a known individual or genetic population",
            "No additional 2000+730 prospective species or future flower-colour outcomes are examined",
            "No climate anomalies, pigments, photographic exposure or plant fitness are inferred",
            "First-page-only search with observer/photo-ID exclusion may miss later eligible records",
        ]
    }
    outdir.mkdir(parents=True,exist_ok=True)
    (outdir/"result.json").write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    (outdir/"species_metadata_status.json").write_text(json.dumps(list(by_species.values()),indent=2)+"\n")
    (outdir/"year_metadata_candidates.json").write_text(json.dumps(technical,indent=2)+"\n")
    print(json.dumps(report,indent=2,sort_keys=True))
    return report


def anchor_context_gap2(source:pd.DataFrame,row)->tuple[float,float,set[str]]:
    """Same frozen site/year as original; unlike gap-one, 0 or 1 observers may be known."""
    from audit_fcp_siteyear_same_month_opportunity_20261008 import gc_matrix_km
    sub=source.loc[source.inat_taxon_id.astype(str)==str(row.inat_taxon_id)]
    anchor=sub.loc[sub.photo_id.astype(str)==str(row.source_anchor_photo_id)]
    if len(anchor)!=1:
        raise RuntimeError("Frozen anchor not found")
    lat=float(anchor.latitude.iloc[0]);lon=float(anchor.longitude.iloc[0])
    old=sub.loc[
        (sub.year==int(row.target_year))&
        (sub.month==int(row.calendar_month))&
        (sub.observer.astype(str)!="")
    ]
    known=set()
    if len(old):
        d=gc_matrix_km(
            pd.concat([pd.Series([lat]),old.latitude],ignore_index=True).to_numpy(float),
            pd.concat([pd.Series([lon]),old.longitude],ignore_index=True).to_numpy(float)
        )[0,1:]
        known=set(old.loc[d<=MAX_SITE_RADIUS_KM+1e-8,"observer"].astype(str))
    needed=int(row.required_additional_distinct_observer_photos)
    if len(known)!=int(row.existing_distinct_observers_in_target_year):
        raise RuntimeError("Frozen year observer count mismatch")
    if needed!=2-len(known) or not 1<=needed<=2:
        raise RuntimeError("Frozen missing-observer gap mismatch")
    return lat,lon,known


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--queue",required=True,type=Path)
    for c in ("discovery","validation","third"):
        p.add_argument("--"+c,required=True,type=Path)
    p.add_argument("--breadth",required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    a=p.parse_args()
    run(a.queue,{c:getattr(a,c) for c in ("discovery","validation","third")},a.breadth,a.outdir)


if __name__=="__main__":
    main()

#!/usr/bin/env python3
"""Independent real FCP U100 4,730-taxon source ID / species-frame audit.

All five *original* cohorts are checked before the 256+256 unchanged ROI
workers finish. Read frozen source photo metadata only; never requery photos,
run pixel inference or calculate biological diversity/statistical significance.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

OLD = {
    "discovery": {
        "n": 50000, "taxa": 500,
        "sha256": "ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4",
    },
    "validation": {
        "n": 50000, "taxa": 500,
        "sha256": "0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6",
    },
    "prospective": {
        "n": 49900, "taxa": 499,
        "sha256": "57630fc9f281bce94a0c40a70aaf7bce879dde93d6154175adcd021e8f5c1186",
    },
}
COHORTS = ("discovery", "validation", "prospective", "wave4", "new_U100")
FULL_U100 = 4730
AUTHORIZED = 4478
UNMET = FULL_U100 - AUTHORIZED
FULL_PHOTOS = AUTHORIZED * 100
ACTUAL = {**OLD, "wave4": {"n": 49800, "taxa": 498},
          "new_U100": {"n": 248100, "taxa": 2481}}


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for data in iter(lambda:handle.read(1<<20), b""):
            h.update(data)
    return h.hexdigest()


def checked_source(table: pd.DataFrame, label: str, *, n: int, taxa: int) -> dict:
    required = {"photo_id", "inat_taxon_id", "species"}
    if required - set(table):
        raise RuntimeError(f"{label}: original photographic source keys missing {required-set(table)}")
    frame=table[list(required)+(["observation_id"] if "observation_id" in table else [])].copy()
    if frame.isna().any().any():
        raise RuntimeError(f"{label}: original photo/taxon/observation identity missing")
    photoids=pd.to_numeric(frame["photo_id"],errors="raise").astype("int64")
    taxonids=pd.to_numeric(frame["inat_taxon_id"],errors="raise").astype("int64")
    if ((photoids<=0).any() or (taxonids<=0).any()
        or frame["species"].astype(str).str.strip().eq("").any()):
        raise RuntimeError(f"{label}: malformed source numeric/photo or species identity")
    frame["photo_id"]=photoids
    frame["inat_taxon_id"]=taxonids
    frame["species"]=frame["species"].astype(str).str.strip()
    if len(frame)!=n or frame.photo_id.nunique()!=n:
        raise RuntimeError(f"{label}: source photographic row count/ID identity drift")
    if frame.inat_taxon_id.nunique()!=taxa or frame.species.nunique()!=taxa:
        raise RuntimeError(f"{label}: original number of species not exact")
    if not frame.groupby("inat_taxon_id").size().eq(100).all():
        raise RuntimeError(f"{label}: source violates unchanged 100 photos/species")
    if frame[["inat_taxon_id","species"]].drop_duplicates().shape[0]!=taxa:
        raise RuntimeError(f"{label}: taxon name changed within original 100")
    obsids=None
    if "observation_id" in frame:
        obsids=pd.to_numeric(frame["observation_id"],errors="raise").astype("int64")
        if (obsids<=0).any() or obsids.nunique()!=n:
            raise RuntimeError(f"{label}: duplicated original iNaturalist observation ID")
    return {"label":label, "taxa":taxa,"photo_rows":n,
            "photographs":set(photoids.tolist()),
            "taxon_ids":set(taxonids.tolist()),
            "species_names":set(frame.species.tolist()),
            "taxon_species":set(zip(taxonids.tolist(),frame.species.tolist())),
            # Temporary forensic lookups only. Do not export all photo identities;
            # the receipt contains at most 20 actual colliding public photo IDs.
            "photo_to_taxon":dict(zip(photoids.tolist(),taxonids.tolist())),
            "taxon_to_species":dict(zip(taxonids.tolist(),frame.species.tolist())),
            "obs_to_taxon":(
                dict(zip(obsids.tolist(),taxonids.tolist())) if obsids is not None else None),
            "observations":set(obsids.tolist()) if obsids is not None else None}


def inspect_all(five: dict[str,dict], original_u100: pd.DataFrame) -> dict:
    if set(five)!=set(COHORTS):
        raise RuntimeError("all five original independent source photo cohorts mandatory")
    if ("inat_taxon_id" not in original_u100
        or "species" not in original_u100
        or "final_status" not in original_u100
        or len(original_u100)!=FULL_U100
        or original_u100.inat_taxon_id.nunique()!=FULL_U100):
        raise RuntimeError("original 4730 taxon frame missing or modified")
    seen_photo:set[int]=set()
    seen_taxa:set[int]=set()
    seen_names:set[str]=set()
    seen_obs:set[int]=set()
    studied:set[tuple[int,str]]=set()
    rows=0
    output={}
    for group in COHORTS:
        item=five[group]
        expected=ACTUAL[group]
        if item["taxa"]!=expected["taxa"] or item["photo_rows"]!=expected["n"]:
            raise RuntimeError(f"{group}: original count contract changed")
        overlap={
            "photos":sorted(item["photographs"] & seen_photo),
            "taxon_ids":sorted(item["taxon_ids"] & seen_taxa),
            "species_names":sorted(item["species_names"] & seen_names),
        }
        obs=item["observations"]
        overlap["observation_ids"]=(
            sorted(obs & seen_obs) if obs is not None else [])
        if any(overlap.values()):
            # Precise biological identity versus numeric photography collisions
            # must never be conflated: a shared binomial can involve DISTINCT
            # iNaturalist taxon IDs; only the observation/photo is the image.
            details={
                "incoming_cohort":group,
                "collision_counts":{k:len(v) for k,v in overlap.items()},
                "collision_examples":{k:v[:20] for k,v in overlap.items()},
                "collision_photo_taxon_records":[
                    {"photo_id":int(pid),
                     "new_taxon_id":int(item["photo_to_taxon"][pid]),
                     "new_species":item["taxon_to_species"][item["photo_to_taxon"][pid]],
                     "old_source_cohorts":[
                         {"cohort":earlier,
                          "taxon_id":int(five[earlier]["photo_to_taxon"][pid]),
                          "species":five[earlier]["taxon_to_species"][
                              five[earlier]["photo_to_taxon"][pid]]}
                         for earlier in output
                         if pid in five[earlier]["photographs"]]}
                    for pid in overlap["photos"][:20]
                ],
                "collision_observation_taxon_records":[
                    {"observation_id":int(oid),
                     "new_taxon_id":int(item["obs_to_taxon"][oid]),
                     "old_source_cohorts":[
                         {"cohort":earlier,"taxon_id":int(five[earlier]["obs_to_taxon"][oid])}
                         for earlier in output
                         if five[earlier]["observations"] is not None
                         and oid in five[earlier]["observations"]]}
                    for oid in overlap["observation_ids"][:20]
                ],
                "previous_cohort_breakdown":{
                    before:{
                        "photos":len(item["photographs"] & five[before]["photographs"]),
                        "taxon_ids":len(item["taxon_ids"] & five[before]["taxon_ids"]),
                        "species_names":len(item["species_names"] & five[before]["species_names"]),
                        "observation_ids":(
                            len(obs & five[before]["observations"])
                            if obs is not None and five[before]["observations"] is not None
                            else None),
                    }
                    for before in output
                },
            }
            print("FCP_ORIGINAL_FIVE_COHORT_REAL_IDENTITY_COLLISION_DIAGNOSIS",
                  json.dumps(details,sort_keys=True,allow_nan=False),flush=True)
            raise RuntimeError(f"{group}: real source identity collision; separate diagnostic counts emitted")
        seen_photo.update(item["photographs"])
        seen_taxa.update(item["taxon_ids"])
        seen_names.update(item["species_names"])
        if obs is not None:
            seen_obs.update(obs)
        studied.update(item["taxon_species"])
        rows+=item["photo_rows"]
        output[group]={
            "real_authentic_original_photo_IDs":item["photo_rows"],
            "species_with_exactly_100_original_source_IDs":item["taxa"],
            "observation_ID_field_available":obs is not None,
            "original_image_pixels_reopened":False,
        }
    if len(seen_photo)!=FULL_PHOTOS or len(seen_taxa)!=AUTHORIZED or rows!=FULL_PHOTOS:
        raise RuntimeError("five-cohort source U100 sample denominator not complete")
    full_taxa=set(original_u100.inat_taxon_id.astype(int))
    if not seen_taxa.issubset(full_taxa):
        raise RuntimeError("source photo taxon ID absent from original U100")
    manifest=set(zip(original_u100.inat_taxon_id.astype(int),
                     original_u100.species.astype(str).str.strip()))
    if studied-manifest:
        raise RuntimeError("source taxon-name relation changed from full original U100")
    counts=original_u100.final_status.value_counts().to_dict()
    baseline={
        "original_149900_photo_measured":sum(ACTUAL[k]["taxa"] for k in COHORTS[:3]),
        "wave4_49800_photo_id_authorized":ACTUAL["wave4"]["taxa"],
        "new_metadata_100_of_100_photo_IDs":ACTUAL["new_U100"]["taxa"],
        "new_H9_metadata_under_100_or_request_error":UNMET,
    }
    if counts!=baseline:
        raise RuntimeError("original 4730 completeness/under-100 species ledger drift")
    if set(original_u100.loc[original_u100.final_status.eq("new_H9_metadata_under_100_or_request_error"),
                            "inat_taxon_id"]) & seen_taxa:
        raise RuntimeError("original less-than-100 photo taxon imputed into measured sample")
    if len(studied)!=AUTHORIZED:
        raise RuntimeError("the same plant taxon carried inconsistent original names")
    return {
        "schema":"fcp_original_full_u100_4730_five_cohort_real_source_id_audit_v1",
        "eligible_original_U100_taxa":FULL_U100,
        "all_five_cohort_original_complete100_photo_taxa":AUTHORIZED,
        "original_eligible_taxa_below_historical_fixed100":UNMET,
        "cross_five_cohort_real_photo_ID_overlaps":0,
        "cross_five_cohort_taxon_ID_overlaps":0,
        "cross_five_cohort_species_name_overlaps":0,
        "cross_five_cohort_observation_ID_overlaps_where_recorded":0,
        "actual_original_source_photo_IDs":rows,
        "all_source_photo_taxon_names_match_original_U100_ledger":True,
        "original_149900_existing_measured_rows_unchanged":True,
        "new_photo_IDs_measured_by_this_audit":False,
        "photographic_colour_assays_or_original_nulls_run":False,
        "not_biological_inference":True,
        "cohorts":output,
    }


def main() -> None:
    ap=argparse.ArgumentParser()
    for key in OLD:
        ap.add_argument("--"+key,required=True,type=Path)
    ap.add_argument("--wave4",required=True,type=Path)
    ap.add_argument("--wave4-metadata-receipt",required=True,type=Path)
    ap.add_argument("--new",required=True,type=Path)
    ap.add_argument("--new-metadata-receipt",required=True,type=Path)
    ap.add_argument("--original-U100",required=True,type=Path)
    ap.add_argument("--out",required=True,type=Path)
    args=ap.parse_args()
    result={}
    for key,item in OLD.items():
        path=getattr(args,key)
        if sha(path)!=item["sha256"]:
            raise RuntimeError(f"{key}: frozen 149900 input bytes modified")
        result[key]=checked_source(pd.read_csv(path,low_memory=False),key,
                                   n=item["n"],taxa=item["taxa"])
    w=json.loads(args.wave4_metadata_receipt.read_text())
    wave=args.wave4
    if (w.get("schema")!="fcp_unchanged_highdepth_wave4_original_metadata_v1"
        or w.get("new_exact_100_species")!=498
        or w.get("new_exact_100_photo_rows")!=49800
        or w.get("outputs",{}).get(wave.name)!=sha(wave)):
        raise RuntimeError("unauthenticated original wave4 photo source metadata")
    result["wave4"]=checked_source(pd.read_csv(wave,low_memory=False),
                                    "wave4",n=49800,taxa=498)
    n=json.loads(args.new_metadata_receipt.read_text())
    new=args.new
    if (n.get("schema")!="fcp_u100_all4730_original_h9_only_v1"
        or n.get("stage")!="all_4730_H9_metadata_terminal"
        or n.get("verified_against_all_4_complete_shards") is not True
        or n.get("new_complete_100_photo_species")!=2481
        or n.get("new_unmet_100_photo_species")!=252
        or n.get("new_original_photo_IDs")!=248100
        or n.get("outputs_sha256",{}).get(new.name)!=sha(new)):
        raise RuntimeError("original 2733 U100 globally reconciled H9 source unverified")
    result["new_U100"]=checked_source(pd.read_csv(new,low_memory=False),
                                        "new_U100",n=248100,taxa=2481)
    full=args.original_U100
    if n["outputs_sha256"].get(full.name)!=sha(full):
        raise RuntimeError("original full U100 species status SHA drift")
    final=inspect_all(result,pd.read_csv(full,low_memory=False))
    final["input_original_sha256"]={k:sha(getattr(args,k)) for k in OLD}
    final["wave4_original_source_sha256"]=sha(wave)
    final["new_U100_original_source_sha256"]=sha(new)
    final["all_U100_species_status_sha256"]=sha(full)
    args.out.mkdir(parents=True,exist_ok=True)
    (args.out/"five_cohort_original_447800_photo_source_identity_audit.json").write_text(
        json.dumps(final,indent=2,sort_keys=True,allow_nan=False)+"\n")
    print("FCP_4730_FIVE_COHORT_REAL_PHOTO_SOURCE_IDS_VERIFIED",
          json.dumps(final,sort_keys=True,allow_nan=False))


if __name__=="__main__":
    main()

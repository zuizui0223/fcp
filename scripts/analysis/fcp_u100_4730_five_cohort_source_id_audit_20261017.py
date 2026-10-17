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
from scripts.analysis.fcp_u100_original_photo_reuse_20261017 import (
    PINNED_OLD_DISCOVERY_PHOTO_COLLISIONS,
    PINNED_OLD_DISCOVERY_OBSERVATION_COLLISIONS,
)

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
FINAL_EFFECTIVE_TAXA = 4475
FINAL_EFFECTIVE_PHOTOS = 447500
EXCLUDED_NEW_TAXA = frozenset({62367, 75633, 1045161})
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
    """Verify all raw photo identities, then freeze only the independent 447500."""
    if set(five)!=set(COHORTS):
        raise RuntimeError("all five original complete-100-photo source cohorts required")
    needed={"inat_taxon_id","species","final_status"}
    if (needed-set(original_u100)
        or len(original_u100)!=FULL_U100
        or original_u100.inat_taxon_id.nunique()!=FULL_U100):
        raise RuntimeError("original U100 4730 taxon frame invalid")
    seen_photo=set()
    seen_obs=set()
    seen_taxa=set()
    seen_species=set()
    studied=set()
    output={}
    all_reused_photo=set()
    all_reused_obs=set()
    raw_rows=0
    for group in COHORTS:
        item=five[group]
        contract=ACTUAL[group]
        if item["taxa"]!=contract["taxa"] or item["photo_rows"]!=contract["n"]:
            raise RuntimeError(f"{group}: original 100-photo cohort counts drift")
        if item["observations"] is None:
            raise RuntimeError(f"{group}: required source observation IDs missing")
        photo_reused=seen_photo & item["photographs"]
        obs_reused=seen_obs & item["observations"]
        if item["taxon_ids"] & seen_taxa or item["species_names"] & seen_species:
            raise RuntimeError(f"{group}: original taxon ID/name repeated across cohorts")
        if group!="new_U100" and (photo_reused or obs_reused):
            raise RuntimeError(f"{group}: earlier frozen original source identity reused")
        all_reused_photo.update(photo_reused)
        all_reused_obs.update(obs_reused)
        seen_photo.update(item["photographs"])
        seen_obs.update(item["observations"])
        seen_taxa.update(item["taxon_ids"])
        seen_species.update(item["species_names"])
        studied.update(item["taxon_species"])
        raw_rows+=item["photo_rows"]
        output[group]={
            "raw_original_source_photo_ID_occurrences":item["photo_rows"],
            "original_taxa_with_100_selected_photo_IDs":item["taxa"],
            "cross_previous_source_photo_IDs":sorted(photo_reused),
            "cross_previous_source_observation_IDs":sorted(obs_reused),
            "photo_pixels_read_or_reclassified":False,
        }
    if (all_reused_photo!=PINNED_OLD_DISCOVERY_PHOTO_COLLISIONS
        or all_reused_obs!=PINNED_OLD_DISCOVERY_OBSERVATION_COLLISIONS):
        raise RuntimeError("real old/new photo/observation reuse differs from frozen evidence")
    # The source-photo identifier collisions belong to three different newly
    # H9-authorized taxa; never select an alternative photo or keep 99 rows.
    new=five["new_U100"]
    affected={int(new["photo_to_taxon"][p]) for p in all_reused_photo}
    affected|={int(new["obs_to_taxon"][oid]) for oid in all_reused_obs}
    if affected!=EXCLUDED_NEW_TAXA:
        raise RuntimeError("three original source-reuse taxa were not correctly identified")
    if (raw_rows!=FULL_PHOTOS or len(seen_photo)!=FULL_PHOTOS-len(all_reused_photo)
        or len(seen_taxa)!=AUTHORIZED or len(studied)!=AUTHORIZED
        or raw_rows-100*len(affected)!=FINAL_EFFECTIVE_PHOTOS
        or len(seen_taxa)-len(affected)!=FINAL_EFFECTIVE_TAXA):
        raise RuntimeError("447500 fixed-100 source analysis denominator differs from reality")
    full_taxa=set(original_u100.inat_taxon_id.astype(int))
    if not seen_taxa.issubset(full_taxa):
        raise RuntimeError("original real source taxon absent from U100 frame")
    manifest=set(zip(original_u100.inat_taxon_id.astype(int),
                     original_u100.species.astype(str).str.strip()))
    if studied-manifest:
        raise RuntimeError("original photographic source taxon names differ from U100")
    baseline={
        "original_149900_photo_measured":sum(ACTUAL[k]["taxa"] for k in COHORTS[:3]),
        "wave4_49800_photo_id_authorized":ACTUAL["wave4"]["taxa"],
        "new_metadata_100_of_100_photo_IDs":ACTUAL["new_U100"]["taxa"],
        "new_H9_metadata_under_100_or_request_error":UNMET,
    }
    if original_u100.final_status.value_counts().to_dict()!=baseline:
        raise RuntimeError("original H9 complete/incomplete species ledger drift")
    incomplete=set(original_u100.loc[
        original_u100.final_status.eq("new_H9_metadata_under_100_or_request_error"),
        "inat_taxon_id"].astype(int))
    if incomplete & seen_taxa or incomplete & affected:
        raise RuntimeError("historical less-than-100 taxa accidentally included")
    return {
        "schema":"fcp_original_full_u100_4730_independent_447500_source_audit_v2",
        "original_eligible_U100_taxa":FULL_U100,
        "raw_H9_complete_100photo_source_taxa":AUTHORIZED,
        "raw_selected_original_source_photo_ID_occurrences":FULL_PHOTOS,
        "globally_distinct_original_source_photo_IDs":len(seen_photo),
        "globally_distinct_original_observation_IDs":len(seen_obs),
        "reused_frozen_discovery_photo_IDs":sorted(all_reused_photo),
        "reused_frozen_discovery_observation_IDs":sorted(all_reused_obs),
        "technically_excluded_entire_new_100photo_taxa":sorted(affected),
        "original_source_100photo_taxa_retain_for_analysis":FINAL_EFFECTIVE_TAXA,
        "original_source_100photo_rows_retain_for_analysis":FINAL_EFFECTIVE_PHOTOS,
        "original_U100_100photo_underfilled_taxa":UNMET,
        "source_reuse_excluded_additional_taxa":len(affected),
        "original_149900_photo_measurements_unchanged":True,
        "no_99_photo_species_reinsertion":True,
        "not_photo_colour_or_ecological_statistical_inference":True,
        "all_source_taxon_names_match_frozen_U100":True,
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

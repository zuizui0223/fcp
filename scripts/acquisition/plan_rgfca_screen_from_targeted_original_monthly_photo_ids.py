#!/usr/bin/env python3
"""Build colour-blind 64+64 RGFCA photo screening from EXACT original dated IDs.

Use after extract_rgfca_1000_original_monthly_photo_ids.py has scanned all
four authentic provider tables, audited global raw photo collisions, and
reconciled EVERY selected taxon count to independently verified 42,111 counts.

Unlike the original full-93-million staged SQLite path, the targeted complete
source photo-ID manifest now supplies ~3.43 million rows of the 1,000 *fixed*
taxa. Neither a live API query nor prior floral colour measurements can be
silently substituted.

Outputs exactly the same selected-source CSV/receipt API consumed by the
existing 128-partition original-image ROI-v4 blind freeze implementation.
Photos excluded by the bounded FIRST screen remain opportunities; do not
call those taxa monomorphic or change the frozen 200x250 species scheduling.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import heapq
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

from fcp_pipeline.photo_polymorphism_screen import sha256, observer_split
from scripts.acquisition.extract_rgfca_1000_original_monthly_photo_ids import (
    ARCHIVE_SHA, COUNT_SHA, CSV_FIELDS, FRAME_SHA, FULL_PHOTO_OPPORTUNITIES,
    RGFCA_ENROLLED_PHOTOS, ROSTER_SHA, SNAPSHOT, load_cohort,
)
from scripts.acquisition.plan_fcp_photo_first_screen import (
    SELECTION, SELECTED_FIELDS, SPECIES_FIELDS, _cell, _choice_hash, _select,
)

URL=re.compile(r"^https://inaturalist-open-data\.s3\.amazonaws\.com/photos/(\d+)/large\.(JPG|JPEG|PNG|jpg|jpeg|png)$")
SCHEMA="fcp_rgfca_targeted_original_20260927_photo_ID_census_v1"
STATUS="COMPLETE_TARGETED_1000_MONTHLY_SOURCE_IDS_COUNT_MATCHED_NO_PHOTO_PIXELS"


def create(manifest:Path, manifest_receipt:Path, frame:Path, counts:Path,
           independent_audit:Path, roster:Path, roster_receipt:Path,
           output:Path, *, photos_per_fold:int=64, max_per_observer:int=8,
           cell_candidate_limit:int=128, synthetic:bool=False)->dict:
    if (photos_per_fold<1 or photos_per_fold>512
        or max_per_observer<1 or max_per_observer>photos_per_fold
        or cell_candidate_limit<photos_per_fold or cell_candidate_limit>2048):
        raise ValueError("Invalid screen compute sample budget, not a biological species cap")
    if output.exists() and any(output.iterdir()):
        raise ValueError("Cannot replace earlier original pre-colour screen selection")
    names,selected_count=load_cohort(
        frame,counts,independent_audit,roster,roster_receipt,synthetic=synthetic)
    targeted=json.loads(manifest_receipt.read_text(encoding="utf-8"))
    if (targeted.get("schema")!=SCHEMA
        or targeted.get("status")!=STATUS
        or targeted.get("source_snapshot_date")!=SNAPSHOT
        or targeted.get("selected_all_photo_IDs_csv_sha256")!=sha256(manifest)
        or targeted.get("full_42111_photo_opportunity_counts_sha256")!=sha256(counts)
        or targeted.get("frozen_1000_source_enrolment_sha256")!=sha256(roster)
        or targeted.get("all_enrolled_species_match_independent_full_source_photo_ID_counts") is not True
        or targeted.get("original_image_pixels_opened") is not False
        or targeted.get("historical_flower_colour_CSV_imported") is not False
        or targeted.get("photos_NOT_classified_or_species_FCP_inferred") is not True
        or targeted.get("not_full_93million_ID_census_for_all_42111_taxa") is not True
        or (not synthetic and (
             targeted.get("source_archive_sha256")!=ARCHIVE_SHA
             or targeted.get("enrolled_outcome_blind_taxa")!=1000
             or targeted.get("full_original_frame_taxa")!=42111
             or targeted.get("n_selected_original_photo_IDs")!=RGFCA_ENROLLED_PHOTOS
             or sha256(frame)!=FRAME_SHA
             or sha256(counts)!=COUNT_SHA
             or sha256(roster)!=ROSTER_SHA))):
        raise ValueError("Targeted RGFCA photo screening lacks full original dated-source ID/collision provenance")

    # A single scan of all 3.43m metadata rows accumulates bounded candidate
    # heaps by 5-degree source cell and observer-heldout split per taxon.
    # Never touch source photo pixels, old colour outcomes, temperature/soil,
    # or iNaturalist API data during inclusion.
    output.mkdir(parents=True,exist_ok=True)
    selected_file=output/"photo_screen_original_metadata_selected.csv"
    opp_file=output/"photo_screen_all_species_opportunities.csv"
    seen_taxon=None
    current_photos=0
    scanned=0
    retained=0
    groups=[defaultdict(list),defaultdict(list)]
    selection_counts={}
    summed_counts=Counter()
    observed_species=set()
    selected_photo_ids=set()

    with selected_file.open("w",encoding="utf-8",newline="") as sink:
        writer=csv.DictWriter(sink,fieldnames=SELECTED_FIELDS,lineterminator="\n")
        writer.writeheader()
        def finish():
            nonlocal groups,retained,current_photos
            if seen_taxon is None:return
            if current_photos!=selected_count[seen_taxon]:
                raise ValueError(f"Source taxon {seen_taxon} exact original photo-ID count not preserved")
            picked=[_select(groups[f],photos_per_fold,max_per_observer)
                    for f in (0,1)]
            selection_counts[seen_taxon]=(len(picked[0]),len(picked[1]))
            for fold,photos in enumerate(picked):
                for score,pid,observation,observer,lat,lon,license_code,extension in photos:
                    if pid in selected_photo_ids:
                        raise ValueError("One global original photo ID was selected twice")
                    selected_photo_ids.add(pid)
                    writer.writerow({
                        "inat_taxon_id":seen_taxon,
                        "species":names[seen_taxon],
                        "photo_id":pid,"observation_id":observation,
                        "observer_id":observer,"latitude":lat,"longitude":lon,
                        "photo_url_large":f"https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/large.{extension}",
                        "photo_license":license_code,
                        "source_snapshot_date":SNAPSHOT,"screen_fold":fold,
                        "selection_hash":f"{score:016x}",
                    })
                    retained+=1
            groups=[defaultdict(list),defaultdict(list)]
            current_photos=0

        with manifest.open("r",encoding="utf-8",newline="") as raw:
            reader=csv.DictReader(raw)
            if tuple(reader.fieldnames or ())!=CSV_FIELDS:
                raise ValueError("RGFCA targeted original four-table metadata columns changed")
            for record in reader:
                tid_str,pid=record["inat_taxon_id"],record["photo_id"]
                if (not tid_str.isdecimal() or not pid.isdecimal()
                    or int(tid_str) not in selected_count
                    or record["species"]!=names[int(tid_str)]
                    or record["source_snapshot_date"]!=SNAPSHOT):
                    raise ValueError("RGFCA source photo taxon/date/ID switched from original frame")
                tid=int(tid_str)
                if seen_taxon!=tid:
                    if seen_taxon is not None and tid<seen_taxon:
                        raise ValueError("Original source photo IDs not sorted by taxon")
                    finish()
                    seen_taxon=tid
                    observed_species.add(tid)
                url=record["photo_url_large"]
                expected_url=URL.fullmatch(url)
                if not expected_url or expected_url.group(1)!=pid:
                    raise ValueError("Original monthly source photo ID/URL mismatch")
                uuid=record["observation_id"].strip()
                observer=record["observer_id"].strip()
                if (not uuid or not observer.isdecimal()
                    or not record["photo_license"] or not record["attribution"]
                    or record["source_snapshot_date"]!=SNAPSHOT):
                    raise ValueError("Original photographic observation/creator or licensing absent")
                lat=float(record["latitude"])
                lon=float(record["longitude"])
                if not (math.isfinite(lat) and math.isfinite(lon)
                        and abs(lat)<=90 and abs(lon)<=180):
                    raise ValueError("Original source photo geographic metadata missing")
                fold=observer_split(tid,observer)
                cell=_cell(lat,lon)
                hashscore=_choice_hash(tid,pid)
                heap=groups[fold][cell]
                candidate=(-hashscore,-int(pid),
                          (hashscore,pid,uuid,observer,lat,lon,
                           record["photo_license"],expected_url.group(2)))
                if len(heap)<cell_candidate_limit:
                    heapq.heappush(heap,candidate)
                elif candidate>heap[0]:
                    heapq.heapreplace(heap,candidate)
                scanned+=1
                current_photos+=1
                summed_counts[tid]+=1
        finish()
    if observed_species!=set(selected_count) or any(
            summed_counts[tid]!=n for tid,n in selected_count.items()):
        raise ValueError("Targeted original monthly photo-ID manifest lost one of the 1000 taxa")
    if scanned!=targeted["n_selected_original_photo_IDs"]:
        raise ValueError("Selected dated-source photo-ID scan differs from sealed source receipt")
    if not synthetic and (scanned!=RGFCA_ENROLLED_PHOTOS or retained>128000):
        raise ValueError("RGFCA new-photo source selection violates original study frame")

    with counts.open("r",encoding="utf-8",newline="") as file:
        source=list(csv.DictReader(file))
    with opp_file.open("w",encoding="utf-8",newline="") as dest:
        writer=csv.DictWriter(dest,fieldnames=SPECIES_FIELDS,lineterminator="\n")
        writer.writeheader()
        for row in source:
            tid=int(row["inat_taxon_id"])
            n=int(row["n_unique_source_photo_ids"])
            fold0,fold1=selection_counts.get(tid,(0,0))
            if fold0+fold1>n:
                raise ValueError("First colour-blind sample exceeds source photo opportunities")
            writer.writerow({
                "inat_taxon_id":tid,"species":row["species"],
                "source_photo_opportunities":n,
                "selected_discovery_photos":fold0,
                "selected_holdout_photos":fold1,
                "total_selected_for_image_screen":fold0+fold1,
                "source_remaining_unmeasured":n-fold0-fold1,
                "screen_unmeasured_is_not_monomorphism":True,
            })
    result={
        "schema":"fcp_photo_first_42111_outcome_blind_sample_schedule_v1",
        "status":"PHOTO_IDS_SELECTED_FROM_FULL_ORIGINAL_SOURCE_NOT_YET_MEASURED",
        "source_snapshot_date":SNAPSHOT,
        "source_species":len(names),
        "original_source_photo_ids_globally_considered":(
            FULL_PHOTO_OPPORTUNITIES if not synthetic
            else sum(int(x["n_unique_source_photo_ids"]) for x in source)),
        "original_RGFCA_1000_full_taxon_photo_IDs_checked":scanned,
        "RGFCA_enrolment_species_count":len(selected_count),
        "RGFCA_enrolment_selected_taxa_sha256":sha256(roster),
        "RGFCA_enrolment_metadata_only_no_literature_or_colour":True,
        "first_screen_photos_per_observer_fold_not_biological_cap":photos_per_fold,
        "per_observer_max_in_one_bounded_screen":max_per_observer,
        "geographic_cell_degrees":5,
        "outcome_blind_image_selection_no_floral_colour_opened":True,
        "source_photo_counts_sha256":sha256(counts),
        "original_source_frame_sha256":sha256(frame),
        "selected_photo_IDs":retained,
        "photo_selection_csv_sha256":sha256(selected_file),
        "species_opportunity_ledger_sha256":sha256(opp_file),
        "original_1000_taxa_targeted_monthly_manifest_sha256":sha256(manifest),
        "original_targeted_photo_ID_receipt_sha256":sha256(manifest_receipt),
        "selected_sample_NOT_full_93million_photo_census":True,
        "selected_sample_photo_pixels_measured":False,
        "photographic_polymorphism_classifications_computed":False,
        "single_photo_state_will_not_be_reported_as_monomorphic":True,
        "literature_polymorphism_labels_used":False,
        "prior_149900_measured_colour_records_used":False,
        "warning":"Photo IDs colour-blind selection is NOT full original image pixels. The 42000 non-enrolled taxa are unresolved not monomorphic; structural flower organs can mimic morphs.",
    }
    (output/"photo_screen_source_selection_receipt.json").write_text(
        json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    return result


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    for flag in ("targeted-dated-photo-ID-CSV","targeted-dated-photo-ID-receipt",
                 "original-42111-frame","original-42111-counts",
                 "independent-source-audit","rgfca-1000-roster","rgfca-roster-receipt",
                 "outdir"):
        p.add_argument("--"+flag,type=Path,required=True)
    p.add_argument("--photos-per-observer-fold",type=int,default=64)
    p.add_argument("--max-per-observer",type=int,default=8)
    a=p.parse_args(argv)
    report=create(a.targeted_dated_photo_ID_CSV,a.targeted_dated_photo_ID_receipt,
                  a.original_42111_frame,a.original_42111_counts,
                  a.independent_source_audit,a.rgfca_1000_roster,
                  a.rgfca_roster_receipt,a.outdir,
                  photos_per_fold=a.photos_per_observer_fold,
                  max_per_observer=a.max_per_observer)
    print("FCP_1000_REAL_ORIGINAL_DATED_SOURCE_PHOTO_SCREEN_PLAN",
          json.dumps({k:report[k] for k in (
              "source_species","original_RGFCA_1000_full_taxon_photo_IDs_checked",
              "selected_photo_IDs","selected_sample_photo_pixels_measured")},sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())

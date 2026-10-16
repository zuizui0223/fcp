#!/usr/bin/env python3
"""Gate the FIRST genuine original-monthly RGFCA 64-photo ROI-v4 compute unit.

Works only AFTER original four-table 2026-09-27 source IDs were globally
deduplicated/reconciled, outcome-blind photos selected, and all 128 worker
partitions SHA256 frozen. The routine does NOT download/open any image, infer
floral colours, or recycle prior live-API pilot photographs.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path

from fcp_pipeline.photo_polymorphism_screen import sha256
from scripts.acquisition.freeze_fcp_photo_first_blind_workers import (
    ACQUISITION, FREEZE_STATUS, WORKER,
)
from scripts.acquisition.run_fcp_photo_first_roi_work_unit import ScreenSourceIndex
from scripts.acquisition.extract_rgfca_1000_original_monthly_photo_ids import (
    ARCHIVE_SHA, COUNT_SHA, FRAME_SHA, RGFCA_ENROLLED_PHOTOS,
    ROSTER_SHA, SNAPSHOT,
)

TARGET_SCHEMA="fcp_rgfca_targeted_original_20260927_photo_ID_census_v1"
TARGET_STATUS="COMPLETE_TARGETED_1000_MONTHLY_SOURCE_IDS_COUNT_MATCHED_NO_PHOTO_PIXELS"
SELECTION_SCHEMA="fcp_photo_first_42111_outcome_blind_sample_schedule_v1"
SELECTION_STATUS="PHOTO_IDS_SELECTED_FROM_FULL_ORIGINAL_SOURCE_NOT_YET_MEASURED"
N_PART=128
BATCH_SIZE=64


def read_doc(path:Path):
    q=json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(q,dict):
        raise ValueError("Original photo provenance receipt is not an object")
    return q


def prepare(source_root:Path,outdir:Path,*,synthetic=False)->dict:
    if outdir.exists() and any(outdir.iterdir()):
        raise ValueError("Never overwrite a preselected original ROI work unit")
    full=source_root/"rgfca_1000_exact_source_photo_IDs"
    manifest=full/"rgfca_1000_complete_original_source_photo_IDs.csv"
    original_receipt=full/"rgfca_1000_targeted_monthly_photo_ID_receipt.json"
    screening=source_root/"rgfca_1000_photo_colour_blind_selection"
    selected=screening/"photo_screen_original_metadata_selected.csv"
    sample_receipt=screening/"photo_screen_source_selection_receipt.json"
    firewall=source_root/"rgfca_1000_original_photo_blind_worker_freeze"
    freeze_receipt=firewall/"photo_screen_freeze.json"
    original=read_doc(original_receipt)
    sample=read_doc(sample_receipt)
    frozen=read_doc(freeze_receipt)
    original_sha=sha256(manifest)
    sample_sha=sha256(selected)
    if (original.get("schema")!=TARGET_SCHEMA
        or original.get("status")!=TARGET_STATUS
        or original.get("source_snapshot_date")!=SNAPSHOT
        or original.get("selected_all_photo_IDs_csv_sha256")!=original_sha
        or original.get("all_enrolled_species_match_independent_full_source_photo_ID_counts") is not True
        or original.get("original_image_pixels_opened") is not False
        or original.get("historical_flower_colour_CSV_imported") is not False
        or original.get("photos_NOT_classified_or_species_FCP_inferred") is not True
        or sample.get("schema")!=SELECTION_SCHEMA
        or sample.get("status")!=SELECTION_STATUS
        or sample.get("source_snapshot_date")!=SNAPSHOT
        or sample.get("original_1000_taxa_targeted_monthly_manifest_sha256")!=original_sha
        or sample.get("original_targeted_photo_ID_receipt_sha256")!=sha256(original_receipt)
        or sample.get("photo_selection_csv_sha256")!=sample_sha
        or sample.get("outcome_blind_image_selection_no_floral_colour_opened") is not True
        or sample.get("selected_sample_photo_pixels_measured") is not False
        or sample.get("prior_149900_measured_colour_records_used") is not False
        or frozen.get("status")!=FREEZE_STATUS
        or frozen.get("source_snapshot_date")!=SNAPSHOT
        or frozen.get("worker_partitions")!=N_PART
        or frozen.get("pre_colour_selection_csv_sha256")!=sample_sha
        or frozen.get("source_photo_screen_selection_receipt_sha256")!=sha256(sample_receipt)
        or frozen.get("image_worker_sees_only_three_fields") is not True
        or frozen.get("no_old_photo_colour_rows_loaded") is not True
        or frozen.get("no_literature_positive_labels_used") is not True
        or frozen.get("original_image_pixels_opened") is not False
        or frozen.get("colour_states_inferred") is not False
        or frozen.get("screen_selected_original_photo_IDs")!=sample.get("selected_photo_IDs")):
        raise ValueError("No SHA-sealed original dated photo-ID→blind-worker source")

    if not synthetic and (
        original.get("n_selected_original_photo_IDs")!=RGFCA_ENROLLED_PHOTOS
        or original.get("enrolled_outcome_blind_taxa")!=1000
        or original.get("full_original_frame_taxa")!=42111
        or original.get("source_archive_sha256")!=ARCHIVE_SHA
        or original.get("frozen_1000_source_enrolment_sha256")!=ROSTER_SHA
        or original.get("full_42111_photo_opportunity_counts_sha256")!=COUNT_SHA
        or sample.get("source_species")!=42111
        or frozen.get("source_frame_sha256")!=FRAME_SHA
        or frozen.get("source_photo_counts_sha256")!=COUNT_SHA
        or not 0<sample.get("selected_photo_IDs",0)<=128000):
        raise ValueError("Original Sept 2026 bulk-photo RGFCA cohort identity drift")

    expected=frozen.get("frozen_file_sha256",{})
    inventory={"metadata_join_key.csv"}
    for n in range(N_PART):
        inventory.update({f"worker_partition_{n:03}.csv",
                          f"acquisition_partition_{n:03}.csv"})
    if set(expected)!=inventory:
        raise ValueError("All 128 frozen blind workers and source keys required")
    for filename in sorted(inventory):
        if expected[filename]!=sha256(firewall/filename):
            raise ValueError("Immutable original photo worker partition/source key changed")

    first=None
    for part in range(N_PART):
        with (firewall/f"worker_partition_{part:03}.csv").open(
             "r",newline="",encoding="utf-8") as stream:
            rows=csv.DictReader(stream)
            if tuple(rows.fieldnames or ())!=WORKER:
                raise ValueError("Original colour-blind worker column contract changed")
            if next(rows,None) is not None:
                first=part
                break
    if first is None:
        raise ValueError("No original source photos in immutable RGFCA worker manifest")
    index=ScreenSourceIndex(firewall,first,synthetic=synthetic)
    worker,acq,n_batches=index.get_batch(first,0,BATCH_SIZE)
    if (not 0<len(worker)<=BATCH_SIZE or tuple(worker.columns)!=WORKER
        or tuple(acq.columns)!=ACQUISITION
        or any(not url.startswith(
            "https://inaturalist-open-data.s3.amazonaws.com/photos/")
            for url in acq.photo_url_large)):
        raise ValueError("ROI worker is not a bounded original source photograph selection")
    ids=set(worker.measurement_id)
    seen={}
    with (firewall/"metadata_join_key.csv").open(
         "r",newline="",encoding="utf-8") as src:
        for row in csv.DictReader(src):
            mid=row["measurement_id"]
            if mid in ids:
                if mid in seen:
                    raise ValueError("Source ROI photo ID keyed more than once")
                seen[mid]=row["photo_id"]
    if set(seen)!=ids or len(set(seen.values()))!=len(ids):
        raise ValueError("Blinded source measurement/photo biological metadata join incomplete")

    doc={
        "schema":"fcp_rgfca_first_original_dated_64photo_roi_unit_v1",
        "status":"ORIGINAL_SOURCE_64_PHOTO_ROI_UNIT_READY_NOT_MEASURED",
        "source_snapshot_date":SNAPSHOT,
        "original_dated_monthly_archive_sha256":original["source_archive_sha256"],
        "n_rgfca_enrolled_source_taxa":original["enrolled_outcome_blind_taxa"],
        "n_rgfca_all_original_source_photo_IDs":original["n_selected_original_photo_IDs"],
        "n_colblind_selected_source_photo_IDs":sample["selected_photo_IDs"],
        "n_photo_IDs_in_this_unit":len(worker),
        "partition":first,"batch":0,"batch_size":BATCH_SIZE,
        "other_compute_batches_in_partition":n_batches-1,
        "original_full_source_manifest_sha256":original_sha,
        "original_monthly_targeted_source_receipt_sha256":sha256(original_receipt),
        "precolour_photo_selection_receipt_sha256":sha256(sample_receipt),
        "frozen_screen_worker_sha256":sha256(freeze_receipt),
        "original_blind_worker_3_fields_only":True,
        "photo_IDs_selected_without_colour_or_environment":True,
        "original_image_pixels_opened":False,
        "original_ROI_v4_completed":False,
        "original_1000_taxon_photo_screen_complete":False,
        "photo_polymorphism_Ppoly_W_H2_climate_or_genetics_computed":False,
        "no_old_colours_or_live_API_photo_IDs_used":True,
    }
    outdir.mkdir(parents=True,exist_ok=True)
    (outdir/"first_original_ROI64_ready_receipt.json").write_text(
        json.dumps(doc,sort_keys=True,indent=2)+"\n",encoding="utf-8")
    return doc


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--admitted-monthly-source-artifact-root",required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    x=p.parse_args(argv)
    print("FCP_PREPARED_EXACT_MONTHLY_FIRST_BLIND_ROI64",
          json.dumps(prepare(x.admitted_monthly_source_artifact_root,x.outdir),
                     sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())

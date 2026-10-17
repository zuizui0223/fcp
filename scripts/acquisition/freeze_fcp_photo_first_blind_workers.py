#!/usr/bin/env python3
"""Freeze ONLY the outcome-blind ORIGINAL FCP photo-screen sample as ROI workers.

Never claims a complete 93,209,351-photo biological measurement. The original
image evaluator sees exactly (measurement_id,image_filename,photo_license).
Taxonomy, geography, photo URLs and C/S evidence remain outside its process.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import sqlite3
from contextlib import ExitStack
from pathlib import Path

from fcp_pipeline.photo_polymorphism_screen import (
    FRAME_SHA, SOURCE_COUNT_SHA, SOURCE_PHOTOS, SNAPSHOT, observer_split, sha256,
)
from scripts.acquisition.plan_fcp_photo_first_screen import (
    SELECTED_FIELDS, SELECTION, _choice_hash,
)
from scripts.acquisition.build_fcp_uncapped_blind_batches import new_id

PARTITIONS = 128
WORKER = ("measurement_id", "image_filename", "photo_license")
ACQUISITION = (*WORKER[:2], "photo_url_large", "photo_license")
METADATA = (
    "measurement_id", "inat_taxon_id", "species", "photo_id",
    "observation_id", "observer_id", "latitude", "longitude", "screen_fold",
)
FREEZE_STATUS = "SEALED_ORIGINAL_PHOTO_SCREEN_SOURCE_BLIND_NOT_FULL_CENSUS"
URL = re.compile(r"^https://inaturalist-open-data\.s3\.amazonaws\.com/photos/(\d+)/large\.(?:jpg|jpeg|png)$", re.I)


def freeze(selected: Path, selection_receipt: Path, source_frame: Path,
           source_counts: Path, dest: Path, *, _synthetic=False) -> dict:
    if dest.exists() and any(dest.iterdir()):
        raise ValueError("Cannot overwrite an original photo-screen worker freeze")
    selection = json.loads(selection_receipt.read_text(encoding="utf-8"))
    if (selection.get("schema") != "fcp_photo_first_42111_outcome_blind_sample_schedule_v1"
        or selection.get("status") != "PHOTO_IDS_SELECTED_FROM_FULL_ORIGINAL_SOURCE_NOT_YET_MEASURED"
        or selection.get("source_snapshot_date") != SNAPSHOT
        or selection.get("outcome_blind_image_selection_no_floral_colour_opened") is not True
        or selection.get("selected_sample_photo_pixels_measured") is not False
        or selection.get("selected_sample_NOT_full_93million_photo_census") is not True
        or selection.get("literature_polymorphism_labels_used") is not False
        or selection.get("prior_149900_measured_colour_records_used") is not False
        or selection.get("photo_selection_csv_sha256") != sha256(selected)
        or selection.get("original_source_frame_sha256") != sha256(source_frame)
        or selection.get("source_photo_counts_sha256") != sha256(source_counts)
        or (not _synthetic and
            (selection.get("source_species") != 42111
             or selection.get("original_source_photo_ids_globally_considered") != SOURCE_PHOTOS
             or sha256(source_frame) != FRAME_SHA
             or sha256(source_counts) != SOURCE_COUNT_SHA))):
        raise ValueError("Photo screen has no valid pre-colour full original-source selection receipt")
    dest.mkdir(parents=True,exist_ok=True)
    protected = []
    count_by_fold = [0,0]
    n_rows = 0
    seen_db = dest / ".screen_identity_tmp.sqlite"
    try:
        with sqlite3.connect(seen_db) as db, ExitStack() as stack:
            db.execute("CREATE TABLE ids (photo_id TEXT PRIMARY KEY, observation TEXT, observer TEXT)")
            writers = []
            keys = []
            for part in range(PARTITIONS):
                wk = f"worker_partition_{part:03}.csv"
                ac = f"acquisition_partition_{part:03}.csv"
                protected.extend((wk, ac))
                handle = stack.enter_context((dest / wk).open("w", encoding="utf-8", newline=""))
                writer = csv.DictWriter(handle, fieldnames=WORKER, lineterminator="\n")
                writer.writeheader()
                writers.append(writer)
                handle = stack.enter_context((dest / ac).open("w", encoding="utf-8", newline=""))
                writer = csv.DictWriter(handle, fieldnames=ACQUISITION, lineterminator="\n")
                writer.writeheader()
                keys.append(writer)
            join_file = "metadata_join_key.csv"
            protected.append(join_file)
            join = stack.enter_context((dest / join_file).open("w", encoding="utf-8", newline=""))
            meta = csv.DictWriter(join, fieldnames=METADATA, lineterminator="\n")
            meta.writeheader()
            with selected.open("r", encoding="utf-8", newline="") as f:
                source = csv.DictReader(f)
                if list(source.fieldnames or []) != SELECTED_FIELDS:
                    raise ValueError("Original photo selection columns are not expected")
                for row in source:
                    tid, pid, obs, observer = (
                        row["inat_taxon_id"], row["photo_id"],
                        row["observation_id"], row["observer_id"])
                    if not tid.isdecimal() or not pid.isdecimal() or not obs or not observer:
                        raise ValueError("Missing original screen photo/observation/observer identity")
                    taxon = int(tid)
                    if (row["source_snapshot_date"] != SNAPSHOT
                        or row["screen_fold"] != str(observer_split(taxon,observer))
                        or row["selection_hash"] != f"{_choice_hash(taxon,pid):016x}"
                        or not row["species"].strip()
                        or not URL.fullmatch(row["photo_url_large"])
                        or URL.fullmatch(row["photo_url_large"]).group(1) != pid
                        or not row["photo_license"]):
                        raise ValueError("Photo-screen source metadata/observer split/URL drift")
                    try:
                        db.execute("INSERT INTO ids VALUES (?,?,?)",(pid,obs,observer))
                    except sqlite3.IntegrityError as exc:
                        raise ValueError("One source photo present twice in first-stage FCP screen") from exc
                    mid = new_id(pid)
                    filename = mid + ".jpg"
                    partition = int(hashlib.sha256(mid.encode()).hexdigest()[:16],16) % PARTITIONS
                    wk = dict(measurement_id=mid,image_filename=filename,
                              photo_license=row["photo_license"])
                    writers[partition].writerow(wk)
                    keys[partition].writerow({
                        "measurement_id":mid,"image_filename":filename,
                        "photo_url_large":row["photo_url_large"],
                        "photo_license":row["photo_license"],
                    })
                    meta.writerow({
                        "measurement_id":mid,"inat_taxon_id":taxon,
                        "species":row["species"],"photo_id":pid,
                        "observation_id":obs,"observer_id":observer,
                        "latitude":row["latitude"],"longitude":row["longitude"],
                        "screen_fold":row["screen_fold"],
                    })
                    count_by_fold[int(row["screen_fold"])] += 1
                    n_rows += 1
            db.commit()
        if n_rows != selection.get("selected_photo_IDs"):
            raise ValueError("New screen worker packets do not preserve the source-selected original photo IDs")
        protected_shas = {p:sha256(dest / p) for p in protected}
        frozen = {
            "schema": "fcp_original_photo_screen_blind_firewall_v1",
            "status": FREEZE_STATUS,
            "source_snapshot_date": SNAPSHOT,
            "source_frame_sha256": sha256(source_frame),
            "source_photo_counts_sha256": sha256(source_counts),
            "pre_colour_selection_csv_sha256": sha256(selected),
            "source_photo_screen_selection_receipt_sha256": sha256(selection_receipt),
            "source_species": selection["source_species"],
            "source_all_provider_photo_IDs": selection["original_source_photo_ids_globally_considered"],
            "screen_selected_original_photo_IDs": n_rows,
            "worker_partitions": PARTITIONS,
            "photo_screen_fold_0": count_by_fold[0],
            "photo_screen_fold_1": count_by_fold[1],
            "frozen_file_sha256": protected_shas,
            "image_worker_sees_only_three_fields": True,
            "source_screen_NOT_complete_photo_census": True,
            "original_image_pixels_opened": False,
            "no_old_photo_colour_rows_loaded": True,
            "no_literature_positive_labels_used": True,
            "colour_states_inferred": False,
        }
        freeze_file = dest / "photo_screen_freeze.json"
        freeze_file.write_text(json.dumps(frozen,ensure_ascii=False,
                                          sort_keys=True,indent=2)+"\n",encoding="utf-8")
        return frozen
    finally:
        for suffix in ("", "-shm", "-wal"):
            path=Path(str(seen_db)+suffix)
            if path.exists():
                path.unlink()


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--selected-metadata",required=True,type=Path)
    p.add_argument("--selection-receipt",required=True,type=Path)
    p.add_argument("--source-frame",required=True,type=Path)
    p.add_argument("--source-counts",required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    args=p.parse_args(argv)
    result=freeze(args.selected_metadata,args.selection_receipt,
                  args.source_frame,args.source_counts,args.outdir)
    print(json.dumps({k:result[k] for k in (
        "status","source_species","screen_selected_original_photo_IDs",
        "photo_screen_fold_0","photo_screen_fold_1")},sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())

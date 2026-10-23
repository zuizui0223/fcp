#!/usr/bin/env python3
"""Complete/reassemble ONLY the original FCP photo-first SCREEN image terminals.

This is a bounded sampled photo screen, not an all-93,209,351 photo census.
Every pre-colour selected screen ID needs one original ROI-v4 terminal status.
Outputs are location-aware only AFTER the independent blind ROI measurements.
The historical 149,900 measured outcomes cannot enter this reassembler.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sqlite3
from pathlib import Path

import pandas as pd

from fcp_pipeline.fcp_fresh_study import (
    BIO, MORPHS, validate_method_attestation,
)
from fcp_pipeline.photo_polymorphism_screen import (
    FRAME_SHA, SOURCE_COUNT_SHA, SOURCE_PHOTOS, SNAPSHOT, PROTOCOL, sha256,
)
from scripts.acquisition.freeze_fcp_photo_first_blind_workers import (
    ACQUISITION, FREEZE_STATUS, METADATA, WORKER,
)

RECEIPT_SCHEMA = "fcp_sealed_photo_screen_terminal_v1"
SCREEN_STATUS = "SEALED_PHOTO_SCREEN_UNBIASED_NOT_FINAL"
BATCH_SIZE = 64


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _check_source(firewall: Path, *, synthetic=False) -> dict:
    doc = _read(firewall/"photo_screen_freeze.json")
    if (doc.get("schema") != "fcp_original_photo_screen_blind_firewall_v1"
        or doc.get("status") != FREEZE_STATUS
        or doc.get("source_snapshot_date") != SNAPSHOT
        or doc.get("worker_partitions") != 128
        or doc.get("image_worker_sees_only_three_fields") is not True
        or doc.get("source_screen_NOT_complete_photo_census") is not True
        or doc.get("original_image_pixels_opened") is not False
        or doc.get("colour_states_inferred") is not False
        or doc.get("no_old_photo_colour_rows_loaded") is not True
        or doc.get("no_literature_positive_labels_used") is not True
        or (not synthetic and
            (doc.get("source_frame_sha256") != FRAME_SHA
             or doc.get("source_photo_counts_sha256") != SOURCE_COUNT_SHA
             or doc.get("source_all_provider_photo_IDs") != SOURCE_PHOTOS
             or doc.get("source_species") != 42111))):
        raise ValueError("Not an original pre-colour photo-screen blind firewall")
    for filename, digest in doc.get("frozen_file_sha256", {}).items():
        if sha256(firewall/filename) != digest:
            raise ValueError("Frozen screen worker or source metadata changed")
    return doc


def join(firewall: Path, sealed: Path, outdir: Path, *,
         _synthetic_test=False) -> dict:
    if outdir.exists() and any(outdir.iterdir()):
        raise ValueError("Do not overwrite original photo-screen reassembly")
    doc = _check_source(firewall, synthetic=_synthetic_test)
    outdir.mkdir(parents=True,exist_ok=True)
    db_file = outdir/".screen_join_work.sqlite"
    destination = outdir/"newly_remeasured_photo_screen.csv"
    receipt_path = outdir/"sealed_photo_screen_terminal.json"
    expected_count = doc["screen_selected_original_photo_IDs"]
    global_seen = 0
    roi_good = 0
    complete = False
    method_fingerprint = None
    try:
        with sqlite3.connect(db_file) as db:
            db.execute("""CREATE TABLE expected(
                measurement_id TEXT PRIMARY KEY,
                photo_id TEXT UNIQUE, taxon INTEGER NOT NULL, species TEXT NOT NULL,
                observation TEXT NOT NULL, observer TEXT NOT NULL,
                latitude TEXT, longitude TEXT, screen_fold INTEGER NOT NULL)""")
            db.execute("CREATE TABLE observed(measurement_id TEXT PRIMARY KEY, row_json TEXT NOT NULL)")
            join_file = firewall/"metadata_join_key.csv"
            with join_file.open("r",newline="",encoding="utf-8") as f:
                reader=csv.DictReader(f)
                if tuple(reader.fieldnames or ()) != METADATA:
                    raise ValueError("Screen selected metadata joining key changed")
                for row in reader:
                    if (not row["inat_taxon_id"].isdecimal()
                        or not row["photo_id"].isdecimal()
                        or row["screen_fold"] not in ("0","1")
                        or not row["species"] or not row["observer_id"]
                        or not row["observation_id"]):
                        raise ValueError("Unknown source screen observation / species")
                    try:
                        db.execute("INSERT INTO expected VALUES(?,?,?,?,?,?,?,?,?)", (
                            row["measurement_id"],row["photo_id"],
                            int(row["inat_taxon_id"]),row["species"],
                            row["observation_id"],row["observer_id"],
                            row["latitude"],row["longitude"],int(row["screen_fold"])))
                    except sqlite3.IntegrityError as exc:
                        raise ValueError("Repeated source-photo or measurement ID in screen join") from exc
            exp_rows=db.execute("SELECT COUNT(*) FROM expected").fetchone()[0]
            if exp_rows != expected_count:
                raise ValueError("Missing selected source photos before ROI reassembly")
            for part in range(128):
                worker_path=firewall/f"worker_partition_{part:03}.csv"
                worker=pd.read_csv(worker_path,dtype=str).fillna("")
                if tuple(worker.columns)!=WORKER or worker.measurement_id.duplicated().any():
                    raise ValueError("Original blind worker packet is invalid")
                partition_dir=sealed/f"partition_{part:03}"
                needed_batches=math.ceil(len(worker)/BATCH_SIZE)
                if needed_batches==0:
                    if partition_dir.is_dir() and list(partition_dir.glob("terminal_p*.csv")):
                        raise ValueError("Unexpected ROI terminal for empty selected photo partition")
                    continue
                if not partition_dir.is_dir():
                    raise ValueError(f"Missing original ROI terminal partition {part}")
                expected_names={f"terminal_p{part:03}_b{batch:07}.csv"
                                for batch in range(needed_batches)}
                actual_names={p.name for p in partition_dir.glob("terminal_p*.csv")}
                if expected_names != actual_names:
                    raise ValueError("Extra/missing screen terminal CSV work unit")
                marker=partition_dir/"frozen_method_attestation.json"
                if not _synthetic_test:
                    fingerprint=validate_method_attestation(marker)
                else:
                    fingerprint={"synthetic_test_only":True,
                                 "no_old_measurement_csv_inputs":True}
                if method_fingerprint is None:
                    method_fingerprint=fingerprint
                elif method_fingerprint != fingerprint:
                    raise ValueError("Different original ROI model/weights mixed across photo-screen batches")
                for batch in range(needed_batches):
                    stem=f"terminal_p{part:03}_b{batch:07}"
                    photo_csv=partition_dir/(stem+".csv")
                    terminal=_read(partition_dir/(stem+".json"))
                    wrapper=_read(partition_dir/(f"screen_p{part:03}_b{batch:07}.json"))
                    sub=worker.iloc[batch*BATCH_SIZE:(batch+1)*BATCH_SIZE].copy()
                    import hashlib
                    expected_subset_sha=hashlib.sha256(sub.to_csv(index=False).encode()).hexdigest()
                    if (terminal.get("status")!="COMPLETE_WORK_UNIT"
                        or terminal.get("partition")!=part or terminal.get("batch")!=batch
                        or terminal.get("n_source_photo_ids")!=len(sub)
                        or terminal.get("n_terminal_photo_ids")!=len(sub)
                        or terminal.get("selected_source_worker_sha256")!=expected_subset_sha
                        or terminal.get("terminal_csv_sha256")!=sha256(photo_csv)
                        or terminal.get("reused_historical_measurement_csv") is not False
                        or terminal.get("image_bytes_persisted") is not False
                        or terminal.get("full_study_complete") is not False
                        or terminal.get("method") != fingerprint
                        or wrapper.get("status")!="ORIGINAL_PHOTO_SCREEN_UNIT_MEASURED_NOT_ALL_FCP"
                        or wrapper.get("partition")!=part or wrapper.get("batch")!=batch
                        or wrapper.get("screen_firewall_sha256")!=sha256(firewall/"photo_screen_freeze.json")
                        or wrapper.get("terminal_sha256")!=sha256(photo_csv)
                        or wrapper.get("full_93million_original_photo_census_complete") is not False
                        or wrapper.get("historical_photo_colours_reused") is not False
                        or (wrapper.get("synthetic_test_only") is True and not _synthetic_test)):
                        raise ValueError("Sealed original ROI-v4 work unit is incomplete or changed")
                    with photo_csv.open("r",newline="",encoding="utf-8") as f:
                        records=csv.DictReader(f)
                        needed={"measurement_id","morph","measurement_status",
                                *(f"flower_fraction_{p}" for p in BIO)}
                        if not needed.issubset(records.fieldnames or []):
                            raise ValueError("Original ROI-v4 nine-colour output missing")
                        if {"inat_taxon_id","species","latitude","longitude",
                            "photo_url_large","photo_id"} & set(records.fieldnames or []):
                            raise ValueError("Original ROI worker was exposed to biology/source coordinates")
                        keys=[]
                        for row in records:
                            good=(row["measurement_status"]=="classified_four_state_morph")
                            if (good!=(row["morph"] in MORPHS)
                                or (not good and row["morph"]!="mixed_uncertain")):
                                raise ValueError("Photo colour state inconsistent with original ROI-v4 QC")
                            vals=[float(row[f"flower_fraction_{p}"]) for p in BIO]
                            if good and (any(not math.isfinite(x) or x<0 or x>1 for x in vals)
                                         or abs(sum(vals)-1)>1e-7):
                                raise ValueError("Nine-colour original ROI-v4 palette invalid")
                            keys.append(row["measurement_id"])
                            try:
                                db.execute("INSERT INTO observed VALUES(?,?)",
                                           (row["measurement_id"],json.dumps(row,sort_keys=True)))
                            except sqlite3.IntegrityError as exc:
                                raise ValueError("Duplicate original source ROI terminal") from exc
                            global_seen+=1
                            roi_good+=int(good)
                    if len(keys)!=len(sub) or set(keys)!=set(sub.measurement_id.astype(str)):
                        raise ValueError("ROI terminal contains wrong original source photo IDs")
            if (global_seen!=expected_count
                or db.execute("SELECT COUNT(*) FROM observed LEFT JOIN expected USING(measurement_id) WHERE expected.measurement_id IS NULL").fetchone()[0]
                or db.execute("SELECT COUNT(*) FROM expected LEFT JOIN observed USING(measurement_id) WHERE observed.measurement_id IS NULL").fetchone()[0]):
                raise ValueError("Screen missing or unexpected original source photograph ROI terminal")
            output_fields=["inat_taxon_id","species","photo_id","observation_id",
                           "observer_id","latitude","longitude","screen_fold",
                           "morph","measurement_status","global_classifiable",
                           *(f"flower_fraction_{c}" for c in BIO)]
            with destination.open("w",newline="",encoding="utf-8") as h:
                w=csv.DictWriter(h,fieldnames=output_fields,lineterminator="\n")
                w.writeheader()
                for tax,sp,pid,obs,observer,lat,lon,fold,raw in db.execute(
                    """SELECT e.taxon,e.species,e.photo_id,e.observation,e.observer,
                              e.latitude,e.longitude,e.screen_fold,o.row_json
                       FROM expected e JOIN observed o USING(measurement_id)
                       ORDER BY e.taxon,CAST(e.photo_id AS INTEGER)"""):
                    photo=json.loads(raw)
                    morph=photo["morph"]
                    w.writerow({
                        "inat_taxon_id":tax,"species":sp,"photo_id":pid,
                        "observation_id":obs,"observer_id":observer,
                        "latitude":lat,"longitude":lon,"screen_fold":fold,
                        "morph":morph,
                        "measurement_status":photo["measurement_status"],
                        "global_classifiable":str(morph in MORPHS).lower(),
                        **{f"flower_fraction_{c}":photo[f"flower_fraction_{c}"] for c in BIO},
                    })
        summary = {
            "schema":RECEIPT_SCHEMA,"status":SCREEN_STATUS,
            "protocol":PROTOCOL,
            "source_date":SNAPSHOT,
            "source_species":doc["source_species"],
            "source_frame_sha256":FRAME_SHA if _synthetic_test else doc["source_frame_sha256"],
            "new_photo_table_sha256":sha256(destination),
            "n_terminal_measurements":global_seen,
            "n_classifiable_original_flower_photos":roi_good,
            "photo_ID_selection_pre_colour":True,
            "source_photo_IDs_sealed":True,
            "original_ROI_v4_method_attested":not _synthetic_test,
            "historical_measured_flower_colours_reused":False,
            "screen_not_entire_original_93million_census":True,
            "full_source_image_measurement_complete":False,
            "source_screen_firewall_sha256":sha256(firewall/"photo_screen_freeze.json"),
            "genetic_FCP_claim":False,
            "W_H2_or_environment_inference":False,
            "synthetic_test_only":_synthetic_test,
        }
        receipt_path.write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n",
                                encoding="utf-8")
        complete=True
        return summary
    finally:
        if db_file.exists():
            db_file.unlink()
        if not complete:
            for p in (destination,receipt_path):
                p.unlink(missing_ok=True)


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--screen-firewall",required=True,type=Path)
    p.add_argument("--screen-terminal-root",required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    a=p.parse_args(argv)
    print(json.dumps(join(a.screen_firewall,a.screen_terminal_root,a.outdir),sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Extract the outcome-blind RGFCA 1,000 taxa from the EXACT monthly photo metadata.

Unlike the original full-42,111 SQL stage this keeps only selected taxon
observations + photo IDs on disk, while still reading/byte-hashing all four
dated source tables. A second independent pass over the WHOLE provider photos
table checks selected IDs for collisions with *any* other provider observation.
The resulting per-taxon IDs must agree EXACTLY with the separately audited
full-universe photo counts or the entire result is rejected.

This trades two streaming reads for dramatically smaller SQLite storage.
Neither source photo counts, an API-created-before query, nor 200 RGFCA
realizations are flower-colour pixel measurements.

The second raw-photo collision pass is deliberately CONSERVATIVE: if a
collision from a non-qualifying external observation could cause the selected
ID universe to differ from the original complete-frame importer, the exact
full-frame taxon-count gate will fail rather than silently inventing photos.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sqlite3
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from fcp_pipeline.fcp_allphotos import opportunity_frame, sha256_file
from scripts.acquisition.import_inat_open_data_snapshot import (
    accepted_extension, allowed_geo, integer, lineage_source_taxon,
    normalize_license, read_tsv, source_crosswalk, verify_source_archive_pin,
)

SNAPSHOT = "2026-09-27"
FRAME_SHA = "5d871bd1f1190c68fd513bb527e6c93de1606dc35ff8cfff5f6187850dc382cc"
COUNT_SHA = "ef26c16252dd785f2d08bde3b18938c7624d73a54386fa3d3974dceb09ed97a9"
ROSTER_SHA = "88aba565e1bdeef48ad3dc5eec0121c52448fda434aca849e542bb615f9d547c"
ARCHIVE_SHA = "bfc1a97c436bf5989fccd43f12392953128489b76b92cea24da335ad73da4d8e"
ORIGINAL_TABLE_HASHES = {
    "taxa": "71c1b2909c3c1c409da71c6f778939b00c5e1cc97a10449b755aa569c700e2d0",
    "observations": "6c8b93cd808bdb0503690a60719a2d0a6478d3494eda7c8623db8c7f454b125b",
    "observers": "735d8112f24c1d1951622bc09fee4a6a07850d87937b3aa77d6c65289e357d63",
    "photos": "0033f67f2362acae0509314c4a7b8429f1f6ff85323cb3c86eeb7cbb83ecc6b6",
}
FULL_TAXA = 42111
ENROLLED = 1000
FULL_PHOTO_OPPORTUNITIES = 93209351
RGFCA_ENROLLED_PHOTOS = 3430502
CSV_FIELDS = (
    "inat_taxon_id", "species", "photo_id", "observation_id", "observer_id",
    "latitude", "longitude", "observed_on", "positional_accuracy_m",
    "photo_license", "photo_url_large", "source_snapshot_date",
)


def load_cohort(frame: Path, count_path: Path, independent_audit: Path,
                roster_path: Path, roster_receipt: Path,
                *, synthetic: bool = False) -> tuple[dict[int,str],dict[int,int]]:
    if not synthetic and any((sha256_file(frame) != FRAME_SHA,
                              sha256_file(count_path) != COUNT_SHA,
                              sha256_file(roster_path) != ROSTER_SHA)):
        raise ValueError("Original outcome-blind species frame/census/RGFCA enrolment SHA mismatch")
    if synthetic:
        with frame.open("r", encoding="utf-8", newline="") as h:
            expected_taxa=sum(1 for _ in csv.DictReader(h))
    else:
        expected_taxa=FULL_TAXA
    all_taxa = opportunity_frame(frame, expected_taxa)
    names = dict(zip(all_taxa.inat_taxon_id.astype(int),all_taxa.species.astype(str)))
    if not synthetic and len(names)!=FULL_TAXA:
        raise ValueError("Original 42111 source taxonomy incomplete")
    with count_path.open("r", encoding="utf-8", newline="") as handle:
        reader=csv.DictReader(handle)
        if reader.fieldnames!=["inat_taxon_id","species","n_unique_source_photo_ids","snapshot_date"]:
            raise ValueError("Original full-frame source photo-count schema changed")
        counts = {}
        for row in reader:
            tid = integer(row["inat_taxon_id"])
            raw = row["n_unique_source_photo_ids"]
            if (tid is None or tid not in names or tid in counts
                or row["species"]!=names[tid] or not raw.isdecimal()
                or row["snapshot_date"]!=SNAPSHOT):
                raise ValueError("Original 42111 audit ledger contains modified taxonomy or photo counts")
            counts[tid]=int(raw)
    if set(counts)!=set(names):
        raise ValueError("Full 42111 taxon opportunity ledger incomplete")
    audit=json.loads(independent_audit.read_text(encoding="utf-8"))
    if (not synthetic and (
        audit.get("schema")!="fcp_42111_original_photo_opportunity_audit_v1"
        or audit.get("status")!="VERIFIED_COMPLETE_COUNTS_ONLY_NOT_COLOUR_MEASUREMENT"
        or audit.get("source_frame_sha256")!=FRAME_SHA
        or audit.get("count_ledger_sha256")!=COUNT_SHA
        or audit.get("all_source_photo_opportunity_ids")!=FULL_PHOTO_OPPORTUNITIES
        or audit.get("original_source_taxon_IDs_and_names_independently_matched") is not True
        or audit.get("all_original_image_pixels_measured") is not False)):
        raise ValueError("Original 42111 taxon census not independently audited")
    receipt=json.loads(roster_receipt.read_text(encoding="utf-8"))
    if (receipt.get("schema")!="fcp_rgfca_photo_polymorphism_balanced_sample_v1"
        or receipt.get("status")!="ORIGINAL_SOURCE_BALANCED_REPETITION_PLAN_NO_COLOUR_MEASURED"
        or receipt.get("selected_species_csv_sha256")!=sha256_file(roster_path)
        or receipt.get("new_flower_images_opened") is not False
        or receipt.get("old_149900_flower_colours_used") is not False
        or receipt.get("photographic_poly_probabilities_measured") is not False
        or (not synthetic and (
            receipt.get("selected_species")!=ENROLLED
            or receipt.get("original_source_count_sha256")!=COUNT_SHA
            or receipt.get("original_frame_sha256")!=FRAME_SHA))):
        raise ValueError("RGFCA pre-colour 200-realization species inclusion changed")
    selected={}
    expected={}
    with roster_path.open("r",encoding="utf-8",newline="") as handle:
        reader=csv.DictReader(handle)
        if reader.fieldnames!=["inat_taxon_id","species","source_photo_opportunities",
                                "opportunity_stratum","enrolment_hash"]:
            raise ValueError("RGFCA 1000 species frozen schema changed")
        for row in reader:
            tid=integer(row["inat_taxon_id"])
            raw=row["source_photo_opportunities"]
            if (tid is None or tid in selected or tid not in names
                or row["species"]!=names[tid] or not raw.isdecimal()
                or int(raw)!=counts[tid] or (not synthetic and int(raw)<128)):
                raise ValueError("Enrolled species not in independent original source photo counts")
            selected[tid]=names[tid]
            expected[tid]=int(raw)
    if (len(selected)!=(len(selected) if synthetic else ENROLLED)
        or len(selected)!=receipt["selected_species"]
        or (not synthetic and sum(expected.values())!=RGFCA_ENROLLED_PHOTOS)):
        raise ValueError("RGFCA 1000-species selected source photo population drifted")
    return names,expected


def audit_tables(source_pin: Path, table_files: dict[str,Path], *,
                 synthetic: bool=False) -> dict:
    hashes = {kind:sha256_file(path) for kind,path in table_files.items()}
    if set(hashes)!=set(ORIGINAL_TABLE_HASHES) or (not synthetic and hashes!=ORIGINAL_TABLE_HASHES):
        raise ValueError("Source four original table SHA256 mismatch")
    receipt=verify_source_archive_pin(source_pin,SNAPSHOT,hashes)
    if not synthetic and receipt["source_archive_sha256"]!=ARCHIVE_SHA:
        raise ValueError("Fixed 2026-09-27 original archive SHA changed")
    return {**receipt, "source_table_sha256":hashes}


def process(frame: Path, count_path: Path, independent_audit: Path,
            roster_path: Path, roster_receipt: Path,
            source_pin: Path, tables: dict[str,Path], outdir: Path,
            *, synthetic: bool=False) -> dict:
    if outdir.exists() and any(outdir.iterdir()):
        raise ValueError("Never overwrite an earlier source frozen RGFCA photo-ID ledger")
    names,expected=load_cohort(frame,count_path,independent_audit,
                               roster_path,roster_receipt,synthetic=synthetic)
    fixed=audit_tables(source_pin,tables,synthetic=synthetic)
    outdir.mkdir(parents=True,exist_ok=True)
    sqlite_path=outdir/"_selected_rgfca_source_photo_IDs.sqlite"
    trace=outdir/"scanning_real_20260927_rgfca_source_progress.json"
    progress={"schema":"fcp_targeted_1000_dated_source_progress_v1",
              "status":"IN_PROGRESS_NOT_ADMITTED_FOR_ROI_OR_FCP",
              "source_snapshot_date":SNAPSHOT,
              "selected_species":len(expected),
              "original_photo_pixels_opened":False,
              "old_flower_colours_used":False}
    def record(phase,scanned):
        progress.update(phase=phase, rows_scanned=scanned,
                        updated_utc=datetime.now(timezone.utc).isoformat())
        pending=trace.with_suffix(".tmp")
        pending.write_text(json.dumps(progress,indent=2,sort_keys=True)+"\n")
        os.replace(pending,trace)

    # Use the SAME full-42111 descendant crosswalk as the complete original
    # source importer. Restricting ancestry to 1,000 before this operation
    # would misassign descendants whose most-specific source species was
    # not enrolled.
    crosswalk=source_crosswalk(tables["taxa"],set(names))
    db=sqlite3.connect(sqlite_path)
    db.execute("PRAGMA journal_mode=DELETE")
    db.execute("PRAGMA synchronous=NORMAL")
    db.execute("PRAGMA cache_size=-100000")
    try:
        db.execute("""CREATE TABLE observations(
           uuid TEXT PRIMARY KEY, taxon INTEGER NOT NULL, observer TEXT NOT NULL,
           observed_on TEXT, lat REAL, lon REAL, accuracy TEXT)""")
        db.execute("""CREATE TABLE photos(
           photo_id TEXT PRIMARY KEY, uuid TEXT NOT NULL,
           license TEXT NOT NULL, extension TEXT NOT NULL,
           photographer TEXT, taxon INTEGER NOT NULL, ambiguous INTEGER NOT NULL DEFAULT 0)""")
        db.execute("""CREATE TEMP TABLE candidate(
           photo_id TEXT, uuid TEXT, license TEXT, extension TEXT, photographer TEXT)""")
        eligible_obs=0
        scanned_obs=0
        batch=[]
        required={"observation_uuid","taxon_id","observer_id","latitude",
                  "longitude","positional_accuracy","quality_grade","observed_on"}
        record("reading_original_282m_observations_for_selected_taxa",0)
        for raw in read_tsv(tables["observations"],required):
            scanned_obs+=1
            if str(raw["quality_grade"]).strip().casefold() not in {"research grade","research"}:
                continue
            uuid=str(raw["observation_uuid"]).strip()
            original_tax=integer(raw["taxon_id"])
            tax=crosswalk.get(original_tax, original_tax if original_tax in names else None)
            geo=allowed_geo(raw["latitude"],raw["longitude"])
            observer=integer(raw["observer_id"])
            if not uuid or tax not in expected or geo is None or observer is None:
                continue
            batch.append((uuid,tax,str(observer),raw["observed_on"],*geo,
                          str(raw["positional_accuracy"] or "")))
            if len(batch)>=10000:
                db.executemany("INSERT INTO observations VALUES (?,?,?,?,?,?,?)",batch)
                db.commit()
                eligible_obs+=len(batch)
                batch.clear()
            if scanned_obs%20_000_000==0:
                record("reading_original_observations",scanned_obs)
        if batch:
            db.executemany("INSERT INTO observations VALUES (?,?,?,?,?,?,?)",batch)
            db.commit()
            eligible_obs+=len(batch)
        del batch
        record("all_original_observations_scanned_only_1000_stored",scanned_obs)

        pass_candidates=0
        source_photo_rows=0
        batch=[]
        required={"photo_id","observation_uuid","observer_id","extension","license"}
        def flush_photo():
            nonlocal pass_candidates
            if not batch:return
            db.executemany("INSERT INTO candidate VALUES (?,?,?,?,?)",batch)
            pass_candidates+=db.execute("""SELECT COUNT(*)
                FROM candidate c JOIN observations o ON c.uuid=o.uuid""").fetchone()[0]
            db.execute("""INSERT INTO photos(photo_id,uuid,license,extension,photographer,taxon)
                SELECT c.photo_id,c.uuid,c.license,c.extension,c.photographer,o.taxon
                FROM candidate c JOIN observations o ON c.uuid=o.uuid WHERE 1
                ON CONFLICT(photo_id) DO UPDATE SET
                   ambiguous=CASE WHEN photos.uuid!=excluded.uuid
                              OR photos.taxon!=excluded.taxon
                              OR photos.license!=excluded.license
                              OR photos.extension!=excluded.extension
                              OR photos.photographer!=excluded.photographer
                           THEN 1 ELSE photos.ambiguous END""")
            db.execute("DELETE FROM candidate")
            db.commit()
            batch.clear()
        record("scanning_original_501m_photos_for_selected_observations",0)
        for raw in read_tsv(tables["photos"],required):
            source_photo_rows+=1
            pid=integer(raw["photo_id"])
            ext=accepted_extension(raw["extension"])
            lic=normalize_license(raw["license"])
            uuid=str(raw["observation_uuid"]).strip()
            if pid is not None and ext is not None and lic is not None and uuid:
                batch.append((str(pid),uuid,lic,ext,str(raw.get("observer_id") or "")))
            if len(batch)>=10000:flush_photo()
            if source_photo_rows%20_000_000==0:
                record("scanning_original_photos_for_selected_observation",source_photo_rows)
        flush_photo()
        del batch
        record("all_original_photos_scanned_for_selected_observations",source_photo_rows)

        # The provider full-frame importer resolves photo-ID collisions across
        # 42,111 species. We cannot pretend within-only dedup equals its
        # globally audited result. Check each selected ID against ALL raw photo
        # records in a SECOND full scan; conservatively mark any valid
        # source-ID reuse with different raw metadata, even when the nonselected
        # observation is not source-eligible. Exact independently audited
        # per-species counts are therefore a mandatory admission gate below.
        candidate_ids={pid for (pid,) in db.execute("SELECT photo_id FROM photos")}
        ambiguous_before=db.execute(
            "SELECT COUNT(*) FROM photos WHERE ambiguous=1").fetchone()[0]
        replay=0
        scan2=0
        to_check=[]
        def flush_replay():
            if not to_check:return
            db.executemany("""UPDATE photos SET ambiguous=1 WHERE photo_id=?
                AND (uuid!=? OR license!=? OR extension!=? OR photographer!=?)""",
                to_check)
            db.commit()
            to_check.clear()
        record("second_global_raw_photo_ID_collision_safety_pass",0)
        for raw in read_tsv(tables["photos"],required):
            scan2+=1
            pid=integer(raw["photo_id"])
            if pid is None or str(pid) not in candidate_ids:
                continue
            ext=accepted_extension(raw["extension"])
            lic=normalize_license(raw["license"])
            uuid=str(raw["observation_uuid"]).strip()
            if ext is None or lic is None or not uuid:
                continue
            to_check.append((str(pid),uuid,lic,ext,str(raw.get("observer_id") or "")))
            replay+=1
            if len(to_check)>=10000:flush_replay()
            if scan2%20_000_000==0:record("global_raw_photo_collision_pass",scan2)
        flush_replay()
        del to_check,candidate_ids
        record("global_raw_photo_collision_pass_complete",scan2)
        if scan2!=source_photo_rows:
            raise ValueError("Second source photos pass missed raw photo rows")

        counts=Counter({int(k):int(n) for k,n in db.execute(
            "SELECT taxon,COUNT(*) FROM photos WHERE ambiguous=0 GROUP BY taxon")})
        discrepancies={str(t):{"source_expected":n,"targeted_conservative":counts[t]}
                       for t,n in expected.items() if counts[t]!=n}
        if discrepancies:
            report=outdir/"FAILED_selected_source_taxon_count_mismatches.json"
            report.write_text(json.dumps({
                "schema":"fcp_targeted_1000_dated_source_discrepancy_v1",
                "status":"REJECTED_NOT_ORIGINAL_PHOTO_ID_FCP",
                "n_mismatched_selected_taxa":len(discrepancies),
                "sample_mismatches":dict(list(sorted(discrepancies.items()))[:30]),
                "explanation":"Targeted globally conservative source ID set must be identical to the independently audited 42111 full frame. Count mismatch forbids ROI admission.",
                "original_image_pixels_opened":False,
            },indent=2,sort_keys=True)+"\n")
            raise ValueError(
                f"Original dated RGFCA photo IDs did not reconcile to the globally audited full universe; mismatching taxa={len(discrepancies)}")
        if not synthetic and sum(counts.values())!=RGFCA_ENROLLED_PHOTOS:
            raise ValueError("Fixed 1000 source original photo ID total mismatch")

        db.execute("""CREATE INDEX selected_rgfca_taxon_photo_numeric_order
            ON photos(taxon,CAST(photo_id AS INTEGER)) WHERE ambiguous=0""")
        db.commit()
        manifest=outdir/"rgfca_1000_complete_original_source_photo_IDs.csv"
        sql="""SELECT p.taxon,p.photo_id,o.uuid,o.observer,o.lat,o.lon,
                      o.observed_on,o.accuracy,p.license,p.extension
               FROM photos p INDEXED BY selected_rgfca_taxon_photo_numeric_order
               JOIN observations o ON p.uuid=o.uuid
               WHERE p.ambiguous=0
               ORDER BY p.taxon,CAST(p.photo_id AS INTEGER)"""
        plan=[r[3] for r in db.execute("EXPLAIN QUERY PLAN "+sql)]
        if any("USE TEMP B-TREE FOR ORDER BY" in z for z in plan):
            raise ValueError("Target source photo manifestation tries to sort into unbounded temp disk")
        written=0
        with manifest.open("w",newline="",encoding="utf-8") as file:
            writer=csv.writer(file,lineterminator="\n")
            writer.writerow(CSV_FIELDS)
            for t,pid,uuid,observer,lat,lon,day,accuracy,lic,ext in db.execute(sql):
                writer.writerow((
                    t,names[t],pid,uuid,observer,lat,lon,day,accuracy,
                    lic,f"https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/large.{ext}",
                    SNAPSHOT))
                written+=1
        if written!=sum(expected.values()):
            raise ValueError("Immutable original source photo manifest lost selected photo IDs")
        counts_path=outdir/"rgfca_1000_original_photo_ID_counts.csv"
        with counts_path.open("w",newline="",encoding="utf-8") as file:
            writer=csv.writer(file,lineterminator="\n")
            writer.writerow(["inat_taxon_id","species","n_original_source_photo_IDs"])
            for t in sorted(expected):writer.writerow((t,names[t],counts[t]))
        doc={
            "schema":"fcp_rgfca_targeted_original_20260927_photo_ID_census_v1",
            "status":"COMPLETE_TARGETED_1000_MONTHLY_SOURCE_IDS_COUNT_MATCHED_NO_PHOTO_PIXELS",
            "source_snapshot_date":SNAPSHOT,
            "full_original_frame_taxa":len(names),
            "enrolled_outcome_blind_taxa":len(expected),
            "full_42111_photo_opportunity_counts_sha256":sha256_file(count_path),
            "frozen_1000_source_enrolment_sha256":sha256_file(roster_path),
            "independent_42111_audit_sha256":sha256_file(independent_audit),
            **fixed,
            "n_original_observations_scanned":scanned_obs,
            "n_selected_research_grade_geo_observations":eligible_obs,
            "n_original_photo_rows_scanned_twice":source_photo_rows,
            "n_source_selected_photo_metadata_rows_matching_target_obs":pass_candidates,
            "n_conservative_within_selected_photo_ID_conflicts_before_global_replay":ambiguous_before,
            "n_global_raw_photo_records_checked_against_selected_IDs":replay,
            "n_selected_original_photo_IDs":written,
            "all_enrolled_species_match_independent_full_source_photo_ID_counts":True,
            "selected_all_photo_IDs_csv_sha256":sha256_file(manifest),
            "selected_taxon_counts_csv_sha256":sha256_file(counts_path),
            "photo_ID_selection_colour_blind":True,
            "original_image_pixels_opened":False,
            "historical_flower_colour_CSV_imported":False,
            "observational_plant_photo_opportunity_is_not_flower_morph":True,
            "not_full_93million_ID_census_for_all_42111_taxa":True,
            "photos_NOT_classified_or_species_FCP_inferred":True,
            "subsequent_ROI_screen_additional_observer_site_independence_needed":True,
        }
        report=outdir/"rgfca_1000_targeted_monthly_photo_ID_receipt.json"
        report.write_text(json.dumps(doc,sort_keys=True,indent=2)+"\n")
        progress.update(status="COMPLETE_TARGETED_SOURCE_PHOTO_ID_LEDGER_ONLY",
                        full_ID_ledger_sha256=doc["selected_all_photo_IDs_csv_sha256"])
        record("targeted_original_photo_ID_CSV_frozen_no_pixels",written)
        return doc
    finally:
        db.close()


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    for flag in ("source-frame","source-counts","independent-audit",
                 "rgfca-cohort","rgfca-receipt","dated-archive-pin",
                 "source-taxa-gz","source-observations-gz","source-observers-gz",
                 "source-photos-gz","outdir"):
        p.add_argument("--"+flag,required=True,type=Path)
    a=p.parse_args(argv)
    files={"taxa":a.source_taxa_gz,"observations":a.source_observations_gz,
           "observers":a.source_observers_gz,"photos":a.source_photos_gz}
    summary=process(a.source_frame,a.source_counts,a.independent_audit,
                    a.rgfca_cohort,a.rgfca_receipt,a.dated_archive_pin,
                    files,a.outdir)
    print("FCP_RGFCA_TARGETED_EXACT_20260927_PHOTO_ID_CENSUS_NO_PIXELS",
          json.dumps({k:summary[k] for k in (
              "enrolled_outcome_blind_taxa","n_selected_original_photo_IDs",
              "all_enrolled_species_match_independent_full_source_photo_ID_counts",
              "original_image_pixels_opened")},sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())

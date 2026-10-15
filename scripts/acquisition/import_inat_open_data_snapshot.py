#!/usr/bin/env python3
"""Make a COMPLETE outcome-blind FCP photo ledger from a frozen iNaturalist Open Data snapshot.

Provider-recommended BULK alternative to crawling tens of thousands of taxa
through the observations API. No historical flower-colour CSV, image pixels,
pigment labels, environment rasters or statistical results are read.

Input: the official tab-delimited gz metadata files taxa, observations, observers, photos.
Output: the same per-species complete manifest/checkpoint contract consumed by
`fcp_pipeline.fcp_fresh_study.freeze`, retaining *all* licensed photographs.

For large sources, run on a disk-backed workstation, not ephemeral CI storage.
This SQLite streaming fallback trades execution speed for strict memory bounds.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import re
import sqlite3
import shutil
from pathlib import Path

from fcp_pipeline.fcp_allphotos import opportunity_frame, opportunity_schedule, sha256_file

LICENSES = {"cc0", "cc-by", "cc-by-sa", "cc-by-nc", "cc-by-nc-sa"}
SNAPSHOT_SOURCE = "inaturalist_open_dataset_monthly_snapshot"
AUTHORITATIVE_42111_SHA256 = "5d871bd1f1190c68fd513bb527e6c93de1606dc35ff8cfff5f6187850dc382cc"
FIELDS = ("inat_taxon_id", "species", "observation_id", "photo_id", "photo_url_large",
          "photo_license", "attribution", "observer_id", "observed_on", "latitude",
          "longitude", "positional_accuracy_m", "source_snapshot_date", "source_query")


def read_tsv(path: Path, required: set[str]):
    with gzip.open(path, "rt", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE)
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path.name}: provider snapshot schema missing {sorted(missing)}")
        yield from reader


def integer(x):
    try:
        v = int(x)
        return v if v > 0 else None
    except (ValueError, TypeError):
        return None


def allowed_geo(lat, lon):
    try:
        a, b = float(lat), float(lon)
    except (ValueError, TypeError):
        return None
    if not math.isfinite(a) or not math.isfinite(b) or abs(a) > 90 or abs(b) > 180:
        return None
    return a, b


def lineage_source_taxon(taxon: int, ancestry: str, included: set[int]) -> int | None:
    """Choose most deeply nested nominated taxon, not an outcome-based choice."""
    raw = [*re.split(r"[/\\]", str(ancestry or "")), str(taxon)]
    hits = [int(x) for x in raw if x.isdigit() and int(x) in included]
    return hits[-1] if hits else None


def source_crosswalk(taxa: Path, included: set[int]) -> dict[int, int]:
    cross = {}
    for rec in read_tsv(taxa, {"taxon_id", "ancestry"}):
        ident = integer(rec["taxon_id"])
        if ident is None:
            continue
        parent = lineage_source_taxon(ident, rec["ancestry"], included)
        if parent is not None:
            if ident in cross and cross[ident] != parent:
                raise ValueError("Provider taxon appears in two different taxonomic lineages")
            cross[ident] = parent
    return cross


def normalize_license(raw: str) -> str | None:
    value = str(raw or "").strip().casefold().replace("_", "-").replace(" ", "-")
    return value if value in LICENSES else None


def accepted_extension(raw: str) -> str | None:
    # S3 object keys are CASE-SENSITIVE. Preserve exactly the provider's extension
    # spelling (e.g. "JPG") rather than silently constructing a broken ".jpg" URL.
    ext = str(raw or "").strip().lstrip(".")
    return ext if ext.casefold() in {"jpg", "jpeg", "png"} else None


def _digest(path: Path):
    return sha256_file(path)



def verify_source_archive_pin(path: Path, snapshot_date: str, hashes: dict[str, str]) -> dict:
    """Pin the original *dated* 36-GB S3 tar before the fresh photo study."""
    receipt = json.loads(path.read_text(encoding="utf-8"))
    if (receipt.get("schema") != "fcp_dated_inaturalist_metadata_bundle_pin_v1"
            or receipt.get("status") != "FOUR_ORIGINAL_METADATA_TABLES_EXTRACTED_AND_PINNED"
            or receipt.get("snapshot_date") != snapshot_date
            or receipt.get("metadata_files_sha256") != hashes
            or receipt.get("source_images_opened") is not False
            or receipt.get("historical_photo_colours_read") is not False):
        raise ValueError("Frozen monthly archive attestation missing or disagrees with inputs")
    key = str(receipt.get("s3_object_key") or "")
    date_compact = snapshot_date.replace("-", "")
    if key != f"metadata/inaturalist-open-data-{date_compact}.tar.gz":
        raise ValueError("Dated immutable snapshot archive key differs from frozen date")
    bundle_sha = str(receipt.get("full_bundle_sha256") or "")
    if len(bundle_sha) != 64 or any(x not in "0123456789abcdef" for x in bundle_sha):
        raise ValueError("Original full archive SHA256 is missing from source pin")
    return {"source_archive_key": key, "source_archive_sha256": bundle_sha,
            "source_archive_pin_sha256": sha256_file(path)}


def ingest(source_frame: Path, taxa: Path, observations: Path, photos: Path,
           outdir: Path, snapshot_date: str, *, observers: Path | None = None,
           source_archive_pin: Path | None = None,
           shards: int = 128, expected_species: int = 42111,
           input_sha256: dict[str, str] | None = None,
           prune_verified_temporary_tables: bool = False) -> dict:
    from datetime import date
    if date.fromisoformat(snapshot_date).isoformat() != snapshot_date:
        raise ValueError("Snapshot must have an explicit immutable YYYY-MM-DD date")
    if observers is None:
        raise ValueError("Official observer table required for copyright/source attribution")
    if shards < 1:
        raise ValueError("shards must be positive")
    if outdir.exists() and any(outdir.iterdir()):
        raise ValueError("Source output must be empty; cannot mix snapshots/API sources")
    outdir.mkdir(parents=True, exist_ok=True)
    if expected_species == 42111 and sha256_file(source_frame) != AUTHORITATIVE_42111_SHA256:
        raise ValueError("Authoritative 42,111-taxon frame byte SHA256 mismatch")
    frame = opportunity_frame(source_frame, expected_species)
    plan = opportunity_schedule(frame, shards)
    names = dict(zip(frame.inat_taxon_id.astype(int), frame.species.astype(str)))
    partition = dict(zip(plan.inat_taxon_id.astype(int), plan.partition.astype(int)))
    hashes = {label: _digest(p) for label, p in [("taxa", taxa), ("observations", observations), ("observers", observers), ("photos", photos)]}
    if input_sha256 is None or set(input_sha256) != set(hashes) or any(
            input_sha256[k] != hashes[k] for k in hashes):
        raise ValueError("Every official monthly snapshot input must have a matching pre-pinned SHA256")
    if expected_species == 42111 and source_archive_pin is None:
        raise ValueError("Full 42,111 species require original dated source archive attestation")
    attestation = (verify_source_archive_pin(source_archive_pin, snapshot_date, hashes)
                   if source_archive_pin is not None else {})
    if prune_verified_temporary_tables:
        # Delete only the verified, dated *scratch* gz table copies; do not
        # touch caller-owned or symlinked source files.
        if source_archive_pin is None:
            raise ValueError('Scratch prune requires immutable dated archive attestation')
        original=(taxa,observations,observers,photos)
        labels=('taxa','observations','observers','photos')
        if (any(x.is_symlink() or x.parent.resolve()!=source_archive_pin.parent.resolve()
                for x in original) or
            any(x.name!=f'{label}.csv.gz' for label,x in zip(labels,original))):
            raise ValueError('Cannot prune external, symlinked or non-original metadata files')
    crosswalk = source_crosswalk(taxa, set(names))
    path = outdir / "_working_fcp_snapshot.sqlite"
    db = sqlite3.connect(path)
    source_obs = 0
    source_photo = 0
    pass_obs = 0
    pass_photos = 0
    try:
        db.execute("PRAGMA synchronous=NORMAL")
        db.execute("PRAGMA journal_mode=DELETE")
        db.execute("PRAGMA cache_size=-128000")
        db.execute("CREATE TABLE observations (uuid TEXT PRIMARY KEY, taxon INTEGER, observer TEXT, observed_on TEXT, lat REAL, lon REAL, accuracy TEXT)")
        db.execute("CREATE TABLE photos (photo_id TEXT PRIMARY KEY, uuid TEXT, license TEXT, extension TEXT, photographer TEXT, taxon INTEGER NOT NULL, ambiguous INTEGER NOT NULL DEFAULT 0)")
        db.execute("CREATE TEMP TABLE photo_candidates (photo_id TEXT, uuid TEXT, license TEXT, extension TEXT, photographer TEXT)")
        db.execute("CREATE TABLE observer_labels (observer_id TEXT PRIMARY KEY, login TEXT, name TEXT)")
        batch = []
        required = {"observation_uuid", "taxon_id", "observer_id", "latitude", "longitude", "positional_accuracy", "quality_grade", "observed_on"}
        for rec in read_tsv(observations, required):
            source_obs += 1
            if str(rec["quality_grade"]).strip().casefold() not in {"research grade", "research"}:
                continue
            uuid = str(rec["observation_uuid"]).strip()
            tid = integer(rec["taxon_id"])
            tax = crosswalk.get(tid, tid if tid in names else None)
            geo = allowed_geo(rec["latitude"], rec["longitude"])
            obsid = integer(rec["observer_id"])
            if not uuid or tax not in names or geo is None or obsid is None:
                continue
            batch.append((uuid, int(tax), str(obsid), rec["observed_on"], *geo,
                          str(rec["positional_accuracy"] or "")))
            if len(batch) >= 10000:
                db.executemany("INSERT INTO observations VALUES (?,?,?,?,?,?,?)", batch)
                db.commit()
                pass_obs += len(batch)
                batch.clear()
        if batch:
            db.executemany("INSERT INTO observations VALUES (?,?,?,?,?,?,?)", batch)
            pass_obs += len(batch)
            db.commit()
        del batch
        observer_batch=[]
        observer_rows=0
        for rec in read_tsv(observers, {"observer_id", "login", "name"}):
            observer_rows += 1
            ident=integer(rec["observer_id"])
            if ident is None:
                continue
            observer_batch.append((str(ident),str(rec["login"] or ""),str(rec["name"] or "")))
            if len(observer_batch)>=10000:
                db.executemany("INSERT INTO observer_labels VALUES (?,?,?)",observer_batch)
                db.commit()
                observer_batch.clear()
        if observer_batch:
            db.executemany("INSERT INTO observer_labels VALUES (?,?,?)",observer_batch)
            db.commit()
        del observer_batch
        photos_batch = []
        required = {"photo_id", "observation_uuid", "observer_id", "extension", "license"}
        def flush_photo_batch():
            nonlocal pass_photos
            if not photos_batch:
                return
            db.executemany("INSERT INTO photo_candidates VALUES (?,?,?,?,?)", photos_batch)
            # JOIN in SQLite instead of 1 Python-SQLite round-trip per source photo;
            # 100+ million original photo rows otherwise become impractical.
            pass_photos += db.execute("""SELECT COUNT(*) FROM photo_candidates AS c
                INNER JOIN observations AS o ON c.uuid=o.uuid""").fetchone()[0]
            db.execute("""INSERT INTO photos(photo_id,uuid,license,extension,photographer,taxon)
                SELECT c.photo_id,c.uuid,c.license,c.extension,c.photographer,o.taxon
                FROM photo_candidates AS c JOIN observations AS o ON c.uuid=o.uuid
                WHERE 1
                ON CONFLICT(photo_id) DO UPDATE SET
                   ambiguous=CASE WHEN photos.uuid!=excluded.uuid
                              OR photos.taxon!=excluded.taxon
                              OR photos.license!=excluded.license
                              OR photos.extension!=excluded.extension
                              OR photos.photographer!=excluded.photographer
                          THEN 1 ELSE photos.ambiguous END
                """)
            db.execute("DELETE FROM photo_candidates")
            db.commit()
            photos_batch.clear()
        for rec in read_tsv(photos, required):
            source_photo += 1
            pid = integer(rec["photo_id"])
            ext = accepted_extension(rec["extension"])
            lic = normalize_license(rec["license"])
            if pid is None or ext is None or lic is None:
                continue
            uuid = str(rec["observation_uuid"]).strip()
            if not uuid:
                continue
            photos_batch.append((str(pid), uuid, lic, ext, str(rec.get("observer_id") or "")))
            if len(photos_batch) >= 10000:
                flush_photo_batch()
        flush_photo_batch()
        del photos_batch
        # All raw source rows and hashes have been read. Preserve original
        # bundle and all four table digests in source_archive_pin.json, while
        # reclaiming disposable compressed scratch copies BEFORE the large
        # per-taxon JSON manifest and blind-worker expansion.
        pruned_table_bytes=0
        if prune_verified_temporary_tables:
            for p in (taxa,observations,observers,photos):
                pruned_table_bytes+=p.stat().st_size
                p.unlink()
            print('FCP_VERIFIED_SCRATCH_SOURCE_TABLES_RELEASED',
                  json.dumps({'freed_bytes':pruned_table_bytes,
                              'tables_sha256':hashes},sort_keys=True),flush=True)
        ambiguous_photos=int(db.execute("SELECT COUNT(*) FROM photos WHERE ambiguous=1").fetchone()[0])
        unique_valid_photos=int(db.execute("SELECT COUNT(*) FROM photos WHERE ambiguous=0").fetchone()[0])
        # Frozen census as complete source manifests, even for zero-photo taxa.
        contracts = {}
        for shard in range(shards):
            directory = outdir / f"shard_{shard:03}_of_{shards:03}"
            directory.mkdir(parents=True, exist_ok=True)
            counts = int((plan.partition == shard).sum())
            contract = {"schema": "fcp_uncapped_shard_source_contract_v1",
                        "source_species_frame_sha256": _digest(source_frame),
                        "expected_source_species": len(frame), "n_source_species_in_shard": counts,
                        "shard": shard, "nshards": shards, "created_d2": snapshot_date,
                        "source_kind": SNAPSHOT_SOURCE,
                        "created_d2_field_semantics": "immutable_source_snapshot_date_not_API_created_d2_query",
                        "input_metadata_sha256": hashes,
                        **attestation,
                        "photo_cap": None, "outcome_blind": True}
            (directory / "source_contract.json").write_text(json.dumps(contract, indent=2) + "\n",encoding="utf-8")
            contracts[shard] = directory
        current_tax = None
        writer = None
        counter = {}
        attributable=0
        unresolved=0
        # Materialise photo.taxon during the streaming UUID join. The previous
        # cross-table ORDER BY observations.taxon required an enormous temporary
        # SQLite B-tree *while* writing the uncapped manifests, exhausting the
        # runner after the 36GB download had already succeeded. A partial
        # covering order index lets SQLite stream in taxon/photo-ID order.
        # Build it only after all source rows have been consumed (and optional
        # SHA-pinned disposable input tables released).
        db.execute("""CREATE INDEX photos_taxon_numeric_order
            ON photos(taxon, CAST(photo_id AS INTEGER))
            WHERE ambiguous=0""")
        db.commit()
        manifest_sql = """SELECT photos.taxon,observations.uuid,observations.observer,
            observations.observed_on,observations.lat,observations.lon,observations.accuracy,
            photos.photo_id,photos.license,photos.extension,photos.photographer,
            observer_labels.login,observer_labels.name
            FROM photos INDEXED BY photos_taxon_numeric_order
            JOIN observations ON photos.uuid=observations.uuid
            LEFT JOIN observer_labels ON photos.photographer=observer_labels.observer_id
            WHERE photos.ambiguous=0
            ORDER BY photos.taxon,CAST(photos.photo_id AS INTEGER)"""
        plan_details = [str(p[3]) for p in db.execute("EXPLAIN QUERY PLAN "+manifest_sql)]
        if any("USE TEMP B-TREE FOR ORDER BY" in p.upper() for p in plan_details):
            raise RuntimeError("FCP manifest query would spill a full-source ORDER BY onto scratch disk")
        if not any("photos_taxon_numeric_order" in p for p in plan_details):
            raise RuntimeError("FCP photo manifest order index was not used")
        print("FCP_MANIFEST_INDEXED_STREAM_PLAN",json.dumps(plan_details),flush=True)
        rows = db.execute(manifest_sql)
        emitted_photo_rows = 0
        for (tax, uuid, observer, observed_on,lat,lon,accuracy,pid,lic,ext,photographer,login,person_name) in rows:
            if tax != current_tax:
                if writer is not None: writer.close()
                current_tax = int(tax)
                directory = contracts[partition[current_tax]] / "species" / str(current_tax)
                directory.mkdir(parents=True, exist_ok=True)
                writer = (directory / "photo_manifest.ndjson").open("w",encoding="utf-8")
            # Attribution uses the PHOTO creator, who may differ from the observation owner.
            photo_author=integer(photographer)
            who=str(person_name or "").strip() or str(login or "").strip()
            if who:
                attributable+=1
            else:
                unresolved+=1
                who=f"iNaturalist photo creator {photo_author}" if photo_author else "iNaturalist photographer not resolved"
            profile=f"https://www.inaturalist.org/users/{photo_author}" if photo_author else "https://www.inaturalist.org/"
            attribution=(f"{who} ({profile}), CC0 public domain" if lic=="cc0"
                         else f"© {who} ({profile}), {lic.upper()}")
            obj = dict(zip(FIELDS,(int(tax),names[tax],str(uuid),str(pid),
                       f"https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/large.{ext}",
                       lic,attribution,observer,str(observed_on or ""),lat,lon,
                       accuracy,snapshot_date,SNAPSHOT_SOURCE)))
            writer.write(json.dumps(obj,sort_keys=True,ensure_ascii=False)+"\n")
            counter[tax] = counter.get(tax,0)+1
            emitted_photo_rows += 1
            # Fail closed with space for a diagnostic artifact, rather than
            # allowing a disk-full OSError after hours of manifest emission.
            if expected_species == 42111 and emitted_photo_rows % 100000 == 0:
                free_bytes = shutil.disk_usage(outdir).free
                if free_bytes < 6_000_000_000:
                    raise RuntimeError(
                        f"FCP_MANIFEST_DISK_RESERVE_EXHAUSTED: free={free_bytes}; "
                        f"emitted={emitted_photo_rows}, total={unique_valid_photos}; "
                        "resume only on a provisioned persistent-volume runner")
        if writer is not None: writer.close()
        if emitted_photo_rows != unique_valid_photos:
            raise ValueError(
                f"Source-photo ID loss/duplication at manifest materialization: "
                f"{emitted_photo_rows} emitted != {unique_valid_photos} unique valid")
        for tid, name in names.items():
            directory = contracts[partition[tid]] / "species" / str(tid)
            directory.mkdir(parents=True, exist_ok=True)
            manifest = directory / "photo_manifest.ndjson"
            manifest.touch(exist_ok=True)
            rec={"schema":"fcp_uncapped_metadata_census_v1", "taxon_id":tid,
                 "species":name,"created_d2":snapshot_date,
                 "source_kind":SNAPSHOT_SOURCE, "snapshot_source_sha256": hashes,
                 "cursor":0,"batches":0,"status":"complete", "page_size":None,
                 "photo_cap":None,"colour_pixels_opened":False,
                 "source_photo_rows":counter.get(tid,0),"metadata_sha256":_digest(manifest)}
            (directory/"checkpoint.json").write_text(json.dumps(rec,indent=2)+"\n",encoding="utf-8")
        receipt={"schema":"fcp_open_data_bulk_complete_source_v1",
                 "status":"COMPLETE_SNAPSHOT_METADATA_ALL_NOMINATED_TAXA",
                 "source_kind":SNAPSHOT_SOURCE,"source_snapshot_date":snapshot_date,
                 "outcome_blind":True,"old_photo_measurements_read":False,
                 "image_pixels_opened":False,
                 "original_species_frame_sha256":_digest(source_frame),
                 "taxa_in_initial_frame":len(names),
                 "n_source_observations_scanned":source_obs,
                 "n_target_taxon_research_grade_observations":pass_obs,
                 "n_source_observers_scanned":observer_rows,
                 "n_source_photos_scanned":source_photo,
                 "n_target_taxon_licensed_photo_metadata_rows":pass_photos,
                 "n_target_taxon_licensed_photos":unique_valid_photos,
                 "n_photo_manifest_rows_written":emitted_photo_rows,
                 "manifest_indexed_stream_without_sqlite_order_temp_btree":True,
                 "n_conflicting_photo_ids_quarantined":ambiguous_photos,
                 "n_target_photos_with_named_creator":attributable,
                 "n_target_photos_creator_unresolved":unresolved,
                 "n_taxa_with_source_photos":len(counter),
                 "n_taxa_without_source_photos":len(names)-len(counter),
                 "metadata_snapshots_sha256":hashes,
                 "verified_original_dated_archive":attestation,
                 "photo_cap":None,"all_available_snapshot_photo_ids_included":True,
                 "source_has_immutable_monthly_cutoff":True,
                 "not_literally_live_API_created_d2_window":True,
                 "unavailable_or_unlicensed_photo_not_called_polymorphic":True,
                 "original_S3_extension_case_preserved":True,
                 "conflicting_reused_photo_IDs_quarantined_without_morph_selection":True,
                 "pruned_pinned_temporary_metadata_after_full_scan":bool(prune_verified_temporary_tables),
                 "pruned_temporary_table_bytes":pruned_table_bytes}
        (outdir/"bulk_source_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n",encoding="utf-8")
        return receipt
    finally:
        db.close()
        if path.exists(): path.unlink()


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--species-frame",type=Path,required=True)
    p.add_argument("--taxa",type=Path,required=True)
    p.add_argument("--observations",type=Path,required=True)
    p.add_argument("--observers",type=Path,required=True)
    p.add_argument("--photos",type=Path,required=True)
    p.add_argument("--snapshot-date",required=True)
    p.add_argument("--snapshot-pin",type=Path,help="Full archive SHA and four extracted table SHA receipt; mandatory for full 42111 study")
    p.add_argument("--shards",type=int,default=128)
    p.add_argument("--expected-species",type=int,default=42111)
    p.add_argument("--output-root",type=Path,required=True)
    p.add_argument("--taxa-sha256",required=True)
    p.add_argument("--observations-sha256",required=True)
    p.add_argument("--observers-sha256",required=True)
    p.add_argument("--photos-sha256",required=True)
    p.add_argument("--prune-verified-temporary-tables",action="store_true",
                   help="Only for byte-verified disposable dated source copies, never caller-owned metadata")
    args=p.parse_args(argv)
    print(json.dumps(ingest(args.species_frame,args.taxa,args.observations,args.photos,
                            args.output_root,args.snapshot_date,shards=args.shards,
                            expected_species=args.expected_species,observers=args.observers,
                            source_archive_pin=args.snapshot_pin,
                            input_sha256={"taxa": args.taxa_sha256,
                                          "observations": args.observations_sha256,
                                          "observers": args.observers_sha256,
                                          "photos": args.photos_sha256},
                            prune_verified_temporary_tables=args.prune_verified_temporary_tables),indent=2))
    return 0


if __name__=="__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Make a COMPLETE outcome-blind FCP photo ledger from a frozen iNaturalist Open Data snapshot.

Provider-recommended BULK alternative to crawling tens of thousands of taxa
through the observations API. No historical flower-colour CSV, image pixels,
pigment labels, environment rasters or statistical results are read.

Input: the official tab-delimited gz metadata files taxa, observations, photos.
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
    ext = str(raw or "").strip().lower().lstrip(".")
    return ext if ext in {"jpg", "jpeg", "png"} else None


def _digest(path: Path):
    return sha256_file(path)


def ingest(source_frame: Path, taxa: Path, observations: Path, photos: Path,
           outdir: Path, snapshot_date: str, *, shards: int = 128,
           expected_species: int = 42111, input_sha256: dict[str, str] | None = None) -> dict:
    from datetime import date
    if date.fromisoformat(snapshot_date).isoformat() != snapshot_date:
        raise ValueError("Snapshot must have an explicit immutable YYYY-MM-DD date")
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
    hashes = {label: _digest(p) for label, p in [("taxa", taxa), ("observations", observations), ("photos", photos)]}
    if input_sha256 is None or set(input_sha256) != set(hashes) or any(
            input_sha256[k] != hashes[k] for k in hashes):
        raise ValueError("Every official monthly snapshot input must have a matching pre-pinned SHA256")
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
        db.execute("CREATE TABLE photos (photo_id TEXT PRIMARY KEY, uuid TEXT, license TEXT, extension TEXT, photographer TEXT)")
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
        photos_batch = []
        required = {"photo_id", "observation_uuid", "observer_id", "extension", "license"}
        for rec in read_tsv(photos, required):
            source_photo += 1
            pid = integer(rec["photo_id"])
            ext = accepted_extension(rec["extension"])
            lic = normalize_license(rec["license"])
            if pid is None or ext is None or lic is None:
                continue
            uuid = str(rec["observation_uuid"]).strip()
            if not uuid or db.execute("SELECT 1 FROM observations WHERE uuid=?",(uuid,)).fetchone() is None:
                continue
            photos_batch.append((str(pid), uuid, lic, ext, str(rec.get("observer_id") or "")))
            if len(photos_batch) >= 10000:
                db.executemany("INSERT INTO photos VALUES (?,?,?,?,?)", photos_batch)
                db.commit()
                pass_photos += len(photos_batch)
                photos_batch.clear()
        if photos_batch:
            db.executemany("INSERT INTO photos VALUES (?,?,?,?,?)", photos_batch)
            pass_photos += len(photos_batch)
            db.commit()
        del photos_batch
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
                        "photo_cap": None, "outcome_blind": True}
            (directory / "source_contract.json").write_text(json.dumps(contract, indent=2) + "\n",encoding="utf-8")
            contracts[shard] = directory
        current_tax = None
        writer = None
        counter = {}
        rows = db.execute("SELECT observations.taxon,observations.uuid,observations.observer,observations.observed_on,observations.lat,observations.lon,observations.accuracy, photos.photo_id,photos.license,photos.extension,photos.photographer FROM photos JOIN observations ON photos.uuid=observations.uuid ORDER BY observations.taxon,CAST(photos.photo_id AS INTEGER)")
        for (tax, uuid, observer, observed_on,lat,lon,accuracy,pid,lic,ext,photographer) in rows:
            if tax != current_tax:
                if writer is not None: writer.close()
                current_tax = int(tax)
                directory = contracts[partition[current_tax]] / "species" / str(current_tax)
                directory.mkdir(parents=True, exist_ok=True)
                writer = (directory / "photo_manifest.ndjson").open("w",encoding="utf-8")
            obj = dict(zip(FIELDS,(int(tax),names[tax],str(uuid),str(pid),
                       f"https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/large.{ext}",
                       lic,"",observer,str(observed_on or ""),lat,lon,
                       accuracy,snapshot_date,SNAPSHOT_SOURCE)))
            writer.write(json.dumps(obj,sort_keys=True,ensure_ascii=False)+"\n")
            counter[tax] = counter.get(tax,0)+1
        if writer is not None: writer.close()
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
                 "n_source_photos_scanned":source_photo,
                 "n_target_taxon_licensed_photos":pass_photos,
                 "n_taxa_with_source_photos":len(counter),
                 "n_taxa_without_source_photos":len(names)-len(counter),
                 "metadata_snapshots_sha256":hashes,
                 "photo_cap":None,"all_available_snapshot_photo_ids_included":True,
                 "source_has_immutable_monthly_cutoff":True,
                 "not_literally_live_API_created_d2_window":True,
                 "unavailable_or_unlicensed_photo_not_called_polymorphic":True}
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
    p.add_argument("--photos",type=Path,required=True)
    p.add_argument("--snapshot-date",required=True)
    p.add_argument("--shards",type=int,default=128)
    p.add_argument("--expected-species",type=int,default=42111)
    p.add_argument("--output-root",type=Path,required=True)
    p.add_argument("--taxa-sha256",required=True)
    p.add_argument("--observations-sha256",required=True)
    p.add_argument("--photos-sha256",required=True)
    args=p.parse_args(argv)
    print(json.dumps(ingest(args.species_frame,args.taxa,args.observations,args.photos,
                            args.output_root,args.snapshot_date,shards=args.shards,
                            expected_species=args.expected_species,
                            input_sha256={"taxa": args.taxa_sha256,
                                          "observations": args.observations_sha256,
                                          "photos": args.photos_sha256}),indent=2))
    return 0


if __name__=="__main__":
    raise SystemExit(main())

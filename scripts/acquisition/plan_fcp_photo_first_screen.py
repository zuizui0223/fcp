#!/usr/bin/env python3
"""Outcome-blind, OBSERVER-disjoint FCP image screening from a dated source SQLite.

Read the complete 93,209,351 original source photo IDs ONCE, but select a
bounded first-stage photo screen per taxon BEFORE downloading image pixels.
The 128-photo first stage is a discovery/held-out screen budget, NOT a species
exclusion or cap on subsequent evidence. All unscreened species remain unknown.

Requires the fully verified 2026-09-27 source SQLite stage; cannot operate on
the 42,111-row metadata-only counts, historical photo colour outcomes, any
literature C/S labels, or an arbitrary undocumented image collection.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import heapq
import json
import math
import itertools
import sqlite3
from collections import Counter, defaultdict
from pathlib import Path

from fcp_pipeline.photo_polymorphism_screen import (
    FRAME_SHA, SOURCE_COUNT_SHA, SOURCE_PHOTOS, SNAPSHOT, load_frame,
    observer_split, sha256,
)

STAGE_STATUS = "SOURCE_SQLITE_STAGED_SHARDS_NOT_YET_EMITTED"
STAGE_SCHEMA = "fcp_staged_sqlite_photo_universe_v1"
SELECTION = "FCP_PHOTO_SCREEN_20261016_SNAPSHOT_COLOUR_BLIND_V1"
SELECTED_FIELDS = [
    "inat_taxon_id", "species", "photo_id", "observation_id", "observer_id",
    "latitude", "longitude", "photo_url_large", "photo_license",
    "source_snapshot_date", "screen_fold", "selection_hash",
]
SPECIES_FIELDS = [
    "inat_taxon_id", "species", "source_photo_opportunities",
    "selected_discovery_photos", "selected_holdout_photos",
    "total_selected_for_image_screen", "source_remaining_unmeasured",
    "screen_unmeasured_is_not_monomorphism",
]


def validate_stage(source_dir: Path, frame: Path, *, synthetic: bool = False):
    receipt_path = source_dir / "fcp_staged_sqlite_source_receipt.json"
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    sqlite_path = source_dir / "_working_fcp_snapshot.sqlite"
    if (receipt.get("schema") != STAGE_SCHEMA
        or receipt.get("status") != STAGE_STATUS
        or receipt.get("source_snapshot_date") != SNAPSHOT
        or receipt.get("source_species_frame_sha256") != (sha256(frame) if synthetic else FRAME_SHA)
        or receipt.get("n_target_taxon_licensed_photos") != (
            receipt.get("n_target_taxon_licensed_photos") if synthetic else SOURCE_PHOTOS)
        or receipt.get("full_source_census_complete") is not False
        or receipt.get("full_study_original_images_measured") is not False
        or receipt.get("historical_flower_colour_csv_used") is not False
        or receipt.get("original_photo_ids_deduplicated_globally") is not True
        or receipt.get("verified_original_dated_archive", {}).get("source_archive_sha256") is None
        or not sqlite_path.is_file()
        or sqlite_path.stat().st_size != receipt.get("sqlite_bytes")
        or receipt.get("sqlite_basename") != sqlite_path.name
        or sha256(sqlite_path) != receipt.get("sqlite_sha256")):
        raise ValueError("Original global 93M dated source photo-ID SQLite not fully frozen")
    db = sqlite3.connect(sqlite_path.as_uri() + "?mode=ro&immutable=1", uri=True)
    db.execute("PRAGMA query_only=ON")
    if db.execute("PRAGMA quick_check").fetchone()[0] != "ok":
        db.close()
        raise ValueError("Original source SQLite integrity failed")
    return receipt, db


def _choice_hash(taxon: int, photo_id: str) -> int:
    data = f"{SELECTION}|{taxon}|{photo_id}".encode("utf-8")
    return int(hashlib.sha256(data).hexdigest()[:16], 16)


def _cell(lat: float, lon: float) -> tuple[int, int]:
    return (int(math.floor((lat + 90) / 5)), int(math.floor((lon + 180) / 5)))


def _select(pool: dict, target: int, max_per_observer: int) -> list[tuple]:
    """Deterministic geographic round robin; no morph, climate or literature data."""
    per_cell = {}
    for cell, heap in pool.items():
        per_cell[cell] = sorted((entry[2] for entry in heap),
                                key=lambda c: (c[0], int(c[1])))
    chosen, used_obs, observer_counts = [], set(), Counter()
    cells = sorted(per_cell)
    while len(chosen) < target:
        progress = False
        for cell in cells:
            remaining = per_cell[cell]
            while remaining:
                item = remaining.pop(0)
                score, pid, obs, observer, lat, lon, license_name, ext = item
                if obs in used_obs or observer_counts[observer] >= max_per_observer:
                    continue
                chosen.append(item)
                used_obs.add(obs)
                observer_counts[observer] += 1
                progress = True
                break
            if len(chosen) >= target:
                break
        if not progress:
            break
    return chosen


def plan(source_dir: Path, frame: Path, source_counts: Path, outdir: Path, *,
         photos_per_fold: int = 64, max_per_observer: int = 8,
         cell_candidate_limit: int = 128, _synthetic: bool = False,
         rgfca_enrolment: Path | None = None,
         rgfca_receipt: Path | None = None) -> dict:
    if (photos_per_fold < 1 or photos_per_fold > 512 or
        max_per_observer < 1 or max_per_observer > photos_per_fold or
        cell_candidate_limit < photos_per_fold or cell_candidate_limit > 2048):
        raise ValueError("Invalid first-stage image SCREEN resource limit")
    if outdir.exists() and any(outdir.iterdir()):
        raise ValueError("Never overwrite an existing pre-colour photo screen")
    taxa = load_frame(frame, source_counts, synthetic=_synthetic)
    if (rgfca_enrolment is None) != (rgfca_receipt is None):
        raise ValueError("RGFCA species pool and independently frozen schedule receipt are inseparable")
    enrolled = None
    enrolment_sha = None
    if rgfca_enrolment is not None:
        doc = json.loads(rgfca_receipt.read_text(encoding="utf-8"))
        enrolment_sha = sha256(rgfca_enrolment)
        if (doc.get("schema") != "fcp_rgfca_photo_polymorphism_balanced_sample_v1"
            or doc.get("status") != "ORIGINAL_SOURCE_BALANCED_REPETITION_PLAN_NO_COLOUR_MEASURED"
            or doc.get("selected_species_csv_sha256") != enrolment_sha
            or doc.get("old_149900_flower_colours_used") is not False
            or doc.get("new_flower_images_opened") is not False
            or doc.get("photographic_poly_probabilities_measured") is not False
            or (not _synthetic and (
                doc.get("selected_species") != 1000
                or doc.get("source_frame_species") != 42111
                or doc.get("original_source_count_sha256") != SOURCE_COUNT_SHA
                or doc.get("original_frame_sha256") != FRAME_SHA
                or doc.get("source_photo_opportunities") != SOURCE_PHOTOS))):
            raise ValueError("RGFCA species subsample must be frozen from the real 42111 source before photo colour opening")
        enrolled = set()
        with rgfca_enrolment.open("r", newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            fields = ["inat_taxon_id", "species", "source_photo_opportunities",
                      "opportunity_stratum", "enrolment_hash"]
            if reader.fieldnames != fields:
                raise ValueError("RGFCA source-only selected taxon file columns changed")
            for row in reader:
                tid = int(row["inat_taxon_id"])
                if (tid in enrolled or tid not in taxa
                    or taxa[tid]["species"] != row["species"]
                    or taxa[tid]["n_source_photo_ids"] != int(row["source_photo_opportunities"])
                    or (not _synthetic and int(row["source_photo_opportunities"]) < 128)):
                    raise ValueError("RGFCA species identity or original source photograph opportunity changed")
                enrolled.add(tid)
        if len(enrolled) != int(doc["selected_species"]):
            raise ValueError("RGFCA enrolment does not include every pre-colour selected taxon")
    stage, db = validate_stage(source_dir, frame, synthetic=_synthetic)
    outdir.mkdir(parents=True, exist_ok=True)
    photo_path = outdir / "photo_screen_original_metadata_selected.csv"
    opportunity_path = outdir / "photo_screen_all_species_opportunities.csv"
    n_scanned = n_selected = n_taxa_with_photos = 0
    count_by_species = {}
    source_counts_by_taxon = Counter()
    # Explicit index use avoids an unbounded temporary source-wide SQLite sort.
    sql = """SELECT photos.taxon, photos.photo_id, photos.uuid, observations.observer,
                    observations.lat, observations.lon, photos.license, photos.extension
             FROM photos INDEXED BY photos_taxon_numeric_order
             JOIN observations ON photos.uuid=observations.uuid
             WHERE photos.ambiguous=0
             ORDER BY photos.taxon, CAST(photos.photo_id AS INTEGER)"""
    if enrolled is not None:
        # Query only 1,000 OUTCOME-BLIND enrolled taxa through the original
        # immutable taxon/photo index. All 93M source IDs were audited earlier.
        # Avoid rescanning all 93M rows or choosing a species based on colour.
        sql = """SELECT photos.taxon, photos.photo_id, photos.uuid, observations.observer,
                        observations.lat, observations.lon, photos.license, photos.extension
                 FROM photos INDEXED BY photos_taxon_numeric_order
                 JOIN observations ON photos.uuid=observations.uuid
                 WHERE photos.ambiguous=0 AND photos.taxon=?
                 ORDER BY CAST(photos.photo_id AS INTEGER)"""
        qp = [x[3] for x in db.execute("EXPLAIN QUERY PLAN " + sql, (min(enrolled),))]
        photo_source = itertools.chain.from_iterable(
            db.execute(sql, (taxon,)) for taxon in sorted(enrolled))
    else:
        qp = [x[3] for x in db.execute("EXPLAIN QUERY PLAN " + sql)]
        photo_source = db.execute(sql)
    if not any("photos_taxon_numeric_order" in x for x in qp):
        raise ValueError("Original globally unique photo-ID index not used")
    if any("USE TEMP B-TREE FOR ORDER BY" in x.upper() for x in qp):
        raise ValueError("Potential source-wide temporary SQLite sort detected")
    cur_taxon = None
    pools = [defaultdict(list), defaultdict(list)]
    photo_cols = SELECTED_FIELDS

    with photo_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=photo_cols, lineterminator="\n")
        writer.writeheader()

        def flush():
            nonlocal n_selected, n_taxa_with_photos, pools
            if cur_taxon is None:
                return
            ident = cur_taxon
            n_source = source_counts_by_taxon[ident]
            if n_source != taxa[ident]["n_source_photo_ids"]:
                raise ValueError(f"Original source photo counts drifted for taxon {ident}")
            n_taxa_with_photos += int(n_source > 0)
            picked = [
                _select(pools[fold], photos_per_fold, max_per_observer)
                for fold in range(2)
            ]
            selected_ids = set()
            for fold, batch in enumerate(picked):
                for score, pid, obs, observer, lat, lon, license_name, extension in batch:
                    if pid in selected_ids:
                        raise ValueError("One source photo ID selected in both independent screen folds")
                    selected_ids.add(pid)
                    writer.writerow({
                        "inat_taxon_id": ident, "species": taxa[ident]["species"],
                        "photo_id": pid, "observation_id": obs,
                        "observer_id": observer, "latitude": lat,
                        "longitude": lon,
                        "photo_url_large": (
                            f"https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/large.{extension}"),
                        "photo_license": license_name,
                        "source_snapshot_date": SNAPSHOT, "screen_fold": fold,
                        "selection_hash": f"{score:016x}",
                    })
            n_chosen = len(selected_ids)
            n_selected += n_chosen
            count_by_species[ident] = (len(picked[0]), len(picked[1]))
            if n_chosen > n_source:
                raise ValueError("Screen included more photos than the source universe")
            pools = [defaultdict(list), defaultdict(list)]

        try:
            for (taxon, pid, obs, observer, lat, lon, license_name, extension) in photo_source:
                taxon = int(taxon)
                if taxon not in taxa:
                    raise ValueError("Unknown source taxon in original SQLite")
                if cur_taxon != taxon:
                    if cur_taxon is not None and taxon < cur_taxon:
                        raise ValueError("Original photos are not in monotonic taxon order")
                    flush()
                    cur_taxon = taxon
                pid = str(pid)
                obs = str(obs)
                observer = str(observer)
                if not pid.isdecimal() or not obs or not observer:
                    raise ValueError("Invalid immutable photo / observation / observer ID")
                lat, lon = float(lat), float(lon)
                if not math.isfinite(lat) or not math.isfinite(lon):
                    raise ValueError("Invalid geographic source coordinates")
                score = _choice_hash(taxon, pid)
                fold = observer_split(taxon, observer)
                cell = _cell(lat, lon)
                heap = pools[fold][cell]
                entry = (-score, -int(pid),
                         (score, pid, obs, observer, lat, lon, license_name, extension))
                if len(heap) < cell_candidate_limit:
                    heapq.heappush(heap, entry)
                elif entry > heap[0]:
                    heapq.heapreplace(heap, entry)
                source_counts_by_taxon[taxon] += 1
                n_scanned += 1
            flush()
        finally:
            db.close()

    checked_taxa = sorted(enrolled) if enrolled is not None else sorted(taxa)
    for taxon in checked_taxa:
        if source_counts_by_taxon[taxon] != taxa[taxon]["n_source_photo_ids"]:
            raise ValueError(f"Original source photo counts differ for taxon {taxon}")
    expected_scanned = (sum(taxa[x]["n_source_photo_ids"] for x in enrolled)
                        if enrolled is not None else stage["n_target_taxon_licensed_photos"])
    if n_scanned != expected_scanned:
        raise ValueError("Selected original source photo IDs were not all checked")
    with opportunity_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=SPECIES_FIELDS, lineterminator="\n")
        writer.writeheader()
        for ident in sorted(taxa):
            n = taxa[ident]["n_source_photo_ids"]
            first, reserve = count_by_species.get(ident, (0, 0))
            writer.writerow({
                "inat_taxon_id": ident, "species": taxa[ident]["species"],
                "source_photo_opportunities": n,
                "selected_discovery_photos": first,
                "selected_holdout_photos": reserve,
                "total_selected_for_image_screen": first + reserve,
                "source_remaining_unmeasured": n - first - reserve,
                "screen_unmeasured_is_not_monomorphism": True,
            })
    receipt = {
        "schema": "fcp_photo_first_42111_outcome_blind_sample_schedule_v1",
        "status": "PHOTO_IDS_SELECTED_FROM_FULL_ORIGINAL_SOURCE_NOT_YET_MEASURED",
        "source_snapshot_date": SNAPSHOT,
        "source_species": len(taxa),
        "original_source_photo_ids_globally_considered": stage["n_target_taxon_licensed_photos"],
        "original_species_with_source_photo": sum(
            v["n_source_photo_ids"] > 0 for v in taxa.values()),
        "original_source_SQLite_photo_rows_examined": n_scanned,
        "RGFCA_enrolment_species_count": len(enrolled) if enrolled is not None else None,
        "RGFCA_enrolment_selected_taxa_sha256": enrolment_sha,
        "RGFCA_enrolment_metadata_only_no_literature_or_colour": enrolled is not None,
        "first_screen_photos_per_observer_fold_not_biological_cap": photos_per_fold,
        "per_observer_max_in_one_bounded_screen": max_per_observer,
        "geographic_cell_degrees": 5,
        "outcome_blind_image_selection_no_floral_colour_opened": True,
        "source_sqlite_sha256": stage["sqlite_sha256"],
        "source_photo_counts_sha256": sha256(source_counts),
        "original_source_frame_sha256": sha256(frame),
        "selected_photo_IDs": n_selected,
        "photo_selection_csv_sha256": sha256(photo_path),
        "species_opportunity_ledger_sha256": sha256(opportunity_path),
        "selected_sample_NOT_full_93million_photo_census": True,
        "selected_sample_photo_pixels_measured": False,
        "photographic_polymorphism_classifications_computed": False,
        "single_photo_state_will_not_be_reported_as_monomorphic": True,
        "literature_polymorphism_labels_used": False,
        "prior_149900_measured_colour_records_used": False,
        "warning": ("A first-stage 64/64 per-observer-fold photo budget is NOT an "
                    "exclusion/absence criterion; unresolved taxa must receive "
                    "predeclared additional independent image sampling before "
                    "any nonpolymorphic claim."),
    }
    if not _synthetic and (
        len(taxa) != 42111
        or stage["n_target_taxon_licensed_photos"] != SOURCE_PHOTOS
        or sum(x["n_source_photo_ids"] > 0 for x in taxa.values()) != 42095
        or n_selected > (len(enrolled) if enrolled is not None else len(taxa))*photos_per_fold*2
    ):
        raise ValueError("Original 42111 source-based RGFCA photo sample did not reconcile")
    (outdir / "photo_screen_source_selection_receipt.json").write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8")
    return receipt


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--staged-source", required=True, type=Path)
    parser.add_argument("--original-frame", required=True, type=Path)
    parser.add_argument("--source-counts", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--photos-per-fold", type=int, default=64)
    parser.add_argument("--max-per-observer", type=int, default=8)
    parser.add_argument("--rgfca-enrolment", type=Path,
                        help="Byte-frozen 1000 species selected independently of flower colours")
    parser.add_argument("--rgfca-receipt", type=Path,
                        help="Immutable original-source 200-realization photo-sampling design receipt")
    args = parser.parse_args(argv)
    got = plan(args.staged_source, args.original_frame, args.source_counts,
               args.outdir, photos_per_fold=args.photos_per_fold,
               max_per_observer=args.max_per_observer,
               rgfca_enrolment=args.rgfca_enrolment,
               rgfca_receipt=args.rgfca_receipt)
    print("FCP_ORIGINAL_PHOTO_SCREEN_METADATA_ONLY", json.dumps(got, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

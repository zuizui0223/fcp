#!/usr/bin/env python3
"""HARD GATE: reconcile a completed 42,111-taxon all-photo-ID source freeze
against the already independently audited, same-date original provider census.

This is intentionally unusable on the count-only receipt alone. It never
opens original flower pixels or historical flower colour outcomes.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

SOURCE_DATE = "2026-09-27"
ARCHIVE_SHA = "bfc1a97c436bf5989fccd43f12392953128489b76b92cea24da335ad73da4d8e"
FRAME_SHA = "5d871bd1f1190c68fd513bb527e6c93de1606dc35ff8cfff5f6187850dc382cc"
REFERENCE_FIELDS = ("inat_taxon_id", "species", "n_unique_source_photo_ids", "snapshot_date")
FULL_FIELDS = ("inat_taxon_id", "species", "n_unique_target_photo_ids",
               "source_status", "source_manifest_sha256")


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def validate(baseline_csv: Path, baseline_audit: Path, complete_source_csv: Path,
             complete_source_receipt: Path, firewall_freeze: Path,
             *, _synthetic_test: bool = False) -> dict:
    population = 4 if _synthetic_test else 42111
    base = json.loads(baseline_audit.read_text(encoding="utf-8"))
    full = json.loads(complete_source_receipt.read_text(encoding="utf-8"))
    freeze = json.loads(firewall_freeze.read_text(encoding="utf-8"))
    if (base.get("schema")!="fcp_42111_original_photo_opportunity_audit_v1"
        or base.get("status")!="VERIFIED_COMPLETE_COUNTS_ONLY_NOT_COLOUR_MEASUREMENT"
        or base.get("source_frame_species")!=population
        or base.get("source_frame_sha256")!=FRAME_SHA
        or base.get("source_archive_sha256")!=ARCHIVE_SHA
        or base.get("snapshot_date")!=SOURCE_DATE
        or base.get("count_ledger_sha256")!=digest(baseline_csv)
        or base.get("original_source_taxon_IDs_and_names_independently_matched") is not True
        or base.get("all_source_photo_id_manifests_frozen") is not False
        or base.get("all_original_image_pixels_measured") is not False
        or base.get("W_D_H2_or_climate_soil_results") is not False):
        raise ValueError("No admitted independently audited original dated full species opportunity baseline")
    if (full.get("status")!="COMPLETE_42111_SOURCE_METADATA_FROZEN_PIXEL_MEASUREMENT_STILL_NOT_DONE"
        or full.get("complete_metadata_census") is not True
        or full.get("target_frame_species")!=population
        or full.get("species_count_checkpoints_audited")!=population
        or full.get("original_archive_sha256")!=ARCHIVE_SHA
        or full.get("immutable_source_date")!=SOURCE_DATE
        or full.get("full_pixel_measurement_complete") is not False
        or full.get("biological_h2_and_environment_inference_ready") is not False
        or full.get("old_measured_flower_colours_used") is not False
        or full.get("no_per_species_photo_cap") is not True
        or full.get("frozen_location_blind_worker_partitions")!=128
        or freeze.get("status")!="FRESH_BLIND_METADATA_FROZEN"
        or freeze.get("source_species")!=population
        or freeze.get("worker_packets")!=128
        or freeze.get("historical_measured_CSV_read") is not False
        or freeze.get("all_available_licensed_photo_ids") !=
             full.get("n_all_available_licensed_original_photo_ids")):
        raise ValueError("The candidate full source is NOT a genuinely complete frozen original photo-ID census")

    baseline: dict[int, tuple[str, int]] = {}
    with baseline_csv.open("r",encoding="utf-8",newline="") as f:
        reader=csv.DictReader(f)
        if tuple(reader.fieldnames or ())!=REFERENCE_FIELDS:
            raise ValueError("Wrong audited source count ledger columns")
        for row in reader:
            tax=int(row["inat_taxon_id"])
            if (tax<=0 or tax in baseline or row["snapshot_date"]!=SOURCE_DATE or
                not str(row["n_unique_source_photo_ids"]).isdecimal()):
                raise ValueError("Corrupted original provider source baseline taxa/photographs")
            baseline[tax]=(row["species"],int(row["n_unique_source_photo_ids"]))
    if len(baseline)!=population or sum(n for _,n in baseline.values())!=base["all_source_photo_opportunity_ids"]:
        raise ValueError("Baseline does not retain the original provider's full photograph counts")

    seen: set[int] = set()
    total=0
    with complete_source_csv.open("r",encoding="utf-8",newline="") as f:
        reader=csv.DictReader(f)
        if tuple(reader.fieldnames or ())!=FULL_FIELDS:
            raise ValueError("Completed full-source photo-ID receipt fields have changed")
        for row in reader:
            tax_txt=row["inat_taxon_id"]
            n_txt=row["n_unique_target_photo_ids"]
            check=str(row["source_manifest_sha256"])
            if (not tax_txt.isdecimal() or not n_txt.isdecimal()
                or row["source_status"]!="complete"
                or len(check)!=64 or any(c not in "0123456789abcdef" for c in check)):
                raise ValueError("Unsealed or malformed per-taxon original photo-ID manifest receipt")
            tax,n=int(tax_txt),int(n_txt)
            if tax in seen or tax not in baseline:
                raise ValueError("Changed or duplicated original provider taxon after full photo freeze")
            seen.add(tax)
            if (row["species"],n)!=baseline[tax]:
                raise ValueError(f"Full original photo-ID count/name differs from independently audited provider universe at taxon {tax}")
            total+=n
    if (seen!=set(baseline) or total!=base["all_source_photo_opportunity_ids"]
        or total!=full["n_all_available_licensed_original_photo_ids"]):
        raise ValueError("All 42,111 original source photos did not reconcile with the frozen source firewall")
    return {
        "schema":"fcp_complete_photo_ID_source_vs_independent_provider_audit_v1",
        "status":"REAL_ALL_42111_ORIGINAL_PHOTO_ID_COUNTS_RECONCILED_NO_PIXELS",
        "snapshot_date":SOURCE_DATE,
        "species":population,
        "distinct_original_photo_ids":total,
        "source_count_ledger_sha256":digest(baseline_csv),
        "fully_frozen_count_ledger_sha256":digest(complete_source_csv),
        "source_status_receipt_sha256":digest(complete_source_receipt),
        "source_firewall_receipt_sha256":digest(firewall_freeze),
        "original_photo_IDs_covered_without_count_loss":True,
        "full_study_original_pixels_classified":False,
        "flower_colour_H2_W_or_environment_results":False,
    }


def main(argv=None) -> int:
    parser=argparse.ArgumentParser(description=__doc__)
    for flag in ("baseline-counts", "independent-audit", "complete-source-counts",
                 "full-source-receipt", "firewall-freeze", "output"):
        parser.add_argument("--"+flag,required=True,type=Path)
    a=parser.parse_args(argv)
    result=validate(a.baseline_counts,a.independent_audit,
                    a.complete_source_counts,a.full_source_receipt,a.firewall_freeze)
    if a.output.exists():
        raise ValueError("Do not overwrite an existing completed source audit")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,ensure_ascii=False,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("FCP_ALL_REAL_SOURCE_PHOTO_ID_COUNTERS_MATCH",json.dumps(result,sort_keys=True))
    return 0


if __name__=="__main__":
    raise SystemExit(main())

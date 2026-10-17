#!/usr/bin/env python3
"""Plan 128 FCP original-image work partitions from independently audited source counts.

This is a *capacity* ledger, NOT a substitute for the immutable all-photo-ID
freeze or any original ROI-v4 pixels. Every source photo remains in scope.
The partition hash exactly matches fcp_pipeline.fcp_allphotos.opportunity_schedule.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

FRAME_SHA = "5d871bd1f1190c68fd513bb527e6c93de1606dc35ff8cfff5f6187850dc382cc"
SOURCE_DATE = "2026-09-27"
FIELDS = ("inat_taxon_id", "species", "n_unique_source_photo_ids", "snapshot_date")
REPORT_STATUS = "VERIFIED_COMPLETE_COUNTS_ONLY_NOT_COLOUR_MEASUREMENT"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as src:
        for chunk in iter(lambda: src.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def partition_for(taxon_id: int, species: str, partitions: int) -> int:
    identity = f"FCP_ALL_OBSERVABLE_20261014|{taxon_id}|{species}"
    return int(hashlib.sha256(identity.encode("utf-8")).hexdigest()[:16], 16) % partitions


def plan(counts_file: Path, independent_audit: Path, *,
         partitions: int = 128, batch_size: int = 64,
         _synthetic_test: bool = False) -> tuple[list[dict], dict]:
    if partitions != 128 and not _synthetic_test:
        raise ValueError("Production schedule MUST preserve the original 128 fixed partitions")
    if not 1 <= partitions <= 128 or not 1 <= batch_size <= 1024:
        raise ValueError("Invalid compute planning dimension (not a biological cap)")
    verification = json.loads(independent_audit.read_text(encoding="utf-8"))
    if (verification.get("schema") != "fcp_42111_original_photo_opportunity_audit_v1"
        or verification.get("status") != REPORT_STATUS
        or verification.get("source_frame_sha256") != FRAME_SHA
        or verification.get("snapshot_date") != SOURCE_DATE
        or verification.get("count_ledger_sha256") != sha256_file(counts_file)
        or verification.get("original_source_taxon_IDs_and_names_independently_matched") is not True
        or verification.get("all_source_photo_id_manifests_frozen") is not False
        or verification.get("all_original_image_pixels_measured") is not False
        or verification.get("W_D_H2_or_climate_soil_results") is not False):
        raise ValueError("Real original photo counts lack valid independent original-taxon audit")
    target = 4 if _synthetic_test else 42111
    if verification.get("source_frame_species") != target:
        raise ValueError("Taxon counts are not the expected original species frame")

    shards = [{
        "partition": k,
        "source_species": 0,
        "source_photo_ids": 0,
        "source_species_ge_40_photos": 0,
        "species_with_no_source_photos": 0,
        "largest_species_source_photo_ids": 0,
    } for k in range(partitions)]
    seen: set[int] = set()
    with counts_file.open("r", encoding="utf-8", newline="") as src:
        reader = csv.DictReader(src)
        if tuple(reader.fieldnames or ()) != FIELDS:
            raise ValueError("Unrecognised original count CSV columns")
        for record in reader:
            taxon = str(record["inat_taxon_id"])
            count = str(record["n_unique_source_photo_ids"])
            species = str(record["species"])
            if not taxon.isdecimal() or not count.isdecimal() or not species.strip() or record["snapshot_date"] != SOURCE_DATE:
                raise ValueError("Malformed original source taxon/photo count")
            ident, n = int(taxon), int(count)
            if ident <= 0 or ident in seen:
                raise ValueError("Duplicate/missing source taxon ID in capacity plan")
            seen.add(ident)
            dest = shards[partition_for(ident, species, partitions)]
            dest["source_species"] += 1
            dest["source_photo_ids"] += n
            dest["source_species_ge_40_photos"] += (n >= 40)
            dest["species_with_no_source_photos"] += (n == 0)
            dest["largest_species_source_photo_ids"] = max(dest["largest_species_source_photo_ids"], n)

    photo_total = sum(s["source_photo_ids"] for s in shards)
    if (len(seen) != target or photo_total != verification.get("all_source_photo_opportunity_ids")
        or sum(s["species_with_no_source_photos"] for s in shards) != verification.get("taxa_with_zero_photos")
        or sum(s["source_species_ge_40_photos"] for s in shards) !=
            verification.get("minimum_source_photo_depth_upper_bound_for_W", {}).get("40")):
        raise ValueError("128-shard original photo planning did not reconcile every counted source photo")
    for s in shards:
        s["estimated_64_photo_work_units"] = math.ceil(s["source_photo_ids"] / batch_size)
    receipt = {
        "schema": "fcp_real_original_source_128_partition_capacity_v1",
        "status": "CAPACITY_ESTIMATE_ONLY_ORIGINAL_PHOTO_IDS_NOT_FROZEN",
        "source_date": SOURCE_DATE,
        "species": target,
        "photos": photo_total,
        "source_photo_count_ledger_sha256": verification["count_ledger_sha256"],
        "source_frame_sha256": FRAME_SHA,
        "independent_audit_sha256": sha256_file(independent_audit),
        "source_photo_partitions": partitions,
        "compute_batch_size_not_photo_cap": batch_size,
        "nominal_photo_work_units": sum(s["estimated_64_photo_work_units"] for s in shards),
        "largest_partition_photo_ids": max(s["source_photo_ids"] for s in shards),
        "smallest_partition_photo_ids": min(s["source_photo_ids"] for s in shards),
        "largest_partition_id": max(shards, key=lambda s: s["source_photo_ids"])["partition"],
        "all_42111_taxa_included_regardless_of_photo_capacity": target == 42111,
        "no_species_or_photo_inclusion_cap": True,
        "source_photo_id_manifests_complete": False,
        "actual_original_ROI_v4_photo_measurements_completed": False,
        "H2_W_999world_or_ecological_results": False,
        "notes": [
            "Every number is provider source PHOTO metadata (not visible flowers or phenotype).",
            "64-photo compute units are planning bins; an actual frozen worker packet may cross species boundaries.",
            "The independent full photo-ID freeze and pixel/morph outcome are still required before biological inference.",
            "Photo availability >=40 is an upper bound on species with >=40 subsequently classifiable original photos."
        ],
    }
    return shards, receipt


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts", required=True, type=Path)
    parser.add_argument("--independent-audit", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args(argv)
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise ValueError("Capacity plan output cannot overwrite an existing source ledger")
    rows, receipt = plan(args.counts, args.independent_audit)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    table = args.output_dir / "source_128_partition_photo_capacity.csv"
    with table.open("w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    receipt["partition_capacity_csv_sha256"] = sha256_file(table)
    summary = args.output_dir / "source_128_partition_capacity_receipt.json"
    summary.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("FCP_42111_PHOTO_PARTITION_CAPACITY_ONLY", json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

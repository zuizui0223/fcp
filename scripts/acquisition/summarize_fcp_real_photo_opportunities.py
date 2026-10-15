#!/usr/bin/env python3
"""Independently audit REAL FCP 42,111-taxon source photo counts (never colours).

Consumes only the 2026-09-27 photo-opportunity CSV + original source receipts.
It cannot claim an all-photo-ID freeze, a flower colour, D, W, or ecology.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

FRAME_SHA = "5d871bd1f1190c68fd513bb527e6c93de1606dc35ff8cfff5f6187850dc382cc"
ARCHIVE_SHA = "bfc1a97c436bf5989fccd43f12392953128489b76b92cea24da335ad73da4d8e"
TABLE_SHA = {
    "taxa": "71c1b2909c3c1c409da71c6f778939b00c5e1cc97a10449b755aa569c700e2d0",
    "observations": "6c8b93cd808bdb0503690a60719a2d0a6478d3494eda7c8623db8c7f454b125b",
    "observers": "735d8112f24c1d1951622bc09fee4a6a07850d87937b3aa77d6c65289e357d63",
    "photos": "0033f67f2362acae0509314c4a7b8429f1f6ff85323cb3c86eeb7cbb83ecc6b6",
}
HEADERS = ["inat_taxon_id", "species", "n_unique_source_photo_ids", "snapshot_date"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def audit(counts_file: Path, source_receipt_file: Path,
          execution_receipt_file: Path, *, species: int = 42111,
          frame_file: Path | None = None) -> dict:
    source = json.loads(source_receipt_file.read_text(encoding="utf-8"))
    execution = json.loads(execution_receipt_file.read_text(encoding="utf-8"))
    if (source.get("schema") != "fcp_full_source_photo_opportunity_counts_v1"
            or source.get("status") != "COMPLETE_REAL_SOURCE_PHOTO_COUNTS_ONLY_NO_PHOTO_ID_FIREWALL"
            or execution.get("status") != "COMPLETE_REAL_42111_PHOTO_OPPORTUNITY_COUNTS_ONLY_NOT_FULL_SOURCE"
            or source.get("source_species") != species
            or execution.get("target_frame_species") != species
            or source.get("source_species_frame_sha256") != FRAME_SHA
            or execution.get("original_archive_sha256") != ARCHIVE_SHA
            or source.get("source_archive_attestation", {}).get("source_archive_sha256") != ARCHIVE_SHA
            or source.get("source_input_sha256") != TABLE_SHA
            or execution.get("all_four_extracted_table_sha256") != TABLE_SHA
            or source.get("snapshot_date") != "2026-09-27"
            or execution.get("immutable_source_date") != "2026-09-27"
            or source.get("counts_cover_every_target_taxon") is not True
            or source.get("photo_cap") is not None
            or execution.get("source_photo_opportunity_counts_completed") is not True
            or source.get("all_original_photo_ID_manifests_frozen") is not False
            or source.get("full_source_census_complete") is not False
            or source.get("full_study_original_images_measured") is not False
            or source.get("historical_flower_colour_csv_used") is not False
            or source.get("source_image_pixels_opened") is not False
            or execution.get("complete_metadata_census") is not False
            or execution.get("complete_original_images_measured") is not False
            or execution.get("new_W_or_environment_results_computed") is not False
            or execution.get("old_measured_flower_colours_used") is not False
            or execution.get("photo_id_manifest_universe_frozen") is not False):
        raise ValueError("No verified REAL full-frame source-photo count receipt, or biological gate drift")

    digest = sha256(counts_file)
    if digest != source.get("counts_file_sha256") or digest != execution.get("photo_counts_csv_sha256"):
        raise ValueError("Source count file fingerprint changed")
    counts: list[int] = []
    names_seen: set[int] = set()
    with counts_file.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames != HEADERS:
            raise ValueError("Unexpected original photo opportunity ledger columns")
        for row in reader:
            tid = str(row["inat_taxon_id"])
            n = str(row["n_unique_source_photo_ids"])
            if (not tid.isdecimal() or not n.isdecimal()
                    or not str(row["species"]).strip()
                    or row["snapshot_date"] != "2026-09-27"):
                raise ValueError("Invalid original photo-count/taxon row")
            key, value = int(tid), int(n)
            if key <= 0 or key in names_seen:
                raise ValueError("Duplicate or invalid source taxon ID")
            names_seen.add(key)
            counts.append(value)

    # A count receipt containing 42,111 plausible IDs is NOT evidence those
    # are the original 42,111 biological taxa. Independently match exact
    # species identities and names against the byte-frozen original frame.
    # Refuse all real full-frame reports if that original frame is unavailable.
    original_frame_matched = False
    if species == 42111 and frame_file is None:
        raise ValueError("Real 42111 audit needs independently verified original taxon frame")
    if frame_file is not None:
        if species == 42111 and sha256(frame_file) != FRAME_SHA:
            raise ValueError("Original 42111 taxon frame SHA256 mismatch")
        expected_names = {}
        with frame_file.open("r", encoding="utf-8", newline="") as f:
            original = csv.DictReader(f)
            if not {"inat_taxon_id", "species"}.issubset(set(original.fieldnames or ())):
                raise ValueError("Original source frame has no taxon ID/species fields")
            for row in original:
                ident = str(row["inat_taxon_id"])
                if not ident.isdecimal() or int(ident) <= 0 or int(ident) in expected_names:
                    raise ValueError("Original source frame contains duplicate or invalid taxon IDs")
                expected_names[int(ident)] = str(row["species"])
        if len(expected_names) != species:
            raise ValueError("Original source frame cardinality mismatches expected 42111")
        # Re-read COUNT ledger rather than trusting only aggregate set/count.
        matched_count = 0
        with counts_file.open("r", encoding="utf-8", newline="") as f:
            for row in csv.DictReader(f):
                ident = int(row["inat_taxon_id"])
                if ident not in expected_names or row["species"] != expected_names[ident]:
                    raise ValueError("Original 42111 source taxon ID or taxonomic name mismatch")
                matched_count += 1
        if matched_count != species or names_seen != set(expected_names):
            raise ValueError("Original 42111 source taxon ID coverage mismatch")
        original_frame_matched = True

    n_taxa = len(counts)
    n_with = sum(n > 0 for n in counts)
    n_zero = n_taxa - n_with
    total = sum(counts)
    if (n_taxa != species
            or n_with != source.get("n_taxa_with_source_photos")
            or n_zero != source.get("n_taxa_without_source_photos")
            or n_with != execution.get("n_taxa_with_source_photos")
            or n_zero != execution.get("n_taxa_without_source_photos")
            or total != source.get("n_target_taxon_licensed_photos")
            or total != execution.get("n_frozen_source_photo_opportunities")):
        raise ValueError("Photo-count aggregation disagrees with original source receipt")

    desc = sorted(counts)
    def percentile(q: float) -> int:
        return desc[int((n_taxa - 1) * q)]
    upper_one_percent = sum(desc[-max(1, (n_taxa + 99) // 100):])
    summary = {
        "schema": "fcp_42111_original_photo_opportunity_audit_v1",
        "status": "VERIFIED_COMPLETE_COUNTS_ONLY_NOT_COLOUR_MEASUREMENT",
        "snapshot_date": "2026-09-27",
        "source_frame_species": species,
        "source_archive_sha256": ARCHIVE_SHA,
        "source_frame_sha256": FRAME_SHA,
        "original_source_taxon_IDs_and_names_independently_matched": original_frame_matched,
        "count_ledger_sha256": digest,
        "all_source_photo_opportunity_ids": total,
        "taxa_with_at_least_one_photo": n_with,
        "taxa_with_zero_photos": n_zero,
        "minimum_source_photo_depth_upper_bound_for_W": {
            str(threshold): sum(n >= threshold for n in counts)
            for threshold in (1, 10, 20, 40, 60, 80, 100)
        },
        "photo_count_per_taxon_quantiles": {
            "p0": desc[0], "p25": percentile(.25), "p50": percentile(.5),
            "p75": percentile(.75), "p90": percentile(.90),
            "p95": percentile(.95), "p99": percentile(.99),
            "max": desc[-1]
        },
        "fraction_of_source_photos_in_top_one_percent_of_taxa": (
            upper_one_percent / total if total else None
        ),
        "all_source_photo_id_manifests_frozen": False,
        "all_original_image_pixels_measured": False,
        "flower_colour_results": False,
        "W_D_H2_or_climate_soil_results": False,
        "note": (
            "All thresholds denote available source photographs, NOT evaluable "
            "flower pixels; >=40 source photos is only an upper bound on species "
            "with >=40 CLASSIFIABLE new flower-colour photographs."
        ),
    }
    return summary


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--counts", type=Path, required=True)
    ap.add_argument("--source-receipt", type=Path, required=True)
    ap.add_argument("--execution-receipt", type=Path, required=True)
    ap.add_argument("--original-frame", type=Path,
                    help="Byte-frozen original 42,111 taxon file; required for full original-source audit")
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args(argv)
    result = audit(args.counts, args.source_receipt, args.execution_receipt,
                   frame_file=args.original_frame)
    if args.output.exists():
        raise ValueError("Audit output already exists: refuse silent overwrite")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("FCP_COMPLETE_REAL_SOURCE_PHOTO_COUNTS_ONLY", json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

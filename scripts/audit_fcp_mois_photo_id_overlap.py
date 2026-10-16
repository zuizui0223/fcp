#!/usr/bin/env python3
"""Audit actual photo-ID overlap: historical MOIS atlas vs frozen FCP high-depth cohorts.

Only original identifiers are compared; no model is fit or reclassified.
Expected SHA256 and denominators are frozen in the October 2026 source receipts.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import zipfile
from collections import Counter
from datetime import date
from pathlib import Path

MOIS_ZIP_SHA = "2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0"
FCP_COHORT_SHA = {
    "discovery": "ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4",
    "validation": "0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6",
    "prospective": "57630fc9f281bce94a0c40a70aaf7bce879dde93d6154175adcd021e8f5c1186",
}
FCP_EXPECTED = {"discovery": 50000, "validation": 50000, "prospective": 49900}
KNOWN_HUES = {"white", "yellow_orange", "red_pink", "blue_purple"}
MOIS_MEMBERS = (
    "fcp42111_real_site_abiotic/breadth_flower_colour_climate_soil.csv.gz",
    "fcp42111_real_site_abiotic/taxon_cell_flower_colour_climate_soil.csv.gz",
)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for piece in iter(lambda: f.read(1 << 20), b""):
            h.update(piece)
    return h.hexdigest()


def id_digest(ids: set[str]) -> str:
    """Digest of sorted ids without publishing individual source identifiers."""
    h = hashlib.sha256()
    for ident in sorted(ids, key=lambda x: (len(x), x)):
        h.update(ident.encode("ascii"))
        h.update(b"\n")
    return h.hexdigest()


def get_mois(path: Path, *, verify_hash: bool = True) -> tuple[dict[str, tuple[str, str, str]], dict]:
    digest = sha256_file(path)
    if verify_hash and digest != MOIS_ZIP_SHA:
        raise ValueError(f"Original MOIS photographic zip SHA mismatch: {digest}")
    records = {}
    source_rows = 0
    duplicates = 0
    conflicts = 0
    with zipfile.ZipFile(path) as zipf:
        for member in MOIS_MEMBERS:
            with gzip.open(io.BytesIO(zipf.read(member)), "rt", encoding="utf-8", newline="") as f:
                reader = csv.DictReader(f)
                if not {"photo_id", "observation_id", "inat_taxon_id", "morph"}.issubset(reader.fieldnames or []):
                    raise ValueError(f"Required original source columns missing from {member}")
                for row in reader:
                    source_rows += 1
                    pid = row["photo_id"].strip()
                    if not pid:
                        raise ValueError("Blank original MOIS photo identifier")
                    item = (row["morph"].strip(), row["inat_taxon_id"].strip(), row["observation_id"].strip())
                    if pid in records:
                        duplicates += 1
                        prev = records[pid]
                        if item != prev:
                            conflicts += 1
                            # Conflicting metadata must never be silently treated as new evidence.
                            # Keep the original breadth-level row, but report all conflicts.
                        continue
                    records[pid] = item
    if verify_hash and (len(records), source_rows, duplicates) != (100543, 127448, 26905):
        raise ValueError("MOIS original source census has changed")
    states = Counter(rec[0] for rec in records.values())
    n_classified = sum(states[k] for k in KNOWN_HUES)
    if verify_hash and (n_classified, len(records) - n_classified) != (45953, 54590):
        raise ValueError(f"MOIS original four-hue denominator drift: {states}")
    return records, {
        "archive_sha256": digest, "original_rows": source_rows,
        "duplicate_photo_rows": duplicates, "duplicate_conflicting_source_rows": conflicts,
        "distinct_photo_ids": len(records),
        "classified_photo_ids": n_classified,
        "unclassified_photo_ids": len(records) - n_classified,
        "photo_id_set_sha256": id_digest(set(records)),
    }


def get_fcp(path: Path, cohort: str, *, verify_hash: bool = True) -> tuple[dict[str, tuple[str, str]], dict]:
    digest = sha256_file(path)
    if verify_hash and digest != FCP_COHORT_SHA[cohort]:
        raise ValueError(f"Original FCP {cohort} measured-table SHA mismatch: {digest}")
    ids = {}
    rows = 0
    duplicate_rows = 0
    with path.open("rt", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        if "photo_id" not in (reader.fieldnames or []):
            raise ValueError(f"Missing FCP photo_id column in {cohort}; fields={reader.fieldnames}")
        for row in reader:
            rows += 1
            pid = row["photo_id"].strip()
            if not pid:
                raise ValueError(f"Blank FCP photo ID, {cohort} row {rows}")
            if pid in ids:
                duplicate_rows += 1
            else:
                ids[pid] = (row.get("morph", "").strip(), row.get("inat_taxon_id", "").strip())
    if verify_hash and (rows, len(ids), duplicate_rows) != (FCP_EXPECTED[cohort], FCP_EXPECTED[cohort], 0):
        raise ValueError(f"Original FCP cohort denominator drift: {cohort} {rows}/{len(ids)}")
    return ids, {
        "source_sha256": digest, "rows": rows, "distinct_photo_ids": len(ids),
        "duplicate_photo_rows": duplicate_rows, "photo_id_set_sha256": id_digest(set(ids)),
    }


def audit(mois: dict, cohorts: dict[str, dict]) -> dict:
    overall = set()
    per = {}
    fcp_cross_duplicates = 0
    for name in ("discovery", "validation", "prospective"):
        fcp = cohorts[name]
        shared = set(fcp) & set(mois)
        fcp_cross_duplicates += len(set(fcp) & overall)
        overall.update(fcp)
        both_class = 0
        exact_agreement = 0
        conflicts_taxon = 0
        for pid in shared:
            m_colour, m_taxon, _ = mois[pid]
            f_colour, f_taxon = fcp[pid]
            if m_colour in KNOWN_HUES and f_colour in KNOWN_HUES:
                both_class += 1
                if m_colour == f_colour:
                    exact_agreement += 1
            if m_taxon and f_taxon and m_taxon != f_taxon:
                conflicts_taxon += 1
        per[name] = {
            "fcp_photo_ids": len(fcp), "shared_with_mois_photo_ids": len(shared),
            "fraction_of_fcp_cohort_shared": len(shared) / len(fcp),
            "fraction_of_mois_source_shared": len(shared) / len(mois),
            "shared_mois_four_hue_classifiable": sum(mois[p][0] in KNOWN_HUES for p in shared),
            "shared_both_four_hue_classifiable": both_class,
            "shared_both_label_exact_agreement": exact_agreement,
            "shared_both_label_disagreement": both_class - exact_agreement,
            "matched_photo_taxon_id_conflicts_if_both_recorded": conflicts_taxon,
            "shared_photo_id_set_sha256": id_digest(shared),
        }
    union_overlap = overall & set(mois)
    return {
        "cohorts": per,
        "fcp_total_distinct_photo_ids": len(overall),
        "fcp_cross_cohort_duplicate_photo_ids": fcp_cross_duplicates,
        "shared_distinct_photo_ids": len(union_overlap),
        "fraction_of_fcp_all_photo_ids_shared": len(union_overlap) / len(overall),
        "fraction_of_mois_all_original_photo_ids_shared": len(union_overlap) / len(mois),
        "shared_photo_id_set_sha256": id_digest(union_overlap),
        "can_claim_independent_source_replication": False,
        "note": "Matching photo IDs quantify literal source reuse; disjoint IDs would still share provider and measurement-system provenance.",
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mois-zip", required=True, type=Path)
    for name in FCP_EXPECTED:
        ap.add_argument(f"--{name}", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--fixture-mode", action="store_true", help="Disable frozen checksum and denominator assertions only for synthetic tests.")
    a = ap.parse_args()
    mois, mm = get_mois(a.mois_zip, verify_hash=not a.fixture_mode)
    cohorts = {}
    meta = {}
    for name in FCP_EXPECTED:
        cohorts[name], meta[name] = get_fcp(getattr(a, name), name, verify_hash=not a.fixture_mode)
    overlap = audit(mois, cohorts)
    if not a.fixture_mode and overlap["fcp_cross_cohort_duplicate_photo_ids"] != 0:
        raise ValueError("FCP three cohort photo IDs unexpectedly intersect")
    result = {
        "schema": "fcp_mois_original_photo_id_overlap_audit_v1",
        "audit_type": "cross_study_original_photo_id_provenance_only_not_biological_model",
        "report_date": "2026-10-16",
        "frozen_mois_source": mm,
        "frozen_fcp_cohorts": meta,
        "comparison": overlap,
        "original_mois_artifact": "zuizui0223/fcp actions artifact 11587115220 from run 37861582446",
        "fcp_input_commits": {"discovery_and_validation": "5142f7951af0dde5364bb047a566d67e8c479e51",
                              "prospective": "7e538e5c51c05a7cc47b2fcf53eea92634c8a863"},
        "scope_warning": "The overlap audit is a technical nonindependence disclosure, NOT evidence of within-species adaptation, monthly climate effects or independent botanical colour replication.",
    }
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("PHOTO_ID_AUDIT_COMPLETE", json.dumps({
        "mois": mm["distinct_photo_ids"], "fcp": overlap["fcp_total_distinct_photo_ids"],
        "shared": overlap["shared_distinct_photo_ids"],
        "by_cohort": {k: v["shared_with_mois_photo_ids"] for k, v in overlap["cohorts"].items()}
    }, sort_keys=True))


if __name__ == "__main__":
    main()

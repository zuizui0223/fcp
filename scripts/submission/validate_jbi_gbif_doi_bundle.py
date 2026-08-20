#!/usr/bin/env python3
"""Validate the exact GBIF occurrence and DOI-preparation bundle."""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent
BUNDLE = ROOT / "docs/supporting/jbi_gbif_doi_bundle"
PROTOCOL = ROOT / "docs/jbi_gbif_doi_protocol.md"
EXPECTED = {
    "rows": 58455,
    "species": 34,
    "parent_datasets": 389,
    "unique_occurrence_keys": 58455,
    "source_csv_sha256": "b0614a729acde5a1daab599d52c39ac4018583e1f73d83e5304ac0afa6f6e7ad",
    "exact_archive_sha256": "f25ae0cf2c84c45ae461a932d6c6063edda64591913a2495e4a3da82d573f094",
}
EXPECTED_BASIS = {
    "PRESERVED_SPECIMEN",
    "MATERIAL_SAMPLE",
    "HUMAN_OBSERVATION",
    "MACHINE_OBSERVATION",
    "OBSERVATION",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def fail(message: str) -> None:
    raise SystemExit(message)


def main() -> None:
    manifest_path = BUNDLE / "jbi_gbif_doi_bundle_manifest.json"
    archive_path = BUNDLE / "jbi_gbif_exact_occurrence_subset.csv.gz"
    parent_path = BUNDLE / "jbi_gbif_parent_dataset_counts.csv"
    species_path = BUNDLE / "jbi_gbif_exact_species_counts.csv"
    request_path = BUNDLE / "jbi_gbif_broad_download_request.json"
    metadata_path = BUNDLE / "jbi_gbif_derived_dataset_metadata_template.json"
    for path in (manifest_path, archive_path, parent_path, species_path, request_path, metadata_path, PROTOCOL):
        if not path.is_file() or path.stat().st_size == 0:
            fail(f"Missing or empty DOI bundle file: {path}")

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    for key, expected in EXPECTED.items():
        if manifest.get(key) != expected:
            fail(f"Manifest mismatch for {key}: {manifest.get(key)!r} != {expected!r}")
    referenced_hashes = {
        archive_path.name: manifest["exact_archive_sha256"],
        parent_path.name: manifest["parent_counts_sha256"],
        species_path.name: manifest["species_counts_sha256"],
        request_path.name: manifest["broad_download_request_sha256"],
        metadata_path.name: manifest["derived_metadata_template_sha256"],
    }
    for name, expected_hash in referenced_hashes.items():
        actual = sha256(BUNDLE / name)
        if actual != expected_hash:
            fail(f"SHA-256 mismatch for {name}: {actual} != {expected_hash}")

    with parent_path.open(newline="", encoding="utf-8") as handle:
        parent_rows = list(csv.DictReader(handle))
    if len(parent_rows) != EXPECTED["parent_datasets"]:
        fail(f"Parent dataset count mismatch: {len(parent_rows)}")
    if len({row["datasetKey"] for row in parent_rows}) != EXPECTED["parent_datasets"]:
        fail("Parent dataset keys are not unique")
    if sum(int(row["record_count"]) for row in parent_rows) != EXPECTED["rows"]:
        fail("Parent dataset record counts do not sum to 58,455")

    with species_path.open(newline="", encoding="utf-8") as handle:
        species_rows = list(csv.DictReader(handle))
    if len(species_rows) != EXPECTED["species"]:
        fail(f"Species count mismatch: {len(species_rows)}")
    if len({row["canonical_name"] for row in species_rows}) != EXPECTED["species"]:
        fail("Species names are not unique")
    if sum(int(row["record_count"]) for row in species_rows) != EXPECTED["rows"]:
        fail("Species record counts do not sum to 58,455")

    dataset_counts: Counter[str] = Counter()
    species_counts: Counter[str] = Counter()
    occurrence_keys: set[str] = set()
    rows = 0
    with gzip.open(archive_path, "rt", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"gbif_key", "datasetKey", "canonical_name", "decimalLatitude", "decimalLongitude"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            fail(f"Exact archive missing columns: {sorted(missing)}")
        for row in reader:
            rows += 1
            occurrence_key = row["gbif_key"].strip()
            dataset_key = row["datasetKey"].strip()
            species = row["canonical_name"].strip()
            if not occurrence_key or not dataset_key or not species:
                fail(f"Blank citation field at row {rows + 1}")
            if occurrence_key in occurrence_keys:
                fail(f"Duplicate occurrence key in exact archive: {occurrence_key}")
            occurrence_keys.add(occurrence_key)
            dataset_counts[dataset_key] += 1
            species_counts[species] += 1
    if rows != EXPECTED["rows"]:
        fail(f"Exact archive row count mismatch: {rows}")
    if len(occurrence_keys) != EXPECTED["unique_occurrence_keys"]:
        fail(f"Unique occurrence key mismatch: {len(occurrence_keys)}")
    if dataset_counts != Counter({row["datasetKey"]: int(row["record_count"]) for row in parent_rows}):
        fail("Parent dataset count file does not match exact archive")
    if species_counts != Counter({row["canonical_name"]: int(row["record_count"]) for row in species_rows}):
        fail("Species count file does not match exact archive")

    request = json.loads(request_path.read_text(encoding="utf-8"))
    if request.get("notificationAddresses") != ["REPLACE_WITH_GBIF_ACCOUNT_EMAIL"]:
        fail("Prepared download request notification placeholder changed")
    if request.get("format") != "SIMPLE_CSV":
        fail("Prepared download format must be SIMPLE_CSV")
    if "checklistKey" in request:
        fail("Prepared request must not change taxonomy using checklistKey")
    predicates = request.get("predicate", {}).get("predicates", [])
    taxon = next((p for p in predicates if p.get("key") == "TAXON_KEY"), None)
    basis = next((p for p in predicates if p.get("key") == "BASIS_OF_RECORD"), None)
    if taxon is None or len(taxon.get("values", [])) != EXPECTED["species"]:
        fail("Prepared request must contain 34 TAXON_KEY values")
    if len(set(taxon["values"])) != EXPECTED["species"]:
        fail("Prepared TAXON_KEY values are not unique")
    if basis is None or set(basis.get("values", [])) != EXPECTED_BASIS:
        fail("Prepared BASIS_OF_RECORD values changed")
    required_equals = {
        ("HAS_COORDINATE", "true"),
        ("HAS_GEOSPATIAL_ISSUE", "false"),
        ("OCCURRENCE_STATUS", "PRESENT"),
    }
    present_equals = {(p.get("key"), p.get("value")) for p in predicates if p.get("type") == "equals"}
    if not required_equals.issubset(present_equals):
        fail(f"Prepared equality predicates changed: {present_equals}")

    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if metadata.get("source_url") != "REPLACE_WITH_PERMANENT_ARCHIVE_URL":
        fail("Derived Dataset source URL placeholder changed before a permanent archive exists")
    if metadata.get("record_count") != EXPECTED["rows"] or metadata.get("parent_dataset_count") != EXPECTED["parent_datasets"]:
        fail("Derived Dataset metadata counts changed")
    if metadata.get("exact_archive_sha256") != EXPECTED["exact_archive_sha256"]:
        fail("Derived Dataset metadata archive SHA changed")

    protocol = PROTOCOL.read_text(encoding="utf-8")
    for phrase in (
        "Why two GBIF citation records are needed",
        "GBIF Derived Dataset",
        "389 parent GBIF datasets",
        EXPECTED["exact_archive_sha256"],
        "Items that still require account-holder action",
    ):
        if phrase not in protocol:
            fail(f"GBIF DOI protocol missing phrase: {phrase}")

    print(
        json.dumps(
            {
                "status": "pass",
                "rows": rows,
                "species": len(species_counts),
                "parent_datasets": len(dataset_counts),
                "unique_occurrence_keys": len(occurrence_keys),
                "exact_archive_sha256": sha256(archive_path),
                "request_taxon_keys": len(taxon["values"]),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

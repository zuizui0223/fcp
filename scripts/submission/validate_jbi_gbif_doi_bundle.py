#!/usr/bin/env python3
"""Validate the frozen 34-species GBIF citation bundle."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_BUNDLE = REPO_ROOT / "docs/supporting/jbi_gbif_doi_bundle"
EXPECTED = {
    "rows": 58455,
    "species": 34,
    "parent_datasets": 389,
    "unique_occurrence_keys": 58455,
    "source_csv_sha256": "b0614a729acde5a1daab599d52c39ac4018583e1f73d83e5304ac0afa6f6e7ad",
    "exact_archive_sha256": "f25ae0cf2c84c45ae461a932d6c6063edda64591913a2495e4a3da82d573f094",
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
    ap = argparse.ArgumentParser()
    ap.add_argument("--bundle-dir", type=Path, default=DEFAULT_BUNDLE)
    args = ap.parse_args()
    bundle = args.bundle_dir

    manifest_path = bundle / "jbi_gbif_doi_bundle_manifest.json"
    archive_path = bundle / "jbi_gbif_exact_occurrence_subset.csv.gz"
    parent_path = bundle / "jbi_gbif_parent_dataset_counts.csv"
    species_path = bundle / "jbi_gbif_exact_species_counts.csv"
    request_path = bundle / "jbi_gbif_broad_download_request.json"
    metadata_path = bundle / "jbi_gbif_derived_dataset_metadata_template.json"
    required = [manifest_path, archive_path, parent_path, species_path, request_path, metadata_path]
    for path in required:
        if not path.is_file() or path.stat().st_size == 0:
            fail(f"Missing or empty DOI bundle file: {path}")

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    for key, expected in EXPECTED.items():
        if manifest.get(key) != expected:
            fail(f"Manifest mismatch for {key}: {manifest.get(key)!r} != {expected!r}")
    if sha256(archive_path) != EXPECTED["exact_archive_sha256"]:
        fail("Exact occurrence archive SHA-256 changed")

    with parent_path.open(newline="", encoding="utf-8") as handle:
        parent_rows = list(csv.DictReader(handle))
    with species_path.open(newline="", encoding="utf-8") as handle:
        species_rows = list(csv.DictReader(handle))
    if len(parent_rows) != EXPECTED["parent_datasets"]:
        fail(f"Expected 389 parent GBIF datasets; found {len(parent_rows)}")
    if len(species_rows) != EXPECTED["species"]:
        fail(f"Expected 34 species; found {len(species_rows)}")

    datasets: Counter[str] = Counter()
    species: Counter[str] = Counter()
    occurrence_keys: set[str] = set()
    rows = 0
    with gzip.open(archive_path, "rt", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required_cols = {"gbif_key", "datasetKey", "canonical_name", "decimalLatitude", "decimalLongitude"}
        missing = required_cols - set(reader.fieldnames or [])
        if missing:
            fail(f"Exact archive missing columns: {sorted(missing)}")
        for row in reader:
            rows += 1
            key = row["gbif_key"].strip()
            if not key or key in occurrence_keys:
                fail(f"Blank or duplicate GBIF occurrence key at data row {rows}")
            occurrence_keys.add(key)
            datasets[row["datasetKey"].strip()] += 1
            species[row["canonical_name"].strip()] += 1

    if rows != EXPECTED["rows"] or len(occurrence_keys) != EXPECTED["unique_occurrence_keys"]:
        fail(f"Exact archive occurrence count mismatch: rows={rows}, unique={len(occurrence_keys)}")
    if len(datasets) != EXPECTED["parent_datasets"] or len(species) != EXPECTED["species"]:
        fail(f"Exact archive dimensions changed: datasets={len(datasets)}, species={len(species)}")
    if datasets != Counter({r["datasetKey"]: int(r["record_count"]) for r in parent_rows}):
        fail("Parent-dataset counts do not match the exact archive")
    if species != Counter({r["canonical_name"]: int(r["record_count"]) for r in species_rows}):
        fail("Species counts do not match the exact archive")

    request = json.loads(request_path.read_text(encoding="utf-8"))
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if request.get("notificationAddresses") != ["REPLACE_WITH_GBIF_ACCOUNT_EMAIL"]:
        fail("Prepared broad-download notification placeholder changed")
    if metadata.get("record_count") != EXPECTED["rows"]:
        fail("Derived Dataset metadata row count changed")
    if metadata.get("exact_archive_sha256") != EXPECTED["exact_archive_sha256"]:
        fail("Derived Dataset metadata archive SHA changed")

    print(json.dumps({
        "status": "pass",
        "bundle": str(bundle),
        "rows": rows,
        "species": len(species),
        "parent_datasets": len(datasets),
        "exact_archive_sha256": sha256(archive_path),
    }, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Prepare the GBIF citation bundle for the JBI manuscript.

Two distinct citation objects are prepared:

1. A regular GBIF occurrence-download request representing the broad pre-filtered
   retrieval. GBIF creates a DOI for this asynchronous download after authenticated
   submission.
2. A GBIF Derived Dataset registration bundle representing the exact locally capped
   and coordinate-deduplicated 58,455-row analytical subset.

The script never embeds credentials. Optional download submission uses the GBIF
username supplied on the command line and the password from GBIF_PASSWORD.
"""
from __future__ import annotations

import argparse
import base64
import csv
import gzip
import hashlib
import json
import os
import shutil
import sys
import urllib.error
import urllib.request
import zipfile
from collections import Counter
from pathlib import Path
from typing import Iterable

EXPECTED_SOURCE_CSV_SHA256 = "b0614a729acde5a1daab599d52c39ac4018583e1f73d83e5304ac0afa6f6e7ad"
EXPECTED_ROWS = 58_455
EXPECTED_SPECIES = 34
EXPECTED_PARENT_DATASETS = 389
EXPECTED_UNIQUE_OCCURRENCES = 58_455
EXPECTED_DETERMINISTIC_GZIP_SHA256 = "f25ae0cf2c84c45ae461a932d6c6063edda64591913a2495e4a3da82d573f094"
SOURCE_WORKFLOW_RUN = 30065308023
SOURCE_ARTIFACT_ID = 8586269620
SOURCE_ARTIFACT_DIGEST = "sha256:4247fe3dae52f7723a25df6cd9956fa8b2f220a3fc101e088a3f0cccd768447a"
ACCEPTED_BASIS = [
    "PRESERVED_SPECIMEN",
    "MATERIAL_SAMPLE",
    "HUMAN_OBSERVATION",
    "MACHINE_OBSERVATION",
    "OBSERVATION",
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def locate_artifact_file(root: Path, filename: str) -> Path:
    candidates = list(root.rglob(filename))
    if len(candidates) != 1:
        raise RuntimeError(f"Expected one {filename!r} below {root}; found {len(candidates)}")
    return candidates[0]


def unpack_if_needed(path: Path, workdir: Path) -> Path:
    if path.is_dir():
        return path
    if not zipfile.is_zipfile(path):
        raise ValueError(f"Artifact must be a directory or ZIP archive: {path}")
    target = workdir / "artifact"
    target.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(path) as archive:
        archive.extractall(target)
    return target


def read_taxon_keys(path: Path) -> list[str]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != EXPECTED_SPECIES:
        raise RuntimeError(f"Taxon audit must contain {EXPECTED_SPECIES} rows; found {len(rows)}")
    keys: list[str] = []
    for row in sorted(rows, key=lambda x: x["canonical_name"]):
        value = str(row.get("taxon_key") or "").strip()
        if not value.isdigit():
            raise RuntimeError(f"Invalid numeric GBIF taxon key for {row.get('canonical_name')}: {value!r}")
        keys.append(value)
    if len(set(keys)) != EXPECTED_SPECIES:
        raise RuntimeError("GBIF taxon keys are not unique")
    return keys


def analyse_exact_csv(path: Path) -> dict[str, object]:
    source_hash = sha256(path)
    if source_hash != EXPECTED_SOURCE_CSV_SHA256:
        raise RuntimeError(f"Source occurrence CSV SHA mismatch: {source_hash}")

    dataset_counts: Counter[str] = Counter()
    species_counts: Counter[str] = Counter()
    occurrence_keys: set[str] = set()
    total = 0
    columns: list[str] | None = None

    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        columns = reader.fieldnames
        required = {"canonical_name", "gbif_key", "datasetKey", "decimalLatitude", "decimalLongitude"}
        missing = required - set(columns or [])
        if missing:
            raise RuntimeError(f"Exact occurrence CSV is missing columns: {sorted(missing)}")
        for row in reader:
            total += 1
            dataset_key = str(row.get("datasetKey") or "").strip()
            occurrence_key = str(row.get("gbif_key") or "").strip()
            species = str(row.get("canonical_name") or "").strip()
            if not dataset_key or not occurrence_key or not species:
                raise RuntimeError(f"Missing citation field at data row {total + 1}")
            if occurrence_key in occurrence_keys:
                raise RuntimeError(f"Duplicate GBIF occurrence key: {occurrence_key}")
            occurrence_keys.add(occurrence_key)
            dataset_counts[dataset_key] += 1
            species_counts[species] += 1

    checks = {
        "rows": total,
        "species": len(species_counts),
        "parent_datasets": len(dataset_counts),
        "unique_occurrence_keys": len(occurrence_keys),
    }
    expected = {
        "rows": EXPECTED_ROWS,
        "species": EXPECTED_SPECIES,
        "parent_datasets": EXPECTED_PARENT_DATASETS,
        "unique_occurrence_keys": EXPECTED_UNIQUE_OCCURRENCES,
    }
    if checks != expected:
        raise RuntimeError(f"Exact occurrence checks failed: {checks} != {expected}")
    return {
        **checks,
        "columns": columns,
        "source_csv_sha256": source_hash,
        "dataset_counts": dataset_counts,
        "species_counts": species_counts,
    }


def write_counts(path: Path, first_column: str, counter: Counter[str], *, alphabetical: bool) -> None:
    if alphabetical:
        rows = sorted(counter.items())
    else:
        rows = sorted(counter.items(), key=lambda item: (-item[1], item[0]))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow([first_column, "record_count"])
        writer.writerows(rows)


def deterministic_gzip(source: Path, target: Path) -> None:
    with source.open("rb") as src, target.open("wb") as raw:
        with gzip.GzipFile(
            filename="jbi_gbif_exact_occurrence_subset.csv",
            mode="wb",
            fileobj=raw,
            compresslevel=9,
            mtime=0,
        ) as compressed:
            shutil.copyfileobj(src, compressed)
    digest = sha256(target)
    if digest != EXPECTED_DETERMINISTIC_GZIP_SHA256:
        raise RuntimeError(f"Deterministic gzip SHA mismatch: {digest}")


def download_request(taxon_keys: Iterable[str], email: str) -> dict[str, object]:
    return {
        "notificationAddresses": [email],
        "sendNotification": True,
        "format": "SIMPLE_CSV",
        # No checklistKey is supplied: the numeric taxon keys were resolved against
        # the same legacy GBIF Backbone used by the original v1 species-match calls.
        "predicate": {
            "type": "and",
            "predicates": [
                {"type": "in", "key": "TAXON_KEY", "values": list(taxon_keys)},
                {"type": "equals", "key": "HAS_COORDINATE", "value": "true"},
                {"type": "equals", "key": "HAS_GEOSPATIAL_ISSUE", "value": "false"},
                {"type": "equals", "key": "OCCURRENCE_STATUS", "value": "PRESENT"},
                {"type": "in", "key": "BASIS_OF_RECORD", "values": ACCEPTED_BASIS},
                {
                    "type": "or",
                    "predicates": [
                        {
                            "type": "lessThanOrEquals",
                            "key": "COORDINATE_UNCERTAINTY_IN_METERS",
                            "value": "20000",
                        },
                        {
                            "type": "isNull",
                            "parameter": "COORDINATE_UNCERTAINTY_IN_METERS",
                        },
                    ],
                },
            ],
        },
    }


def submit_download(request: dict[str, object], username: str) -> str:
    password = os.environ.get("GBIF_PASSWORD")
    if not password:
        raise RuntimeError("GBIF_PASSWORD is required for --submit-download")
    credentials = base64.b64encode(f"{username}:{password}".encode()).decode()
    payload = json.dumps(request).encode("utf-8")
    http_request = urllib.request.Request(
        "https://api.gbif.org/v1/occurrence/download/request",
        data=payload,
        method="POST",
        headers={
            "Authorization": f"Basic {credentials}",
            "Content-Type": "application/json",
            "Accept": "text/plain, application/json",
            "User-Agent": "fcp-jbi-gbif-doi-preparation/1.0",
        },
    )
    try:
        with urllib.request.urlopen(http_request, timeout=120) as response:
            body = response.read().decode("utf-8").strip()
    except urllib.error.HTTPError as error:
        detail = error.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"GBIF download request failed ({error.code}): {detail}") from error
    if not body:
        raise RuntimeError("GBIF returned an empty download key")
    return body


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact", required=True, type=Path)
    parser.add_argument("--taxon-audit", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--notification-email", default="REPLACE_WITH_GBIF_ACCOUNT_EMAIL")
    parser.add_argument("--gbif-username")
    parser.add_argument("--submit-download", action="store_true")
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    workdir = args.outdir / ".work"
    workdir.mkdir(exist_ok=True)
    artifact_root = unpack_if_needed(args.artifact, workdir)
    source_csv = locate_artifact_file(artifact_root, "gbif_occurrence_sample_paginated.csv")
    analysis = analyse_exact_csv(source_csv)

    exact_archive = args.outdir / "jbi_gbif_exact_occurrence_subset.csv.gz"
    parent_counts = args.outdir / "jbi_gbif_parent_dataset_counts.csv"
    species_counts = args.outdir / "jbi_gbif_exact_species_counts.csv"
    request_path = args.outdir / "jbi_gbif_broad_download_request.json"
    metadata_path = args.outdir / "jbi_gbif_derived_dataset_metadata_template.json"
    manifest_path = args.outdir / "jbi_gbif_doi_bundle_manifest.json"

    deterministic_gzip(source_csv, exact_archive)
    write_counts(parent_counts, "datasetKey", analysis["dataset_counts"], alphabetical=False)
    write_counts(species_counts, "canonical_name", analysis["species_counts"], alphabetical=True)
    taxon_keys = read_taxon_keys(args.taxon_audit)
    request = download_request(taxon_keys, args.notification_email)
    request_path.write_text(json.dumps(request, indent=2) + "\n", encoding="utf-8")

    metadata = {
        "title": "Quality-filtered GBIF occurrence subset for spatial organization of flower-colour variation",
        "source_url": "REPLACE_WITH_PERMANENT_ARCHIVE_URL",
        "description": (
            "Exact 58,455-record GBIF-mediated occurrence subset used for the paginated "
            "sensitivity analysis in the Journal of Biogeography manuscript. Records were "
            "retrieved for 34 audited plant species, capped at 3,000 records per species, "
            "restricted to accepted basis-of-record categories, filtered for geospatial "
            "issues and coordinate uncertainty, and deduplicated at 0.001-degree coordinate "
            "precision. The accompanying parent-dataset count file credits 389 contributing "
            "GBIF datasets."
        ),
        "record_count": EXPECTED_ROWS,
        "parent_dataset_count": EXPECTED_PARENT_DATASETS,
        "parent_dataset_counts_file": parent_counts.name,
        "exact_archive_file": exact_archive.name,
        "exact_archive_sha256": sha256(exact_archive),
        "source_csv_sha256": analysis["source_csv_sha256"],
        "optional_original_download_doi": None,
        "registration_date": None,
    }
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")

    manifest = {
        "source_workflow_run": SOURCE_WORKFLOW_RUN,
        "source_artifact_id": SOURCE_ARTIFACT_ID,
        "source_artifact_digest": SOURCE_ARTIFACT_DIGEST,
        "source_csv_filename": source_csv.name,
        "source_csv_sha256": analysis["source_csv_sha256"],
        "source_csv_bytes": source_csv.stat().st_size,
        "rows": analysis["rows"],
        "species": analysis["species"],
        "parent_datasets": analysis["parent_datasets"],
        "unique_occurrence_keys": analysis["unique_occurrence_keys"],
        "exact_archive": exact_archive.name,
        "exact_archive_bytes": exact_archive.stat().st_size,
        "exact_archive_sha256": sha256(exact_archive),
        "parent_counts_file": parent_counts.name,
        "parent_counts_sha256": sha256(parent_counts),
        "species_counts_file": species_counts.name,
        "species_counts_sha256": sha256(species_counts),
        "broad_download_request_file": request_path.name,
        "broad_download_request_sha256": sha256(request_path),
        "derived_metadata_template_file": metadata_path.name,
        "derived_metadata_template_sha256": sha256(metadata_path),
        "taxon_key_count": len(taxon_keys),
        "download_taxonomy_note": (
            "No checklistKey is supplied because the original numeric taxon keys were "
            "resolved with the legacy GBIF Backbone through the v1 species-match API."
        ),
        "citation_boundary": (
            "The regular GBIF download DOI represents the broad predicate-based retrieval. "
            "The Derived Dataset DOI must represent the exact locally capped and "
            "coordinate-deduplicated archive."
        ),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    shutil.rmtree(workdir)

    if args.submit_download:
        if not args.gbif_username:
            raise RuntimeError("--gbif-username is required with --submit-download")
        key = submit_download(request, args.gbif_username)
        (args.outdir / "jbi_gbif_download_key.txt").write_text(key + "\n", encoding="utf-8")
        print(f"GBIF download key: {key}")

    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

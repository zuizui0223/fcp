#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

U0_SPECIES = 42111
U100_SPECIES = 4730
LEGACY_EXCLUDED_SPECIES = 1000
P100_SPECIES = 3730
P500_SELECTED_SPECIES = 500
THIRD_CANDIDATE_SPECIES = 3230
U100_MIN = 100
PARENT_WORKFLOW_RUN = 34708044962
PARENT_ARTIFACT_ID = 10302477571
PARENT_ARTIFACT_DIGEST = "sha256:7e43763e3f1d9fa3ff108154dc6ef9443430b70fdb62532e2f6845582d607abb"
P100_SHA256 = "1473aad680fe2fa84903c5e11ee104fd0eceef828957f16c2ef1a35f9dd6993c"
P500_SHA256 = "f54d07fb2338a20a3f92808b58f0ebc948ab0d5af3cb01fb26702f861184b7a4"
EXPECTED_CANDIDATE_SHA256 = "7fc0337074b55a91c0b50773a8d5c3ef82e877074c8b3cacc21f9af58e0ba77e"
ALLOWED_SOURCE_COLUMNS = {
    "inat_taxon_id",
    "species",
    "after_observer_cap",
    "selection_hash",
    "prospective_rank",
}
OUTPUT_COLUMNS = ["inat_taxon_id", "species", "after_observer_cap"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_outcome_blind_source(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        fields = set(reader.fieldnames or [])
        required = set(OUTPUT_COLUMNS)
        if not required.issubset(fields):
            raise RuntimeError(f"missing required source columns: {sorted(required - fields)}")
        unexpected = fields - ALLOWED_SOURCE_COLUMNS
        if unexpected:
            raise RuntimeError(f"forbidden or unexpected source columns: {sorted(unexpected)}")
        rows = []
        for row in reader:
            taxon_id = int(row["inat_taxon_id"])
            species = row["species"].strip()
            depth = int(row["after_observer_cap"])
            if depth < U100_MIN:
                raise RuntimeError(f"species {species} below frozen U100 gate: {depth}")
            rows.append(
                {
                    "inat_taxon_id": taxon_id,
                    "species": species,
                    "after_observer_cap": depth,
                }
            )
    return rows


def _validate_unique(rows: list[dict], label: str) -> None:
    taxon_ids = [r["inat_taxon_id"] for r in rows]
    species = [r["species"] for r in rows]
    if len(taxon_ids) != len(set(taxon_ids)):
        raise RuntimeError(f"duplicate inat_taxon_id in {label}")
    if len(species) != len(set(species)):
        raise RuntimeError(f"duplicate species in {label}")


def derive_candidate_rows(p100_rows: list[dict], p500_rows: list[dict]) -> list[dict]:
    _validate_unique(p100_rows, "P100")
    _validate_unique(p500_rows, "P500")
    p100_by_id = {r["inat_taxon_id"]: r for r in p100_rows}
    p500_ids = {r["inat_taxon_id"] for r in p500_rows}
    missing = sorted(p500_ids - set(p100_by_id))
    if missing:
        raise RuntimeError(f"P500 contains taxon ids outside P100: {missing[:10]}")
    for r in p500_rows:
        parent = p100_by_id[r["inat_taxon_id"]]
        if parent["species"] != r["species"]:
            raise RuntimeError(f"species identity mismatch for taxon {r['inat_taxon_id']}")
    out = [
        {k: r[k] for k in OUTPUT_COLUMNS}
        for r in p100_rows
        if r["inat_taxon_id"] not in p500_ids
    ]
    out.sort(key=lambda r: (r["inat_taxon_id"], r["species"]))
    _validate_unique(out, "third cohort candidate frame")
    return out


def write_candidate_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=OUTPUT_COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def build(p100_path: Path, p500_path: Path, out_dir: Path, enforce_fingerprints: bool = True) -> dict:
    if enforce_fingerprints:
        p100_sha = sha256(p100_path)
        p500_sha = sha256(p500_path)
        if p100_sha != P100_SHA256:
            raise RuntimeError(f"P100 SHA256 drift: {p100_sha}")
        if p500_sha != P500_SHA256:
            raise RuntimeError(f"P500 SHA256 drift: {p500_sha}")
    else:
        p100_sha = sha256(p100_path)
        p500_sha = sha256(p500_path)

    p100 = read_outcome_blind_source(p100_path)
    p500 = read_outcome_blind_source(p500_path)
    if enforce_fingerprints and len(p100) != P100_SPECIES:
        raise RuntimeError(f"P100 count drift: {len(p100)}")
    if enforce_fingerprints and len(p500) != P500_SELECTED_SPECIES:
        raise RuntimeError(f"P500 count drift: {len(p500)}")

    candidates = derive_candidate_rows(p100, p500)
    if enforce_fingerprints and len(candidates) != THIRD_CANDIDATE_SPECIES:
        raise RuntimeError(f"third candidate count drift: {len(candidates)}")

    out_dir.mkdir(parents=True, exist_ok=True)
    candidate_path = out_dir / "candidate_pool.csv"
    write_candidate_csv(candidate_path, candidates)
    candidate_sha = sha256(candidate_path)
    if enforce_fingerprints and candidate_sha != EXPECTED_CANDIDATE_SHA256:
        raise RuntimeError(f"candidate SHA256 drift: {candidate_sha}")

    result = {
        "analysis": "polymorphism_h2_third_cohort_preopening_candidate_frame",
        "date_jst": "2026-09-16",
        "status": "PREOUTCOME_ELIGIBILITY_FROZEN",
        "biological_outcomes_opened": False,
        "lineage": {
            "U0_species": U0_SPECIES,
            "U100_species": U100_SPECIES,
            "legacy_high_depth_species_excluded_before_P100": LEGACY_EXCLUDED_SPECIES,
            "P100_species": len(p100),
            "P500_selected_species_excluded_here": len(p500),
            "third_candidate_species": len(candidates),
            "parent_selection_workflow_run": PARENT_WORKFLOW_RUN,
            "parent_artifact_id": PARENT_ARTIFACT_ID,
            "parent_artifact_digest": PARENT_ARTIFACT_DIGEST,
        },
        "source_files": {
            "p100": {"sha256": p100_sha, "expected_sha256": P100_SHA256},
            "p500": {"sha256": p500_sha, "expected_sha256": P500_SHA256},
        },
        "candidate_pool": {
            "file": "candidate_pool.csv",
            "columns": OUTPUT_COLUMNS,
            "sha256": candidate_sha,
            "species": len(candidates),
            "minimum_after_observer_cap": min(r["after_observer_cap"] for r in candidates),
            "p500_taxon_overlap": len({r["inat_taxon_id"] for r in candidates} & {r["inat_taxon_id"] for r in p500}),
            "p500_species_overlap": len({r["species"] for r in candidates} & {r["species"] for r in p500}),
        },
        "outcome_firewall": {
            "allowed_source_columns": sorted(ALLOWED_SOURCE_COLUMNS),
            "output_columns": OUTPUT_COLUMNS,
            "flower_colour_opened": False,
            "morph_opened": False,
            "D_opened": False,
            "palette_opened": False,
            "H2_W_opened": False,
            "P500_measurement_table_read": False,
            "P500_recovered_H2_read": False,
            "H3_predictors_read": False,
        },
        "selection_status": "candidate_universe_only_no_third_cohort_sample_selected_yet",
        "next_gate": "freeze a new outcome-blind third-cohort sample from this 3230-species universe, then pass synthetic end-to-end H2 serialization/artifact qualification before biological opening",
    }
    (out_dir / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--p100", type=Path, required=True)
    p.add_argument("--p500", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--no-fingerprint-check", action="store_true")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    result = build(args.p100, args.p500, args.out_dir, enforce_fingerprints=not args.no_fingerprint_check)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

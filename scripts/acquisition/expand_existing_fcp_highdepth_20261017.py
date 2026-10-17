#!/usr/bin/env python3
"""Expand frozen FCP 149,900-photo study with NEW species, using the SAME methods.

Stage 1: SHA-authenticated outcome-blind selection from the original U100 frame.
Stage 2: identical H9 one-page iNaturalist acquisition and species x 100 gate.
No image colours, ROI pixels, genetic interpretation or original manuscript results
are opened or changed by this script.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

SOURCE_COMMIT = "ca69930986e39e3ea9b2d2f97ab247045ba9db0b"
P100_SHA = "1473aad680fe2fa84903c5e11ee104fd0eceef828957f16c2ef1a35f9dd6993c"
P500_SHA = "f54d07fb2338a20a3f92808b58f0ebc948ab0d5af3cb01fb26702f861184b7a4"
CANDIDATE_3230_SHA = "7fc0337074b55a91c0b50773a8d5c3ef82e877074c8b3cacc21f9af58e0ba77e"
THIRD_SELECTED_SHA = "16db6a2fc265f0ab20e7dae4e6f9058b907f5c19af3ddab5b6077797dd1def59"
SALT = "FCP_HIGHDEPTH_WAVE4_20261017_FIXED_SELECTION_V1"
SELECTION_N = 500
PER_SPECIES = 100
EXCLUSION_EARLIER_N = 178_462
THIRD_CANDIDATE_METADATA_N = 49_999
EXCLUSION_TOTAL_N = EXCLUSION_EARLIER_N + THIRD_CANDIDATE_METADATA_N
EXCLUSION_KEYS = ("observation_id", "photo_id")
EXCLUSION_FILES = (
    "random_photo_first_h9_exclusion_ledger_v1.csv",
    "random_photo_first_h9_fresh_metadata_v1.csv",
    "rgfca_42111_species_breadth_measured.csv.gz",
    "polymorphism_h2_p500_candidate_metadata_v1.csv.gz",
    "polymorphism_h2_third_cohort_candidate_metadata_v1.csv.gz",
)
SOURCE_COLUMNS = ("inat_taxon_id", "species", "after_observer_cap")
ALLOWED_SOURCE_COLUMNS = set(SOURCE_COLUMNS) | {"selection_hash", "prospective_rank"}
ALLOWED_LICENSES = ("cc0", "cc-by", "cc-by-sa", "cc-by-nc", "cc-by-nc-sa")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def assert_sha(path: Path, expected: str) -> None:
    actual = sha256(path)
    if actual != expected:
        raise RuntimeError(f"immutable SHA drift for {path}: {actual} != {expected}")


def source_rows(path: Path, expected_sha: str, n: int) -> list[dict]:
    assert_sha(path, expected_sha)
    with path.open(newline="", encoding="utf8") as handle:
        reader = csv.DictReader(handle)
        names = set(reader.fieldnames or ())
        if not set(SOURCE_COLUMNS) <= names or names - ALLOWED_SOURCE_COLUMNS:
            raise RuntimeError("source contains missing/forbidden columns")
        rows = []
        for row in reader:
            item = {
                "inat_taxon_id": int(row["inat_taxon_id"]),
                "species": row["species"].strip(),
                "after_observer_cap": int(row["after_observer_cap"]),
            }
            if not item["species"] or item["after_observer_cap"] < PER_SPECIES:
                raise RuntimeError("historical capacity/identity rule changed")
            rows.append(item)
    if len(rows) != n or len({r["inat_taxon_id"] for r in rows}) != n or len({r["species"] for r in rows}) != n:
        raise RuntimeError("historical unique species/frame size changed")
    return rows


def freeze_selection(p100_path: Path, p500_path: Path, third_path: Path, out_dir: Path) -> dict:
    p100 = source_rows(p100_path, P100_SHA, 3730)
    p500 = source_rows(p500_path, P500_SHA, 500)
    all_p100 = {r["inat_taxon_id"]: r for r in p100}
    by_id = {r["inat_taxon_id"]: r for r in p500}
    for tid, r in by_id.items():
        if tid not in all_p100 or all_p100[tid]["species"] != r["species"]:
            raise RuntimeError("parent selected species are not within original P100")
    candidates = sorted((r for r in p100 if r["inat_taxon_id"] not in by_id), key=lambda r: (r["inat_taxon_id"], r["species"]))
    if len(candidates) != 3230:
        raise RuntimeError("3230 original pre-colour candidates not reproduced")
    # Verify byte identity of the ORIGINAL candidate pool, not merely membership.
    source_blob = "inat_taxon_id,species,after_observer_cap\n" + "".join(
        f'{r["inat_taxon_id"]},{r["species"]},{r["after_observer_cap"]}\n' for r in candidates
    )
    if hashlib.sha256(source_blob.encode("utf8")).hexdigest() != CANDIDATE_3230_SHA:
        raise RuntimeError("frozen original 3230 candidate CSV differs")
    assert_sha(third_path, THIRD_SELECTED_SHA)
    with third_path.open(newline="", encoding="utf8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames != ["prospective_rank", "inat_taxon_id", "species"]:
            raise RuntimeError("original 500 selected manifest schema differs")
        third = list(reader)
    if (len(third) != 500 or [int(r["prospective_rank"]) for r in third] != list(range(1, 501))
            or len({r["inat_taxon_id"] for r in third}) != 500):
        raise RuntimeError("third cohort 500 selected identities have drifted")
    idx = {r["inat_taxon_id"]: r for r in candidates}
    for old in third:
        tid = int(old["inat_taxon_id"])
        if tid not in idx or idx[tid]["species"] != old["species"]:
            raise RuntimeError("third cohort identity drift")
    used = {int(r["inat_taxon_id"]) for r in third}
    available = [r for r in candidates if r["inat_taxon_id"] not in used]
    if len(available) != 2730:
        raise RuntimeError("unallocated 2730 taxa not reproduced")
    ranked = []
    for r in available:
        item = dict(r)
        item["selection_hash"] = hashlib.sha256(
            f'{SALT}|{item["inat_taxon_id"]}|{item["species"]}'.encode("utf8")
        ).hexdigest()
        ranked.append(item)
    ranked.sort(key=lambda r: (r["selection_hash"], r["inat_taxon_id"]))
    selected = ranked[:SELECTION_N]
    for index, item in enumerate(selected, 1):
        item["prospective_rank"] = index
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest = out_dir / "selected_species_manifest.tsv"
    with manifest.open("w", newline="", encoding="utf8") as out:
        writer = csv.DictWriter(out, fieldnames=["prospective_rank", "inat_taxon_id", "species"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({k: row[k] for k in writer.fieldnames} for row in selected)
    richer = out_dir / "selected_species_selection_audit.csv"
    with richer.open("w", newline="", encoding="utf8") as out:
        writer = csv.DictWriter(out, fieldnames=["prospective_rank", "inat_taxon_id", "species", "after_observer_cap", "selection_hash"], lineterminator="\n")
        writer.writeheader()
        writer.writerows({k: row[k] for k in writer.fieldnames} for row in selected)
    record = {
        "schema": "fcp_unchanged_highdepth_wave4_outcome_blind_selection_v1",
        "status": "NEW_500_SPECIES_SELECTED_BEFORE_PHOTO_COLOUR",
        "historical_source_commit": SOURCE_COMMIT,
        "original_species_sampled": 1499,
        "original_photo_rows": 149900,
        "original_U100": 4730,
        "original_legacy_species_excluded": 1000,
        "original_P500_selected_excluded": 500,
        "original_third_selected_excluded": 500,
        "remaining_candidate_species": len(available),
        "wave4_selected_species": SELECTION_N,
        "target_photo_rows_per_species": PER_SPECIES,
        "target_additional_photo_rows": SELECTION_N * PER_SPECIES,
        "selection_salt": SALT,
        "sampling_rule": "exact historical taxon-identity hash ordering, not phenotype/environment/region aware",
        "old_parent_p100_sha256": P100_SHA,
        "old_parent_p500_sha256": P500_SHA,
        "old_3230_candidate_sha256": CANDIDATE_3230_SHA,
        "old_third_selection_sha256": THIRD_SELECTED_SHA,
        "new_selection_manifest_sha256": sha256(manifest),
        "new_selection_audit_sha256": sha256(richer),
        "overlap_with_old_500": 0,
        "overlap_with_old_1000": 0,
        "no_image_pixels_opened": True,
        "no_flower_colour_used": True,
        "original_frozen_results_unchanged": True,
    }
    (out_dir / "selection_result.json").write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    return record


def load_exclusions(paths: list[Path]) -> tuple[set[int], set[int], list[dict]]:
    import pandas as pd

    if len(paths) != len(EXCLUSION_FILES):
        raise RuntimeError("all five old exclusion ledgers are required")
    obs: set[int] = set()
    photo: set[int] = set()
    audit = []
    for path, basename in zip(paths, EXCLUSION_FILES):
        if path.name != basename:
            raise RuntimeError(f"wrong original exclusion ledger: {path}")
        df = pd.read_csv(path, usecols=list(EXCLUSION_KEYS), low_memory=False)
        o = set(pd.to_numeric(df["observation_id"], errors="raise").dropna().astype("int64"))
        p = set(pd.to_numeric(df["photo_id"], errors="raise").dropna().astype("int64"))
        if len(o) != len(df) or len(p) != len(df):
            raise RuntimeError("missing or duplicated original source IDs")
        obs |= o
        photo |= p
        audit.append({"name": basename, "sha256": sha256(path), "rows": len(df), "unique_photo_ids": len(p)})
    if len(obs) != EXCLUSION_TOTAL_N or len(photo) != EXCLUSION_TOTAL_N:
        raise RuntimeError(f"frozen exclusion union changed, obs={len(obs)}, photo={len(photo)}, expected={EXCLUSION_TOTAL_N}")
    return obs, photo, audit


def acquire_metadata(selection_dir: Path, original_exclusions: list[Path], out_dir: Path) -> dict:
    import pandas as pd
    from fcp_pipeline.random_photo_h9_pool import freeze_h9_metadata
    from fcp_pipeline.random_photo_pool import InaturalistObservationClient
    from fcp_pipeline.third_cohort_metadata_gate import (
        metadata_gate_verdict, seal_exact_denominator, validate_fresh_rows, validate_selected_manifest,
    )
    selection = json.loads((selection_dir / "selection_result.json").read_text())
    if selection["status"] != "NEW_500_SPECIES_SELECTED_BEFORE_PHOTO_COLOUR":
        raise RuntimeError("no signed wave4 precolour selection")
    manifest = selection_dir / "selected_species_manifest.tsv"
    assert_sha(manifest, selection["new_selection_manifest_sha256"])
    df = pd.read_csv(manifest, sep="\t")
    validate_selected_manifest(df, expected_n=500)
    exc_obs, exc_photo, exc_audit = load_exclusions(original_exclusions)

    # Original third-cohort acquisition method, arguments and callable unchanged.
    client = InaturalistObservationClient(
        user_agent="zuizui0223-fcp-expanded-highdepth/1.0 (github.com/zuizui0223/fcp)",
        request_interval_seconds=1.05, timeout_seconds=45.0, max_retries=0,
    )
    frozen = freeze_h9_metadata(
        client=client,
        species_frame=df[["inat_taxon_id", "species"]].copy(),
        exclusion_observation_ids=exc_obs, exclusion_photo_ids=exc_photo,
        per_page=200, observer_cap_n=2, fixed_raw_photos=100,
        maximum_positional_accuracy_m=5000, allowed_photo_licenses=ALLOWED_LICENSES,
    )
    species = frozen.species_audit.merge(
        df[["inat_taxon_id", "prospective_rank"]],
        on="inat_taxon_id", how="left", validate="one_to_one",
    ).sort_values("prospective_rank", kind="mergesort").reset_index(drop=True)
    observations = frozen.observations
    if observations.empty:
        raise RuntimeError("no new original photo metadata returned")
    observations = observations.merge(
        df[["inat_taxon_id", "prospective_rank"]],
        on="inat_taxon_id", how="left", validate="many_to_one",
    ).sort_values(["prospective_rank", "h9_selection_order", "observation_id", "photo_id"], kind="mergesort")
    validate_fresh_rows(observations, exc_obs, exc_photo)
    full_species, authorized_rows = seal_exact_denominator(species, observations, target_n=100)
    n_error = int(species["request_error"].fillna("").astype(str).str.len().gt(0).sum())
    verdict, authorized = metadata_gate_verdict(n_error, 500, len(full_species), min_full_species=300, error_ceiling=.05)
    # A metadata pass is NOT permission to run ROI image measurement.
    out_dir.mkdir(parents=True, exist_ok=True)
    out_metadata = out_dir / "wave4_candidate_metadata.csv.gz"
    out_authorized = out_dir / "wave4_authorized_metadata.csv.gz"
    out_audit = out_dir / "wave4_species_audit.csv"
    for table, path in [(observations, out_metadata), (authorized_rows, out_authorized)]:
        table.to_csv(path, index=False, lineterminator="\n", compression={"method": "gzip", "compresslevel": 9, "mtime": 0})
    species.to_csv(out_audit, index=False, lineterminator="\n")
    record = {
        "schema": "fcp_unchanged_highdepth_wave4_original_metadata_v1",
        "status": "METADATA_ONLY_NO_PIXELS",
        "source_commit": SOURCE_COMMIT,
        "original_photos": 149900,
        "new_selected_species": 500,
        "new_exact_100_species": len(full_species),
        "new_exact_100_photo_rows": len(authorized_rows),
        "new_all_candidate_metadata_rows": len(observations),
        "request_errors": n_error,
        "metadata_gate_verdict": verdict,
        "metadata_gate_pass": authorized,
        "photographic_measurement_authorized_by_metadata_alone": False,
        "selection_sha256": sha256(manifest),
        "previous_source_exclusion_photo_ids": len(exc_photo),
        "previous_source_exclusion_observation_ids": len(exc_obs),
        "exclusion_sources": exc_audit,
        "outputs": {x.name: sha256(x) for x in (out_metadata, out_authorized, out_audit)},
        "sampling_method": {
            "one_random_page_per_species": True,
            "quality_grade": "research", "flowering_term_id": 12, "flowering_value_id": 13,
            "rank": "species", "georeferenced": True, "unobscured": True,
            "max_positional_accuracy_m": 5000, "per_page": 200,
            "observer_cap": 2, "geographic_maximin_100": True,
            "licenses": list(ALLOWED_LICENSES), "same_h9_code_commit": SOURCE_COMMIT,
        },
        "original_manuscript_numbers_changed": False,
        "image_pixels_opened": False,
        "no_colour_or_environment_screening": True,
    }
    (out_dir / "metadata_result.json").write_text(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    return record


def main() -> None:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="mode", required=True)
    first = sub.add_parser("select")
    first.add_argument("--p100", type=Path, required=True)
    first.add_argument("--p500", type=Path, required=True)
    first.add_argument("--third-selected", type=Path, required=True)
    first.add_argument("--out", type=Path, required=True)
    second = sub.add_parser("metadata")
    second.add_argument("--selected", type=Path, required=True)
    second.add_argument("--exclusions", nargs=5, type=Path, required=True)
    second.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    record = (freeze_selection(args.p100, args.p500, args.third_selected, args.out)
              if args.mode == "select" else acquire_metadata(args.selected, args.exclusions, args.out))
    print("FCP_WAVE4_" + args.mode.upper() + "_COMPLETE")
    print(json.dumps(record, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()

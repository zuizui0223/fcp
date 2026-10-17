#!/usr/bin/env python3
"""Append independently measured FCP wave4 photos to the original 149,900.

Only adds rows and source/provenance columns. Existing measurements, morph labels
and species-level analysis methods are NEVER recomputed, relabelled or replaced.
The result is an expanded input file suitable for exactly the original analyses.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

OLD_SOURCE = {
    "discovery": {
        "rows": 50000,
        "species": 500,
        "sha256": "ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4",
    },
    "validation": {
        "rows": 50000,
        "species": 500,
        "sha256": "0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6",
    },
    "prospective": {
        "rows": 49900,
        "species": 499,
        "sha256": "57630fc9f281bce94a0c40a70aaf7bce879dde93d6154175adcd021e8f5c1186",
    },
}
FOUR = frozenset({"white", "yellow_orange", "red_pink", "blue_purple"})
REQUIRED = frozenset({"photo_id", "inat_taxon_id", "species", "global_classifiable", "morph"})
TARGET_PHOTOS_PER_SPECIES = 100


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def bool_series(x: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(x):
        return x.fillna(False).astype(bool)
    return x.fillna("").astype(str).str.lower().str.strip().isin({"true", "1", "yes"})


def enforce_and_tag(df: pd.DataFrame, *, label: str, species: int, rows: int, target: int) -> tuple[pd.DataFrame, dict]:
    missing = REQUIRED - set(df.columns)
    if missing:
        raise ValueError(f"{label}: original analysis columns missing: {sorted(missing)}")
    x = df.copy()
    x["photo_id"] = x["photo_id"].astype("string").str.strip()
    x["inat_taxon_id"] = pd.to_numeric(x["inat_taxon_id"], errors="raise").astype("int64")
    x["species"] = x["species"].astype("string").str.strip()
    if x.photo_id.isna().any() or x.photo_id.eq("").any() or x.photo_id.duplicated().any():
        raise ValueError(f"{label}: missing/duplicated original photo IDs")
    if x.inat_taxon_id.isna().any() or x.species.isna().any() or x.species.eq("").any():
        raise ValueError(f"{label}: missing original taxon identities")
    if len(x) != rows:
        raise ValueError(f"{label}: expected {rows} actual original photographs, found {len(x)}")
    if x.inat_taxon_id.nunique() != species or x.species.nunique() != species:
        raise ValueError(f"{label}: original taxon count/identity drift")
    if not x.groupby("inat_taxon_id", sort=True).size().eq(target).all():
        raise ValueError(f"{label}: original 100-photo per-species denomination changed")
    if x[["inat_taxon_id", "species"]].drop_duplicates().groupby("inat_taxon_id").size().max() != 1:
        raise ValueError(f"{label}: inconsistent species names per taxon")
    valid = bool_series(x["global_classifiable"])
    if (valid & ~x["morph"].astype(str).isin(FOUR)).any():
        raise ValueError(f"{label}: four-state original classifier integrity violated")
    x.insert(0, "fcp_source_cohort", label)
    return x, {
        "cohort": label,
        "photo_rows": int(len(x)),
        "named_species": int(x.species.nunique()),
        "four_state_classifiable_rows": int((valid & x.morph.astype(str).isin(FOUR)).sum()),
        "nonclassifiable_rows": int((~valid).sum()),
    }


def combine_old_and_wave4(
    paths: dict[str, Path],
    new_path: Path,
    terminal_path: Path,
    *,
    old_frame_contract: dict = OLD_SOURCE,
    target: int = TARGET_PHOTOS_PER_SPECIES,
) -> tuple[pd.DataFrame, dict]:
    """Pure testable assembly; caller chooses old contract, production uses immutable constants."""
    if set(paths) != set(old_frame_contract):
        raise ValueError("all three separately frozen original cohorts are required")
    out = []
    audits = []
    fingerprints = {}
    seen_photos: set[str] = set()
    seen_taxa: set[int] = set()
    seen_species: set[str] = set()
    for name in ("discovery", "validation", "prospective"):
        path = paths[name]
        contract = old_frame_contract[name]
        fingerprint = sha256_file(path)
        if fingerprint != contract["sha256"]:
            raise ValueError(f"{name}: immutable original source SHA256 mismatch")
        frame, a = enforce_and_tag(pd.read_csv(path, low_memory=False), label=name,
                                   rows=contract["rows"], species=contract["species"], target=target)
        _check_separately_sampled(frame, seen_photos, seen_taxa, seen_species, name)
        out.append(frame)
        audits.append(a)
        fingerprints[name] = fingerprint
    old_rows = sum(a["photo_rows"] for a in audits)
    if old_frame_contract is OLD_SOURCE and old_rows != 149900:
        raise ValueError("original high-depth 149900-photo denominator changed")
    receipt = json.loads(terminal_path.read_text())
    new_sha = sha256_file(new_path)
    if (receipt.get("schema") != "fcp_same_149900_method_wave4_real_photo_terminals_v1"
        or receipt.get("original_results_changed") is not False
        or receipt.get("measured_sha256") != new_sha
        or receipt.get("terminal_partitions") != 256
        or receipt.get("measurement_original_fixed_colour_worker") is not True):
        raise ValueError("extra original photo measurements lack completed unchanged-method terminal")
    n = int(receipt["wave4_100photo_species"])
    num = int(receipt["wave4_new_photo_ids"])
    if not 300 <= n <= 500 or num != n * target:
        raise ValueError("extra wave original 100-photo sampling denominator changed")
    new, audit = enforce_and_tag(pd.read_csv(new_path, low_memory=False), label="wave4",
                                 rows=num, species=n, target=target)
    if int(audit["four_state_classifiable_rows"]) != receipt["wave4_classifiable_photo_rows"]:
        raise ValueError("extra original classified photo row count changed")
    _check_separately_sampled(new, seen_photos, seen_taxa, seen_species, "wave4")
    out.append(new)
    audits.append(audit)
    fingerprints["wave4"] = new_sha
    combined = pd.concat(out, axis=0, ignore_index=True, sort=False)
    expected = old_rows + num
    if (len(combined) != expected or combined.photo_id.nunique() != expected
        or combined.inat_taxon_id.nunique() != sum(a["named_species"] for a in audits)):
        raise ValueError("expanded study unexpectedly duplicated or lost original photos")
    result = {
        "schema": "fcp_original_149900_plus_same_method_wave4_merged_v1",
        "append_only": True,
        "original_study_unchanged": True,
        "original_source_rows": old_rows,
        "new_original_photo_rows": num,
        "all_original_photo_rows": expected,
        "original_named_species": sum(a["named_species"] for a in audits[:-1]),
        "new_named_species": n,
        "all_named_species": int(combined.species.nunique()),
        "cohorts": audits,
        "input_sha256": fingerprints,
        "only_new_column": "fcp_source_cohort",
        "four_state_method_relabelled": False,
        "existing_50km_null_reweighted": False,
        "external_iNaturalist_source_independent": False,
    }
    return combined, result


def _check_separately_sampled(x: pd.DataFrame, old_photos: set[str], old_taxa: set[int],
                               old_species: set[str], label: str) -> None:
    photos = set(x.photo_id.astype(str))
    taxa = set(x.inat_taxon_id.astype(int))
    species = set(x.species.astype(str))
    if photos & old_photos or taxa & old_taxa or species & old_species:
        raise ValueError(f"{label}: original-photo or species identity reused across cohorts")
    old_photos.update(photos)
    old_taxa.update(taxa)
    old_species.update(species)


def main() -> None:
    ap = argparse.ArgumentParser()
    for name in OLD_SOURCE:
        ap.add_argument("--" + name, required=True, type=Path)
    ap.add_argument("--wave4", required=True, type=Path)
    ap.add_argument("--wave4-terminal", required=True, type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    args = ap.parse_args()
    original = {name: getattr(args, name) for name in OLD_SOURCE}
    combined, receipt = combine_old_and_wave4(original, args.wave4, args.wave4_terminal)
    args.outdir.mkdir(parents=True, exist_ok=True)
    out = args.outdir / "fcp_original149900_plus_wave4_same_method_photos.csv.gz"
    combined.to_csv(out, index=False, lineterminator="\n",
                    compression={"method": "gzip", "compresslevel": 9, "mtime": 0})
    receipt["merged_csv_sha256"] = sha256_file(out)
    (args.outdir / "expanded_data_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    print("FCP_149900_APPEND_ONLY_SAME_METHOD_COMPLETE", json.dumps({
        k: receipt[k] for k in ("original_source_rows", "new_original_photo_rows",
                               "all_original_photo_rows", "all_named_species")
    }, sort_keys=True))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Repeat the EXACT FCP 50-km spatial analysis on 149,900 old + genuine new photos.

No ecological estimator, flower ROI, palette, radius, sample gate or null design
changes. This merely applies the existing historical species-conditioned
statistics to a larger, photo-ID-disjoint pool, preserving the original results
and showing old, new and combined outputs side by side.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis.run_polymorphism_distributed_polymorphism_posthoc_20261007 import (
    MORPHS,
    RADII_KM,
    PERMUTATIONS,
    MIN_LOCAL_PAIRS,
    PRIMARY_RADIUS,
    analyse_species,
    as_bool,
    cohort_summary,
)

SOURCE_KEY = "fcp_source_cohort"
GROUPS = ("discovery", "validation", "prospective", "wave4")
ORIGINAL_EXPECTED = {
    "discovery": (166, 0.020529254583812922, 0.005),
    "validation": (181, 0.018672971642749295, 0.005),
    "prospective": (204, 0.01468491968437185, 0.005),
}
FROZEN_ORIGINAL_P = 0.005
N_ORIGINAL_PHOTOS = 149900


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def import_exact_measured_rows(merged: Path, receipt: dict) -> pd.DataFrame:
    if receipt.get("schema") != "fcp_original_149900_plus_same_method_wave4_merged_v1":
        raise ValueError("not the verified 149900 plus original-method wave4 input")
    if receipt.get("original_study_unchanged") is not True or receipt.get("append_only") is not True:
        raise ValueError("the original 149900 must not be changed")
    if receipt.get("four_state_method_relabelled") is not False or receipt.get("existing_50km_null_reweighted") is not False:
        raise ValueError("a scientific method change was detected")
    if sha256(merged) != receipt.get("merged_csv_sha256"):
        raise ValueError("expanded original photo CSV checksum changed")
    old_sha = {
        "discovery": "ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4",
        "validation": "0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6",
        "prospective": "57630fc9f281bce94a0c40a70aaf7bce879dde93d6154175adcd021e8f5c1186",
    }
    if any(receipt.get("input_sha256", {}).get(k) != v for k, v in old_sha.items()):
        raise ValueError("the prior immutable 149900 source fingerprints changed")
    df = pd.read_csv(merged, low_memory=False)
    if len(df) != receipt.get("all_original_photo_rows") or df.photo_id.nunique() != len(df):
        raise ValueError("merged photo ID denominator differs from verified append")
    if df[SOURCE_KEY].isna().any() or set(df[SOURCE_KEY].astype(str)) != set(GROUPS):
        raise ValueError("source cohort identity is missing or changed")
    if len(df) - N_ORIGINAL_PHOTOS != receipt["new_original_photo_rows"]:
        raise ValueError("new photo rows are not a pure extension")
    n = df.groupby(SOURCE_KEY).size()
    if any(int(n[k]) != expected for k, expected in (
        ("discovery", 50000), ("validation", 50000),
        ("prospective", 49900), ("wave4", int(receipt["new_original_photo_rows"])),
    )):
        raise ValueError("old or new cohort sample sizes differ from frozen inputs")
    # Old loader also sorts photo IDs numerically; coercion is harmless only
    # because iNaturalist photo IDs are already unique unsigned integers.
    numeric_photo_ids = pd.to_numeric(df["photo_id"], errors="raise")
    if numeric_photo_ids.isna().any() or np.any(numeric_photo_ids < 0):
        raise ValueError("photo identifier cannot be matched to historical source")
    df["photo_id"] = numeric_photo_ids.astype(np.int64)
    if df["inat_taxon_id"].isna().any():
        raise ValueError("taxon identity missing")
    return df


def original_eligible(df: pd.DataFrame) -> pd.DataFrame:
    """Exactly original load() filtering, minus its fixed-initial-file SHA gate."""
    keep = as_bool(df["global_classifiable"]) & df["morph"].astype(str).isin(MORPHS)
    d = df.loc[keep].copy()
    counts = d.groupby("inat_taxon_id").size()
    d = d.loc[d["inat_taxon_id"].isin(counts[counts >= 40].index)].copy()
    return d.sort_values(["inat_taxon_id", "photo_id"], kind="stable").reset_index(drop=True)


def run_four_cohorts(df: pd.DataFrame) -> tuple[dict, pd.DataFrame]:
    if set(df[SOURCE_KEY].dropna().astype(str)) != set(GROUPS):
        raise ValueError("four original method cohorts required")
    all_rows: list[dict] = []
    all_nulls = {r: [] for r in RADII_KM}
    per_cohort: dict[str, dict] = {}
    for cohort in GROUPS:
        d = original_eligible(df.loc[df[SOURCE_KEY].eq(cohort)])
        # Historical third-cohort original seed namespace is preserved for
        # the prospective group. New wave4 uses the same frozen algorithm.
        old_seed_label = "third" if cohort in ("prospective", "wave4") else cohort
        rows: list[dict] = []
        nulls: dict[float, list[np.ndarray]] = {r: [] for r in RADII_KM}
        for _, g in d.groupby("inat_taxon_id", sort=True):
            g = g.dropna(subset=["latitude", "longitude"])
            if len(g) < 40:
                continue
            row, per = analyse_species(g, old_seed_label)
            row["report_cohort"] = cohort
            rows.append(row)
            for radius in RADII_KM:
                nulls[radius].append(per[radius])
                all_nulls[radius].append(per[radius])
        if not rows:
            raise ValueError(f"{cohort} has no eligible original FCP photographic taxa")
        species = pd.DataFrame(rows)
        all_rows.extend(rows)
        per_cohort[cohort] = {
            "n_measured_photos": int(df[SOURCE_KEY].eq(cohort).sum()),
            "n_classifiable_species": int(d.inat_taxon_id.nunique()),
            "radii": {str(int(rad)): result for rad, result in cohort_summary(species, nulls).items()},
        }
    combined = pd.DataFrame(all_rows)
    total = cohort_summary(combined, all_nulls)
    return {
        "cohorts": per_cohort,
        "pooled": {str(int(rad)): result for rad, result in total.items()},
        "pooled_n_spatial_eligible_species": int(combined.inat_taxon_id.nunique()),
    }, combined


def require_reproduced_original_three(result: dict) -> None:
    for group, (n, dep, p) in ORIGINAL_EXPECTED.items():
        actual = result["cohorts"][group]["radii"]["50"]
        if (
            actual["n_species"] != n
            or abs(actual["mean_depletion"] - dep) > 1e-10
            or actual["mean_depletion_p_upper"] != p
        ):
            raise ValueError(
                f"FROZEN_ORIGINAL_50KM_REPLAY_MISMATCH:{group}:"
                f"n={actual['n_species']},dep={actual['mean_depletion']},p={actual['mean_depletion_p_upper']}"
            )


def verify_existing_wave4_inference(result: dict, receipt: dict, original_wave4_result: dict) -> dict:
    if (
        original_wave4_result.get("schema") != "fcp_same_149900_method_wave4_statistical_replication_v1"
        or original_wave4_result.get("original_frozen_results_unchanged") is not True
        or original_wave4_result.get("wave4_measured_table_sha256") != receipt["input_sha256"]["wave4"]
    ):
        raise ValueError("extra same-method original photo inference receipt not authenticated")
    historic = original_wave4_result.get("spatial_original_method", {}).get("50", {})
    new = result["cohorts"]["wave4"]["radii"]["50"]
    if (
        int(historic.get("n_species", -1)) != new["n_species"]
        or abs(float(historic.get("mean_depletion", np.nan)) - new["mean_depletion"]) > 1e-10
        or historic.get("mean_depletion_p_upper") != new["mean_depletion_p_upper"]
    ):
        raise ValueError("added original method 50km effect changed during 149900+ append/replay")
    return {
        "new_wave4_original_50km_verified": True,
        "original_W_H1_methods_kept_from_wave4": True,
        "wave4_h1": original_wave4_result.get("h1_original_method"),
        "wave4_H2_W": original_wave4_result.get("h2_original_W_method"),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--combined", required=True, type=Path)
    ap.add_argument("--receipt", required=True, type=Path)
    ap.add_argument("--wave4-original-inference", required=True, type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    args = ap.parse_args()
    rc = json.loads(args.receipt.read_text())
    combined = import_exact_measured_rows(args.combined, rc)
    stat, species = run_four_cohorts(combined)
    require_reproduced_original_three(stat)
    new_stat = json.loads(args.wave4_original_inference.read_text())
    linked = verify_existing_wave4_inference(stat, rc, new_stat)
    out = {
        "schema": "fcp_149900_plus_wave4_original_50km_replayed_all_species_v1",
        "role": "same_spatial_statistic_more_genuine_original_photos",
        "frozen_original_source_rows": 149900,
        "new_photo_rows": int(rc["new_original_photo_rows"]),
        "total_photo_rows": len(combined),
        "original_method": {
            "radii_km": RADII_KM, "primary_50km": PRIMARY_RADIUS,
            "min_local_pairs": MIN_LOCAL_PAIRS, "permutations": PERMUTATIONS,
            "vertex_colour_multiset_preserved": True,
            "equal_species_weights": True,
        },
        "all_original_frozen_cohort_50km_values_reproduced": True,
        "baseline_W_verdict_not_changed": True,
        "not_independent_photo_provider_replication": True,
        "posthoc_primary_50km_remains_posthoc": True,
        "results": stat,
        "wave4_linked_old_H1_W": linked,
        "source_sha256": {
            "merged": rc["merged_csv_sha256"],
            "append_receipt": sha256(args.receipt),
            "wave4_original_method_inference": sha256(args.wave4_original_inference),
        },
        "hard_nonclaims": [
            "no different spatial estimator or tuned distance parameter",
            "no new adaptive, genetic or selection claim",
            "this is phenotypic photographic variation, not botanical population frequency",
        ],
    }
    args.outdir.mkdir(parents=True, exist_ok=True)
    species.to_csv(args.outdir / "original_149900_plus_wave4_species_50km.csv", index=False)
    path = args.outdir / "expanded_fcp_original_50km_result.json"
    path.write_text(json.dumps(out, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("FCP_EXPANDED_ALL_ORIGINAL_PHOTOS_SAME_50KM_ANALYSIS_COMPLETE")
    print(json.dumps({
        "old_three_reproduced": True,
        "total_photos": len(combined),
        "new_photos": out["new_photo_rows"],
        "original_50km_n": {c: stat["cohorts"][c]["radii"]["50"]["n_species"] for c in GROUPS},
        "pooled_50km": stat["pooled"]["50"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()

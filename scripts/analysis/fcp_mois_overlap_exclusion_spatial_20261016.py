#!/usr/bin/env python3
"""Post hoc FCP 50-km spatial analysis after excluding MOIS-shared original photo IDs.

Both original photographic archives must pass immutable source SHA checks. This
does not change any frozen FCP inference, and shared measurement-system ancestry
remains even when exact photo IDs are disjoint.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(SCRIPTS / "analysis"))
from audit_fcp_mois_photo_id_overlap import get_mois, get_fcp  # noqa: E402
from run_polymorphism_distributed_polymorphism_posthoc_20261007 import (  # noqa: E402
    RADII_KM, analyse_species, cohort_summary, load,
)

COHORTS = {"discovery": "discovery", "validation": "validation", "prospective": "third"}
ORIGINAL_FROZEN = {
    "discovery": (166, 0.020529254583812922),
    "validation": (181, 0.018672971642749295),
    "prospective": (204, 0.01468491968437185),
}
PHOTO_INTERSECTIONS = {"discovery": 1248, "validation": 1326, "prospective": 837}
PRIMARY = "50"


def run_on_df(df: pd.DataFrame, original_cohort: str) -> tuple[pd.DataFrame, dict]:
    rows = []
    by_radius = {radius: [] for radius in RADII_KM}
    for _, g in df.groupby("inat_taxon_id", sort=True):
        # After exclusion, originally eligible taxa may lose the 40-photo gate.
        # Never count a newly unevaluable taxon as a biological single-morph species.
        g = g.dropna(subset=["latitude", "longitude"])
        if len(g) < 40:
            continue
        row, nulls = analyse_species(g, original_cohort)
        rows.append(row)
        for r in RADII_KM:
            by_radius[r].append(nulls[r])
    if not rows:
        raise ValueError("No taxa with valid original photographic spatial support")
    species = pd.DataFrame(rows)
    res = cohort_summary(species, by_radius)[PRIMARY]
    if not res["evaluable"]:
        raise ValueError("No valid 50-km conspecific local colour pairs")
    return species, res


def summarize_valid_50km_species(original_species: pd.DataFrame,
                                 reduced_species: pd.DataFrame, cohort: str) -> dict:
    """Paired effects use species *evaluable at 50 km in both analyses*.

    Some >=40-photo species have <30 local pairs, hence NaN depletion. They
    must not enter paired means/bootstraps, even when D itself is defined.
    """
    a = original_species.set_index("inat_taxon_id")
    b = reduced_species.set_index("inat_taxon_id")
    same = a.index.intersection(b.index)
    before = a.loc[same]
    after = b.loc[same]
    good = np.isfinite(before["depletion_50km"].to_numpy(float)) & np.isfinite(
        after["depletion_50km"].to_numpy(float))
    before = before.iloc[np.flatnonzero(good)]
    after = after.iloc[np.flatnonzero(good)]
    if not len(before):
        raise ValueError("No matched species with finite original and excluded 50-km depletion")
    delta = (after["depletion_50km"] - before["depletion_50km"]).to_numpy(float)
    rng = np.random.default_rng(2026101668 + list(COHORTS).index(cohort))
    sims = delta[rng.integers(0, len(delta), size=(999, len(delta)))].mean(axis=1)
    D0, D1 = before["D_pair"], after["D_pair"]
    out = {
        "n_species_evaluable_both": int(len(before)),
        "original_depletion_mean": float(before["depletion_50km"].mean()),
        "excluded_depletion_mean": float(after["depletion_50km"].mean()),
        "mean_depletion_change": float(delta.mean()),
        "species_bootstrap_95_change_interval": [float(v) for v in np.quantile(sims, [.025, .975])],
        "median_absolute_species_D_change": float(np.median(np.abs(D1 - D0))),
        "spearman_D_before_vs_after": float(D0.corr(D1, method="spearman")),
    }
    if not np.isfinite([out["mean_depletion_change"], *out["species_bootstrap_95_change_interval"],
                        out["median_absolute_species_D_change"]]).all():
        raise ValueError("Nonfinite paired effect or CI")
    return out


def compare(df: pd.DataFrame, cohort: str, mois_photo_ids: set[str]) -> dict:
    original_label = COHORTS[cohort]
    original_species, full = run_on_df(df, original_label)
    expected_n, expected_mean = ORIGINAL_FROZEN[cohort]
    if full["n_species"] != expected_n or abs(full["mean_depletion"] - expected_mean) > 1e-10:
        raise ValueError(f"Frozen original 50-km result could not be reproduced: {cohort}")
    shared = df["photo_id"].astype(str).isin(mois_photo_ids)
    truncated = df.loc[~shared].copy()
    filtered_species, reduced = run_on_df(truncated, original_label)
    matched = summarize_valid_50km_species(original_species, filtered_species, cohort)
    if matched["n_species_evaluable_both"] != reduced["n_species"]:
        raise ValueError("Exclusion must not create new 50-km-evaluable taxa")
    return {
        "frozen_before_exclusion": full,
        "excluded_classifiable_fcp_photo_rows": int(shared.sum()),
        "remaining_50km_species_count": int(reduced["n_species"]),
        "after_overlap_exclusion": reduced,
        "matched_taxa": matched,
        "caveat": "Original D/50km null analysis held fixed; this is a newly performed post hoc photo-overlap-exclusion sensitivity, not prospective biological confirmation.",
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mois-zip", required=True, type=Path)
    for cohort in COHORTS:
        ap.add_argument("--" + cohort, required=True, type=Path)
    ap.add_argument("--overlap-receipt", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()
    known = json.loads(args.overlap_receipt.read_text(encoding="utf-8"))
    assert known["measured_overlap"]["shared_photo_ids"] == 3411
    mois, mm = get_mois(args.mois_zip)
    mois_ids = set(mois)
    results = {}
    fcp_hashes = {}
    for name, model_name in COHORTS.items():
        path = getattr(args, name)
        # get_fcp validates all original photo IDs and source bytes before
        # cohort eligibility screening; no colours are used to choose exclusion.
        all_photos, meta = get_fcp(path, name)
        intersection = len(set(all_photos) & mois_ids)
        if intersection != PHOTO_INTERSECTIONS[name]:
            raise ValueError(f"Photo-ID overlap changed before colour filtering: {name}")
        if meta["source_sha256"] != known["provenance"]["fcp_csv_sha256"][name]:
            raise ValueError(f"Frozen input fingerprint differs from audit: {name}")
        df = load(path, model_name)  # original paper's unmodified eligibility rules
        results[name] = compare(df, name, mois_ids)
        results[name]["original_all_photo_ID_intersection"] = intersection
        fcp_hashes[name] = meta["source_sha256"]
    result = {
        "schema": "fcp_mois_photoid_disjoint_spatial_sensitivity_v1",
        "role": "posthoc_fcp_50km_source_reuse_falsification",
        "historical_FCP_inference_changed": False,
        "independent_photo_provider_replication": False,
        "original_mois_archive_sha256": mm["archive_sha256"],
        "original_fcp_csv_sha256": fcp_hashes,
        "excluded_source_rule": "exclude original photo IDs belonging to MOIS 100543 originals without using phenotype, place or season",
        "original_cross_paper_shared_photo_ids": 3411,
        "cohorts": results,
        "scientific_ceiling": "Assess robustness of within-species photographed local colour depletion to literal MOIS photo reuse only. Not a climate-selection or genetically maintained polymorphism test.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("FCP_MOIS_DISJOINT_50KM_COMPLETE",
          json.dumps({name: {"n_species": res["remaining_50km_species_count"],
                             "mean_depletion": res["after_overlap_exclusion"]["mean_depletion"],
                             "p": res["after_overlap_exclusion"]["mean_depletion_p_upper"]}
                      for name, res in results.items()}, sort_keys=True))


if __name__ == "__main__":
    main()

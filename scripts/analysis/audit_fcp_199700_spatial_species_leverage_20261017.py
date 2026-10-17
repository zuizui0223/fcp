#!/usr/bin/env python3
"""Post hoc species-leverage audit of the permanently released FCP 199,700 real photos.

Does NOT rerun nulls or improve a p value. Uses the archived, unchanged per-species
50-km depletion and checks exact historical cohort/pooled means before describing
how broadly the positive mean is distributed. Four cohorts, NOT pending fifth.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

RELEASE = "fcp-real-original-199700-wave4-50km-H1-W-20261017"
ARCHIVE_SHA = "844c255f28406b5158a8ac3e27d94b32aad48d785efd3c4ae512fb30ee1871ee"
RECEIPT_SHA = "26c93dc7226bd0b2a96edbd626d143f26a91fc9cf1ed4b78e28bf276dbf78eeb"
COHORTS = ("discovery", "validation", "prospective", "wave4")
EXPECTED = {
    "discovery": (166, 0.020529254583812922, 0.005),
    "validation": (181, 0.018672971642749295, 0.005),
    "prospective": (204, 0.01468491968437185, 0.005),
    "wave4": (190, 0.017293596294403593, 0.005),
    "pooled": (741, 0.01763720788097046, 0.005),
}


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fd:
        for data in iter(lambda: fd.read(1 << 20), b""):
            h.update(data)
    return h.hexdigest()


def demand(ok: bool, why: str) -> None:
    if not ok:
        raise ValueError("FCP_SPATIAL_LEVERAGE_NOT_AUTHENTICATED: " + why)


def summarise(series: pd.Series, taxon: pd.Series) -> dict:
    """Independent species = descriptive units, never photo-pair pseudoreplicates."""
    effects = np.asarray(series, dtype=float)
    IDs = np.asarray(taxon, dtype=np.int64)
    demand(len(effects) > 0 and len(np.unique(IDs)) == len(effects),
           "no eligible taxa or repeated taxon identity")
    demand(np.all(np.isfinite(effects)) and np.all(np.abs(effects) <= 1.0),
           "invalid per-species depletion")
    n = len(effects)
    order = np.lexsort((IDs, -effects))
    leave_top = {}
    for prop in (0.05, 0.10, 0.20):
        cutoff = max(1, int(math.ceil(n * prop)))
        remaining = effects[order[cutoff:]]
        leave_top[f"{int(prop*100)}pct"] = {
            "n_removed": cutoff,
            "remaining_species": int(len(remaining)),
            "remaining_mean_depletion": float(remaining.mean()) if len(remaining) else None,
        }
    return {
        "n_species": n,
        "n_positive": int((effects > 0).sum()),
        "n_zero": int((effects == 0).sum()),
        "n_negative": int((effects < 0).sum()),
        "fraction_positive": float(np.mean(effects > 0)),
        "fraction_zero": float(np.mean(effects == 0)),
        "mean_depletion": float(effects.mean()),
        "median_depletion": float(np.median(effects)),
        "q25_depletion": float(np.quantile(effects, 0.25)),
        "q75_depletion": float(np.quantile(effects, 0.75)),
        "leave_largest_positive_effects_out": leave_top,
    }


def analyse(species: pd.DataFrame, statistical: dict, *, expected: dict = EXPECTED) -> dict:
    names = {
        "report_cohort", "inat_taxon_id", "n_classifiable", "local_pairs_50km",
        "D_pair", "D_local_50km", "depletion_50km"
    }
    demand(names.issubset(species.columns), "archived old per-species source schema changed")
    demand(species["inat_taxon_id"].notna().all()
           and not species["inat_taxon_id"].duplicated().any(),
           "same taxon repeated across four source cohorts")
    demand(set(species.report_cohort.astype(str)) == set(COHORTS),
           "missing or added cohort identities")
    demand(len(species) == statistical["results"]["pooled_n_spatial_eligible_species"],
           "original complete per-species source count changed")
    rows = {}
    eligible = []
    for label in COHORTS:
        d = species.loc[species.report_cohort == label].copy()
        x = d.loc[pd.to_numeric(d.depletion_50km, errors="coerce").notna()].copy()
        demand((x.n_classifiable >= 40).all() and (x.local_pairs_50km >= 30).all(),
               f"{label}: photo or 50km local-pair gate drift")
        demand(np.allclose(x.D_pair - x.D_local_50km, x.depletion_50km, atol=1e-10),
               f"{label}: original Dpair-Dlocal contrast identity changed")
        result = summarise(x.depletion_50km, x.inat_taxon_id)
        n, effect, p = expected[label]
        old = statistical["results"]["cohorts"][label]["radii"]["50"]
        demand(result["n_species"] == n and old["n_species"] == n
               and abs(result["mean_depletion"] - effect) < 1e-10
               and abs(old["mean_depletion"] - effect) < 1e-10
               and old["mean_depletion_p_upper"] == p,
               f"{label}: archived/frozen 50km original statistic inconsistent")
        result["frozen_original_50km_upper_p"] = p
        rows[label] = result
        eligible.append(x)
    pooled = pd.concat(eligible, ignore_index=True)
    summary = summarise(pooled.depletion_50km, pooled.inat_taxon_id)
    n, effect, p = expected["pooled"]
    original = statistical["results"]["pooled"]["50"]
    demand(summary["n_species"] == n and original["n_species"] == n
           and abs(summary["mean_depletion"] - effect) < 1e-10
           and abs(original["mean_depletion"] - effect) < 1e-10
           and original["mean_depletion_p_upper"] == p
           and sum(rows[k]["n_species"] for k in COHORTS) == n,
           "pooled old four-cohort effect/eligibility changed")
    summary["frozen_original_50km_upper_p"] = p
    return {
        "schema": "fcp_199700_archived_50km_species_leverage_audit_v1",
        "status": "descriptive_exploratory_post_hoc_only",
        "real_measured_photos": 199700,
        "four_original_source_species": 1997,
        "spatially_eligible_50km_species": n,
        "primary_radius_km": 50,
        "radius_is_post_hoc": True,
        "reported_original_vertex_null_replicates": 199,
        "independent_species_units_not_photo_pairs": True,
        "not_new_statistical_confirmation": True,
        "no_independent_photo_provider": True,
        "not_a_population_morph_frequency_census": True,
        "no_causal_local_adaptation_inference": True,
        "pooled": summary,
        "cohorts": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--external-receipt", type=Path, required=True)
    parser.add_argument("--extracted-root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    a = parser.parse_args()
    demand(all(p.is_file() for p in (a.archive, a.external_receipt)), "release files missing")
    demand(sha(a.archive) == ARCHIVE_SHA and sha(a.external_receipt) == RECEIPT_SHA,
           "immutable public 199700 release bytes do not match pinned published digest")
    receipt = json.loads(a.external_receipt.read_text())
    demand(receipt.get("archive_sha256") == ARCHIVE_SHA
           and receipt.get("archive_bytes") == a.archive.stat().st_size
           and receipt.get("all_original_measured_photo_rows") == 199700
           and receipt.get("all_original_100photo_taxa") == 1997
           and receipt.get("pooled_original_50km_evaluable_species") == 741
           and receipt.get("original_50km_radius_posthoc") is True
           and receipt.get("original_199_permutation_50km_null_unchanged") is True
           and receipt.get("pre_reporting_fix_50km_zero_field_is_not_the_actual_spatial_result") is True,
           "published original four-cohort empirical receipt inconsistent")
    species_path = a.extracted_root / "combined/original_149900_plus_wave4_species_50km.csv"
    result_path = a.extracted_root / "combined/expanded_fcp_original_50km_result.json"
    for p in (species_path, result_path):
        demand(p.is_file(), "archived real per-species or original null result missing")
        rel = str(p.relative_to(a.extracted_root))
        demand(receipt.get("files_sha256", {}).get(rel) == sha(p),
               f"{rel}: per-file digest failed")
    statistical = json.loads(result_path.read_text())
    demand(statistical.get("schema") == "fcp_149900_plus_wave4_original_50km_replayed_all_species_v1"
           and statistical.get("all_original_frozen_cohort_50km_values_reproduced") is True
           and statistical.get("total_photo_rows") == 199700,
           "archived original source null output incompatible")
    result = analyse(pd.read_csv(species_path), statistical)
    result["public_release_tag"] = RELEASE
    result["published_archive_sha256"] = ARCHIVE_SHA
    result["published_receipt_sha256"] = RECEIPT_SHA
    result["per_species_table_sha256"] = sha(species_path)
    result["original_four_cohort_result_sha256"] = sha(result_path)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("FCP_199700_AUTHENTIC_REAL_SPECIES_LEVERAGE", json.dumps({
        "pooled": result["pooled"],
        "cohorts": result["cohorts"],
        "source_sha": result["published_archive_sha256"],
    }, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()

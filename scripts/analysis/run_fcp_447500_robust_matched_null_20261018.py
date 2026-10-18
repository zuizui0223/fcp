#!/usr/bin/env python3
"""Post hoc distributional robustness of the authentic five-cohort 50-km result.

Replays the original vertex null once per species, never changes the classifier,
spatial gate, labels, radius, seeds or original confirmatory decisions. All
source identities, original means/nulls and the complete old four-cohort robust
result must match before any result files are written.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis.run_fcp_199700_robust_matched_null_20261018 import statistics
from scripts.analysis.run_polymorphism_distributed_polymorphism_posthoc_20261007 import (
    MORPHS, PERMUTATIONS, PRIMARY_RADIUS, MIN_LOCAL_PAIRS, analyse_species, as_bool,
)

COHORTS = ("discovery", "validation", "prospective", "wave4", "new_U100")
SEEDS = dict(zip(COHORTS, ("discovery", "validation", "third", "third", "third")))
COUNTS = dict(zip(COHORTS, (500, 500, 499, 498, 2478)))
CLASSIFIED = dict(zip(COHORTS, (25377, 24885, 25788, 24585, 126472)))
LOCAL_COUNTS = dict(zip(COHORTS, (166, 181, 204, 190, 1020)))
ARCHIVE_SHA = "99f7ee68c3fd009e8ab60fd0cfcdf1fb407e3a8de07edd38341b27889310bc25"
HISTORICAL_SHA = "543e4a544310086b800d31e410eed176240c88eb7ce2544a25eaef324488547c"
INPUTS = {
    "photos": ("measured/fcp_all4730_original_method_measured_photos.csv.gz",
               "8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6"),
    "receipt": ("measured/fcp_all4730_original_photo_combine_receipt.json",
                "63403d03222f5485ab3dc935bc08081742c4c984666b7ae85e24d78800f1ac3c"),
    "species": ("spatial/fcp_all4730_unchanged_original_50km_species.csv",
                "960e8c4bae8a9ecb01fefbef06b6b88cd8fdc63d86af3e4c354cab4eb1aedd35"),
    "spatial_result": ("spatial/fcp_all4730_unchanged_original_50km_results.json",
                       "79aeab0555ec2efd865606a9bdc3e06574a5168c38029ee16e22b945e1a9eacf"),
}
ROBUST = ("median", "fraction_positive", "drop_top_10pct", "drop_top_20pct")
TOL = 1e-10


def demand(ok, message):
    if not ok:
        raise ValueError("FCP_447500_ROBUST_SOURCE_OR_NULL_INVALID: " + message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def close(a, b, context):
    demand(np.isfinite(a) and np.isfinite(b) and abs(float(a) - float(b)) <= TOL,
           f"{context}: numeric mismatch {a!r} != {b!r}")


def canonical_photos(photo):
    """Keep the original stable taxon/photo order, with no ambiguous identity."""
    demand({"inat_taxon_id", "photo_id", "fcp_source_cohort"}.issubset(photo.columns),
           "photo identity columns missing")
    demand(not photo[["inat_taxon_id", "photo_id", "fcp_source_cohort"]].isna().any().any(),
           "missing photo/taxon/cohort identity")
    demand(photo.photo_id.is_unique, "photo IDs are not globally unique")
    demand(set(photo.fcp_source_cohort).issubset(COHORTS), "unknown photo cohort")
    demand((photo.groupby("inat_taxon_id").fcp_source_cohort.nunique() == 1).all(),
           "taxon spans multiple cohorts")
    return photo.sort_values(["inat_taxon_id", "photo_id"], kind="stable").reset_index(drop=True)


def verify_species(row, old, context):
    for key in ("n_classifiable", "D_pair", "local_pairs_50km", "D_local_50km",
                "depletion_50km", "null_mean_depletion_50km"):
        a, b = row.get(key, np.nan), old.get(key, np.nan)
        if pd.isna(a) and pd.isna(b):
            continue
        close(a, b, f"{context}/{key}")


def verify_mean(values, nulls, frozen, context):
    stat = statistics(values, nulls)["mean"]
    close(frozen["n_species"], len(values), context + ": n_species changed")
    for key, old_key in (("observed", "mean_depletion"),
                         ("null_mean", "mean_depletion_null_mean"),
                         ("null_q025", "mean_depletion_null_q025"),
                         ("null_q975", "mean_depletion_null_q975"),
                         ("upper_p_exploratory", "mean_depletion_p_upper")):
        close(stat[key], frozen[old_key], context + "/" + old_key)


def compare_tree(actual, expected, path="historical"):
    """Check every archived numeric/statistical field, including raw and Holm p."""
    if isinstance(expected, dict):
        demand(isinstance(actual, dict) and set(actual) == set(expected), path + ": keys differ")
        for key, value in expected.items():
            compare_tree(actual[key], value, path + "/" + key)
    elif isinstance(expected, bool) or expected is None or isinstance(expected, str):
        demand(actual == expected, path + ": identity differs")
    elif isinstance(expected, (int, float)):
        close(actual, expected, path)
    else:
        demand(actual == expected, path + ": value differs")


def verify_historical(computed, historical):
    demand(set(computed) == set(COHORTS[:4]) | {"pooled"}, "historical cohort set changed")
    compare_tree(computed, historical["tests"])


def authenticate(root, historical_path):
    paths = {}
    for key, (name, digest) in INPUTS.items():
        path = root / name
        demand(path.is_file() and sha(path) == digest, key + ": immutable source SHA256 mismatch")
        paths[key] = path
    demand(historical_path.is_file() and sha(historical_path) == HISTORICAL_SHA,
           "complete original 199700 robust result SHA256 mismatch")
    receipt = json.loads(paths["receipt"].read_text())
    demand(receipt["total_original_pixel_terminal_photo_rows"] == 447500
           and receipt["total_original_100photo_measured_species"] == 4475
           and receipt["merged_photo_table_sha256"] == INPUTS["photos"][1]
           and receipt["original_149900_photo_rows_unchanged"] is True
           and receipt["old_50km_and_W_estimators_unchanged"] is True,
           "combined original measurement receipt mismatch")
    return paths


def validate_frames(photo, spatial):
    photo = canonical_photos(photo)
    demand(len(photo) == 447500 and photo.inat_taxon_id.nunique() == 4475,
           "actual measured photo/species denominator mismatch")
    demand((photo.groupby("inat_taxon_id").size() == 100).all(), "taxon does not have exactly 100 photos")
    counts = photo.groupby("fcp_source_cohort").inat_taxon_id.nunique().to_dict()
    demand(counts == COUNTS, "original cohort taxon counts mismatch")
    classified = as_bool(photo.global_classifiable) & photo.morph.isin(MORPHS)
    demand(classified.equals(as_bool(photo.global_classifiable)), "classified state label invalid")
    demand(photo.loc[classified].groupby("fcp_source_cohort").size().to_dict() == CLASSIFIED,
           "227107 classified photo cohort totals mismatch")
    eligible = photo.loc[classified].copy()
    n = eligible.groupby("inat_taxon_id").size()
    eligible = eligible.loc[eligible.inat_taxon_id.isin(n[n >= 40].index)]
    demand(eligible.inat_taxon_id.nunique() == 3297, "D-eligible denominator mismatch")
    demand(spatial.inat_taxon_id.is_unique and len(spatial) == 3297,
           "saved spatial species identities/count mismatch")
    demand(set(eligible.inat_taxon_id) == set(spatial.inat_taxon_id), "photo/spatial taxon set differs")
    mapping = photo.groupby("inat_taxon_id").fcp_source_cohort.first()
    demand((spatial.source_cohort.to_numpy() == spatial.inat_taxon_id.map(mapping).to_numpy()).all(),
           "photo/spatial source-cohort mapping differs")
    demand(spatial.loc[spatial.depletion_50km.notna()].groupby("source_cohort").size().to_dict()
           == LOCAL_COUNTS, "saved 1761-species eligibility mismatch")
    return eligible


def compute(photo, spatial, original, historical):
    demand(PERMUTATIONS == 199 and PRIMARY_RADIUS == 50 and MIN_LOCAL_PAIRS == 30,
           "original spatial constants changed")
    demand(original["original_199_vertex_null_changed"] is False
           and original["original_50km_radius_changed"] is False
           and original["original_species_equal_weighting_changed"] is False
           and original["old_three_original_149900_effects_exactly_reproduced"] is True,
           "original null-method identity mismatch")
    eligible = validate_frames(photo, spatial)
    table = spatial.set_index("inat_taxon_id", verify_integrity=True)
    observed, nulls, rows, cohort_labels = [], [], [], []
    for index, (taxon, g) in enumerate(eligible.groupby("inat_taxon_id", sort=True), 1):
        source = str(g.fcp_source_cohort.iloc[0])
        demand(g.photo_id.is_monotonic_increasing, "original photo-ID ordering drift")
        row, worlds = analyse_species(g, SEEDS[source])
        verify_species(row, table.loc[taxon], f"taxon/{taxon}")
        if np.isfinite(row["depletion_50km"]):
            z = worlds[50.]
            demand(z.shape == (199,) and np.isfinite(z).all(), "original vertex-null array invalid")
            observed.append(row["depletion_50km"])
            nulls.append(z)
            cohort_labels.append(source)
            rows.append({**{k: row[k] for k in ("inat_taxon_id", "species", "n_classifiable", "D_pair",
                                                "local_pairs_50km", "D_local_50km", "depletion_50km",
                                                "null_mean_depletion_50km")},
                         "source_cohort": source})
        if index % 100 == 0 or index == 3297:
            print(f"FCP_447500_ROBUST_PROGRESS taxa={index}/3297 local={len(rows)}/1761", flush=True)
    d, z, c = np.array(observed), np.stack(nulls), np.array(cohort_labels)
    demand(len(d) == 1761, "actual local-pair denominator drift")
    groups, old_groups = {}, {}
    for cohort in COHORTS:
        use = c == cohort
        demand(int(use.sum()) == LOCAL_COUNTS[cohort], cohort + ": local eligibility drift")
        verify_mean(d[use], z[use], original["all_sources"][cohort]["50"], cohort)
        groups[cohort] = {"eligible_species": int(use.sum()),
                          "frozen_50km_mean_upper_p": original["all_sources"][cohort]["50"]["mean_depletion_p_upper"],
                          "exploratory_matched_null": statistics(d[use], z[use])}
        if cohort in COHORTS[:4]:
            old_groups[cohort] = groups[cohort]
    verify_mean(d, z, original["pooled_original_50km"], "pooled")
    groups["pooled"] = {"eligible_species": len(d), "frozen_50km_mean_upper_p": .005,
                        "exploratory_matched_null": statistics(d, z)}
    old = c != "new_U100"
    demand(int(old.sum()) == 741, "historical four-cohort denominator drift")
    old_groups["pooled"] = {"eligible_species": 741,
                             "frozen_50km_mean_upper_p": .005,
                             "exploratory_matched_null": statistics(d[old], z[old])}
    verify_historical(old_groups, historical)
    result = {
        "schema": "fcp_real_447500_five_cohort_robust_50km_matched_null_v1",
        "date": "2026-10-18", "status": "complete_post_hoc_source_authenticated",
        "actual_measured_photos": 447500, "all_measured_species": 4475,
        "four_state_classifiable_photos": 227107, "D_eligible_species": 3297,
        "spatially_eligible_species": 1761,
        "method": {"radius_km": 50, "radius_post_hoc": True, "vertex_permutations": 199,
                   "equal_species_weighting": True, "seed_mapping": SEEDS,
                   "top_tail_selection_recomputed_in_every_null_world": True,
                   "primary_new_family": {"group": "pooled", "statistics": list(ROBUST),
                                          "adjustment": "Holm within four post hoc robust tests"},
                   "cohort_families": "separate four-test Holm per cohort; exploratory transport, not independent confirmation",
                   "global_project_multiplicity_adjusted": False},
        "verification": {"each_original_taxon_checked": 3297,
                         "original_five_cohort_means_and_null_quantiles_checked": True,
                         "complete_historical_four_cohort_robust_result_checked": True,
                         "comparison_absolute_tolerance": TOL},
        "tests": groups,
        "historical_four_cohort_reproduced": old_groups,
        "nonclaims": ["Not genetically verified polymorphism or field-population morph frequencies",
                      "Not evidence identifying selection, adaptation, gene flow or plasticity",
                      "Same iNaturalist provider and original image pipeline",
                      "Post hoc analyses do not change original prospective H2 decisions",
                      "Trimmed observed effect can be negative despite exceeding its matched trimmed null"],
    }
    return result, pd.DataFrame(rows), d, z, c


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extracted-root", type=Path, required=True)
    parser.add_argument("--historical-result", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    paths = authenticate(args.extracted_root, args.historical_result)
    print("FCP_447500_ROBUST_ALL_SOURCE_HASHES_AUTHENTICATED", flush=True)
    columns = ["inat_taxon_id", "species", "photo_id", "fcp_source_cohort", "morph",
               "global_classifiable", "latitude", "longitude"]
    photo = pd.read_csv(paths["photos"], usecols=columns, low_memory=False)
    result, species, observed, nulls, cohorts = compute(
        photo, pd.read_csv(paths["species"]), json.loads(paths["spatial_result"].read_text()),
        json.loads(args.historical_result.read_text()))
    result["source"] = {"release_tag": "fcp-all4730-original-method-actual-colour-results-20261017",
                        "archive_sha256": ARCHIVE_SHA,
                        "input_files": {name: {"path": val[0], "sha256": val[1]} for name, val in INPUTS.items()},
                        "historical_full_result_sha256": HISTORICAL_SHA}
    # No completion artifact is written before every historical and source gate passes.
    args.outdir.mkdir(parents=True, exist_ok=True)
    summary = [{"group": group, "statistic": name, **stat}
               for group, data in result["tests"].items()
               for name, stat in data["exploratory_matched_null"].items()]
    species.to_csv(args.outdir / "species_50km_observed.csv", index=False)
    pd.DataFrame(summary).to_csv(args.outdir / "group_statistic_summary.csv", index=False)
    np.savez_compressed(args.outdir / "matched_vertex_nulls.npz",
                        taxon_ids=species.inat_taxon_id.to_numpy(np.int64),
                        cohort_codes=np.array([COHORTS.index(x) for x in cohorts], dtype=np.int8),
                        observed=observed, nulls=nulls)
    result["output_sha256"] = {name: sha(args.outdir / name) for name in
                               ("species_50km_observed.csv", "group_statistic_summary.csv", "matched_vertex_nulls.npz")}
    (args.outdir / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("FCP_447500_ROBUST_MATCHED_NULL_COMPLETE " + json.dumps(result["tests"]["pooled"], sort_keys=True), flush=True)


if __name__ == "__main__":
    main()

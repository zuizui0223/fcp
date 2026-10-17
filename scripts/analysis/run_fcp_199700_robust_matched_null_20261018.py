#!/usr/bin/env python3
"""Post hoc robust-statistic null replay for genuine four-cohort FCP flower photographs.

The prior original hypothesis was the species-equal *mean* of per-species
50-km global-minus-local dissimilarity. It was supported in 199 original
composition-preserving vertex permutations. Earlier outlier deletion was
DESCRIPTIVE ONLY. This audit recomputes the *same original null assignments*
from independently archived, SHA-authenticated measured rows and applies the
same order statistic/positive-fraction/median functional within each null
world. No new pixels, methods, cohort claims, radius selection, or reweighting.

All additional tests are exploratory and post-outcome; never market them as
prospective validation, multiplicity-adjusted discovery, or ecological cause.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis.audit_fcp_199700_spatial_species_leverage_20261017 import (
    ARCHIVE_SHA, RECEIPT_SHA, COHORTS, EXPECTED, sha, demand,
)
from scripts.analysis.replay_fcp_149900_plus_wave4_same_spatial_20261017 import (
    import_exact_measured_rows, original_eligible,
)
from scripts.analysis.run_polymorphism_distributed_polymorphism_posthoc_20261007 import (
    PERMUTATIONS, PRIMARY_RADIUS, MIN_LOCAL_PAIRS, analyse_species,
)

SAMPLE_SIZE = 199700
RELEASE = "fcp-real-original-199700-wave4-50km-H1-W-20261017"
EXPECTED_SOURCE_NAMES = {
    "photos": "combined/fcp_original149900_plus_wave4_same_method_photos.csv.gz",
    "source_receipt": "combined/expanded_data_receipt.json",
    "species": "combined/original_149900_plus_wave4_species_50km.csv",
    "statistical": "combined/expanded_fcp_original_50km_result.json",
}
SEED_LABELS = {"discovery": "discovery", "validation": "validation",
               "prospective": "third", "wave4": "third"}


def validate_release(archive: Path, external_receipt: Path,
                     extracted_root: Path) -> tuple[dict, dict[str, Path]]:
    demand(archive.is_file() and external_receipt.is_file(),
           "original 199700 measured-photo release and external receipt mandatory")
    demand(sha(archive) == ARCHIVE_SHA and sha(external_receipt) == RECEIPT_SHA,
           "SHA drift in immutable original release")
    receipt = json.loads(external_receipt.read_text(encoding="utf-8"))
    demand(receipt.get("archive_sha256") == ARCHIVE_SHA
           and receipt.get("archive_bytes") == archive.stat().st_size
           and receipt.get("all_original_measured_photo_rows") == SAMPLE_SIZE
           and receipt.get("all_original_100photo_taxa") == 1997
           and receipt.get("pooled_original_50km_evaluable_species") == 741
           and receipt.get("original_199_permutation_50km_null_unchanged") is True
           and receipt.get("original_50km_radius_posthoc") is True,
           "published original source/scientific identity mismatch")
    paths = {key: extracted_root / val
             for key, val in EXPECTED_SOURCE_NAMES.items()}
    for key, path in paths.items():
        demand(path.is_file(), f"{key}: missing original released archive member")
        demand(receipt.get("files_sha256", {}).get(EXPECTED_SOURCE_NAMES[key]) == sha(path),
               f"{key}: SHA256 of archived original measured/scientific bytes drift")
    return receipt, paths


def statistics(values: np.ndarray, nulls: np.ndarray) -> dict:
    """Apply each statistic identically to observations and 199 matched null worlds.

    Null columns are paired across all independently permuted species exactly
    as in the original pooled cohort mean. The top 10/20% is removed AFRESH
    in each null world, never the originally observed selected species IDs.
    """
    d = np.asarray(values, dtype=float)
    z = np.asarray(nulls, dtype=float)
    demand(d.ndim == 1 and len(d) >= 2 and z.shape == (len(d), PERMUTATIONS),
           "species/null size mismatch")
    demand(np.isfinite(d).all() and np.isfinite(z).all()
           and np.abs(d).max() <= 1 and np.abs(z).max() <= 1,
           "missing, nonfinite or impossible matched 50-km species nulls")
    ordered = np.sort(d)
    order_null = np.sort(z, axis=0)
    methods = {
        "mean": (float(d.mean()), z.mean(axis=0)),
        "median": (float(np.median(d)), np.median(z, axis=0)),
        "fraction_positive": (float(np.mean(d > 0)), (z > 0).mean(axis=0)),
    }
    for fraction in (0.10, 0.20):
        removed = int(math.ceil(len(d) * fraction))
        key = f"drop_top_{int(100 * fraction)}pct"
        methods[key] = (float(ordered[:-removed].mean()),
                        order_null[:-removed].mean(axis=0))
    out = {}
    for method, (observed, null_replicates) in methods.items():
        p = float((1 + np.count_nonzero(null_replicates >= observed))
                  / (PERMUTATIONS + 1))
        removed = int(math.ceil(len(d) * int(method.split("_")[2][:2]) / 100)) \
            if method.startswith("drop_top_") else 0
        out[method] = {
            "observed": observed,
            "null_mean": float(np.mean(null_replicates)),
            "null_median": float(np.median(null_replicates)),
            "null_q025": float(np.quantile(null_replicates, 0.025)),
            "null_q975": float(np.quantile(null_replicates, 0.975)),
            "upper_p_exploratory": p,
            "n_species": int(len(d)),
            "n_removed_highest_effect": removed,
            "null_selection_recomputed_in_every_permutation": True,
        }
    # Family-wise Holm step-down adjustment within the four newly explored
    # robust statistics for each cohort. The original mean is a separately
    # frozen test, not a fifth newly chosen robust statistic.
    family = ("median", "fraction_positive", "drop_top_10pct", "drop_top_20pct")
    sorted_family = sorted(family, key=lambda k: (out[k]["upper_p_exploratory"], k))
    floor = 0.0
    for rank, name in enumerate(sorted_family):
        raw = out[name]["upper_p_exploratory"]
        floor = max(floor, min(1.0, (len(family) - rank) * raw))
        out[name]["holm_p_within_four_robust_tests_exploratory"] = floor
        out[name]["multiplicity_family_size"] = len(family)
    out["mean"]["holm_p_within_four_robust_tests_exploratory"] = None
    out["mean"]["multiplicity_family_size"] = 0
    return out


def matched_from_photos(measured: pd.DataFrame,
                        frozen_individual: pd.DataFrame,
                        frozen_aggregate: dict,
                        *, expected: dict = EXPECTED) -> dict:
    demand(PERMUTATIONS == 199 and PRIMARY_RADIUS == 50
           and MIN_LOCAL_PAIRS == 30, "original spatial method constants changed")
    require = {"inat_taxon_id", "report_cohort", "depletion_50km", "D_pair",
               "local_pairs_50km", "n_classifiable", "D_local_50km"}
    demand(require.issubset(frozen_individual.columns),
           "original released species table incomplete")
    demand(not frozen_individual["inat_taxon_id"].duplicated().any()
           and set(frozen_individual["report_cohort"]) == set(COHORTS),
           "frozen 4-cohort per-taxon identities drift")
    demand(frozen_aggregate.get("schema") ==
           "fcp_149900_plus_wave4_original_50km_replayed_all_species_v1"
           and frozen_aggregate.get("total_photo_rows") == SAMPLE_SIZE
           and frozen_aggregate.get("all_original_frozen_cohort_50km_values_reproduced") is True
           and frozen_aggregate.get("posthoc_primary_50km_remains_posthoc") is True
           and frozen_aggregate.get("original_method", {}).get("permutations") == 199
           and frozen_aggregate.get("original_method", {}).get("equal_species_weights") is True
           and frozen_aggregate.get("original_method", {}).get("min_local_pairs") == 30,
           "archived historical inference protocol cannot be authenticated")
    table = frozen_individual.set_index("inat_taxon_id", verify_integrity=True)
    n_full = 0
    observed_by = {}
    null_by = {}
    summaries = {}
    for cohort in COHORTS:
        rows = []
        nullrows = []
        subset = original_eligible(measured.loc[measured.fcp_source_cohort.eq(cohort)])
        original_rows = frozen_individual.loc[frozen_individual.report_cohort.eq(cohort)]
        counted = set()
        for taxon, g in subset.groupby("inat_taxon_id", sort=True):
            g = g.dropna(subset=["latitude", "longitude"])
            if len(g) < 40:
                continue
            old, perm = analyse_species(g, SEED_LABELS[cohort])
            taxon = int(taxon)
            demand(taxon in table.index and table.loc[taxon, "report_cohort"] == cohort,
                   f"{cohort}: previously frozen individual taxon absent")
            current = table.loc[taxon]
            counted.add(taxon)
            demand(int(current["n_classifiable"]) == old["n_classifiable"]
                   and int(current["local_pairs_50km"]) == old["local_pairs_50km"]
                   and abs(float(current["D_pair"]) - old["D_pair"]) < 1e-10,
                   f"{cohort}: classifiable count, local-pair graph or colour labels changed")
            old_dep = float(old["depletion_50km"])
            frozen_dep = float(current["depletion_50km"])
            if np.isnan(old_dep):
                demand(np.isnan(frozen_dep),
                       f"{cohort}: frozen species falsely became evaluable")
                continue
            demand(np.isfinite(frozen_dep) and abs(frozen_dep - old_dep) < 1e-10
                   and abs(float(current["D_local_50km"]) - old["D_local_50km"]) < 1e-10,
                   f"{cohort}: frozen observed local-depletion statistic changed")
            z = perm[PRIMARY_RADIUS]
            demand(z.shape == (PERMUTATIONS,) and np.isfinite(z).all(),
                   f"{cohort}: original matched permutation missing")
            rows.append(old_dep)
            nullrows.append(z)
        demand(counted == set(original_rows["inat_taxon_id"].astype(int)),
               f"{cohort}: source-to-archived species count/identity mismatch")
        n_full += len(counted)
        n, target, p_original = expected[cohort]
        demand(len(rows) == n and len(nullrows) == n,
               f"{cohort}: actual 50km local-pair eligibility differs from original")
        obs = np.array(rows, float)
        world = np.stack(nullrows, axis=0)
        original_stat = frozen_aggregate["results"]["cohorts"][cohort]["radii"]["50"]
        original_p_replayed = (1 + (world.mean(axis=0) >= obs.mean()).sum()) / 200
        demand(abs(obs.mean() - target) < 1e-10
               and original_stat["n_species"] == n
               and abs(original_stat["mean_depletion"] - target) < 1e-10
               and original_stat["mean_depletion_p_upper"] == p_original
               and float(original_p_replayed) == p_original,
               f"{cohort}: historical 50km mean null not reproduced exactly")
        observed_by[cohort], null_by[cohort] = obs, world
        summaries[cohort] = {
            "eligible_species": int(n),
            "frozen_50km_mean_upper_p": p_original,
            "exploratory_matched_null": statistics(obs, world),
        }
    demand(n_full == frozen_aggregate["results"]["pooled_n_spatial_eligible_species"],
           "previously frozen all-eligible species table size changed")
    pooled_obs = np.concatenate([observed_by[k] for k in COHORTS])
    pooled_null = np.vstack([null_by[k] for k in COHORTS])
    n, target, p_original = expected["pooled"]
    pooled_stat = frozen_aggregate["results"]["pooled"]["50"]
    replayed_p = (1 + (pooled_null.mean(axis=0) >= pooled_obs.mean()).sum()) / 200
    demand(len(pooled_obs) == n
           and abs(pooled_obs.mean() - target) < 1e-10
           and pooled_stat["n_species"] == n
           and abs(pooled_stat["mean_depletion"] - target) < 1e-10
           and pooled_stat["mean_depletion_p_upper"] == p_original
           and float(replayed_p) == p_original,
           "pooled original 50km null/mean could not be reconstructed")
    summaries["pooled"] = {
        "eligible_species": int(n),
        "frozen_50km_mean_upper_p": p_original,
        "exploratory_matched_null": statistics(pooled_obs, pooled_null),
    }
    return {
        "schema": "fcp_real_199700_four_cohort_robust_50km_matched_null_v1",
        "status": "post_hoc_exploratory_not_multiple_test_corrected",
        "actual_measured_photos": 199700,
        "all_measured_species": 1997,
        "frozen_original_spatial_evaluable_species": 741,
        "same_original_199_vertex_permutations": True,
        "all_cohorts_old_mean_nulls_replayed": True,
        "radius_50km_post_hoc": True,
        "no_new_photo_classification": True,
        "not_independent_source_replication": True,
        "not_field_population_morph_frequency": True,
        "no_local_adaptation_or_selection_claim": True,
        "trim_rule_applied_again_in_each_permutation": True,
        "tests": summaries,
    }


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--archive", type=Path, required=True)
    p.add_argument("--external-receipt", type=Path, required=True)
    p.add_argument("--extracted-root", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    _, paths = validate_release(a.archive, a.external_receipt, a.extracted_root)
    data_receipt = json.loads(paths["source_receipt"].read_text(encoding="utf-8"))
    measured = import_exact_measured_rows(paths["photos"], data_receipt)
    output = matched_from_photos(
        measured,
        pd.read_csv(paths["species"]),
        json.loads(paths["statistical"].read_text(encoding="utf-8")),
    )
    output["source_public_release"] = RELEASE
    output["source_archive_sha256"] = ARCHIVE_SHA
    output["source_external_receipt_sha256"] = RECEIPT_SHA
    output["source_files_sha256"] = {k: sha(v) for k, v in paths.items()}
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(output, indent=2, sort_keys=True, allow_nan=False) + "\n",
                     encoding="utf-8")
    focal = output["tests"]["pooled"]["exploratory_matched_null"]
    print("FCP_ORIGINAL_199700_ROBUST_MATCHED_NULL_COMPLETE",
          json.dumps({
              "pooled_median": focal["median"],
              "pooled_fraction_positive": focal["fraction_positive"],
              "pooled_drop_top_10pct": focal["drop_top_10pct"],
              "pooled_drop_top_20pct": focal["drop_top_20pct"],
              "frozen_mean_original_p":
                  output["tests"]["pooled"]["frozen_50km_mean_upper_p"],
          }, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()

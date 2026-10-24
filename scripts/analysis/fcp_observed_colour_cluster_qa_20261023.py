#!/usr/bin/env python3
"""Audit apparent 1–4 photographic colour categories versus source-support and geography.

One immutable row per original iNaturalist species. NOT a model identifying biological
morphs, genetic clustering, or 1–4 continuous colour-space mixture components.
Do not drop low-quality rows from the original source population or rewrite P1–P4.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "results/fcp_np_all4475_unified_species_ledger_20261021.csv"
SOURCE_SHA = "23d3daadc9799c2e7b5867d54ccc1ed98369ba79a81b72dc7f831e5197e63c49"
PHOTO_SHA = "8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6"
PHOTO_CLASSES = (
    "miss_coarse_photo_white",
    "miss_coarse_photo_yellow_orange",
    "miss_coarse_photo_red_pink",
    "miss_coarse_photo_blue_purple",
)
EXPECTED_N = 4475
EXPECTED_PHOTOS = 447500
EXPECTED_CLASSIFIED = 227107


def check(ok, why):
    if not ok:
        raise ValueError("FCP_CLUSTER_QA_SOURCE_ERROR: " + why)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as src:
        for block in iter(lambda: src.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def rows_from_csv(path: Path, *, verify=True):
    if verify:
        check(sha256(path) == SOURCE_SHA, "4,475-species ledger SHA-256 mismatch")
    result = []
    with path.open(newline="", encoding="utf-8") as src:
        reader = csv.DictReader(src)
        required = {*PHOTO_CLASSES, "species", "inat_taxon_id", "source_cohort",
                    "n_classified_photo", "n_local50km_colour_photo_pairs",
                    "global_vs_local_50km_colour_mismatch_delta",
                    "global_valid_geo_mismatch_D", "n_missing_colour_label",
                    "n_distinct_observers_source"}
        check(required <= set(reader.fieldnames or []), "required frozen source columns missing")
        ids = set()
        for r in reader:
            taxon = int(r["inat_taxon_id"])
            check(taxon not in ids, "duplicate species ID")
            ids.add(taxon)
            counts = tuple(int(r[c]) for c in PHOTO_CLASSES)
            check(all(c >= 0 for c in counts), "negative photo-colour count")
            n = int(r["n_classified_photo"])
            missing = int(r["n_missing_colour_label"])
            check(sum(counts) == n, "photo labels do not sum to classified records")
            check(n + missing == 100, "measured species does not have 100 original photos")
            pairs = int(r["n_local50km_colour_photo_pairs"])
            value = r["global_vs_local_50km_colour_mismatch_delta"]
            delta = float(value) if value not in ("", None) else None
            if pairs == 0:
                check(delta is None, "source-no-pair NA was converted to a false 0")
            else:
                check(delta is not None and math.isfinite(delta), "missing valid photo-pair mismatch")
            richness = sum(c > 0 for c in counts)
            check(1 <= richness <= 4, "invalid coarse-state photo richness")
            result.append({
                "taxon_id": taxon, "genus": r["species"].split()[0],
                "source_cohort": r["source_cohort"],
                "counts": counts, "k": richness, "n": n, "missing": missing,
                "n_observers": int(r["n_distinct_observers_source"]),
                "n_classes_ge5": sum(c >= 5 for c in counts),
                "n_classes_ge10": sum(c >= 10 for c in counts),
                "min_observed_class_n": min(c for c in counts if c > 0),
                "pairs": pairs, "delta": delta,
            })
    if verify:
        check(len(result) == EXPECTED_N, "taxon denominator drift")
        check(sum(x["n"] for x in result) == EXPECTED_CLASSIFIED, "classified source count drift")
        check(sum(x["n"] + x["missing"] for x in result) == EXPECTED_PHOTOS,
              "original complete source count drift")
    return result


def mean(values):
    return sum(values) / len(values) if values else None


def stats(rows):
    geo = [r for r in rows if r["pairs"] >= 30 and r["n"] >= 40]
    deltas = sorted(r["delta"] for r in geo)
    return {
        "species": len(rows),
        "classifiable_photo_mean": mean([r["n"] for r in rows]),
        "species_ge40_classified": sum(r["n"] >= 40 for r in rows),
        "species_any_50km_pair": sum(r["pairs"] > 0 for r in rows),
        "species_strict_50km_ge40_ge30pairs": len(geo),
        "strict_50km_mean_global_minus_local_mismatch": mean(deltas),
        "strict_50km_median_global_minus_local_mismatch": (
            deltas[len(deltas)//2] if deltas else None
        ),
        "strict_50km_species_with_positive_depletion": sum(d > 0 for d in deltas),
    }


def analyse(rows):
    richness = {str(k): stats([r for r in rows if r["k"] == k]) for k in (1,2,3,4)}
    supported = {}
    for k in (1,2,3,4):
        subset = [r for r in rows if r["k"] == k and r["n"] >= 40
                  and r["min_observed_class_n"] >= 5]
        supported[str(k)] = stats(subset)
    flags = {
        "at_least_two_labels": lambda r: r["k"] >= 2,
        "at_least_two_photo_classes_ge5_and_nclass_ge40":
            lambda r: r["n"] >= 40 and r["n_classes_ge5"] >= 2,
        "at_least_three_photo_classes_ge5_and_nclass_ge40":
            lambda r: r["n"] >= 40 and r["n_classes_ge5"] >= 3,
        "four_photo_classes_ge5_and_nclass_ge40":
            lambda r: r["n"] >= 40 and r["n_classes_ge5"] == 4,
        "at_least_two_photo_classes_ge10_and_nclass_ge40":
            lambda r: r["n"] >= 40 and r["n_classes_ge10"] >= 2,
        "at_least_three_photo_classes_ge10_and_nclass_ge40":
            lambda r: r["n"] >= 40 and r["n_classes_ge10"] >= 3,
        "four_photo_classes_ge10_and_nclass_ge40":
            lambda r: r["n"] >= 40 and r["n_classes_ge10"] == 4,
    }
    support_layers = {name: stats([r for r in rows if predicate(r)])
                      for name,predicate in flags.items()}
    data = {
        "schema": "fcp_source4475_photographic_label_richness_and_geo_quality_audit_v1",
        "analysis_status": "exploratory post hoc retrospective source quality/sampling audit",
        "measured_photo_sha256": PHOTO_SHA,
        "source_4475_species_ledger_sha256": SOURCE_SHA,
        "n_original_photos": EXPECTED_PHOTOS,
        "n_classified_photos": EXPECTED_CLASSIFIED,
        "n_species": EXPECTED_N,
        "observed_colour_label_richness_1to4": richness,
        "each_detected_colour_ge5_and_atleast40photos": supported,
        "nested_photo_support_sensitivity": support_layers,
        "all_independent_observers_per_colour_unknown_in_aggregate": True,
        "minor_colour_threshold_is_exploratory_not_biological_validation": True,
        "photo_colour_classes_are_not_unsupervised_9_palette_modes": True,
        "two_means_Hellinger_estimator_does_not_estimate_mode_count": True,
        "only_observed_geographic_delta_no_new_stratum_conditioned_spatial_null": True,
        "original_matched_source_50km_null_not_retested_by_this_summary": True,
        "matched_genotype_or_petal_stage_fitness_evidence": False,
        "no_exact_coordinates_observation_or_photo_ids_exported": True,
        "strict_5km_month_matched_positive_in_original_test": False,
        "strict_5km_source_null_p_ge5pairs": 0.300,
    }
    return data


def self_test():
    synthetic=[
        dict(taxon_id=1,genus="A",source_cohort="x",counts=(40,0,0,0),k=1,
             n=40,missing=60,n_observers=40,n_classes_ge5=1,n_classes_ge10=1,
             min_observed_class_n=40,pairs=35,delta=0.0),
        dict(taxon_id=2,genus="B",source_cohort="x",counts=(35,5,0,0),k=2,
             n=40,missing=60,n_observers=40,n_classes_ge5=2,n_classes_ge10=1,
             min_observed_class_n=5,pairs=30,delta=0.07),
        dict(taxon_id=3,genus="C",source_cohort="x",counts=(19,10,10,1),k=4,
             n=40,missing=60,n_observers=40,n_classes_ge5=3,n_classes_ge10=3,
             min_observed_class_n=1,pairs=0,delta=None)
    ]
    r=analyse(synthetic)
    check(r["observed_colour_label_richness_1to4"]["4"]["species"]==1,"synthetic k4")
    check(r["each_detected_colour_ge5_and_atleast40photos"]["4"]["species"]==0,
          "synthetic one singleton is not robust 4")
    check(r["nested_photo_support_sensitivity"][
          "at_least_three_photo_classes_ge5_and_nclass_ge40"]["species"]==1,
          "synthetic three well-supported modes")
    check(r["nested_photo_support_sensitivity"][
          "at_least_two_photo_classes_ge5_and_nclass_ge40"][
          "strict_50km_species_with_positive_depletion"]==1,"synthetic geo")
    print("FCP_COLOUR_CLUSTER_QA_SYNTHETIC_TEST_PASS")


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--source",type=Path,default=SOURCE)
    parser.add_argument("--output",type=Path,default=None)
    parser.add_argument("--self-test",action="store_true")
    args=parser.parse_args()
    if args.self_test:
        self_test()
        return
    result=analyse(rows_from_csv(args.source))
    # Frozen ledger-level arithmetic checks; NOT a new matched permutation p-value.
    check([result["observed_colour_label_richness_1to4"][str(i)]["species"]
           for i in range(1,5)]==[286,1743,1632,814],
          "source coarse-label richness distribution drift")
    check(result["each_detected_colour_ge5_and_atleast40photos"]["4"]["species"]==11,
          "four-label ≥5-per-colour source denominator drift")
    check(result["nested_photo_support_sensitivity"][
          "at_least_two_photo_classes_ge5_and_nclass_ge40"]["species"]==1592,
          "≥2 evidenced photo colours denominator drift")
    check(result["nested_photo_support_sensitivity"][
          "at_least_three_photo_classes_ge5_and_nclass_ge40"]["species"]==309,
          "≥3 evidenced photo colours denominator drift")
    check(result["nested_photo_support_sensitivity"][
          "four_photo_classes_ge10_and_nclass_ge40"]["species"]==0,
          "four ≥10-each colours unexpected source drift")
    check(result["observed_colour_label_richness_1to4"]["1"]["strict_50km_mean_global_minus_local_mismatch"]==0,
          "monolabel taxon cannot have positive categorical spatial depletion")
    allstrict=[r for r in rows_from_csv(args.source) if r["n"]>=40 and r["pairs"]>=30]
    check(len(allstrict)==1761,"original strict support drift")
    check(abs(mean([r["delta"] for r in allstrict])-0.0176149021685411)<1e-11,
          "original frozen total geographic depletion drift")
    output=json.dumps(result,indent=2,ensure_ascii=False,sort_keys=True)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(output,encoding="utf-8")
    print("FCP_COLOUR_CLUSTER_QA_ALL4475_SOURCE_VALIDATED")
    print(output)


if __name__ == "__main__":
    main()

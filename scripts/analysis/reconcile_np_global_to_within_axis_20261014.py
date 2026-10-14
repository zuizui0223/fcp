#!/usr/bin/env python3
"""Reconcile already measured global floral colours and 3 high-depth white-axis cohorts.

This is a reporting-only exact-arithmetic audit. It does not refit the white-axis
clustering, recompute H2 permutation p-values, reopen photo pixels or convert
historical candidate metadata into a newly measured independent cohort.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
NP = ROOT / "results" / "np_unified_20261014"
ORIG = ROOT / "results" / "polymorphism_white_axis_targeted_test_20260912" / "result.json"
PROS = ROOT / "results" / "polymorphism_h2_third_cohort_prospective_white_axis_20260917" / "result.json"
MEAS = ROOT / "results" / "polymorphism_h2_third_cohort_prospective_measurement_20260917" / "result.json"
GLOBAL = NP / "source_pr135_global_colour_and_crosscell.json"
FUTURE = NP / "source_pr133_next2000_metadata_only.json"
BALANCE = NP / "source_pr158_white_balance.json"
DEST = NP / "global_to_within_white_axis_coverage_and_reconciliation.json"


def calculate():
    old = json.loads(ORIG.read_text(encoding="utf8"))
    future = json.loads(FUTURE.read_text(encoding="utf8"))
    global_data = json.loads(GLOBAL.read_text(encoding="utf8"))
    photo = json.loads(BALANCE.read_text(encoding="utf8"))
    pros = json.loads(PROS.read_text(encoding="utf8"))
    measured = json.loads(MEAS.read_text(encoding="utf8"))
    assert old["status"] == "targeted_post_audit_decomposition_not_untouched_confirmatory_test"
    assert pros["status"] == "untouched_prospective_test_of_previously_frozen_axis"
    assert future["biological_outcome_state"]["new_image_measurement_performed"] is False
    assert future["new_prospective_sampled_species"] == 2000
    assert future["permanently_unopened_reserve_species"] == 730
    assert future["unused_capacity_species_after_all_allocations"] == 2730
    assert measured["measurement_evaluable_species"] == 377
    assert old["fixed_axis_palette_loadings"] == pros["fixed_axis_palette_loadings"]
    gd = global_data["breadth_single_original_photo_per_species"]
    cx = global_data["observer_disjoint_crosscell_morph_comparison"]
    rp = next(row for row in cx["regional_and_global_descriptive"] if row["region"] == "ALL")
    n_pairs = rp["n_fixed_pairs"]
    n_classified_pairs = rp["both_photos_classified"]
    n_disc = rp["discordant_photo_colour"]
    n_unknown = rp["one_or_both_unclassifiable"]
    assert n_pairs == 13416 and n_classified_pairs == 3196
    assert n_disc == 797 and n_unknown == 10220
    assert n_classified_pairs + n_unknown == n_pairs
    assert n_disc + rp["same_photo_colour"] == n_classified_pairs
    assert math.isclose(n_disc / n_classified_pairs, rp["discordance_fraction_among_both_classified"], abs_tol=1e-15)
    assert math.isclose(n_disc / n_pairs, rp["no_assumption_full_fixed_pair_discordance_bounds"][0], abs_tol=1e-15)
    assert math.isclose((n_disc + n_unknown) / n_pairs, rp["no_assumption_full_fixed_pair_discordance_bounds"][1], abs_tol=1e-15)
    assert gd["species_source_rows"] == 42111 and gd["species_four_colour_classifiable"] == 18457
    assert gd["species_four_colour_classifiable"] + gd["unclassifiable"] == 42111
    assert sum(gd["colour_counts"].values()) == gd["species_four_colour_classifiable"]
    assert photo["source"]["unique_original_photos"] == 100543
    assert photo["source"]["unique_original_4colour_classifiable_photos"] == 45953
    assert photo["source"]["original_unclassifiable_photos"] == 54590

    counts = []
    for label, row, tier in (
        ("discovery", old["thresholds"], "discovery"),
        ("validation", old["thresholds"], "reserve"),
        ("prospective", pros["thresholds"], None),
    ):
        record = {"cohort": label, "inferential_role": "untouched prospective confirmation" if label == "prospective" else "retrospective target localization"}
        for suffix, short in (("primary_0_10", "primary"), ("strict_0_20", "strict")):
            t = row[suffix] if tier is None else row[suffix][tier]
            n = int(t["species"])
            w = float(t.get("observed_W", t.get("observed_mean_squared_white_axis_alignment", float("nan"))))
            d = t["structured_null_summary"]
            null_mean = float(d["mean"])
            p = t.get("structured_null_upper_p")
            assert 0 < n <= 500 and 0 <= w <= 1 and 0 <= null_mean <= 1
            assert p is not None and 0 < p < 1
            record[short] = {
                "eligible_species": n, "observed_W": w,
                "source_null_mean_W": null_mean, "source_null_median_W": float(d["q50"]),
                "source_original_null_p": p,
                "source_null_replicates": 999
            }
        counts.append(record)

    combined = {}
    for tier in ("primary", "strict"):
        total_n = sum(x[tier]["eligible_species"] for x in counts)
        weighted_w = sum(x[tier]["eligible_species"] * x[tier]["observed_W"] for x in counts)/total_n
        weighted_null_mean = sum(x[tier]["eligible_species"] * x[tier]["source_null_mean_W"] for x in counts)/total_n
        assert total_n == (439 if tier == "primary" else 226)
        combined[tier] = {
            "eligible_species_total": total_n,
            "eligible_species_by_cohort": {x["cohort"]: x[tier]["eligible_species"] for x in counts},
            "species_equal_weighted_observed_W": weighted_w,
            "species_equal_weighted_cohort_null_MEANS": weighted_null_mean,
            "observed_minus_weighted_null_means": weighted_w-weighted_null_mean,
            "no_pooled_null_worlds_constructed": True,
            "no_pooled_p_value_computed": True,
            "do_not_relabel_post_audit_cohorts_as_prospective": True
        }
    assert measured["species"] == 499 and sum((369, 363, measured["measurement_evaluable_species"])) == 1109
    return {
        "schema": "fcp_np_global_to_within_white_axis_reporting_20261014_v1",
        "status": "EXISTING_SOURCE_ARITHMETIC_ONLY_NO_NEW_SOURCE_PHOTOS_OR_PERMUTATIONS",
        "source": {
            "global_measured_PR": 135,
            "global_source_image_commit": global_data["origin_commit"],
            "global_provenance_success_run": global_data["execution_receipt"]["workflow_run_id"],
            "old_H2_target_status": old["status"],
            "new_H2_confirm_status": pros["status"],
            "new_2000_species_metadata_PR": 133
        },
        "worldwide_four_colour_coverage": {
            "one_photo_nominal_species": 42111,
            "one_photo_colour_classifiable_species": gd["species_four_colour_classifiable"],
            "one_photo_unknown_species": gd["unclassifiable"],
            "one_photo_original_four_colour_counts": gd["colour_counts"],
            "taxon_geographic_cell_photo_records": global_data["regional_taxon_cell_measurement"]["all_taxon_cell_rows"],
            "taxon_cell_colour_classifiable": global_data["regional_taxon_cell_measurement"]["classified_taxon_cell_rows"],
            "merged_unique_photo_ids": photo["source"]["unique_original_photos"],
            "merged_photo_four_colour_classifiable": photo["source"]["unique_original_4colour_classifiable_photos"],
            "merged_photo_unknown": photo["source"]["original_unclassifiable_photos"],
            "global_colour_polymorphism_prevalence_estimable": False
        },
        "cross_geographic_cell_photo_pairs": {
            "n_fixed_observer_disjoint_pairs": n_pairs,
            "both_classifiable_pairs": n_classified_pairs,
            "discordant_four_colour_pairs": n_disc,
            "unknown_hue_in_either": n_unknown,
            "discordance_among_both_classified": n_disc / n_classified_pairs,
            "no_assumption_full_frame_lower": n_disc / n_pairs,
            "no_assumption_full_frame_upper": (n_disc+n_unknown)/n_pairs,
            "distinct_individual_plants_or_genotypes_identified": False,
            "geographic_variation_not_within_population_coexistence": True
        },
        "high_depth_sample_coverage": {
            "n_historical_high_depth_photographed_species": 1499,
            "n_retained_high_depth_photo_ids": 149900,
            "n_at_least_40_classifiable_species_all_cohorts": 1109,
            "n_classifiable_by_cohort": {"discovery": 369, "validation": 363, "prospective": 377},
            "n_prospective_sampled_species": measured["species"],
            "n_prospective_photo_rows": measured["rows"],
            "n_prospective_classifiable_rows": measured["classifiable_rows"]
        },
        "white_axis_already_evaluable_three_cohorts": counts,
        "white_axis_descriptive_all_cohort_summary": combined,
        "new_high_depth_species_not_yet_measured": {
            "historical_u100_metadata_opportunity": 4730,
            "old_selected_taxon_ids_excluded": future["historically_allocated_capacity_species"],
            "never_allocated_capacity_species": future["unused_capacity_species_after_all_allocations"],
            "preselected_2000_species_pixels_unopened": future["new_prospective_sampled_species"],
            "reserve_730_species_unopened": future["permanently_unopened_reserve_species"],
            "max_possible_if_all_new_2000_measured": future["potential_total_photographed_species_if_all_new_sample_successful"],
            "new_colour_axis_confirmation_available": False
        },
        "boundaries": [
            "42,111 single-image or 85,337 taxon-cell records are NOT high-depth within-species colour axes.",
            "439 is the union of separately fitted original 152+129 and prospectively held-out 158 eligible axes, not 439 untouched confirmatory taxa.",
            "Weighted average of null cohort MEANS is descriptive expectation, not a reconstructed joint pooled permutation null or p value.",
            "The historical H2 white state can be affected by camera near-clipping and does not prove pigment absence.",
            "No universal white/coloured cost/benefit fitness balance is identified, regardless of global observed 46.95 percent photo-white share."
        ]
    }


def main():
    cli = argparse.ArgumentParser()
    cli.add_argument("--check", action="store_true")
    cli.add_argument("--write", action="store_true")
    opt = cli.parse_args()
    output=calculate()
    if opt.write:
        DEST.parent.mkdir(parents=True, exist_ok=True)
        DEST.write_text(json.dumps(output,ensure_ascii=False,sort_keys=True,indent=2)+"\n", encoding="utf8")
    if opt.check:
        frozen=json.loads(DEST.read_text(encoding="utf8"))
        assert output == frozen, "Global-to-axis source reconciliation differs from committed reporting-only receipt"
    print("world classified species =",output["worldwide_four_colour_coverage"]["one_photo_colour_classifiable_species"])
    print("cross-cell discordant pairs =",output["cross_geographic_cell_photo_pairs"]["discordant_four_colour_pairs"])
    for tier, d in output["white_axis_descriptive_all_cohort_summary"].items():
        print(tier,d["eligible_species_total"],d["species_equal_weighted_observed_W"],d["species_equal_weighted_cohort_null_MEANS"])


if __name__ == "__main__":
    main()

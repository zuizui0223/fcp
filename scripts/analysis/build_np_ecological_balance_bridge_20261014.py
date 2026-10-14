#!/usr/bin/env python3
"""Reporting-only cross-scale ecological evidence reconciliation for FCP New Phytologist.

This script NEVER re-estimates photo models, selection, fitness, or a null p-value.
It independently reconstructs arithmetic from immutable historical PR receipts.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = ROOT / "results" / "np_unified_20261014"
SOURCES = {
    "balance": "source_pr158_white_balance.json",
    "latitude": "source_pr158_latitude_null.json",
    "missingness": "source_pr158_mnar_sensitivity.json",
    "soil_matched": "source_pr168_soil_matched.json",
    "midrain": "source_pr160_midrain_comixing.json",
}
ORIGINAL_PHOTO_SHA = "2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0"


def compute() -> dict:
    s = {name: json.loads((HERE / filename).read_text(encoding="utf-8"))
         for name, filename in SOURCES.items()}
    b, lat, mn, soil, mid = (s[k] for k in SOURCES)
    photo = b["source"]
    total = photo["unique_original_photos"]
    classified = photo["unique_original_4colour_classifiable_photos"]
    missing = photo["original_unclassifiable_photos"]
    white = b["world_photo_white"]["original_white_photo_count"]
    assert total == 100543 and classified == 45953 and missing == 54590
    assert classified + missing == total and white == 21574
    assert all(x == ORIGINAL_PHOTO_SHA for x in (
        photo["original_artifact_zip_sha256"], lat["source_sha256"],
        mn["source_original_archive_sha256"], mid["source"]["archive_sha256"]))
    assert soil["status"] == "SOURCE_REPLAY_VERIFIED"

    n_bands = lat["original_latitude_band_n_photos"]
    vals = lat["photo_white_fractions_by_original_true_site_latitude_band"]
    assert n_bands == [19642, 25033, 1204]
    assert sum(n_bands) == lat["n_original_classified_geotagged_photos"]
    original_bands = b["world_photo_white"]["true_photo_site_latitude"]
    for i, row in enumerate(original_bands):
        assert row["n"] == n_bands[i]
        assert math.isclose(row["white_n"]/row["n"], vals[i], abs_tol=1e-14)
    range_pp = 100 * (max(vals) - min(vals))
    assert math.isclose(range_pp, lat["observed_range_white_percentage_points"], abs_tol=1e-10)
    within = lat["conditional_colour_photo_label_nulls"]["within_species"]

    co = b["ecological_compensation"]
    observed_var = co["photo_weighted_cell_whiteness_variance"]
    recon = (co["species_compositional_component_variance"]
             + co["within_species_local_component_variance"]
             + 2*co["component_covariance"])
    assert abs(recon-observed_var) < 1e-12
    assert co["component_covariance"] > 0
    assert co["negative_compensation_one_sided_p"] > .95

    wetdry = mn["original_photo_wet_dry_unweighted"]
    raw_wet_dry = wetdry["9"]["white_fraction"] - wetdry["0"]["white_fraction"]
    assert math.isclose(raw_wet_dry, mn["original_photo_wet_minus_dry_unweighted"], abs_tol=1e-12)
    tipping = mn["latent_hue_missing_not_at_random_tipping_curve"]["0.5"]
    assert tipping["feasible_value_between_zero_and_one"]
    matched = soil["soil_addition_on_identical_20786_species"]
    assert matched["CI95"][0] < 0 < matched["CI95"][1]
    p_mid = [v["p_one_sided"] for v in mid["results"].values()
             if v.get("p_one_sided") is not None]
    assert len(p_mid) == 5 and min(p_mid) > .05
    assert mid["study"]["species_eligible_in_all_three_rain_photo_pair_bins_any_study"] == 0

    return {
        "schema": "fcp_np_reporting_only_water_axis_balance_bridge_20261014_v1",
        "status": "RECONCILED_SOURCE_RECEIPTS_NO_NEW_BIOLOGICAL_INFERENCE",
        "frozen_sources": {
            "original_photo_zip_sha256": ORIGINAL_PHOTO_SHA,
            "source_pr158_commit": "0bc53089dc9d1baed207c2ea898b130cb7b07566",
            "source_pr168_commit": "30a855c28ac060bed084337a85cd09ce1092ee56",
            "source_pr160_commit": "2601f1ed689cd96678019f331457178ca1150e89",
            "source_file_names": SOURCES,
        },
        "photographic_composition": {
            "unique_photo_ids": total,
            "classifiable_photos": classified,
            "unclassifiable_photos": missing,
            "white_classifiable_photos": white,
            "fraction_white_among_classified": white/classified,
            "fraction_classification_missing": missing/total,
            "no_assumption_all_photo_white_lower": white/total,
            "no_assumption_all_photo_white_upper": (white+missing)/total,
            "not_population_morph_prevalence": True,
        },
        "latitude_near_half_null": {
            "classifiable_geotagged_photos": sum(n_bands),
            "three_latitude_band_white_fractions": vals,
            "observed_band_range_pp": range_pp,
            "within_species_colour_preserving_null_mean_range_pp": within["null_mean_range_percentage_points"],
            "within_species_colour_preserving_null_p_unusually_even": within["one_sided_p_observed_unusually_even_or_lower_range"],
            "extraordinary_half_white_homeostasis_not_supported": True,
        },
        "geographic_compensation_null": {
            "repeat_species": co["repeat_species"],
            "observed_photo_cell_variance": observed_var,
            "reconstructed_photo_cell_variance": recon,
            "reconstruction_absolute_error": abs(recon-observed_var),
            "composition_vs_local_covariance": co["component_covariance"],
            "negative_compensation_one_sided_p": co["negative_compensation_one_sided_p"],
            "hypothesized_geographic_cancellation_not_supported": True,
        },
        "moisture_opportunity_and_missingness": {
            "unadjusted_wet_minus_dry_photo_white_pp": 100*raw_wet_dry,
            "driest_classified_white_fraction": wetdry["0"]["white_fraction"],
            "wettest_classified_white_fraction": wetdry["9"]["white_fraction"],
            "assumed_dry_unclassified_white_fraction": 0.5,
            "wet_unclassified_white_fraction_needed_to_erase_contrast": tipping["required_unclassified_wet_white_fraction_for_EQUAL_true_dry_wet_white"],
            "unobserved_outcome_colour_mnar_remains_identification_problem": True,
        },
        "soil_added_on_exact_same_species_unseen_genus_test": {
            "n_species": soil["primary36"]["absolute_species_colour"]["species"],
            "environmental_brier_increment": matched["gain"],
            "original_95_interval": matched["CI95"],
            "supports_incremental_soil_information": False,
            "applies_only_to_original_unseen_genus_design": True,
        },
        "intermediate_rainfall_local_cocomposition": {
            "reported_nonhold_tests": len(p_mid),
            "smallest_unadjusted_one_sided_p": min(p_mid),
            "n_species_in_all_three_rainfall_pair_bins": 0,
            "supported_local_midrain_mixing_excess": False,
        },
        "white_axis_separate_frozen_high_depth_source": {
            "n_high_depth_species": 1499,
            "prospective_admissible_species": 158,
            "observed_W": 0.5172457461053418,
            "null_median_W": 0.45714281500506854,
            "p_plus_one": 0.001,
            "source_independent_of_sparse_global_photo_colour_system": False,
            "cross_scale_morph_fitness_mechanism_tested": False,
        },
        "scientific_decision": {
            "source_observational_moisture_correspondence": "BROAD_SAMPLE_ASSOCIATION_ONLY",
            "prospective_white_coloured_direction": "SUPPORTED_CONDITIONAL_ON_WHITE_CLASSIFIER",
            "global_near_half_white": "DESCRIPTIVE_NOT_AN_EXCEPTIONAL_EQUILIBRIUM",
            "geographic_compensation": "NOT_SUPPORTED",
            "midrain_local_colour_mixing": "NOT_SUPPORTED_OR_INADEQUATELY_IDENTIFIED",
            "pigmentation_benefit_cost_balance": "TESTABLE_HYPOTHESIS_NO_PIGMENT_OR_RELATIVE_FITNESS",
        },
    }


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--check", action="store_true")
    p.add_argument("--write", action="store_true")
    a = p.parse_args()
    result = compute()
    dest = HERE / "ecological_balance_reporting_reconciliation.json"
    if a.write:
        dest.write_text(json.dumps(result, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")
    if a.check:
        frozen = json.loads(dest.read_text(encoding="utf-8"))
        assert frozen == result, "Ecological reporting receipt differs from source rederivation"
    print("REPORTING ONLY", result["photographic_composition"]["fraction_white_among_classified"],
          result["latitude_near_half_null"]["within_species_colour_preserving_null_p_unusually_even"],
          result["geographic_compensation_null"]["negative_compensation_one_sided_p"])


if __name__ == "__main__":
    main()

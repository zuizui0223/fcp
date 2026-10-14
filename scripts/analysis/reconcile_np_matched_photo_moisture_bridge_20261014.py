#!/usr/bin/env python3
"""Reconcile historical, identical-original-photo between- and within-species moisture tests.

This is a source-preserving reporting-only bridge. The original PR #158 and
PR #159 regression models, paired bootstrap intervals and H2 results are not
refitted or selected by this script; no photographic outcome is reopened.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
HERE=ROOT/"results/np_unified_20261014"
SOURCE_A=HERE/"source_pr158_matched_species_colour_rainfall.json"
SOURCE_B=HERE/"source_pr159_white_area_colour_rainfall.json"
DEST=HERE/"cross_scale_matched_photo_moisture_bridge.json"
ORIGINAL_ZIP="2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0"


def compute() -> dict:
    a=json.loads(SOURCE_A.read_text(encoding="utf8"))
    b=json.loads(SOURCE_B.read_text(encoding="utf8"))
    assert a["source"]["original_archive_sha256"]==ORIGINAL_ZIP
    assert b["source_sha256"]==ORIGINAL_ZIP
    assert a["source"]["original_phenotype_labels_changed"] is False
    assert a["matched_source_photo_and_species"]["same_photo_species_name_and_colour_outcomes_for_both_scales"]
    matched=a["matched_source_photo_and_species"]
    cont=b["models"]["min2_photos_two_cells"]
    q=matched["between_congener_BIO12_photo_white_slope"]
    r=matched["within_same_species_BIO12_photo_white_slope"]
    assert (matched["n_common_original_photos"],matched["n_common_species"],matched["n_source_genera"])==(21448,6143,787)
    assert math.isclose(q-r,matched["between_minus_within_slope"],abs_tol=1e-14)
    assert matched["paired_genus_bootstrap95_within"][0] < 0 < matched["paired_genus_bootstrap95_within"][1]
    assert matched["paired_genus_bootstrap95_difference"][0] > 0
    assert (cont["n_photos"],cont["n_species"],cont["n_genera"]) == (14667,3351,489)
    q2=cont["between_genus"]["rain_slope"]
    r2=cont["within_species"]["rain_slope"]
    assert math.isclose(q2-r2,cont["paired_between_minus_within"]["effect"],abs_tol=1e-14)
    qual=b["models"]["quality_lenient_1000px_nuisance_max50pct"]["between_genus"]
    assert qual["ci95"][0] < 0 < qual["ci95"][1]
    assert b["source"]["source_observation_date_column_present"] is False
    return {
        "schema":"fcp_np_source_matched_photo_cross_scale_bridge_20261014_v1",
        "status":"SOURCE_ONLY_ARITHMETIC_NO_NEW_PHOTO_CLASSIFICATION_OR_MODEL_REFIT",
        "source_authority":{
            "original_photo_zip_sha256":ORIGINAL_ZIP,
            "categorical_original_PR":158,
            "categorical_original_commit":"0bc53089dc9d1baed207c2ea898b130cb7b07566",
            "continuous_original_PR":159,
            "continuous_original_commit":"a26b826ae13727053b95413ca78a3713286a603e",
            "categorical_source_file":SOURCE_A.name,
            "continuous_source_file":SOURCE_B.name
        },
        "matched_photo_categorical_white":{
            "same_original_photo_count":matched["n_common_original_photos"],
            "same_named_species_count":matched["n_common_species"],
            "same_named_genera_count":matched["n_source_genera"],
            "rainfall_source_SD_mm":a["environmental_covariates"]["standardized_photo_site_BIO12_sd_mm"],
            "covariate_names":a["environmental_covariates"]["adjustment"],
            "outcome":"original same white/yellow-orange/red-pink/blue-purple photo classification in both level-specific fits",
            "between_congener_white_slope_pp_per_SD":100*q,
            "between_95CI_pp":[100*x for x in matched["paired_genus_bootstrap95_between"]],
            "within_species_white_slope_pp_per_SD":100*r,
            "within_95CI_pp":[100*x for x in matched["paired_genus_bootstrap95_within"]],
            "between_minus_within_pp_per_SD":100*(q-r),
            "paired_diff_95CI_pp":[100*x for x in matched["paired_genus_bootstrap95_difference"]],
            "original_bootstrap_replicates":matched["source_paired_genus_bootstrap_repeats"],
            "original_exploratory_one_sided_bootstrap_p":matched["paired_bootstrap_one_sided_p_difference_nonpositive"],
            "positive_within_species_white_rainfall_gradient_confirmed":False,
            "scales_not_same_causal_estimand":True
        },
        "matched_photo_continuous_white_roi_fraction":{
            "same_original_photo_count":cont["n_photos"],
            "same_named_species_count":cont["n_species"],
            "same_named_genera_count":cont["n_genera"],
            "rainfall_source_SD_mm":b["source_global_BIO12_sd_mm"],
            "outcome":"photographic segmented flower ROI fraction of 9 visible palette pixels called white",
            "between_congener_white_area_slope_pp_per_SD":100*q2,
            "between_95CI_pp":[100*x for x in cont["between_genus"]["ci95"]],
            "within_species_white_area_slope_pp_per_SD":100*r2,
            "within_95CI_pp":[100*x for x in cont["within_species"]["ci95"]],
            "between_minus_within_pp_per_SD":100*(q2-r2),
            "paired_diff_95CI_pp":[100*x for x in cont["paired_between_minus_within"]["genus_bootstrap_ci95"]],
            "original_exploratory_one_sided_bootstrap_p":cont["paired_between_minus_within"]["one_sided_p_nonpositive"],
            "quality_restricted_between_slope_pp_per_SD":100*qual["rain_slope"],
            "quality_restricted_between_95CI_pp":[100*x for x in qual["ci95"]],
            "strict_photo_quality_signal_not_confirmed":True,
            "not_spectrophotometry_or_pigment_concentration":True
        },
        "interpretive_decision":{
            "global_average_hue_association_is_not_a_common_within_species_photo_response":True,
            "same_photo_source_reduces_but_does_not_remove_measurement_geography_species_history_confounding":True,
            "source_selected_posthoc_not_independent_prospective_environment_confirm":True,
            "source_missingness_and_exposure_not_fully_corrected":True,
            "not_a_linked_mechanistic_fitness_test":True,
            "do_not_pool_different_SD_different_source_populations":True,
            "high_depth_158_species_W_remains_independent_of_this_source_model_fit":True
        }
    }


def main() -> None:
    p=argparse.ArgumentParser()
    p.add_argument("--write",action="store_true")
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    result=compute()
    if args.write:
        DEST.write_text(json.dumps(result,indent=2,sort_keys=True,ensure_ascii=False)+"\n",encoding="utf8")
    if args.check:
        assert result==json.loads(DEST.read_text(encoding="utf8")), "Cross-scale reporting ledger no longer matches original photo receipts"
    a=result["matched_photo_categorical_white"]
    print("SAME_PHOTO_CROSS_SCALE_SOURCE_AUDIT",
          a["same_original_photo_count"],a["same_named_species_count"],
          a["between_minus_within_pp_per_SD"])


if __name__=="__main__":
    main()

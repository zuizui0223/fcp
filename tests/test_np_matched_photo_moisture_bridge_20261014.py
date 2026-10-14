"""Reporting-only matched-photo global-to-within species moisture bridge tests."""
from __future__ import annotations
import importlib.util
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/"scripts/analysis/reconcile_np_matched_photo_moisture_bridge_20261014.py"
DEST=ROOT/"results/np_unified_20261014/cross_scale_matched_photo_moisture_bridge.json"
spec=importlib.util.spec_from_file_location("np_photo_bridge",SCRIPT)
bridge=importlib.util.module_from_spec(spec)
spec.loader.exec_module(bridge)

def result():
    actual=bridge.compute()
    committed=json.loads(DEST.read_text(encoding="utf8"))
    assert actual==committed
    return actual

def test_original_same_photo_between_within_bio12_contrast_and_source_scope():
    s=result()
    r=s["matched_photo_categorical_white"]
    assert (r["same_original_photo_count"],r["same_named_species_count"],r["same_named_genera_count"])==(21448,6143,787)
    assert r["rainfall_source_SD_mm"]>700
    assert math.isclose(r["between_congener_white_slope_pp_per_SD"],2.722186894224273,abs_tol=1e-11)
    assert math.isclose(r["within_species_white_slope_pp_per_SD"],-0.5106759199980051,abs_tol=1e-11)
    assert math.isclose(r["between_congener_white_slope_pp_per_SD"]-r["within_species_white_slope_pp_per_SD"],
                        r["between_minus_within_pp_per_SD"],abs_tol=1e-10)
    assert r["within_95CI_pp"][0]<0<r["within_95CI_pp"][1]
    assert r["paired_diff_95CI_pp"][0]>0
    assert r["scales_not_same_causal_estimand"] is True

def test_same_picture_fraction_white_not_calibrated_spectrum_and_quality_sensitivity():
    r=result()["matched_photo_continuous_white_roi_fraction"]
    assert (r["same_original_photo_count"],r["same_named_species_count"],r["same_named_genera_count"])==(14667,3351,489)
    assert r["rainfall_source_SD_mm"]>700
    assert r["between_congener_white_area_slope_pp_per_SD"]>0
    assert r["within_95CI_pp"][0]<0<r["within_95CI_pp"][1]
    assert r["quality_restricted_between_95CI_pp"][0]<0<r["quality_restricted_between_95CI_pp"][1]
    assert r["not_spectrophotometry_or_pigment_concentration"] is True

def test_cross_scale_photo_bridge_does_not_upgrade_prospective_selection():
    s=result()
    assert s["status"]=="SOURCE_ONLY_ARITHMETIC_NO_NEW_PHOTO_CLASSIFICATION_OR_MODEL_REFIT"
    d=s["interpretive_decision"]
    assert d["global_average_hue_association_is_not_a_common_within_species_photo_response"] is True
    assert d["source_selected_posthoc_not_independent_prospective_environment_confirm"] is True
    assert d["not_a_linked_mechanistic_fitness_test"] is True
    assert d["high_depth_158_species_W_remains_independent_of_this_source_model_fit"] is True

def test_manuscript_measurement_and_bridge_are_explicit():
    manuscript=(ROOT/"docs/POLYMORPHISM_MANUSCRIPT_NEW_PHYTOLOGIST_UNIFIED_20261014.md").read_text(encoding="utf8")
    intro=manuscript[manuscript.index("## Introduction"):manuscript.index("## Materials and Methods")]
    methods=manuscript[manuscript.index("## Materials and Methods"):manuscript.index("## Table 1.")]
    results=manuscript[manuscript.index("## Results"):manuscript.index("## Discussion")]
    for phrase in ("ROI-v4","EfficientSAM","Hellinger","white","bronze","classifiable","not"):
        assert phrase in methods,phrase
    assert "spectrophotometry" in methods.lower()
    assert "flower-colour polymorphism" in intro.lower()
    assert results.index("### Global photographed flower-colour diversity")<results.index("### Pattern: achromatic–chromatic")
    assert results.index("### Matched original photographs")<results.index("### Pattern: achromatic–chromatic")
    assert "21,448" in results and "6,143" in results and "3.233" in results
    assert "not a shared within-species response" in results.lower() or "within-species" in results.lower()

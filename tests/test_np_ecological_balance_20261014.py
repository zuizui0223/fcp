"""Audit reporting-only ecology-to-white-axis calculations from committed PR receipts."""
import importlib.util
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "scripts/analysis/build_np_ecological_balance_bridge_20261014.py"
OUTPUT = ROOT / "results/np_unified_20261014/ecological_balance_reporting_reconciliation.json"
spec = importlib.util.spec_from_file_location("ecological_balance_bridge", SOURCE)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_historical_source_receipts_are_reconciled_without_refitting():
    derived = module.compute()
    committed = json.loads(OUTPUT.read_text(encoding="utf8"))
    assert derived.keys() == committed.keys()
    for key, value in derived.items():
        assert committed[key].keys() == value.keys() if isinstance(value, dict) else committed[key] == value
    assert math.isclose(
        derived["photographic_composition"]["fraction_white_among_classified"],
        committed["photographic_composition"]["fraction_white_among_classified"],
        abs_tol=1e-12
    )
    assert derived["scientific_decision"] == committed["scientific_decision"]


def test_worldwide_near_half_is_not_fitness_equilibrium():
    r=module.compute()
    p=r["photographic_composition"]
    assert p["classifiable_photos"]==45953
    assert p["unclassifiable_photos"]==54590
    assert p["white_classifiable_photos"]==21574
    assert math.isclose(p["fraction_white_among_classified"],21574/45953,abs_tol=1e-14)
    assert p["no_assumption_all_photo_white_lower"] < 0.5 < p["no_assumption_all_photo_white_upper"]
    assert p["not_population_morph_prevalence"] is True
    l=r["latitude_near_half_null"]
    assert l["within_species_colour_preserving_null_p_unusually_even"]==0.214
    assert l["extraordinary_half_white_homeostasis_not_supported"]


def test_geographic_cancellation_mathematically_rejected_on_source_receipt():
    r=module.compute()["geographic_compensation_null"]
    assert r["composition_vs_local_covariance"] > 0
    assert r["negative_compensation_one_sided_p"] == 0.997
    assert r["reconstruction_absolute_error"] < 1e-12
    assert r["hypothesized_geographic_cancellation_not_supported"]


def test_rainfall_mnar_and_soil_on_identical_species_limit_mechanistic_claim():
    r=module.compute()
    m=r["moisture_opportunity_and_missingness"]
    assert 11 < m["unadjusted_wet_minus_dry_photo_white_pp"] < 12
    assert 0.40 < m["wet_unclassified_white_fraction_needed_to_erase_contrast"] < 0.42
    assert m["unobserved_outcome_colour_mnar_remains_identification_problem"]
    s=r["soil_added_on_exact_same_species_unseen_genus_test"]
    assert s["n_species"]==20786
    assert s["original_95_interval"][0]<0<s["original_95_interval"][1]
    assert s["supports_incremental_soil_information"] is False


def test_local_mixedness_and_prospective_w_remain_distinct():
    r=module.compute()
    q=r["intermediate_rainfall_local_cocomposition"]
    assert q["n_species_in_all_three_rainfall_pair_bins"]==0
    assert q["reported_nonhold_tests"]==5
    assert q["smallest_unadjusted_one_sided_p"]>0.05
    h=r["white_axis_separate_frozen_high_depth_source"]
    assert h["prospective_admissible_species"]==158
    assert h["observed_W"]>h["null_median_W"]
    assert h["cross_scale_morph_fitness_mechanism_tested"] is False

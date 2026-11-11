"""Synthetic seasonal-positive control; not a test of morph fitness."""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT=Path(__file__).resolve().parents[1]/"scripts/analysis/run_polymorphism_moricandia_seasonal_control_20261008.py"
spec=importlib.util.spec_from_file_location("moricandia_seasonal",SCRIPT)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def fixture(*, colour_reversal=False, no_summer=False):
    n=32
    season=["2020-04-15"]*8+["2021-05-15"]*8+["2020-07-15"]*8+["2021-08-15"]*8
    morph=(["blue_purple"]*16+["white"]*16)
    if colour_reversal:
        morph=(["white"]*16+["blue_purple"]*16)
    if no_summer:
        season=["2020-04-15"]*n
    return pd.DataFrame({
        "inat_taxon_id":[m.TARGET_TAXON_ID]*n,
        "species":[m.TARGET_SPECIES]*n,
        "photo_id":np.arange(n),"morph":morph,
        "global_classifiable":[True]*n,
        "latitude":10+np.linspace(0,0.01,n),
        "longitude":np.zeros(n),
        "observer_id":["observer_a","observer_b"]*16,
        "observed_on":season,
    })


def test_published_seasonal_direction_recovers_from_synthetic_photos():
    r=m.analyze(fixture())
    assert r["result_status"]=="positive_control_evaluable"
    assert r["colour_season_test"]["spring_white_fraction"]==0
    assert r["colour_season_test"]["summer_white_fraction"]==1
    assert r["colour_season_test"]["fisher_one_sided_white_enriched_in_summer_p_exploratory"]<0.001
    assert r["site_opportunity_northern_hemisphere"]["n_anchors_with_expected_colour_season_pair"]>0
    assert r["site_opportunity_northern_hemisphere"]["n_expected_pair_anchors_spanning_two_observation_years"]>0


def test_reverse_seasonal_direction_not_positive():
    r=m.analyze(fixture(colour_reversal=True))
    assert r["colour_season_test"]["summer_minus_spring_white_fraction"]<0
    assert r["colour_season_test"]["fisher_one_sided_white_enriched_in_summer_p_exploratory"]>0.5


def test_lack_of_summer_is_unestimable_not_negative_evidence():
    r=m.analyze(fixture(no_summer=True))
    assert r["result_status"]=="not_estimable_positive_control"
    assert not r["colour_season_test"]["estimable"]


def test_other_taxon_does_not_leak_into_defined_positive_control():
    d=fixture()
    d["inat_taxon_id"]=123
    r=m.analyze(d)
    assert r["n_classifiable_target_photos"]==0
    assert r["result_status"]=="not_estimable_positive_control"

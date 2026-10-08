"""Synthetic checks for multi-year floral morph evidence (not selection inference)."""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

S = Path(__file__).resolve().parents[1]/"scripts/analysis/run_polymorphism_local_multiyear_morphs_20261008.py"
sp=importlib.util.spec_from_file_location("multiyear",S)
m=importlib.util.module_from_spec(sp)
sp.loader.exec_module(m)


def fixture(*, spatially_separated=False, same_observer=False, single_year=False):
    # 12 W, 12 red, 12 blue, 12 yellow; every morph repeated each year
    # and across multiple independent observers.
    labels=(["white"]*12 + ["red_pink"]*12
            + ["blue_purple"]*12 + ["yellow_orange"]*12)
    n=len(labels)
    if spatially_separated:
        lat=np.r_[np.zeros(12)+10, np.zeros(36)+15]
    else:
        lat=10.0+np.linspace(0,0.005,n)
    obs=(["person_a","person_b"]*6)*4
    if same_observer:
        obs=["single_user"]*n
    years=([2018]*6+[2020]*6)*4
    if single_year:
        years=[2020]*n
    return pd.DataFrame({
        "latitude":lat, "longitude":np.zeros(n),
        "morph":labels,"year":np.array(years,dtype=float),
        "observer":obs,"inat_taxon_id":[1]*n,
        "species":["Synthetic example"]*n,
        "photo_id":np.arange(n),
        "cohort":["synthetic"]*n,
    })


def test_multiyear_multiobserver_colocalization():
    d=fixture()
    n=m.strict_anchor_neighborhoods(d,10.0)
    assert n["white_nonwhite_multi_year_two_observers"] is True
    assert n["two_nonwhite_hues_multi_year_two_observers"] is True
    assert n["n_white_nonwhite_anchors"]>0
    assert n["n_nonwhite_hue_anchors"]>0
    row=m.summarize_species(d)
    assert row["both_comparisons_available"]
    assert row["diameter_10km_white_nonwhite_multi_year_two_observers"]


def test_distant_sites_cannot_count_as_local_white_colour():
    d=fixture(spatially_separated=True)
    n=m.strict_anchor_neighborhoods(d,50.0)
    assert not n["white_nonwhite_same_neighborhood"]
    assert not n["white_nonwhite_multi_year_two_observers"]
    assert n["two_nonwhite_hues_multi_year_two_observers"]


def test_one_year_and_single_observer_fail_strict_gate():
    for kwargs in ({"single_year":True},{"same_observer":True}):
        d=fixture(**kwargs)
        res=m.strict_anchor_neighborhoods(d,50)
        assert res["white_nonwhite_same_neighborhood"]
        assert not res["white_nonwhite_multi_year_two_observers"]
        assert not res["two_nonwhite_hues_multi_year_two_observers"]


def test_temporal_pair_gate_fail_closed_and_true_when_supported():
    d=fixture()
    assert m.temporal_colour_discordance(d)["temporal_pair_comparison_evaluable"]
    dd=d.copy()
    dd["year"]=2020
    assert not m.temporal_colour_discordance(dd)["temporal_pair_comparison_evaluable"]
    assert m.pair_distance_km(np.array([0, 0]),np.array([0, 0]))[0,1]==0


def test_exact_matched_species_control_and_bins():
    d=m.summarize_species(fixture())
    frame=pd.DataFrame([d])
    summary=m.describe(frame,"synthetic")
    assert summary["n_white_plus_nonwhite_5_each"]==1
    assert summary["n_both_comparisons_available"]==1
    for distance in ("10","25","50"):
        gate=summary["diameters_km"][distance]["joint_opportunity_exact_mcnemar"]
        assert gate["n_matched_species"]==1
        assert sum(gate[k] for k in (
            "both_types_recur","only_white_plus_nonwhite_recur",
            "only_two_nonwhite_hues_recur","neither_recur"))==1


def test_invalid_years_and_observer_missingness_do_not_pass():
    d=fixture()
    d["observer"]=""
    assert not m.strict_anchor_neighborhoods(d,25)["white_nonwhite_multi_year_two_observers"]
    d=fixture()
    d["year"]=np.nan
    assert not m.strict_anchor_neighborhoods(d,25)["white_nonwhite_multi_year_two_observers"]

"""Synthetic fail-closed tests: contemporaneous photographic states, not fitness."""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

SOURCE = Path(__file__).resolve().parents[1] / "scripts/analysis/run_polymorphism_local_contemporaneity_20261008.py"
spec = importlib.util.spec_from_file_location("fcp_contemporaneous", SOURCE)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def sample(
    *,
    shift_white_days=0,
    white_different_year=False,
    split_sites=False,
    one_observer_white=False,
    same_observer_all=False,
):
    # Four categories with 12 photographs each; labels and dates fixed
    # without any field observation or real biological genotype.
    hues = (["white"] * 12 + ["red_pink"] * 12
            + ["blue_purple"] * 12 + ["yellow_orange"] * 12)
    n = len(hues)
    dates = pd.to_datetime(["2020-05-10"]*12 + ["2020-05-12"]*12 +
                           ["2020-05-13"]*12 + ["2020-05-14"]*12)
    dates = pd.Series(dates)
    if shift_white_days:
        dates.loc[:11] += pd.Timedelta(days=shift_white_days)
    if white_different_year:
        dates.loc[:11] -= pd.DateOffset(years=1)
    lat = 10.0 + np.linspace(0, 0.01, n)
    if split_sites:
        lat[:12] += 3.0
    observers = np.asarray(["w1", "w2"]*6 + ["r1", "r2"]*6
                           + ["b1", "b2"]*6 + ["y1", "y2"]*6, dtype=object)
    if one_observer_white:
        observers[:12] = "w_only"
    if same_observer_all:
        observers[:] = "same_photographer"
    return pd.DataFrame({
        "inat_taxon_id": np.repeat(42, n),
        "species": ["Synthetic flower"] * n,
        "photo_id": np.arange(n, dtype=int),
        "morph": hues,
        "latitude": lat,
        "longitude": np.zeros(n),
        "observer": observers,
        "year": dates.dt.year.to_numpy(float),
        "day": dates.dt.dayofyear.to_numpy(float),
        "cohort": ["synthetic"] * n,
    })


def test_close_same_year_mixed_colours_strict_both_categories():
    r = m.one_species_evidence(sample())
    assert r["eligible_white_nonwhite"]
    assert r["eligible_two_nonwhite_hues"]
    for scenario in ("d10_t14", "d10_t30", "d50_t14"):
        assert r[f"{scenario}_has_colour_blind_four_photo_local_year_window"]
        assert r[f"{scenario}_has_strict_white_colour_contemporaneity"]
        assert r[f"{scenario}_has_strict_two_nonwhite_contemporaneity"]
        assert r[f"{scenario}_has_independent_observer_white_colour_pair"]
    assert r["d10_t14_n_years_strict_white_colour"] == 1


def test_photographic_seasonal_replacement_not_contemporaneous():
    r = m.one_species_evidence(sample(shift_white_days=80))
    assert r["eligible_white_nonwhite"]
    assert r["d10_t14_has_colour_blind_four_photo_local_year_window"]
    assert not r["d10_t14_has_independent_observer_white_colour_pair"]
    assert not r["d10_t30_has_strict_white_colour_contemporaneity"]
    assert r["d10_t14_has_strict_two_nonwhite_contemporaneity"]


def test_year_disjunction_not_mistaken_for_same_year_coexistence():
    r = m.one_species_evidence(sample(white_different_year=True))
    assert not r["d50_t30_has_strict_white_colour_contemporaneity"]
    assert r["d10_t14_has_strict_two_nonwhite_contemporaneity"]


def test_geographic_separation_not_called_local_overlap():
    r = m.one_species_evidence(sample(split_sites=True))
    assert not r["d50_t30_has_independent_observer_white_colour_pair"]
    assert not r["d50_t30_has_strict_white_colour_contemporaneity"]


def test_one_white_observer_can_make_weak_not_strict_pair():
    r = m.one_species_evidence(sample(one_observer_white=True))
    assert r["d10_t14_has_independent_observer_white_colour_pair"]
    assert not r["d10_t14_has_strict_white_colour_contemporaneity"]
    assert r["d10_t14_has_strict_two_nonwhite_contemporaneity"]


def test_single_observer_everywhere_fails_both_gates():
    r = m.one_species_evidence(sample(same_observer_all=True))
    assert not r["d10_t14_has_colour_blind_four_photo_local_year_window"]
    assert not r["d10_t14_has_independent_observer_white_colour_pair"]
    assert not r["d10_t14_has_strict_white_colour_contemporaneity"]


def test_prior_ledger_provenance_join_and_opportunity_denominators(tmp_path):
    s = m.one_species_evidence(sample())
    frame = pd.DataFrame([s])
    prior = pd.DataFrame([{
        "cohort": "synthetic", "species": "Synthetic flower",
        "inat_taxon_id": 42,
        "diameter_10km_white_nonwhite_multi_year_two_observers": True,
        "diameter_25km_white_nonwhite_multi_year_two_observers": True,
        "diameter_50km_white_nonwhite_multi_year_two_observers": True,
    }])
    path = tmp_path/"prior.csv"
    prior.to_csv(path, index=False)
    out = m.summary_one_cohort(frame, m.read_prior(path))
    primary = out["scales"]["d10_t14"]
    assert primary["n_species_wc_globally_eligible_and_four_photo_opportunity"] == 1
    assert primary["n_strict_wc_two_photos_two_observers_each"] == 1
    assert primary["prior_multiyear_with_strict_contemporaneous_wc"] == 1
    assert primary["paired_eligible_both_contemporaneous"] == 1
    assert sum(primary[x] for x in (
        "paired_eligible_both_contemporaneous",
        "paired_eligible_only_wc_contemporaneous",
        "paired_eligible_only_hue_contemporaneous",
        "paired_eligible_neither_contemporaneous"
    )) == primary["n_species_both_comparisons_and_opportunity"]


def test_missing_prior_species_must_fail():
    s = pd.DataFrame([m.one_species_evidence(sample())])
    wrong_prior = pd.DataFrame([{
        "cohort": "synthetic", "inat_taxon_id": 999,
        **{f"diameter_{int(d)}km_white_nonwhite_multi_year_two_observers": False
           for d in m.DIAMETERS_KM}
    }])
    with pytest.raises(ValueError, match="missing a species"):
        m.summary_one_cohort(s, wrong_prior)


def test_bounded_neighborhood_both_axes_fail_outside_thresholds():
    r = sample()
    r.loc[12:, "latitude"] += 0.08   # > half 10 km diameter
    r.loc[:11, "day"] += 8            # > half 14-day window
    output = m.one_species_evidence(r)
    assert not output["d10_t14_has_strict_white_colour_contemporaneity"]


def test_independence_logic_handles_missing_observers():
    codes = np.array(["", "w2", "r1", "r2"], dtype=str)
    maskw = np.array([True, True, False, False])
    maskc = ~maskw
    assert not m.two_independent_observers(maskw, codes)
    assert m.two_independent_observers(maskc, codes)
    assert m.pair_has_independent_observers(maskw, maskc, codes)

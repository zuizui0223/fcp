"""Synthetic non-phenotype, stable-site, and repeated-year audit guards."""
from __future__ import annotations

import sys
from pathlib import Path
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "analysis"))
from audit_fcp_siteyear_same_month_opportunity_20261008 import best_site_month


def fixture_rows():
    # Same approximate location, same May flowering month, two distinct years,
    # at least two separate observers per year; one far photo is irrelevant.
    return pd.DataFrame([
        {"photo_id": "1", "latitude": 35.000, "longitude": 139.000, "year": 2021, "month": 5, "observer": "a"},
        {"photo_id": "2", "latitude": 35.001, "longitude": 139.001, "year": 2021, "month": 5, "observer": "b"},
        {"photo_id": "3", "latitude": 35.002, "longitude": 139.002, "year": 2023, "month": 5, "observer": "c"},
        {"photo_id": "4", "latitude": 35.003, "longitude": 139.003, "year": 2023, "month": 5, "observer": "d"},
        {"photo_id": "5", "latitude": 36.000, "longitude": 140.000, "year": 2021, "month": 5, "observer": "e"},
    ])


def test_same_month_different_year_two_observers():
    r = best_site_month(fixture_rows(), 10.0)
    assert r is not None
    assert r["month"] == 5
    assert r["years"] == [2021, 2023]
    assert r["n_photos_qualified_site_month_years"] == 4
    assert r["maximum_within_site_pair_distance_km"] <= 10.0


def test_same_year_only_is_insufficient():
    d = fixture_rows()
    d["year"] = 2023
    assert best_site_month(d, 10.0) is None


def test_different_month_not_wrongly_paired():
    d = fixture_rows()
    d.loc[d.photo_id.isin(["3", "4"]), "month"] = 6
    assert best_site_month(d, 10.0) is None


def test_same_observer_repeated_is_insufficient():
    d = fixture_rows()
    d.loc[d.photo_id.isin(["3", "4"]), "observer"] = "same"
    assert best_site_month(d, 10.0) is None


def test_no_morph_input_or_colour_dependent_site_selection():
    d = fixture_rows()
    d["morph"] = ["white", "white", "red_pink", "red_pink", "white"]
    a = best_site_month(d, 10.0)
    d["morph"] = ["blue_purple"] * len(d)
    b = best_site_month(d, 10.0)
    assert a == b


def test_long_chain_cannot_fake_10km_diameter():
    # 2 photos each year in a 10km-spread pair; cross-year sites >10km apart.
    d = fixture_rows().iloc[:4].copy()
    d.loc[d.photo_id.isin(["3", "4"]), "longitude"] += 0.2
    assert best_site_month(d, 10.0) is None


def test_missing_observers_do_not_count_as_independent():
    d = fixture_rows()
    d.loc[d.photo_id == "4", "observer"] = ""
    assert best_site_month(d, 10.0) is None


def test_bad_radius_rejected():
    with pytest.raises(ValueError, match="Radius"):
        best_site_month(fixture_rows(), 0)

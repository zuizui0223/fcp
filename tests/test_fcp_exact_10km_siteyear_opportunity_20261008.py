"""Exact 10km site-year-month minimal witness tests (retrospective geometry only)."""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"/"analysis"))
from audit_fcp_exact_10km_siteyear_opportunity_20261008 import exact_four_photo_witness
from audit_fcp_siteyear_same_month_opportunity_20261008 import best_site_month


def square():
    # Four corners: horizontal/vertical ~6.67km, diagonal ~9.44km.
    # No sample-photo anchor within 5km of all 4 corners.
    return pd.DataFrame([
        {"photo_id":"a","latitude":0.00,"longitude":0.00,"year":2021,"month":5,"observer":"p"},
        {"photo_id":"b","latitude":0.00,"longitude":0.06,"year":2021,"month":5,"observer":"q"},
        {"photo_id":"c","latitude":0.06,"longitude":0.00,"year":2023,"month":5,"observer":"r"},
        {"photo_id":"d","latitude":0.06,"longitude":0.06,"year":2023,"month":5,"observer":"s"},
    ])


def test_exact_recovers_square_anchored_method_misses():
    a=square()
    assert best_site_month(a,10) is None
    found=exact_four_photo_witness(a,10)
    assert found is not None
    assert found["years"]==[2021,2023]
    assert found["max_diameter_km"]<=10
    assert len(set(found["photo_ids"]))==4


def test_exact_does_not_allow_far_diagonal():
    d=square()
    d.loc[d.photo_id=="d","longitude"]=0.20
    assert exact_four_photo_witness(d,10) is None


def test_exact_does_not_pool_different_months():
    d=square()
    d.loc[d.year==2023,"month"]=6
    assert exact_four_photo_witness(d,10) is None


def test_exact_requires_two_distinct_observers_per_year():
    d=square()
    d.loc[d.year==2023,"observer"]="one"
    assert exact_four_photo_witness(d,10) is None


def test_extra_third_year_is_allowed_but_not_necessary():
    d=square()
    d=pd.concat([d,d.assign(photo_id=lambda x:x.photo_id+"2",year=2024)],ignore_index=True)
    found=exact_four_photo_witness(d,10)
    assert found and len(found["years"])==2


def test_morph_changes_do_not_affect_exact_geometry():
    d=square()
    d["morph"]=["white","white","red_pink","red_pink"]
    first=exact_four_photo_witness(d)
    d["morph"]=["blue_purple"]*len(d)
    second=exact_four_photo_witness(d)
    assert first==second


def test_bad_radius_rejected():
    with pytest.raises(ValueError,match="Diameter"):
        exact_four_photo_witness(square(),0)

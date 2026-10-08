"""Observer-slot shortfall counts are outcome-blind and constrained to observed years."""
import sys
from pathlib import Path
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"/"analysis"))
from audit_fcp_siteyear_replenishment_gap_20261008 import min_observer_gap


def points(a="a", b="b", c="c", d="d"):
    return pd.DataFrame([
        {"photo_id":"11","latitude":35.0,"longitude":139.0,"year":2020,"month":5,"observer":a},
        {"photo_id":"12","latitude":35.001,"longitude":139.0,"year":2020,"month":5,"observer":b},
        {"photo_id":"21","latitude":35.0,"longitude":139.001,"year":2022,"month":5,"observer":c},
        {"photo_id":"22","latitude":35.0,"longitude":139.002,"year":2022,"month":5,"observer":d},
    ])


def test_complete_zero_slots_and_original_eligibility():
    assert min_observer_gap(points())["gap_min_additional_distinct_observer_photos"]==0


def test_one_missing_distinct_observer_not_one_missing_photo():
    x=points()
    x.loc[x.photo_id=="22","observer"]="c"
    r=min_observer_gap(x)
    assert r["gap_min_additional_distinct_observer_photos"]==1
    assert r["existing_photos_per_year"]==[2,2]
    assert r["existing_distinct_observers_per_year"]==[2,1]


def test_two_years_one_observer_each_requires_two_slots():
    r=min_observer_gap(points("a","a","c","c"))
    assert r["gap_min_additional_distinct_observer_photos"]==2


def test_two_years_without_known_observer_requires_four_slots():
    r=min_observer_gap(points("","","",""))
    assert r["gap_min_additional_distinct_observer_photos"]==4


def test_no_two_years_of_same_month_is_unobserved_not_zero_gap():
    x=points()
    x.loc[x.year==2022,"month"]=6
    assert min_observer_gap(x) is None


def test_no_long_distance_links():
    x=points()
    x.loc[x.year==2022,"longitude"]=141.0
    assert min_observer_gap(x) is None


def test_outcome_morph_does_not_change_selection():
    x=points()
    x["morph"]=["white","white","red_pink","red_pink"]
    a=min_observer_gap(x)
    x["morph"]=["blue_purple"]*len(x)
    assert min_observer_gap(x)==a


def test_minimum_across_three_existing_years():
    x=points("a","a","b","b")
    new=points("a","a","c","d").iloc[2:].copy()
    new["year"]=2024
    new["photo_id"]=["31","32"]
    x=pd.concat([x,new],ignore_index=True)
    r=min_observer_gap(x)
    assert r["gap_min_additional_distinct_observer_photos"]==1
    assert r["years"]==[2020,2024] or r["years"]==[2022,2024]

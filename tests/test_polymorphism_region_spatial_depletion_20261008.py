"""Synthetic controls for region-specific FCP spatial structure and generality gating."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

SOURCE=Path(__file__).resolve().parents[1]/"scripts/analysis/run_polymorphism_region_spatial_depletion_20261008.py"
sp=importlib.util.spec_from_file_location("fcp_regional_geographic_itv",SOURCE)
m=importlib.util.module_from_spec(sp)
sp.loader.exec_module(m)


def make_species(taxon=101, pattern="geographic", region_lat=10.0,
                 same_observer=False, unique_years=False):
    # 80 photographs across two sites separated by >50km;
    # both months occur at each site, unlike season-only sorting.
    n=80
    lat=np.r_[np.full(40,region_lat),np.full(40,region_lat+1.8)]
    mon=np.tile([4]*20+[7]*20,2)
    yr=np.full(n,2020,int)
    if unique_years:
        yr=np.arange(n)+1900
    if pattern=="geographic":
        labels=np.r_[np.full(40,"white"),np.full(40,"red_pink")]
    elif pattern=="month":
        labels=np.where(mon==4,"white","red_pink")
    elif pattern=="no_difference":
        labels=np.full(n,"white")
    elif pattern=="no_sort":
        labels=np.tile(["white","red_pink"],n//2)
    else:
        raise ValueError(pattern)
    if same_observer:
        users=["same_observer"]*n
    else:
        users=["user_"+str(i%10) for i in range(n)]
    return pd.DataFrame({
        "inat_taxon_id":np.repeat(taxon,n),"photo_id":np.arange(n)+taxon*1000,
        "species":[f"Genus{taxon} example"]*n,
        "latitude":lat,"longitude":np.zeros(n),
        "abs_latitude":np.abs(lat),
        "morph":labels,
        "month":pd.Series(mon,dtype="Int64"),
        "year":pd.Series(yr,dtype="Int64"),
        "observer":users,
        "cohort":["synthetic"]*n
    })


def test_all_latitude_edges_and_hemisphere_absolute():
    a=np.array([-90,-60,-30,-23.5,-10,0,15,23.5,30,60,90],float)
    p=m.bands_for(a,"primary")
    assert p.tolist()==[2,2,1,0,0,0,0,0,1,2,2]
    q=m.bands_for(a,"geographic_tropics_sensitivity")
    assert q.tolist()==[2,2,1,1,0,0,0,1,1,2,2]


def test_geographic_partition_survives_fixed_month_photo_composition():
    d=make_species(pattern="geographic")
    for pol in m.PAIR_POLICIES:
        a=m.region_species_test(d,"synthetic","primary","low_0_30","month",pol)
        assert a is not None
        row,null=a
        assert row["identifiable_after_calendar"]
        assert row["regional_observed_local_depletion"]>0.48
        assert row["regional_residual_depletion"]>0.30
        assert row["n_local_pairs"]>=m.MIN_NEAR_PAIRS
        assert row["n_distant_pairs"]>=m.MIN_FAR_PAIRS
        assert null.shape==(m.PERMUTATIONS,)
        assert np.std(null)>0


def test_purely_seasonal_sorting_does_not_fake_residual_geography():
    d=make_species(pattern="month")
    for mode in m.CONDITIONAL_TIME:
        a=m.region_species_test(d,"synthetic","primary","low_0_30",mode,"all")
        assert a is not None
        row,null=a
        assert not row["identifiable_after_calendar"]
        assert row["regional_residual_depletion"]==0.0
        assert np.allclose(null,row["regional_observed_local_depletion"])


def test_single_colour_species_stays_in_regional_evaluable_denominator_as_zero():
    d=make_species(pattern="no_difference")
    row,null=m.region_species_test(
        d,"synthetic","primary","low_0_30","month","all")
    assert row["n_photo_nonwhite"]==0
    assert not row["identifiable_after_calendar"]
    assert row["regional_residual_depletion"]==0
    assert len(null)==199


def test_one_observer_cannot_pass_distinct_observer_pairs():
    d=make_species(same_observer=True)
    assert m.region_species_test(d,"synthetic","primary","low_0_30",
                                 "month","different_observer") is None
    assert m.region_species_test(d,"synthetic","primary","low_0_30",
                                 "month","all") is not None


def test_yearmonth_can_be_unidentifiable_but_month_conditioned_informative():
    d=make_species(unique_years=True)
    a=m.region_species_test(d,"synthetic","primary","low_0_30","month","all")
    b=m.region_species_test(d,"synthetic","primary","low_0_30","year_month","all")
    assert a[0]["identifiable_after_calendar"] and a[0]["regional_residual_depletion"]>0
    assert not b[0]["identifiable_after_calendar"]
    assert b[0]["regional_residual_depletion"]==0


def test_same_species_across_region_boundary_does_not_fake_regional_structure():
    d=make_species(pattern="geographic",region_lat=29.0)
    # Each region contains a single (different) site; with no distant same-
    # region pairs the regional metric is non-estimable rather than positive.
    d["abs_latitude"]=np.abs(d.latitude)
    region=m.bands_for(d.latitude.to_numpy(float),"primary")
    low=d.loc[region==0].copy()
    mid=d.loc[region==1].copy()
    assert len(low)==40 and len(mid)==40
    assert m.region_species_test(low,"synthetic","primary","low_0_30",
                                 "month","all") is None
    assert m.region_species_test(mid,"synthetic","primary","middle_30_60",
                                 "month","all") is None


def test_species_level_equal_weights_and_zero_nonexchangeable():
    high=m.region_species_test(make_species(101),"synthetic",
        "primary","low_0_30","month","all")
    zero=m.region_species_test(make_species(102,pattern="month"),"synthetic",
        "primary","low_0_30","month","all")
    out=m.summarize([high[0],zero[0]],[high[1],zero[1]],
                    "synthetic","primary","low_0_30","month","all")
    assert out["n_geographically_evaluable_species"]==2
    assert out["n_calendar_identifiable_species"]==1
    assert out["n_nonidentifiable_kept"]==1
    assert np.isclose(2*out["mean_residual_depletion"],
                      out["mean_identifiable_subset_residual"])
    assert not out["coverage_pass"]


def test_holm_and_three_cohort_fail_closed():
    assert np.allclose(m.holm([0.01,0.03,0.20]),[0.03,0.06,0.20])
    c=make_species(123)
    r,rows=m.analyze_cohort(c,"synthetic")
    low=r["schemes"]["primary"]["regions"]["low_0_30"]
    assert low["opportunity"]["n_species_with_any_region_photos"]==1
    assert low["scenarios"]["month__all"]["n_geographically_evaluable_species"]==1
    assert low["scenarios"]["month__all"]["status"]=="HOLD_LIMITED_REGION_OR_CONDITIONAL_SUPPORT"
    assert low["scenarios"]["month__all"]["holm_p_across_3regions_x_2observers"]==1.0
    assert len(rows)>0
    output={"cohorts":{k:r for k in m.SHA256}}
    decision=m.cross_cohort_decision(output)
    assert decision["primary"]["low_0_30"]["status"]=="HOLD_REGION_COVERAGE"
    assert decision["primary"]["high_60_90"]["status"]=="HOLD_REGION_COVERAGE"


def test_determinism_of_perm_and_spatial_nulls():
    d=make_species()
    a=m.region_species_test(d,"synthetic","primary","low_0_30","month","all")
    b=m.region_species_test(d,"synthetic","primary","low_0_30","month","all")
    assert a[0]==b[0] and np.array_equal(a[1],b[1])


def test_no_geographical_photo_opportunity_returns_hold():
    z=m.summarize([],[],"synthetic","primary","high_60_90","month","all")
    assert not z["estimable"]
    assert z["status"]=="HOLD_NO_REGIONAL_GEOGRAPHIC_COVERAGE"
    assert not z["coverage_pass"]

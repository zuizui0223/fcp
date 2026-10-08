"""Species-wide space-after-season inference: falsification and reproducibility tests."""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/analysis/run_polymorphism_specieswide_space_season_20261008.py"
sp = importlib.util.spec_from_file_location("space_season", SCRIPT)
mod = importlib.util.module_from_spec(sp)
sp.loader.exec_module(mod)


def sample(*, mechanism="space", same_observer=False, n=80):
    assert n == 80
    # 40 photos per site, 2 sites separated by ~333km. Both sites have
    # repeated April and July photographs to avoid false phenology.
    lat = np.r_[np.full(40, 10.0), np.full(40, 13.0)]
    month = np.tile(np.r_[np.full(20, 4), np.full(20, 7)], 2)
    year = np.full(n, 2020)
    if mechanism == "space":
        morph = np.array(["white"]*40 + ["red_pink"]*40)
    elif mechanism == "time":
        morph = np.where(month == 4, "white", "red_pink")
    elif mechanism == "constant":
        morph = np.repeat("white", n)
    else:
        raise ValueError(mechanism)
    observers = ["obs_" + str(i % 12) for i in range(n)]
    if same_observer:
        observers = ["the_same_observer"] * n
    return pd.DataFrame({
        "inat_taxon_id": np.repeat(44, n),
        "species": ["Synthetic angiosperm"] * n,
        "observer": observers,
        "latitude": lat,
        "longitude": np.zeros(n),
        "morph": morph,
        "month": pd.Series(month,dtype="Int64"),
        "quarter": pd.Series((month-1)//3+1,dtype="Int64"),
        "year": pd.Series(year,dtype="Int64"),
        "photo_id": np.arange(n,dtype=int),
    })


def test_synthetic_true_space_partition_survives_month_conditioning():
    g=sample(mechanism="space")
    for mode in mod.STRATA:
        output = mod.species_test(g,"synthetic",mode,"all")
        assert output is not None, mode
        row,null=output
        assert row["observed_local_depletion"]>0.45
        assert row["excess_over_stratified_null"]>0.25
        assert np.std(null)>0
        assert row["n_geographic_local_pairs"]>=mod.MIN_LOCAL_PAIRS


def test_season_only_signal_not_called_residual_space_structure():
    g=sample(mechanism="time")
    a=mod.species_test(g,"synthetic","unconditional","all")
    assert a is not None
    b=mod.species_test(g,"synthetic","month","all")
    # Month-specific photographs carry a single fixed colour per stratum:
    # no label exchange opportunity; result must be HOLD, not "null proved".
    assert b is not None
    assert not b[0]["conditional_identifiable"]
    assert b[0]["excess_over_stratified_null"] == 0
    assert a[0]["observed_local_depletion"]<0.05


def test_colour_constant_null_not_informative():
    g=sample(mechanism="constant")
    for scenario, policy in [("month","all"),("year_month","different_observer")]:
        x=mod.species_test(g,"synthetic",scenario,policy)
        assert x is not None
        assert not x[0]["conditional_identifiable"]
        assert x[0]["excess_over_stratified_null"] == 0


def test_observer_filter_and_partial_dates():
    g=sample(same_observer=True)
    assert mod.species_test(g,"synthetic","month","all") is not None
    assert mod.species_test(g,"synthetic","month","different_observer") is None
    g=sample()
    g.loc[:45,"month"]=pd.NA
    # Date complete-case falls below 40 and must fail.
    assert mod.species_test(g,"synthetic","month","all") is None


def test_composition_and_month_group_fixed_under_every_permutation():
    g=sample(mechanism="space")
    labels=pd.Categorical(g.morph,categories=mod.COLOURS).codes.astype(np.int8)
    indices=mod.groups_for(g,"month")
    rng=np.random.default_rng(123)
    for _ in range(20):
        p=mod.permute_strata(labels,indices,rng)
        assert sorted(p)==sorted(labels)
        for idx in indices:
            assert sorted(p[idx])==sorted(labels[idx])


def test_month_year_can_be_separately_unidentified():
    g=sample()
    # Every photograph has a different year-month. Month-of-year across
    # years retains 20+20 matches; strict month×year cannot permute labels.
    g["year"]=pd.Series(np.arange(80)+1900,dtype="Int64")
    assert mod.species_test(g,"synthetic","month","all") is not None
    x=mod.species_test(g,"synthetic","year_month","all")
    assert x is not None
    assert not x[0]["conditional_identifiable"]
    assert x[0]["excess_over_stratified_null"] == 0


def test_reproducibility_and_cohort_mean():
    g=sample()
    a=mod.species_test(g,"synthetic","month","different_observer")
    b=mod.species_test(g,"synthetic","month","different_observer")
    assert a is not None and b is not None
    assert a[0]==b[0]
    assert np.array_equal(a[1],b[1])
    s=mod.summarize([a[0]],[a[1]],"synthetic","month","different_observer")
    assert s["n_species"]==1
    assert s["n_identifiable_species"]==1
    assert s["n_nonidentified_species_kept"]==0
    assert s["status"]=="coverage_limited_diagnostic_only"
    assert s["mean_excess_over_season_stratified_null"]>0
    assert s["permutation_p_upper"]<=0.05


def test_empty_mode_is_not_promoted_to_result():
    x=mod.summarize([],[],"synthetic","month","all")
    assert not x["estimable"]
    assert not x["meets_frozen_80_species_coverage_gate"]
    assert x["status"]=="not_estimable_no_geographic_photo_coverage"


def test_nonexchangeable_month_is_retained_in_inclusive_mean():
    # Two species: one geographically separated by morph regardless of month,
    # the other has morph solely controlled by month. The latter must contribute
    # exactly zero additional geographic signal, not be dropped from denominator.
    g1=sample(mechanism="space")
    g2=sample(mechanism="time").copy()
    g2["inat_taxon_id"]=45
    yes=mod.species_test(g1,"synthetic","month","all")
    no=mod.species_test(g2,"synthetic","month","all")
    assert yes is not None and no is not None
    assert yes[0]["conditional_identifiable"]
    assert not no[0]["conditional_identifiable"]
    d=mod.summarize([yes[0],no[0]],[yes[1],no[1]],"synthetic","month","all")
    assert d["n_species"]==2
    assert d["n_identifiable_species"]==1
    assert d["n_nonidentified_species_kept"]==1
    assert d["mean_excess_over_season_stratified_null"]>0
    assert abs(2*d["mean_excess_over_season_stratified_null"]-
               d["mean_identifiable_subset_excess_sensitivity"])<1e-10
    assert not d["meets_frozen_80_species_coverage_gate"]

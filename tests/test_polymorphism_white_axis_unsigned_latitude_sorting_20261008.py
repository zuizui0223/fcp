"""Distinguish chromatic-pigment candidate latitude sorting from generic coloured-hue sorting."""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]/"scripts/analysis"
sys.path.insert(0,str(ROOT))
spec=importlib.util.spec_from_file_location("fcp_white_unsigned_sort",ROOT/"run_polymorphism_white_axis_unsigned_latitude_sorting_20261008.py")
a=importlib.util.module_from_spec(spec)
spec.loader.exec_module(a)


def make_plant(pattern="white_axis",same_lat=False,one_month=False):
    rng=np.random.default_rng(42)
    n=120
    lats=np.repeat([10.,30.,50.],40)
    if same_lat: lats[:]=30.
    if pattern=="white_axis":
        blocks=[
            ["white"]*30+["red_pink"]*5+["blue_purple"]*5,
            ["white"]*5+["red_pink"]*17+["blue_purple"]*18,
            ["white"]*5+["red_pink"]*18+["blue_purple"]*17
        ]
    elif pattern=="hue_axis":
        blocks=[
            ["white"]*13+["red_pink"]*25+["blue_purple"]*2,
            ["white"]*14+["red_pink"]*13+["blue_purple"]*13,
            ["white"]*13+["red_pink"]*2+["blue_purple"]*25
        ]
    elif pattern=="month_only":
        blocks=[
            ["white"]*40,
            ["red_pink"]*40,
            ["blue_purple"]*40,
        ]
    else: raise ValueError(pattern)
    labels=[]; months=[]
    for k,grp in enumerate(blocks):
        grp=np.array(grp,object)
        if pattern!="month_only":
            rng.shuffle(grp)
            months.extend(([4,7]*20))
        else:
            months.extend(([4]*40) if k==0 else ([7]*40))
        labels.extend(grp.tolist())
    return pd.DataFrame({
        "inat_taxon_id":np.repeat(101,n),
        "photo_id":np.arange(n),
        "species":["Synthetic plant"]*n,
        "cohort":["synthetic"]*n,
        "latitude":lats,
        "morph":labels,
        "month":pd.Series(months,dtype="Int64"),
    })


def test_bin_edges_no_fake_separation_of_tied_localities():
    x=make_plant()
    b,meta=a.latitude_bins(x)
    assert meta["valid"]
    assert len(np.unique(b))==3
    assert all(len(set(b[x.latitude==z]))==1 for z in [10.,30.,50.])
    x=make_plant(same_lat=True)
    b,meta=a.latitude_bins(x)
    assert not meta["valid"] and np.all(b==-1)


def test_white_axis_geographic_effect_where_chromatic_hues_mix():
    df=make_plant("white_axis")
    row,null=a.species_analysis(df,"synthetic")
    assert row["W_all__eligible"]
    assert row["W_lead__eligible"]
    assert row["hue_lead_next__eligible"]
    assert row["W_all__identifiable"]
    assert row["hue_lead_next__identifiable"]
    assert row["W_all__excess_V2"]>0.05
    assert row["matched_z_diff_W_minus_hue"]>0
    assert all(x.shape==(a.N_PERM,) for x in null.values())


def test_white_specific_hypothesis_not_forced_when_nonwhite_sorting_stronger():
    row,_=a.species_analysis(make_plant("hue_axis"),"synthetic")
    assert row["W_all__eligible"]
    assert row["hue_lead_next__eligible"]
    assert row["hue_lead_next__excess_V2"]>0
    assert row["matched_z_diff_W_minus_hue"]<0


def test_nonexchangeable_month_only_is_not_claimed_as_benefit_cost():
    row,_=a.species_analysis(make_plant("month_only"),"synthetic")
    assert row["W_all__eligible"]
    assert not row["W_all__identifiable"]
    assert row["W_all__excess_V2"]==0
    assert row["hue_lead_next__eligible"]
    assert not row["hue_lead_next__identifiable"]
    assert row["hue_lead_next__excess_V2"]==0


def test_cramer_variance_normalizes_sample_frequency():
    target=np.r_[np.ones(20),np.zeros(20)]
    bins=np.r_[np.zeros(20),np.ones(20)]
    assert a.bin_association(target,bins)==1
    mixed=np.r_[np.ones(10),np.zeros(10),np.ones(10),np.zeros(10)]
    assert a.bin_association(mixed,bins)==0
    assert a.bin_association(np.zeros(40),bins)==0


def test_same_photo_morph_frequency_and_month_labels_remain_fixed():
    row,null=a.species_analysis(make_plant("white_axis"),"synthetic")
    assert row["n_white"]==40
    assert row["n_lead_hue"]==40
    assert row["n_next_hue"]==40
    a1,n1=a.species_analysis(make_plant("white_axis"),"synthetic")
    assert a1==row
    assert all(np.array_equal(n1[k],v) for k,v in null.items())


def test_result_summary_requires_matched_opportunity_and_credible_ci():
    rows=[]
    for p in ["white_axis","hue_axis","month_only"]:
        row,_=a.species_analysis(make_plant(p),"synthetic")
        row["inat_taxon_id"]+=len(rows)
        rows.append(row)
    df=pd.DataFrame(rows)
    ax=a.summarise_axis(df,"synthetic","W_all")
    assert ax["n_geographically_evaluable"]==3
    assert ax["n_conditional_identifiable"]==2
    paired=a.summary_match(df,"synthetic")
    assert paired["n_matched_species"]==2
    assert paired["status"]=="HOLD_SPARSE_PAIRED_AXIS"
    assert 0<=paired["paired_signflip_two_sided_p"]<=1
    assert len(paired["species_bootstrap_95CI"])==2

"""Synthetic falsification tests for the geographic species-equal FCP atlas."""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SRC=Path(__file__).resolve().parents[1]/"scripts/analysis"
sys.path.insert(0,str(SRC))
spec=importlib.util.spec_from_file_location(
    "fcp_latitude_elevation_atlas", SRC/"run_polymorphism_global_geography_bands_20261008.py")
a=importlib.util.module_from_spec(spec)
spec.loader.exec_module(a)


def make_species(taxon, n=80, *, pattern="cline"):
    half=n//2
    lat=np.r_[np.full(half,10.0),np.full(n-half,40.0)]
    elev=np.r_[np.full(half,100.0),np.full(n-half,1800.0)]
    if pattern=="cline":
        labels=["white"]*half+["red_pink"]*(n-half)
    elif pattern=="mixed":
        labels=(["white","red_pink"]*(half//2))+(
            ["white","red_pink"]*((n-half)//2))
    elif pattern=="allwhite":
        labels=["white"]*n
    elif pattern=="allcolour":
        labels=["red_pink"]*n
    else:
        raise ValueError(pattern)
    assert len(labels)==n
    return pd.DataFrame({
        "inat_taxon_id": [taxon]*n,
        "species": [f"Genus{taxon} species"]*n,
        "latitude": lat, "longitude":np.zeros(n),
        "elevation_m":elev, "morph":labels,
    })


def test_fourstate_probability_is_pairwise_unbiased():
    assert a.pair_diversity(np.array([8,0,0,0]))==0
    assert a.pair_diversity(np.array([4,4,0,0]))==32/56
    assert np.isnan(a.pair_diversity(np.array([1,0,0,0])))


def test_exact_abs_lat_and_elevation_band_semantics():
    b=a.region_assignments(np.array([0,14.99,15,29.99,30,45,60,90]),"absolute_latitude")
    assert b.tolist()==[0,0,1,1,2,3,4,4]
    z=a.region_assignments(np.array([-20,249.9,250,999,1000,1999,2000,3200,np.nan]),"elevation")
    assert z.tolist()==[0,0,1,1,2,2,3,3,-1]


def test_strong_signed_latitude_and_elevation_cline_within_species():
    g=make_species(100)
    for axis in ("absolute_latitude","elevation"):
        x=np.abs(g.latitude.to_numpy()) if axis=="absolute_latitude" else g.elevation_m.to_numpy()
        out,nul,mi=a.species_band_metrics(
            x,g.morph.to_numpy(),axis,nperm=199,seed=171)
        assert len(out)==2
        low,high=out
        assert low["n_photos_in_band"]==40
        assert high["n_photos_in_band"]==40
        assert low["delta_white_vs_same_species"]>0.49
        assert high["delta_white_vs_same_species"]< -0.49
        assert low["delta_diversity_vs_same_species"]< -0.48
        assert high["delta_diversity_vs_same_species"]< -0.48
        assert nul.shape==(len(a.NAMES[axis]),199,2)
        assert abs(nul[low["band_index"],:,0].mean())<0.10
        assert mi[low["band_index"]].mean()>0.5


def test_no_within_species_gradients_from_only_species_turnover():
    # Two taxa each sampled in one band, allwhite vs allcolour, are NOT
    # a within-species latitudinal transition, despite a huge pooled contrast.
    x=make_species(1,pattern="allwhite")
    x.loc[x.index,"latitude"]=10
    y=make_species(2,pattern="allcolour")
    y.loc[y.index,"latitude"]=40
    d=pd.concat([x,y],ignore_index=True)
    c,detail=a.analyse(d,"synthetic","absolute_latitude")
    assert detail.empty
    assert c["bands"]["0–15°"]["n_species_with_min8_photos_in_bin_even_without_outside_opportunity"]==1
    assert c["bands"]["30–45°"]["n_species_with_min8_photos_in_bin_even_without_outside_opportunity"]==1
    assert c["bands"]["0–15°"]["n_informative_species"]==0
    assert c["bands"]["30–45°"]["status"]=="NO_WITHIN_SPECIES_BAND_OPPORTUNITY"
    # A pure species-turnover map *does* show white vs coloured geography.
    # That descriptive contrast cannot be misreported as within-species FCP.
    assert c["bands"]["0–15°"]["descriptive_species_equal_white_fraction_all_band_species"]==1
    assert c["bands"]["30–45°"]["descriptive_species_equal_white_fraction_all_band_species"]==0
    assert c["bands"]["0–15°"]["descriptive_between_band_species_turnover_confounding"] is True


def test_noninformative_species_excluded_from_effect_but_retained_as_coverage():
    d=pd.concat([make_species(1,80,pattern="cline"),
                 make_species(2,80,pattern="allwhite")],ignore_index=True)
    c,detail=a.analyse(d,"synthetic","absolute_latitude")
    assert len(detail)==4
    low=c["bands"]["0–15°"]
    assert low["n_informative_species"]==2
    assert np.isclose(low["mean_white_deviation_vs_specieswide"],0.25,atol=1e-9)
    assert low["status"]=="HOLD_SPARSE_WITHIN_SPECIES_BAND_OPPORTUNITY"


def test_photo_sample_scale_does_not_weight_species_by_number_of_photos():
    d=pd.concat([make_species(1,80,pattern="cline"),
                 make_species(2,100,pattern="allwhite")],ignore_index=True)
    c,_=a.analyse(d,"synthetic","absolute_latitude")
    assert np.isclose(c["bands"]["0–15°"]["mean_white_deviation_vs_specieswide"],0.25)


def test_species_assignment_keeps_photo_counts_and_seed_determinism():
    g=make_species(10,80,pattern="cline")
    o1,n1,m1=a.species_band_metrics(np.abs(g.latitude.to_numpy()),g.morph.to_numpy(),
                                    "absolute_latitude",nperm=25,seed=42)
    o2,n2,m2=a.species_band_metrics(np.abs(g.latitude.to_numpy()),g.morph.to_numpy(),
                                    "absolute_latitude",nperm=25,seed=42)
    assert o1==o2
    assert np.array_equal(n1,n2)
    assert np.array_equal(m1,m2)

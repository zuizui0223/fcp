"""Synthetic negative and positive source-photo mixed phenotype candidate gate."""
from pathlib import Path
import importlib.util
import sys
import numpy as np
import pandas as pd
import pytest

p=Path(__file__).resolve().parents[1]/"scripts/analysis/audit_fcp_white_colour_contemporaneous_photo_candidates_20261010.py"
spec=importlib.util.spec_from_file_location("mixed_fcp_site",p)
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

def make_species(*,two_colours=True,same_year=True,observers_unique=True,far=False):
    n=50
    g=pd.DataFrame({
        "inat_taxon_id":[18]*n,"species":["Campanula example"]*n,
        "photo_id":np.arange(n).astype(str),
        "latitude":[35.]*n,
        "longitude":[140. if not far else 140.+i*.2 for i in range(n)],
        "morph":["white"]*25+["blue_purple" if two_colours else "white"]*25,
        "observer":[str(i) if observers_unique else "person1" for i in range(n)],
        "year":[2025 if same_year or i<25 else 2024 for i in range(n)],
        "month":[7]*n,"has_valid_date":[True]*n,
        "is_white":[True]*25+[False if two_colours else True]*25
    })
    return g

def test_exact_10km_same_year_two_photo_and_four_photo_witness():
    g=make_species()
    d=mod.geo_km(g.latitude,g.longitude)
    assert mod.find_mixed_pair(g,d,10.,"same_month_year")["witness"]
    x=mod.four_distinct_observer_anchor(g,d,10.,"same_month_year")
    assert x["witness"] is True
    assert x["maximum_distance_km"]<.001

def test_one_photographer_does_not_become_four_independent_source_observers():
    g=make_species(observers_unique=False)
    d=mod.geo_km(g.latitude,g.longitude)
    assert not mod.find_mixed_pair(g,d,10.,"same_month_year")["witness"]
    assert not mod.four_distinct_observer_anchor(g,d,10.,"same_month_year")["witness"]

def test_single_colour_photo_species_is_not_photo_polymorphism():
    g=make_species(two_colours=False)
    d=mod.geo_km(g.latitude,g.longitude)
    assert not mod.find_mixed_pair(g,d,10.,"same_month_year")["witness"]
    assert not mod.four_distinct_observer_anchor(g,d,10.,"same_month_year")["witness"]

def test_calendar_month_across_years_is_separate_from_same_month_and_year():
    g=make_species(same_year=False)
    d=mod.geo_km(g.latitude,g.longitude)
    assert not mod.find_mixed_pair(g,d,10.,"same_month_year")["witness"]
    assert mod.find_mixed_pair(g,d,10.,"same_calendar_month")["witness"]
    assert not mod.four_distinct_observer_anchor(g,d,10.,"same_month_year")["witness"]
    assert mod.four_distinct_observer_anchor(g,d,10.,"same_calendar_month")["witness"]

def test_far_photographs_fail_site_diameter_even_when_colours_and_dates_match():
    g=make_species(far=True)
    d=mod.geo_km(g.latitude,g.longitude)
    assert not mod.find_mixed_pair(g,d,10.,"same_month_year")["witness"]
    assert not mod.four_distinct_observer_anchor(g,d,10.,"same_month_year")["witness"]

def test_four_distinct_photographers_even_if_shared_names_across_colours():
    g=make_species()
    g.loc[:24,"observer"]=[f"s{i}" for i in range(25)]
    g.loc[25:,"observer"]=[f"s{i}" for i in range(25)]
    d=mod.geo_km(g.latitude,g.longitude)
    # Four on each side is sufficient to construct two-and-two disjoint observers.
    assert mod.four_distinct_observer_anchor(g,d,10.,"same_month_year")["witness"]

def test_original_source_hash_gate_rejects_tampering(tmp_path):
    p=tmp_path/"mutated.csv";p.write_text("ina_photo,colour\n1,white\n")
    with pytest.raises(RuntimeError,match="SHA mismatch"):
        mod.load(p,"discovery")

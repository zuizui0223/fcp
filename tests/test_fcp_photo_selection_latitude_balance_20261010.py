"""Synthetic original-photo missingness and conditional pigment-balance null guards."""
from pathlib import Path
import importlib.util
import sys
import numpy as np
import pandas as pd
import pytest
BASE=Path(__file__).resolve().parents[1]/"scripts/analysis"
sys.path.insert(0,str(BASE))
from audit_fcp_global_photo_rainfall_classifiability_IPW_MNAR_20261010 import (
    load_original,features,fractions,tipping,SOURCE_SHA
)
from audit_fcp_global_photo_latitude_balance_null_20261010 import uniformity_test

def source_frame():
    rows=[]
    for band,lat,rain in [(0,10.,250.),(1,45.,950.),(2,75.,3000.)]:
        for genus in range(5):
            for i in range(42):
                rows.append({
                  "inat_taxon_id":genus,"species":f"G{genus} s","genus":f"G{genus}",
                  "latitude":lat,"longitude":20.,"morph":("white" if i%2 else "blue_purple"),
                  "classified":True,"wc_bio1":10., "wc_bio5":25.,
                  "wc_bio12":rain,"wc_bio15":40.,"wc_elevation_m":100.,
                  "source_panel":"cell","flower_effective_pixels":1000.,
                  "photo_id":len(rows),"rain_decile":0 if band==0 else 9})
    out=pd.DataFrame(rows)
    out["white"]=out.morph.eq("white").astype(int)
    return out

def test_source_checksum_fails_closed(tmp_path):
    p=tmp_path/"fake.zip";p.write_bytes(b"wrong archive")
    with pytest.raises(RuntimeError,match="SHA256"):load_original(p)

def test_explicit_four_hues_not_unverified_genotypes():
    a=source_frame();z=uniformity_test(a,permutations=19)
    assert z["n_original_classified_geotagged_photos"]==630
    assert z["latent_adaptive_selection_equilibrium_unidentified"]
    assert z["conditional_colour_photo_label_nulls"]["within_species"]["group_white_totals_preserved_exactly"]
    assert z["conditional_colour_photo_label_nulls"]["within_genus"]["group_white_totals_preserved_exactly"]

def test_colour_label_missing_in_source_remains_unknown():
    a=source_frame();a.loc[:12,"classified"]=False
    a.loc[:12,"morph"]="mixed_uncertain"
    z=uniformity_test(a,permutations=9)
    assert z["alternative_photo_colour_unknown_count"]==13
    assert z["n_original_classified_geotagged_photos"]==617

def test_propensity_predictor_uses_only_source_available_fields():
    d=source_frame()
    assert features(d,"geo_climate").shape==(len(d),9)
    assert features(d,"geo_climate_plus_flower_area").shape==(len(d),11)
    with pytest.raises(ValueError):features(d,"hypothesis_chasing")

def test_ipw_constant_propensity_preserves_source_photo_white_fraction():
    x=source_frame()
    a,da=fractions(x,np.ones(len(x)))
    b,db=fractions(x,np.ones(len(x))*3.)
    assert np.isclose(da,db)
    assert np.isclose(a["0"]["white_fraction"],b["0"]["white_fraction"])

def test_MNAR_tipping_curve_solves_exact_missing_white_equality():
    d=source_frame()
    d.loc[0:29,"classified"]=False
    d.loc[210:269,"classified"]=False
    c=tipping(d)
    for dry,k in c.items():
        if not k["feasible_value_between_zero_and_one"]:continue
        wet=k["required_unclassified_wet_white_fraction_for_EQUAL_true_dry_wet_white"]
        lo=d[d.rain_decile==0];hi=d[d.rain_decile==9]
        pl=(lo.loc[lo.classified,"morph"].eq("white").sum()+
            len(lo.loc[~lo.classified])*float(dry))/len(lo)
        ph=(hi.loc[hi.classified,"morph"].eq("white").sum()+
            len(hi.loc[~hi.classified])*wet)/len(hi)
        assert np.isclose(pl,ph)

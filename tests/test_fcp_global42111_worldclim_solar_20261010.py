"""Synthetic source-only WorldClim monthly srad, BIO and missingness guards."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"/"analysis"))
import extend_fcp_global42111_worldclim_solar_20261010 as M


@pytest.fixture
def photos(monkeypatch):
    monkeypatch.setattr(M,"ORIGINAL_N",4)
    monkeypatch.setattr(M,"CLASSIFIED_N",3)
    d=pd.DataFrame({
        "inat_taxon_id":[1,2,3,4],"observation_id":[11,12,13,14],"photo_id":[111,112,113,114],
        "species":["G a","G b","H a","K a"],
        "morph":["white","red_pink","blue_purple","UNKNOWN"],
        "measurement_status":["classified_four_state_morph"]*3+["not_classified"],
        "latitude":[12.,13.,14.,np.nan],"longitude":[10.,11.,12.,np.nan],
        "site_geo_status":[M.SOURCE_GEO]*3+["NO_UNOBSCURED_PUBLIC_COORDINATES"],
    })
    for c in ["wc_bio1","wc_bio5","wc_bio12","wc_bio15","wc_elevation_m",
              "soil_pH","soil_SOC","soil_N","soil_clay","soil_available_water_proxy"]:
        d[c]=[1.,2.,3.,np.nan]
    return d


def sampler(path,lon,lat):
    v=np.where(np.isfinite(lon),4.0,np.nan)
    if "srad_" in str(path):
        month=int(Path(path).stem.split("_")[-1])
        return np.where(np.isfinite(lon),1000+100*month,np.nan)
    return v


def test_annual_12month_source_srad_cv_and_original_photo_labels(photos):
    original=photos.copy()
    out,receipt=M.extract(photos,Path("bio"),Path("srad"),sampler=sampler)
    assert len(out)==4
    assert receipt["original_photo_unclassified"]==1
    assert receipt["n_solar_complete_source_taxa"]==3
    assert receipt["n_additional_sun_bio_complete_classified"]==3
    assert out.loc[0,"wc_srad_annual_kj_m2_day"]==pytest.approx(1650)
    assert out.loc[0,"wc_srad_monthly_cv"]>0
    assert out.loc[3,list(M.NEW_FEATURES)].isna().all()
    pd.testing.assert_series_equal(out.morph,original.morph)


def test_any_month_missing_sun_means_missing_annual_exposure(photos):
    def missing(path,lon,lat):
        a=sampler(path,lon,lat)
        if "srad_6.tif" in str(path):a[0]=np.nan
        return a
    out,report=M.extract(photos,Path("bio"),Path("srad"),sampler=missing)
    assert pd.isna(out.loc[0,"wc_srad_annual_kj_m2_day"])
    assert report["n_solar_complete_source_taxa"]==2


def test_negative_sun_radiation_cannot_be_silently_used(photos):
    def negative(path,lon,lat):
        a=sampler(path,lon,lat)
        if "srad_1.tif" in str(path):a[0]=-10
        return a
    with pytest.raises(ValueError,match="Negative"):
        M.extract(photos,Path("bio"),Path("srad"),sampler=negative)


def test_ungeolocated_original_photo_does_not_get_fake_radiation(photos):
    def bad(path,lon,lat):
        return np.array([1.,2.,3.,4.])
    with pytest.raises(ValueError,match="artificial"):
        M.extract(photos,Path("bio"),Path("srad"),sampler=bad)


def test_original_fourstate_status_must_be_preserved(photos):
    photos.loc[3,"measurement_status"]="classified_four_state_morph"
    photos.loc[3,"morph"]="white"
    with pytest.raises(ValueError,match="classified"):
        M.extract(photos,Path("bio"),Path("srad"),sampler=sampler)

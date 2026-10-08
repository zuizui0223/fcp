"""Synthetic source-preserving FCP real-site climate/soil raster attachment."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"/"analysis"))
import attach_fcp_global42111_true_site_climate_soil_20261009 as M


@pytest.fixture
def small_fixed_population(monkeypatch):
    # Synthetic denominator substitution applies only to these unit tests.
    monkeypatch.setattr(M,"BREADTH_ORIGINAL",3)
    monkeypatch.setattr(M,"CELL_ORIGINAL",4)
    points=pd.DataFrame([
        (1,100,111,True,True,12.,10.,M.PHOTO_STATUS),
        (2,200,222,True,True,13.,11.,M.PHOTO_STATUS),
        (3,300,333,True,False,14.,12.,M.PHOTO_STATUS),
        (1,400,444,False,True,15.,13.,M.PHOTO_STATUS),
        (2,500,555,False,True,np.nan,np.nan,"NO_UNOBSCURED_PUBLIC_COORDINATES"),
    ],columns=["inat_taxon_id","observation_id","photo_id","present_in_breadth",
               "present_in_taxon_cell","latitude","longitude","site_geo_status"])
    b=pd.DataFrame({"inat_taxon_id":[1,2,3],"observation_id":[100,200,300],
                    "photo_id":[111,222,333],"morph":["white","red_pink","unclassified"],
                    "measurement_status":["classified_four_state_morph","classified_four_state_morph","not_classified"]})
    c=pd.DataFrame({"inat_taxon_id":[1,2,1,2],"observation_id":[100,200,400,500],
                    "photo_id":[111,222,444,555],"cell_id":[99,99,99,99],
                    "morph":["white","red_pink","yellow_orange","unclassified"],
                    "measurement_status":["classified_four_state_morph"]*3+["not_classified"]})
    return points,b,c


def known_fake_raster(path,lon,lat):
    valid=np.isfinite(lon)&np.isfinite(lat)
    value=1.0 if "wv1500" in str(path) else 3.0 if "wv0033" in str(path) else 20.0
    return np.where(valid,value,np.nan)


def test_real_photo_positions_only_and_soil_depth_weights(small_fixed_population):
    points,b,c=small_fixed_population
    env=M.attach(points,Path("/tmp/wc"),Path("/tmp/elev"),Path("/tmp/soil"),sampler=known_fake_raster)
    assert len(env)==5
    assert int(env.environment_all_complete.sum())==4
    assert pd.isna(env.loc[env.photo_id==555,"soil_pH"].iloc[0])
    assert env.loc[env.photo_id==111,"soil_available_water_proxy"].iloc[0]==pytest.approx(0.2)
    assert not any(s.startswith("cell_centroid") for s in env.columns)


def test_preserve_every_species_and_taxon_cell_even_if_environment_missing(small_fixed_population):
    points,b,c=small_fixed_population
    env=M.attach(points,Path("wc"),Path("elev"),Path("soil"),sampler=known_fake_raster)
    bb=M.left_join_measurement(b,env,"breadth")
    cc=M.left_join_measurement(c,env,"taxon_cell")
    report=M.aggregate(bb,cc,env)
    assert len(bb)==3 and len(cc)==4
    assert report["n_unique_original_photo_IDs"]==5
    assert report["breadth_classifiable_by_covariate"]["environment_all_complete"]==2
    assert report["taxon_cell_classifiable_by_covariate"]["environment_all_complete"]==3
    assert report["geographic_cell_centroid_imputation"] is False
    assert report["new_colour_label_or_photo_pixel_opened"] is False


def test_same_photo_id_with_different_taxon_identity_cannot_duplicate(small_fixed_population):
    points,b,c=small_fixed_population
    bad=pd.concat([points,points.iloc[[0]]],ignore_index=True)
    with pytest.raises(ValueError,match="duplicated"):
        M.check_point_identity(bad)


def test_invalid_geo_has_no_fake_filled_coordinate(small_fixed_population):
    points,b,c=small_fixed_population
    points.loc[points.photo_id==555,"latitude"]=20.
    with pytest.raises(RuntimeError,match="fabricated"):
        M.check_point_identity(points)


def test_source_colour_never_uses_soil_missingness_to_filter(small_fixed_population):
    points,b,c=small_fixed_population
    env=M.attach(points,Path("wc"),Path("elev"),Path("soil"),sampler=known_fake_raster)
    out=M.left_join_measurement(c,env,"taxon_cell")
    assert len(out)==4
    assert len(out.loc[out.morph=="unclassified"])==1
    assert int(out.environment_all_complete.sum())==3

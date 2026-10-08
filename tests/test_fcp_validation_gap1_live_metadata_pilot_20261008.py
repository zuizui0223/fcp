"""Synthetic safety tests: Validation metadata-only pilot, never call live API."""
import sys
from pathlib import Path
import pandas as pd
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"/"acquisition"))
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"/"analysis"))
from audit_fcp_validation_gap1_live_metadata_pilot_20261008 import (
    qualify, anchor_context, public_coordinates
)


def response(**overrides):
    item={
        "id":999901,"taxon":{"id":118778},"observed_on":"2021-08-12",
        "user":{"id":491},"coordinates_obscured":False,
        "geoprivacy":None,"geojson":{"type":"Point","coordinates":[-105.0,40.0]},
        "positional_accuracy":300,"captive":False,
        "photos":[{"id":80002,"license_code":"cc-by"}],
    }
    item.update(overrides)
    return item

def pick(item,**params):
    args=dict(taxon_id=118778,year=2021,month=8,anchor_lat=40,anchor_lon=-105,
              existing_observers={"490"},banned_obs={"999900"},banned_photos={"80001"})
    args.update(params)
    return qualify(item,**args)

def test_accepts_open_licenced_different_observer_only():
    o=pick(response())
    assert len(o)==1 and o[0]["photo_id"]=="80002"
    assert o[0]["public_distance_km"] < 0.1

def test_rejects_obscured_or_private_geography():
    assert public_coordinates(response(geoprivacy="private")) is None
    assert pick(response(coordinates_obscured=True))==[]

def test_rejects_same_observer():
    assert pick(response(user={"id":490}))==[]

def test_rejects_used_photo_and_observation_ids():
    assert pick(response(photos=[{"id":80001,"license_code":"cc-by"}]))==[]
    assert pick(response(id=999900))==[]

def test_rejects_wrong_taxon_or_wrong_month():
    assert pick(response(taxon={"id":999}))==[]
    assert pick(response(observed_on="2021-09-01"))==[]

def test_rejects_nonopen_license_and_inaccurate_geography():
    assert pick(response(photos=[{"id":80002,"license_code":None}]))==[]
    assert pick(response(positional_accuracy=6000))==[]
    assert pick(response(geojson={"coordinates":[-106,40]}))==[]

def test_rejects_captive_and_handles_null_coordinates():
    assert pick(response(captive=True))==[]
    assert pick(response(geojson=None,location=None))==[]

def test_only_existing_year_month_with_one_observer_passes():
    src=pd.DataFrame([
        {"inat_taxon_id":118778,"photo_id":90011478,"latitude":40.0,"longitude":-105.0,
         "month":8,"year":2021,"observer":"490","species":"Gentiana parryi"},
        {"inat_taxon_id":118778,"photo_id":90011479,"latitude":40.0001,"longitude":-105.0001,
         "month":8,"year":2024,"observer":"999","species":"Gentiana parryi"},
    ])
    row=pd.Series({"inat_taxon_id":118778,"source_anchor_photo_id":90011478,
                   "calendar_month":8,"target_year":2021,
                   "existing_distinct_observers_in_target_year":1,
                   "required_additional_distinct_observer_photos":1})
    lat,lon,obs=anchor_context(src,row)
    assert lat==40 and lon==-105 and obs=={"490"}
    row["existing_distinct_observers_in_target_year"]=2
    with pytest.raises(RuntimeError,match="counts drifted"):
        anchor_context(src,row)

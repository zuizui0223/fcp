"""Mock only the original public observation metadata fetch; no network or pixels."""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"/"acquisition"))
from recover_fcp_global42111_all_original_photo_public_geocoordinates_20261009 import (
    make_manifest, unpack, run_shard, SCHEMA, SHARDS
)


def source_tables():
    species=np.arange(1,42112)
    b=pd.DataFrame({"inat_taxon_id":species,"observation_id":100000+species,
                    "photo_id":200000+species})
    # 20,000 source photo IDs are reused unchanged by 85,337 cell records.
    reused=b.iloc[:20000].copy()
    new_n=85337-len(reused)
    ids=np.arange(new_n)
    new=pd.DataFrame({"inat_taxon_id":1+ids%42111,
                      "observation_id":500000+ids,
                      "photo_id":600000+ids})
    c=pd.concat([reused,new],ignore_index=True)
    c["cell_id"]=0
    assert len(c)==85337
    return b,c


def test_exact_all_original_photo_manifest_and_reuse():
    b,c=source_tables()
    m,receipt=make_manifest(b,c)
    assert receipt["original_breadth_source_rows"]==42111
    assert receipt["original_taxon_cell_source_rows"]==85337
    assert receipt["original_taxon_species"]==42111
    assert len(m)==42111+85337-20000
    assert m.observation_id.nunique()==len(m)
    assert set(m.shard_index).issubset(set(range(SHARDS)))
    assert int(m.present_in_breadth.sum())==42111
    assert int(m.present_in_taxon_cell.sum())==85337


def test_photo_id_not_reused_across_different_observation():
    b,c=source_tables()
    c.loc[20000,"photo_id"]=b.photo_id.iloc[0]
    with pytest.raises(RuntimeError,match="photo ID"):
        make_manifest(b,c)


def observation(lat=12.,lon=10.,taxon=3,photo=1234):
    return {"id":340, "taxon":{"id":taxon},"photos":[{"id":photo}],
            "geoprivacy":None,"coordinates_obscured":False,
            "geojson":{"coordinates":[lon,lat]}, "positional_accuracy":50,"captive":False}


def record(cell=99):
    return pd.DataFrame([{"inat_taxon_id":3,"observation_id":340,"photo_id":1234,
                          "original_cell_id":cell,"present_in_breadth":True,
                          "present_in_taxon_cell":True}])


def test_original_photo_observation_and_cell_are_consistent():
    out=unpack(record(),observation())[0]
    assert out["site_geo_status"]=="VALID_PUBLIC_ORIGINAL_PHOTO_POINT"
    assert (out["latitude"],out["longitude"])==(12.,10.)


def test_nonmatching_current_location_never_fabricated():
    assert unpack(record(cell=100),observation())[0]["site_geo_status"]=="CURRENT_PUBLIC_POSITION_OUTSIDE_HISTORICAL_EQUAL_AREA_CELL"
    assert unpack(record(),observation(taxon=55))[0]["site_geo_status"]=="SOURCE_TAXON_IDENTITY_MISMATCH"
    assert unpack(record(),observation(photo=1))[0]["site_geo_status"]=="EXACT_ORIGINAL_PHOTO_MISSING"


def test_private_geography_and_cultivated_not_native_points():
    q=observation()
    q["geoprivacy"]="obscured"
    assert unpack(record(),q)[0]["site_geo_status"]=="NO_UNOBSCURED_PUBLIC_COORDINATES"
    assert unpack(record(),{**observation(),"captive":True})[0]["site_geo_status"]=="CULTIVATED_NOT_WILD"


def test_missing_api_observation_is_unknown_not_biological_absence():
    assert unpack(record(),None)[0]["site_geo_status"]=="SOURCE_OBSERVATION_NOT_RETURNED"


def test_mocked_shard_writes_exact_photo_and_no_colour_outcome(tmp_path):
    m=record().assign(shard_index=0)
    a=run_shard(m,0,tmp_path,client=lambda ids:{"results":[observation()]},sleep=lambda _:None)
    assert a["schema"]==SCHEMA
    assert a["n_source_original_photos"]==1
    assert a["n_exact_photo_public_geocoordinates"]==1
    assert (tmp_path/"original_geo_records_shard_00.csv.gz").exists()
    assert a["images_unopened"] and a["colour_results_unopened"]


def test_failed_api_batch_stays_unresolved(tmp_path):
    def bad(ids):raise TimeoutError("simulated API outage")
    m=record().assign(shard_index=0)
    a=run_shard(m,0,tmp_path,client=bad,sleep=lambda _:None)
    assert a["n_batch_errors"]==1
    assert a["source_state_counts"]=={"API_ERROR_UNRESOLVED":1}

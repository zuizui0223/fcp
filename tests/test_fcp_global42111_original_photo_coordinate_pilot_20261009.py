"""Mocked, rate-bounded global original ID geocoordinate pilot."""
from __future__ import annotations
import sys
from pathlib import Path
import pandas as pd
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"/"acquisition"))
from audit_fcp_global42111_original_photo_coordinate_pilot_20261009 import (
    cell_id_for,public_position,summarize_one,pilot,fetch_batch
)


def photo(sample_cell=99):
    return {"inat_taxon_id":1,"observation_id":123,"photo_id":456,
            "cell_id":sample_cell,"origin":"taxon_cell"}


def observed(**kw):
    o={"id":123,"taxon":{"id":1},"photos":[{"id":456}],
       "geojson":{"type":"Point","coordinates":[10.0,12.0]},
       "positional_accuracy":200,"coordinates_obscured":False,
       "geoprivacy":None,"captive":False}
    o.update(kw)
    return o


def test_exact_true_point_photo_match():
    assert cell_id_for(12,10)==99
    out=summarize_one(photo(),observed())
    assert out["status"]=="EXACT_SOURCE_PHOTO_PUBLIC_COORDINATE_AVAILABLE"
    assert (out["latitude"],out["longitude"])==(12.0,10.0)


def test_private_obscured_is_not_geolocated():
    assert public_position(observed(geoprivacy="obscured")) is None
    assert summarize_one(photo(),observed(coordinates_obscured=True))["status"]=="NO_UNOBSCURED_PUBLIC_COORDINATES"


def test_wrong_cell_and_current_taxon_not_borrowed():
    assert summarize_one(photo(sample_cell=100),observed())["status"]=="LOCATION_CHANGED_OR_SOURCE_CELL_MISMATCH"
    assert summarize_one(photo(),observed(taxon={"id":2}))["status"]=="TAXON_IDENTITY_CHANGED"


def test_original_photo_id_required():
    assert summarize_one(photo(),observed(photos=[{"id":457}]))["status"]=="ORIGINAL_PHOTO_ID_MISSING"


def test_one_batch_mock_and_unresolved_id():
    rows=pd.DataFrame([photo(),{**photo(),"observation_id":124,"photo_id":457}])
    report,out=pilot(rows,client=lambda ids:{"results":[observed()]},pause=lambda t:None)
    assert report["n_photo_metadata_targets"]==2
    assert report["n_exact_source_photo_coordinates"]==1
    assert report["n_requests"]==1
    assert out.status.tolist()==["EXACT_SOURCE_PHOTO_PUBLIC_COORDINATE_AVAILABLE","OBSERVATION_ID_NOT_RETURNED"]


def test_api_error_unknown_not_no_coordinates():
    rows=pd.DataFrame([photo()])
    def fail(ids):raise TimeoutError("synthetic")
    report,out=pilot(rows,client=fail,pause=lambda t:None)
    assert report["n_failed_request_batches"]==1
    assert report["n_exact_source_photo_coordinates"]==0
    assert out.status.iloc[0]=="API_ERROR_UNRESOLVED"


def test_rate_guard_rejects_duplicate_or_overlarge():
    with pytest.raises(ValueError,match="Duplicate"):
        fetch_batch([12,12])
    with pytest.raises(ValueError,match="Duplicate"):
        fetch_batch(list(range(101)))

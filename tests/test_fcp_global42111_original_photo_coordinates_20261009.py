"""Synthetic original-photo-ID geolocation recovery and ocean-centre safety tests."""
from __future__ import annotations
import sys
from pathlib import Path

import pandas as pd
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"/"analysis"))
from recover_fcp_global42111_original_photo_coordinates_20261009 import (
    source_index, join_immutable, grid_cell_for, actual_coordinate_pair
)


def frame(index_coord=True):
    x=pd.DataFrame({
        "inat_taxon_id":[1,2],
        "observation_id":[101,202],
        "photo_id":[1011,2022],
        "cell_id":[99,100],
        "species":["Species one","Species two"],
        "morph":["white","blue_purple"],
        "measurement_status":["classified_four_state_morph"]*2
    })
    # row 5 -> sin(lat) band 4, lon index 9 => 18*4 + 9 = 81 (not 99)
    # row index 5 -> sin(lat) band 5, lon index 9 => 99
    # band 5 is 0.111.. to 0.333.. sin(lat); latitude 12 deg.
    if index_coord:
        x["latitude"]=[12.0,12.0]
        x["longitude"]=[10.0,30.0]
    return x


def test_observed_geometry_not_false_centroid():
    assert grid_cell_for(np.array([12.,12.]),np.array([10.,30.])).tolist()==[99,100]


def test_exact_matching_source_photo_ids_and_positions():
    measured=frame(index_coord=False)
    index=frame(index_coord=True)
    idx,receipt=source_index(index,index.iloc[0:0].copy())
    matched,stats=join_immutable(measured,idx,"taxon_cell")
    assert stats["n_matched_index_id_triples"]==2
    assert stats["n_recorded_coordinate_valid_and_original_cell_consistent"]==2
    assert matched.latitude.tolist()==[12.0,12.0]


def test_unrelated_obs_photo_not_borrowed():
    measured=frame(index_coord=False)
    measured.loc[1,"photo_id"]=333333
    idx,_=source_index(frame(index_coord=True),frame(index_coord=True).iloc[0:0])
    matched,stats=join_immutable(measured,idx,"taxon_cell")
    assert stats["n_recorded_coordinate_valid_and_original_cell_consistent"]==1
    assert pd.isna(matched.latitude.iloc[1])


def test_photo_metadata_without_true_coordinates_means_missing():
    measured=frame(index_coord=False)
    idx,receipt=source_index(measured,measured.iloc[0:0].copy())
    matched,stats=join_immutable(measured,idx,"taxon_cell")
    assert stats["n_recorded_coordinate_valid_and_original_cell_consistent"]==0
    assert receipt["v1"]["coordinate_source_fields"] is None
    assert matched.latitude.isna().all()


def test_wrong_geocell_is_rejected_not_a_fake_soil_coordinate():
    measured=frame(index_coord=False)
    index=frame(index_coord=True)
    index.loc[0,"latitude"]=-50.0
    idx,_=source_index(index,index.iloc[0:0].copy())
    matched,stats=join_immutable(measured,idx,"taxon_cell")
    assert stats["n_recorded_coordinate_valid_and_original_cell_consistent"]==1
    assert stats["n_source_cell_coordinate_mismatch"]==1
    assert pd.isna(matched.latitude.iloc[0])


def test_conflicted_duplicate_origins_are_discarded():
    measured=frame(index_coord=False)
    index=frame(index_coord=True)
    alternate=index.copy()
    alternate.loc[0,"cell_id"]=100
    idx,receipt=source_index(index,alternate)
    matched,stats=join_immutable(measured,idx,"taxon_cell")
    assert receipt["combined"]["n_ambiguous_rows_unusable"]>=2
    assert pd.isna(matched.latitude.iloc[0])


def test_missing_geoidentity_never_invented():
    assert actual_coordinate_pair(pd.DataFrame({"cell_id":[99]})) is None

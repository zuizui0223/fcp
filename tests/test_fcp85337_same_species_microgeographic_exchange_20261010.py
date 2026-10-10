"""Original photographed within-species nearby flower colour exchangeability guards."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"/"analysis"))
import audit_fcp85337_same_species_microgeographic_exchange_20261010 as M


@pytest.fixture
def originals():
    # Species 1: 3 different original source cells photographed <12km apart;
    # Species 2: 2 cells ~167km apart, only same photo-colour;
    # Species 3: two original photos far apart.
    return pd.DataFrame({
        "inat_taxon_id":[1,1,1,2,2,3,3],
        "cell_id":[1,2,3,1,2,1,2],
        "photo_id":[101,102,103,201,202,301,302],
        "latitude":[0.,0.,.1,0.,0.,0.,0.],
        "longitude":[0.,.05,.10,10.,11.5,20.,28.],
        "morph":["white","red_pink","white","white","white","blue_purple","white"],
        "site_geo_status":["VALID_PUBLIC_ORIGINAL_PHOTO_POINT"]*7,
        "eligible_heldout_original_species":[True]*7,
    })


def test_haversine_photo_location_not_geographic_cell_center(originals):
    d=M.geometry(np.array([0,0,0]),np.array([0,.05,.10]))
    assert d[0,1]==pytest.approx(5.55975,abs=.02)
    assert d[0,2]==pytest.approx(11.1195,abs=.03)
    with pytest.raises(ValueError,match="missing or invalid"):
        M.geometry(np.array([0,float("nan")]),np.array([0,1]))


def test_complete_link_species_and_region_sample_no_photo_reuse(originals):
    micro=M.completed_source_microgroups(originals,50)
    assert sorted(map(len,micro))==[1,1,2,3] or sorted(map(len,micro))==[1,1,3,2]
    allidx=np.concatenate(micro)
    assert len(allidx)==len(originals) and len(set(allidx))==len(originals)
    assert all(originals.iloc[x].inat_taxon_id.nunique()==1 for x in micro)


def test_four_thresholds_and_photo_colour_blind_group_selection(originals):
    a=M.audit(originals)
    assert set(a)=={"50","100","250","500"}
    assert a["50"]["n_original_photos_in_two_plus_geographically_nearby_same_species_groups"]==3
    assert a["50"]["n_original_pairs_of_source_photos_within_complete_link_neighbourhoods"]==3
    assert a["50"]["n_source_photo_records_in_colour_exchangeable_nearby_groups"]==3
    assert a["250"]["n_original_photos_in_two_plus_geographically_nearby_same_species_groups"]==5
    assert a["250"]["n_source_photo_records_in_colour_exchangeable_nearby_groups"]==3
    assert a["500"]["readiness"]=="HOLD_INSUFFICIENT_NEARBY_INTRASPECIFIC_COLOUR_EXCHANGE"
    other=originals.copy()
    other["morph"]="white"
    b=M.audit(other)
    assert [a[k]["n_same_species_photo_groups_including_singletons"] for k in a]==[
        b[k]["n_same_species_photo_groups_including_singletons"] for k in b]
    assert all(x["group_membership_photo_colour_blind"] for x in a.values())
    assert b["50"]["n_source_photo_records_in_colour_exchangeable_nearby_groups"]==0


def test_source_group_ignores_outcome_but_fails_if_same_species_cell_duplicated(originals):
    bad=originals.copy()
    bad.loc[1,"cell_id"]=1
    with pytest.raises(ValueError,match="one original photograph"):
        M.completed_source_microgroups(bad,100)


def test_holdout_same_species_training_only_respects_original_162_cell_groups(originals):
    data=originals.copy()
    # Five distinct original source cells necessary for GroupKFold; the actual
    # 85337 source has 128 occupied cells.
    data.loc[4,"cell_id"]=4
    data.loc[6,"cell_id"]=5
    mask=M.source_heldout_mask(data)
    assert len(mask)==len(data)
    assert mask.dtype==bool
    assert bool(mask.sum())
    assert data.loc[mask,"inat_taxon_id"].nunique()>0


def test_original_source_ungeolocated_fails_without_imputation(originals):
    data=originals.copy()
    data.loc[0,"site_geo_status"]="NO_PUBLIC_COORDINATE"
    with pytest.raises(ValueError,match="Ungeolocated"):
        M.audit(data)

"""Synthetic source-photo constrained spatial permutation audit, no external data."""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"/"analysis"))
import audit_fcp_global42111_microspatial_photo_shuffle_20261009 as M


@pytest.fixture
def local_photos():
    rng=np.random.default_rng(12026)
    rows=[]
    for genus in range(22):
        for j in range(15):
            # All source photos within 5km of this original genus-cell anchor.
            lat=10+rng.uniform(-.012,.012)
            lon=-130+10*(genus%9)+rng.uniform(-.012,.012)
            rows.append({
                "inat_taxon_id":genus*15+j+1,
                "genus":f"Genus{genus}",
                "genus_cell_id":f"Genus{genus}|{genus%9}",
                "photo_cell_162":genus%9,
                "latitude":lat,"longitude":lon,
                "morph":("white","red_pink","yellow_orange","blue_purple")[j%4],
                "abs_latitude":abs(lat),
                "lon_sin":np.sin(np.deg2rad(lon)),
                "lon_cos":np.cos(np.deg2rad(lon)),
                "wc_elevation_m":rng.uniform(5,100),
                "wc_bio1":rng.uniform(10,22),
                "wc_bio5":rng.uniform(18,38),
                "wc_bio12":rng.uniform(200,2000),
                "wc_bio15":rng.uniform(10,90)
            })
    return pd.DataFrame(rows)


def test_geodesic_distance_and_complete_link_chain_not_single_link():
    d=M.great_circle_matrix_km([0,0,0],[0,.3,.6])
    assert d[0,1]==pytest.approx(33.36,abs=.1)
    assert d[0,2]>60
    photos=pd.DataFrame({
        "genus_cell_id":["A|5"]*3,"inat_taxon_id":[1,2,3],
        "latitude":[0,0,0],"longitude":[0,.3,.6]})
    clusters=M.microgroups(photos,40.0)
    assert sorted(map(len,clusters))==[1,2]


def test_source_grouping_is_colour_and_row_order_independent(local_photos):
    a=M.microgroups(local_photos,50.)
    randomized=local_photos.copy()
    randomized["morph"]="white"
    b=M.microgroups(randomized,50.)
    assert [list(z) for z in a]==[list(z) for z in b]


def test_same_original_genus_cell_colours_preserved_under_microshuffles(local_photos):
    clusters=M.microgroups(local_photos,50.)
    labels=local_photos.morph.to_numpy(str)
    rng=np.random.default_rng(77)
    shuffled=labels.copy()
    for cluster in clusters:
        if len(cluster)>1:
            shuffled[cluster]=rng.permutation(labels[cluster])
        assert sorted(shuffled[cluster])==sorted(labels[cluster])
    for _,original in local_photos.groupby("genus_cell_id"):
        ids=original.index.to_numpy(int)
        assert sorted(shuffled[ids])==sorted(labels[ids])


def test_unsupported_tight_microneighborhood_has_explicit_hold(local_photos):
    # No two original photos are within essentially zero km apart.
    res=M.assess(local_photos,1e-5,nperm=3,seed=7)
    assert res["status"]=="HOLD_INSUFFICIENT_MICROSPATIAL_EXCHANGEABILITY"
    assert res["p_value"] is None
    assert res["n_refitted_null_models"]==0
    assert res["n_source_species_in_informatively_shufflable_microgroups"]==0


def test_estimable_null_refits_same_method_without_new_flower_labels(local_photos):
    result=M.assess(local_photos,50,nperm=2,seed=2026)
    assert result["status"]=="READY_FOR_CONDITIONAL_SPATIAL_NULL"
    assert result["n_refitted_null_models"]==2
    assert result["n_source_species_in_informatively_shufflable_microgroups"]==330
    assert result["group_membership_uses_only_original_species_ID_genus_and_public_photo_coordinates"]
    assert result["original_genus_cell_colour_counts_always_preserved"] and result["microgroup_colour_counts_preserved_on_each_null"]
    assert 1/3<=result["p_value"]<=1


def test_no_imputed_original_geocoordinates_accepted(local_photos):
    bad=local_photos.copy()
    bad.loc[0,"latitude"]=np.nan
    with pytest.raises(ValueError,match="missing or invalid"):
        M.microgroups(bad,50)

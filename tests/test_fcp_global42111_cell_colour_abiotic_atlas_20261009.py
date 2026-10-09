"""162 original equal-area cell atlas preserves source photo denominators and honest soil gaps."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"/"analysis"))
import build_fcp_global42111_cell_colour_abiotic_atlas_20261009 as M


@pytest.fixture
def world(monkeypatch):
    monkeypatch.setattr(M,"N_SOURCE",6)
    monkeypatch.setattr(M,"N_CLASSIFIED",3)
    records=[]
    for i,(cid,t,colour,status) in enumerate([
        (0,1,"white","classified_four_state_morph"),
        (0,2,"red_pink","classified_four_state_morph"),
        (30,1,"white","classified_four_state_morph"),
        (30,3,"unclassified","unclassified"),
        (79,4,"","unclassified"),
        (79,5,"","unclassified"),
    ]):
        good=i!=5
        row={"cell_id":cid,"inat_taxon_id":t,"morph":colour,
             "measurement_status":status,
             "site_geo_status":"VALID_PUBLIC_ORIGINAL_PHOTO_POINT" if good else "NO_PUBLIC_GEO",
             "environment_climate_complete":good,
             "environment_soil_complete":good and i!=4,
             "environment_all_complete":good and i!=4}
        for v in M.FEATURES:
            row[v]=float(i+1) if row["environment_all_complete"] else np.nan
        records.append(row)
    return pd.DataFrame(records)


def test_original_denominators_map_to_162_cells_and_all_photo_colours(world):
    report,cells,bands,geojson=M.color_and_environment_atlas(world)
    assert len(cells)==162
    assert int(cells.source_taxon_cell_rows.sum())==6
    assert int(cells.n_four_state_classified.sum())==3
    assert report["n_original_source_occupied_cells"]==3
    assert report["original_taxon_cell_photo_colour_counts"]=={
        "white":2,"yellow_orange":0,"red_pink":1,"blue_purple":0}
    assert len(geojson["features"])==162
    assert all(f["geometry"]["type"]=="Polygon" for f in geojson["features"])


def test_empty_cells_unclassified_and_conditional_denominators(world):
    _,cells,_,_=M.color_and_environment_atlas(world)
    row=cells.loc[cells.cell_id==0].iloc[0]
    assert row["fraction_white_given_classified"]==pytest.approx(.5)
    assert row["count_red_pink"]==1
    empty=cells.loc[cells.cell_id==50].iloc[0]
    assert empty.source_taxon_cell_rows==0
    assert pd.isna(empty["fraction_white_given_classified"])
    assert pd.isna(empty["median_real_photo_soil_pH"])


def test_site_environment_extracted_before_cell_display_not_at_centres(world):
    _,cells,_,_=M.color_and_environment_atlas(world)
    row=cells.loc[cells.cell_id==0].iloc[0]
    assert row["median_real_photo_soil_pH"]==pytest.approx(1.5)
    assert row["display_only_cell_centroid_lon"]==-170
    assert row["display_only_cell_centroid_lat"]<0


def test_unlocated_photo_must_not_receive_fake_soil(world):
    world.loc[5,"soil_pH"]=6
    with pytest.raises(ValueError,match="artificial climate or soil"):
        M.color_and_environment_atlas(world)


def test_never_turn_unclassified_into_white_and_never_duplicate_species_cell(world):
    world.loc[4,"morph"]="white"
    _,cells,_,_=M.color_and_environment_atlas(world)
    assert int(cells.n_four_state_classified.sum())==3
    changed=world.copy()
    changed.loc[1,"inat_taxon_id"]=1
    with pytest.raises(ValueError,match="unique identity"):
        M.color_and_environment_atlas(changed)

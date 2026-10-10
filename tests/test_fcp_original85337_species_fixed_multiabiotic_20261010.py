"""Synthetic repeat photos per plant species across independently heldout cells."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"/"analysis"))
import compare_fcp_original85337_species_fixed_multiabiotic_20261010 as M
from extend_fcp_global42111_worldclim_solar_20261010 import NEW_FEATURES


@pytest.fixture
def source(monkeypatch):
    monkeypatch.setattr(M,"N_CELL",800)
    monkeypatch.setattr(M,"N_CLASSIFIED",600)
    monkeypatch.setattr(M,"N_SOURCE_SITES",800)
    monkeypatch.setattr(M,"MIN_TEST",100)
    rng=np.random.default_rng(20261010)
    # Every one of 120 taxa has multiple SOURCE photographed regions; 800 original photo cells.
    i=np.arange(800)
    species=i%120+1
    old=pd.DataFrame({
        "inat_taxon_id":species,
        "species":[f"Genus{t%30} species{t}" for t in species],
        "observation_id":30000+i,"photo_id":40000+i,
        "cell_id":(i//120)%8+45,
        "morph":np.where(i<600,np.array(M.CLASSES)[i%4],"UNCLASSIFIED"),
        "measurement_status":np.where(i<600,"classified_four_state_morph","roi_unavailable"),
        "site_geo_status":"VALID_PUBLIC_ORIGINAL_PHOTO_POINT",
        "latitude":rng.uniform(-30,30,800),
        "longitude":rng.uniform(50,100,800),
    })
    for p in tuple(f for k,v in M.BLOCKS.items() if k!="soil" for f in v):
        if p not in old:old[p]=rng.uniform(0,20,800)
    for p in M.BLOCKS["soil"]:
        if p not in old:old[p]=rng.uniform(0,20,800)
    sites=old[["inat_taxon_id","observation_id","photo_id"]].copy()
    for key in NEW_FEATURES:
        sites[key]=old[key].values
        old=old.drop(columns=[key])
    sites["present_in_breadth"]=True
    sites["present_in_taxon_cell"]=True
    return old,sites


def test_all_photo_ID_triples_match_and_unclassifiable_preserved(source):
    old,sites=source
    merged,rec=M.match_original(old,sites,strict=True)
    assert len(merged)==800
    assert rec["source_four_state_colour_classifiable"]==600
    assert rec["source_unclassifiable_photo_cell_records"]==200
    assert rec["all_original_photo_ID_triples_matched_exactly"]
    assert set(NEW_FEATURES).issubset(merged)


def test_original_photo_id_taxon_swap_fails_closed(source):
    old,sites=source
    sites.loc[0,"photo_id"]=99999
    with pytest.raises(RuntimeError,match="key lost"):
        M.match_original(old,sites)


def test_species_photo_region_folds_never_reuse_same_observation_in_test_train(source):
    old,sites=source
    merged,_=M.match_original(old,sites)
    pop,coverage=M.populations(merged)
    assert coverage["n_original_fourstate_classified_photo_cells"]==600
    assert len(pop["CLIMATE_SOURCE_ONLY"])==600
    assert len(pop["CLIMATE_AND_SOIL"])==600
    out=M.species_out_of_cell(pop["CLIMATE_SOURCE_ONLY"],with_soil=False,nboot=19)
    assert out["status"]=="SOURCE_ORIGINAL_SPECIES_INTERCEPT_OUTOFCELL_ABIOTIC_DIAGNOSTIC"
    assert out["n_heldout_original_photos_with_species_training_intercept"]>=100
    assert out["training_species_colour_means_only_from_other_cells"]
    assert out["species_fixed_intercept_eliminates_unidentified_lineage_static_baseline"]
    assert set(out["source_block_incremental_prediction"])==(
        {"all_abiotic_beyond_species_geography"}|
        {"conditional_"+x for x in M.BLOCKS if x!="soil"})
    assert all(v["same_heldout_source_photos_all_models"] for v in out["source_block_incremental_prediction"].values())


def test_original_photo_colour_not_needed_to_select_taxon_cell_identical_source_sites(source):
    old,sites=source
    merge,_=M.match_original(old,sites)
    changed=old.copy()
    changed["morph"]="white"
    changed["measurement_status"]="classified_four_state_morph"
    alt,_=M.match_original(changed,sites,strict=False)
    assert len(alt)==len(merge)
    pd.testing.assert_series_equal(alt.photo_id,merge.photo_id)
    assert alt[M.BLOCKS["soil"][0]].equals(merge[M.BLOCKS["soil"][0]])


def test_species_intercept_train_only_no_heldout_photo_label_leakage():
    train=pd.DataFrame({
        "inat_taxon_id":[1,1,2,2],"cell_id":[10,11,10,11],
        "morph":["white","white","blue_purple","blue_purple"],
        "wc_bio12":[10.,20.,30.,40.]})
    test=pd.DataFrame({"inat_taxon_id":[1,2],"cell_id":[12,12],
                       "morph":["red_pink","yellow_orange"],"wc_bio12":[15.,35.]})
    a=M.train_predict(train,test,("wc_bio12",),np.array([0,1]))
    test["morph"]="white"
    b=M.train_predict(train,test,("wc_bio12",),np.array([0,1]))
    assert np.allclose(a,b)
    assert a.shape==(2,4)
    assert np.allclose(a.sum(axis=1),1)


def test_no_outofsample_species_heldout_baseline_invented():
    train=pd.DataFrame({"inat_taxon_id":[1,1],"cell_id":[10,11],
                        "morph":["white","white"],"wc_bio12":[10.,20.]})
    test=pd.DataFrame({"inat_taxon_id":[2],"cell_id":[12],
                       "morph":["red_pink"],"wc_bio12":[15.]})
    with pytest.raises(RuntimeError,match="Unknown"):
        M.train_predict(train,test,("wc_bio12",),np.array([0]))

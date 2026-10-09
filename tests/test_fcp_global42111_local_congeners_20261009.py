"""Test local same-genus×same-cell sourced photographs, no outcome-selected groups."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"/"analysis"))
import audit_fcp_global42111_local_congeners_20261009 as M


@pytest.fixture
def original():
    rng=np.random.default_rng(92026)
    n=1000
    i=np.arange(n)
    genus=i%50
    # 50 independent genera, exactly one coarse photographed cell per genus,
    # 20 different named species per local genus-cell group.
    lat=14.0+rng.uniform(-1.2,1.2,n)
    lon=-170+20*(genus%10)+rng.uniform(-2.0,2.0,n)
    d=pd.DataFrame({
        "inat_taxon_id":i+1,
        "species":[f"Genus{g} species{j}" for j,g in enumerate(genus)],
        "morph":[M.CLASSES[(j//50+j%7)%4] for j in i],
        "measurement_status":"classified_four_state_morph",
        "latitude":lat,"longitude":lon,
        "site_geo_status":"VALID_PUBLIC_ORIGINAL_PHOTO_POINT",
        "environment_climate_complete":True,
        "environment_soil_complete":i%8!=0,
        "wc_elevation_m":rng.uniform(10,900,n),
    })
    for f in M.TEMP+M.RAIN+M.SOIL:
        d[f]=rng.normal(10,3,n)
    d.loc[i%8==0,list(M.SOIL)]=np.nan
    return d


def test_sample_source_denominators_and_same_cell_support(original):
    pops,meta=M.photo_populations(original,strict=False)
    assert meta["source_taxa"]==1000
    assert len(pops["CLIMATE_ALL"])==1000
    assert len(pops["CLIMATE_SOIL_COMPLETE"])==875
    audit=M.source_group_coverage(pops["CLIMATE_ALL"])
    assert audit["n_groups_with_three_or_more_original_species"]>=45
    assert audit["n_photo_species_in_groups_with_three_or_more"]==1000


def test_train_baseline_uses_other_species_not_test_photo_label():
    train=pd.DataFrame({
        "inat_taxon_id":[1,2,3],
        "genus_cell_id":["GenusA|50"]*3,
        "morph":["white","white","red_pink"],
        "wc_bio12":[100.,120.,140.],
    })
    test=pd.DataFrame({
        "inat_taxon_id":[4],"genus_cell_id":["GenusA|50"],
        "morph":["blue_purple"],"wc_bio12":[127.],
    })
    base=M.train_only_group_predict(train,test,())
    moist=M.train_only_group_predict(train,test,("wc_bio12",))
    assert base[0,0]==pytest.approx(2/3)
    assert np.isclose(moist.sum(),1.0)
    test["morph"]="white"
    assert np.array_equal(moist,M.train_only_group_predict(train,test,("wc_bio12",)))


def test_taxon_hash_fold_not_outcome_dependent(original):
    d,_=M.photo_populations(original,strict=False)
    f=M.fold_plan(d["CLIMATE_ALL"].reset_index(drop=True))
    changed=d["CLIMATE_ALL"].reset_index(drop=True).copy()
    changed["morph"]="white"
    g=M.fold_plan(changed)
    assert np.array_equal(f,g)
    assert len(np.unique(f))==5


def test_geographic_group_baseline_and_moisture_model_same_taxa(original):
    d,_=M.photo_populations(original,strict=False)
    z=M.test_local_congeners(d["CLIMATE_ALL"],include_soil=False)
    assert z["status"]=="SOURCE_LOCAL_GENUS_CELL_HELDOUT_PHOTO_SPECIES_ANALYSIS"
    assert z["n_original_photo_species_evaluated"]==1000
    assert z["n_original_taxonomic_genera_evaluated"]==50
    assert set(z["heldout_4class_photo_brier"])==set(M.FAMILIES)-{"GENUS_CELL_GEO_ALL_CLIMATE_SOIL"}
    assert {x["source_test_photo_species"] for x in z["heldout_4class_photo_brier"].values()}=={1000}
    for delta in z["incremental_source_photo_prediction"].values():
        assert len(delta["original_cell_cluster_bootstrap_95CI"])==2
        assert len(delta["genus_cluster_bootstrap_95CI"])==2
        assert len(delta["local_genus_cell_cluster_bootstrap_95CI"])==2


def test_soil_completeness_changes_sample_but_never_changes_photos(original):
    d,_=M.photo_populations(original,strict=False)
    full=M.test_local_congeners(d["CLIMATE_SOIL_COMPLETE"],include_soil=True)
    assert full["status"]=="SOURCE_LOCAL_GENUS_CELL_HELDOUT_PHOTO_SPECIES_ANALYSIS"
    assert full["n_original_photo_species_evaluated"]==875
    assert "soil_beyond_local_geography_climate" in full["incremental_source_photo_prediction"]


def test_groups_with_only_two_species_are_held_not_reclassified(original):
    z=original.groupby(original.species.str.split().str[0],sort=True).head(2)
    src,_=M.photo_populations(z,strict=False)
    out=M.test_local_congeners(src["CLIMATE_ALL"],include_soil=False)
    assert out["status"]=="HOLD_INSUFFICIENT_SAME_GENUS_SAME_CELL_SPECIES"
    assert out["original_photo_support"]["n_photo_species_in_groups_with_three_or_more"]==0


def test_missing_photo_status_and_duplicates_never_synthetic_labels(original):
    original.loc[0,"measurement_status"]="UNCLASSIFIABLE"
    pool,coverage=M.photo_populations(original,strict=False)
    assert coverage["source_unclassified_photo_taxa"]==1
    assert len(pool["CLIMATE_ALL"])==999
    duplicate=original.copy()
    duplicate.loc[1,"inat_taxon_id"]=1
    with pytest.raises(ValueError,match="duplicated"):
        M.photo_populations(duplicate,strict=False)


def test_inference_explicitly_not_within_species(original):
    out=M.run(original,strict=False)
    assert out["original_42111_denominator"]==1000
    assert out["no_new_original_photo_colour_or_raster_read"] is True
    assert out["climate_only_local_congeners"]["within_same_genus_and_same_source_162_cell_only"]
    assert out["soil_complete_local_congeners"]["one_photo_per_source_species"] is True

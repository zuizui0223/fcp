"""Synthetic source-photo genotype-free between-species within-genus tests."""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"/"analysis"))
# Load the analysis file explicitly: this test module shares a test_ prefix.
import importlib.util
SCRIPT=Path(__file__).resolve().parents[1]/"scripts"/"analysis"/"test_fcp_global42111_within_genus_environment_20261009.py"
SPEC=importlib.util.spec_from_file_location("fcp_within_genus_analysis_source",SCRIPT)
assert SPEC is not None and SPEC.loader is not None
M=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)


@pytest.fixture
def source():
    n=1200
    ix=np.arange(n)
    rng=np.random.default_rng(19210)
    g=ix%100
    lon=(ix%27-13)*12.7
    lat=(ix%17-8)*6.3
    data=pd.DataFrame({
        "inat_taxon_id":ix+1,
        "species":[f"Genus{k} species{i}" for i,k in enumerate(g)],
        "morph":[M.CLASSES[(i//100+i%7)%4] for i in ix],
        "measurement_status":"classified_four_state_morph",
        "site_geo_status":"VALID_PUBLIC_ORIGINAL_PHOTO_POINT",
        "latitude":lat,"longitude":lon,
        "environment_climate_complete":True,
        "environment_soil_complete":ix%9!=0,
        "wc_elevation_m":rng.uniform(20,2200,n),
    })
    for name in M.CLIMATE+M.SOIL:
        data[name]=rng.normal(loc=10,scale=3,size=n)
    data.loc[ix%9==0,list(M.SOIL)]=np.nan
    return data


def test_source_frame_preserves_all_species_and_soil_missing(source):
    pop,stats=M.source_populations(source,strict=False)
    assert stats["source_single_photo_species"]==1200
    assert len(pop["CLIMATE_ALL"])==1200
    assert len(pop["CLIMATE_SOIL_COMPLETE"])<1200
    assert pop["CLIMATE_ALL"].inat_taxon_id.nunique()==1200
    assert len(M.cell_index(np.array([12.,float("nan")]),np.array([10.,10.])))==2


def test_genus_intercepts_are_train_only_test_morph_cannot_leak():
    train=pd.DataFrame({
        "inat_taxon_id":[1,2,3,4],
        "genus":["A","A","A","B"],
        "morph":["white","white","red_pink","blue_purple"],
        "abs_latitude":[20.,22.,24.,30.],
        "lon_sin":[.1,.2,.3,.4],
        "lon_cos":[.9,.8,.7,.6],
        "wc_elevation_m":[50.,60.,70.,80.],
        "wc_bio1":[10.,12.,16.,20.],
    })
    test=pd.DataFrame({
        "inat_taxon_id":[5,6],
        "genus":["A","B"],
        "morph":["white","yellow_orange"],
        "abs_latitude":[23.,30.],
        "lon_sin":[.25,.4],"lon_cos":[.75,.6],
        "wc_elevation_m":[65.,85.],
        "wc_bio1":[14.,22.],
    })
    idx0,base=M.training_genus_predictors(train,test,())
    idx1,other=M.training_genus_predictors(train,test,("wc_bio1",))
    assert idx0.tolist()==idx1.tolist()==[0]
    assert base.shape==(1,4)
    assert other.shape==(1,4)
    assert base[0,0]==pytest.approx(2/3)
    assert np.isclose(other.sum(),1.0)
    test.loc[0,"morph"]="blue_purple"
    _,changed=M.training_genus_predictors(train,test,("wc_bio1",))
    assert np.array_equal(other,changed)


def test_heldout_cell_and_genus_sample_is_identical_between_models(source):
    populations,_=M.source_populations(source,strict=False)
    out=M.fivefold_within_genus_cv(populations["CLIMATE_SOIL_COMPLETE"],with_soil=True)
    assert out["same_exact_species_and_cell_folds_for_all_feature_families"]
    n=out["source_species_with_train_observed_genus_and_min_two_training_species"]
    assert n>=100
    assert sum(f["n_genus_estimable_test_species"] for f in out["heldout_5fold_original_cell_status"])==n
    assert all(q["n_same_source_test_species"]==n for q in out["brier_score_models"].values())
    assert out["brier_score_models"]["GENUS_BASELINE"]["heldout_multiclass_brier"]>=0
    assert set(("GENUS_GEO_TEMPERATURE","GENUS_GEO_MOISTURE","GENUS_GEO_CLIMATE")).issubset(out["brier_score_models"])
    for term in ("unique_temperature_block_beyond_geo_moisture",
                 "unique_moisture_block_beyond_geo_temperature"):
        item=out["fixed_fold_incremental_gains"][term]
        assert np.isfinite(item["mean_heldout_brier_reduction"])
        assert len(item["genus_equal_genus_cluster_bootstrap_95CI"])==2

    assert len(out["fixed_fold_incremental_gains"]["climate_beyond_genus_geography"]["source_cell_block_bootstrap_95CI"])==2
    for k in out["fixed_fold_incremental_gains"]:
        row=out["fixed_fold_incremental_gains"][k]
        assert len(row["genus_equal_source_cell_bootstrap_95CI"])==2
        assert len(row["genus_equal_genus_cluster_bootstrap_95CI"])==2
        assert row["n_genera_in_genus_equal_sensitivity"]==out["n_distinct_source_genera_in_evaluation"]
        assert 0<=row["fraction_evaluated_genera_with_positive_increment"]<=1
    for m in out["brier_score_models"].values():
        assert m["heldout_genus_equal_multiclass_brier"]>=0


def test_mismatched_soil_never_imputed_into_model(source):
    _,c=M.source_populations(source,strict=False)
    assert c["climate_only_additional_species_not_available_in_soil_complete"]>0
    source.loc[0,"site_geo_status"]="NO_UNOBSCURED_PUBLIC_COORDINATES"
    pool,stats=M.source_populations(source,strict=False)
    assert not 1 in set(pool["CLIMATE_ALL"].inat_taxon_id)
    assert stats["original_photograph_colour_classifiable"]==1200


def test_repeated_taxon_identity_fails_closed(source):
    source.loc[1,"inat_taxon_id"]=1
    with pytest.raises(ValueError,match="Repeated"):
        M.source_populations(source,strict=False)


def test_outcome_status_missing_does_not_create_color_state(source):
    source.loc[0,"measurement_status"]="unclassifiable"
    pop,stats=M.source_populations(source,strict=False)
    assert len(pop["CLIMATE_ALL"])==1199
    assert stats["unclassified_photo_species_preserved_outside_model"]==1


def test_end_to_end_heldout_estimand_not_within_species(source):
    result=M.run(source,strict=False)
    assert result["schema"]=="fcp_global42111_within_genus_compositional_vs_climate_soil_v1"
    assert result["climate_only_population"]["source_species_in_pool"]==1200
    assert result["soil_complete_population"]["source_species_in_pool"]<1200
    assert result["no_photo_pixels_or_new_colour_labels_opened"] is True
    assert "between different species" in result["climate_only_population"]["inference_boundary"].lower()

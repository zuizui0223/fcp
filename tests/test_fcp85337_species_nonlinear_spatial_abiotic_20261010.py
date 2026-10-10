"""Original FCP 85,337 repeated photographed colour species + spherical spatial controls."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"/"analysis"))
import audit_fcp85337_species_nonlinear_spatial_abiotic_20261010 as M
import compare_fcp_original85337_species_fixed_multiabiotic_20261010 as SOURCE
from compare_fcp_global42111_full_climate_sun_soil_blocks_20261010 import BLOCKS,CLASSES
from extend_fcp_global42111_worldclim_solar_20261010 import NEW_FEATURES


@pytest.fixture
def original(monkeypatch):
    monkeypatch.setattr(SOURCE,"N_CELL",240)
    monkeypatch.setattr(SOURCE,"N_CLASSIFIED",200)
    monkeypatch.setattr(SOURCE,"N_SOURCE_SITES",240)
    monkeypatch.setattr(M,"MIN_TEST",30)
    i=np.arange(240);rng=np.random.default_rng(22026)
    orig=pd.DataFrame({
        "inat_taxon_id":i%40+1,
        "species":[f"G{j%20} species{j}" for j in i%40+1],
        "observation_id":3000+i,"photo_id":4000+i,
        "cell_id":i//40+45,
        "morph":np.where(i<200,np.array(CLASSES)[(i//40+i)%4],"not_classified"),
        "measurement_status":np.where(i<200,"classified_four_state_morph","unclassifiable"),
        "latitude":rng.uniform(-25,25,240),"longitude":rng.uniform(-160,160,240),
        "site_geo_status":"VALID_PUBLIC_ORIGINAL_PHOTO_POINT",
    })
    for k in {t for v in BLOCKS.values() for t in v}:
        orig[k]=rng.normal(5,2,240)
    site=orig[["inat_taxon_id","observation_id","photo_id"]].copy()
    for k in NEW_FEATURES:
        site[k]=orig[k].to_numpy(float)
        orig=orig.drop(columns=[k])
    return orig,site


def test_spherical_polynomial_real_photo_geo_exactly_19_basis_columns(original):
    old,site=original
    d=M.add_spatial_terms(old)
    assert len(M.SPATIAL)==19
    assert d[list(M.SPATIAL)].notna().all().all()
    assert d["sphere_monomial_1"].between(-1,1).all()
    fake=old.copy()
    fake.loc[0,"site_geo_status"]="NO_ORIGINAL_PUBLIC_PHOTO_LOCATION"
    out=M.add_spatial_terms(fake)
    assert out.loc[0,list(M.SPATIAL)].isna().all()


def test_all_original_species_photo_abiotic_denominators_kept(original):
    old,site=original
    combined,proof=SOURCE.match_original(old,site)
    pools,receipt=SOURCE.populations(combined)
    assert proof["source_original_taxon_cell_photo_rows"]==240
    assert proof["source_four_state_colour_classifiable"]==200
    assert receipt["n_original_fourstate_unclassified_photo_cells"]==40
    assert len(pools["CLIMATE_SOURCE_ONLY"])==200


def test_fixed_source_region_and_species_photo_support_same_for_all_models(original):
    old,site=original
    frame,_=SOURCE.match_original(old,site)
    pools,_=SOURCE.populations(frame)
    q=M.spatial_species_heldout(pools["CLIMATE_SOURCE_ONLY"],soil=False,nboot=15)
    assert q["status"]=="SOURCE_SPECIES_FIXED_WITH_NONLINEAR_SPATIAL_BASIS_ABIOTIC_EXPLORATION"
    assert q["n_original_outofcell_photo_species_heldout"]==200
    assert q["n_distinct_original_species_heldout"]==40
    assert q["all_model_source_test_photos_and_fold_assignments_identical"]
    assert len(q["model_four_colour_brier_scores"])==3+6
    assert "all_environment_beyond_species_plus_nonlinear_space" in q["incremental_predictive_gains"]
    assert all(len(z["source_cell_cluster_95CI"])==2 for z in q["incremental_predictive_gains"].values())


def test_soil_complete_same_species_all_seven_environment_blocks(original):
    old,site=original
    frame,_=SOURCE.match_original(old,site)
    pools,_=SOURCE.populations(frame)
    v=M.spatial_species_heldout(pools["CLIMATE_AND_SOIL"],soil=True,nboot=13)
    assert v["status"]=="SOURCE_SPECIES_FIXED_WITH_NONLINEAR_SPATIAL_BASIS_ABIOTIC_EXPLORATION"
    assert len(v["model_four_colour_brier_scores"])==10
    assert "unique_soil_beyond_species_nonlinear_space_other_environment" in v["incremental_predictive_gains"]


def test_unclassified_source_photo_never_replaced_with_white(original):
    old,site=original
    old.loc[0,"morph"]="white"
    old.loc[0,"measurement_status"]="ROI_NOT_CLASSIFIED"
    _,photo=SOURCE.match_original(old,site,strict=False)
    assert photo["source_four_state_colour_classifiable"]==199


def test_noncausal_full_response_retains_85337_contract(original):
    old,site=original
    report=M.run(old,site,strict=False,nboot=10)
    assert report["all_original_photo_ids_and_flower_colour_labels_preserved"]
    assert report["source_phylogenetic_species_invariant_component_absorbed_by_species_intercepts"]
    assert report["genuine_spatial_autocorrelation_random_field_not_estimable_from_this_basis_alone"]

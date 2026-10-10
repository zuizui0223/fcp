"""Original within-species photo prediction after TRAIN-only Nyström geodesic covariance."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"/"analysis"))
import audit_fcp85337_true_spatial_kernel_allenv_20261010 as M


def test_genuine_photo_geodesic_kernel_responds_to_actual_site_distance():
    a=M.xyz(np.array([0.,0.,0.]),np.array([0.,1.,20.]))
    K=M.geodesic_exponential_kernel(a,a,100)
    assert K.shape==(3,3)
    assert np.allclose(np.diag(K),1.,atol=1e-5)
    assert 0<K[0,2]<K[0,1]<1
    assert np.linalg.eigvalsh(K).min()>-1e-6
    with pytest.raises(ValueError,match="missing or impossible"):
        M.xyz([np.nan,0],[0,2])


def test_source_landmarks_train_only_and_valid_features():
    xyz=M.xyz([0,0,0,0,1,1,1,1,2,2,2,2],[0,1,2,3,0,1,2,3,0,1,2,3])
    tr,te,r=M.nystrom_train_test(xyz[:10],xyz[10:],scale_km=500,nlandmarks=6)
    assert tr.shape==(10,6) and te.shape==(2,6)
    assert r["landmark_selection_colour_outcome_blind"]
    assert np.isfinite(tr).all() and np.isfinite(te).all()
    tr2,te2,_=M.nystrom_train_test(xyz[:10],xyz[10:],scale_km=500,nlandmarks=6)
    assert np.allclose(tr,tr2) and np.allclose(te,te2)


@pytest.fixture
def repeated_species_cells():
    rng=np.random.default_rng(20261010)
    n=360
    i=np.arange(n)
    cell=i//72+65
    species=i%72+1
    d=pd.DataFrame({
        "inat_taxon_id":species,
        "species":[f"G{s%15} sp{s}" for s in species],
        "cell_id":cell,
        "morph":np.array(M.CLASSES)[(species+i//72)%4],
        "latitude":rng.uniform(-20,20,n),
        "longitude":rng.uniform(-80,80,n),
    })
    d["abs_latitude"]=d.latitude.abs()
    d["lon_sin"]=np.sin(np.deg2rad(d.longitude))
    d["lon_cos"]=np.cos(np.deg2rad(d.longitude))
    for fields in M.BLOCKS.values():
        for col in fields:
            d[col]=rng.uniform(.5,20,n)
    return d


def test_repeated_species_outofregion_predictors_same_rows_and_folds(repeated_species_cells):
    d=repeated_species_cells
    z=M.fit_source(d,soil=False,scale_km=500,nboot=9,nlandmarks=12)
    assert z["status"]=="SOURCE_SPECIES_INTERCEPT_GEODESIC_NYSTROM_SPATIAL_KERNEL_EXPLORATION"
    assert z["n_original_heldout_source_photos_with_training_species"]==360
    assert z["n_distinct_original_source_species_in_heldout"]==72
    assert z["source_colour_means_by_species_training_cells_only"]
    assert z["all_models_same_test_photo_ids_and_geocell_folds"]
    assert z["n_source_training_only_spatial_landmarks"]==12
    assert "all_abiotic_beyond_species_real_geodesic_spatial_covariance" in z["environment_increment"]
    assert set(z["fixed_prediction_nearby_residual_photo_spatial_Moran_diagnostic"])=={"100","500"}
    assert all("n_nearby_directed_edges" in q for q in
               z["fixed_prediction_nearby_residual_photo_spatial_Moran_diagnostic"].values())
    assert set(z["environment_increment"])=={
        "spatial_kernel_beyond_linear_geography",
        "all_abiotic_beyond_species_real_geodesic_spatial_covariance",
        *["unique_"+k+"_beyond_species_spatial_kernel_and_other_abiotic"
          for k in M.BLOCKS if k!="soil"],
    }
    assert all(len(v["species_block_bootstrap_95CI"])==2 and
               len(v["geographical_cell_block_bootstrap_95CI"])==2
               for v in z["environment_increment"].values())


def test_species_site_coverage_and_soil_complete_same_test_ids(repeated_species_cells):
    d=repeated_species_cells
    z=M.fit_source(d,soil=True,scale_km=250,nboot=7,nlandmarks=12)
    assert z["n_original_heldout_source_photos_with_training_species"]==360
    assert "unique_soil_beyond_species_spatial_kernel_and_other_abiotic" in z["environment_increment"]
    assert all(q["n_same_original_heldout_photo_rows"]==360 for q in z["fourstate_model_scores"].values())


def test_insufficient_geographic_or_original_species_fails_closed(repeated_species_cells):
    d=repeated_species_cells.iloc[:50]
    z=M.fit_source(d,soil=False,scale_km=100,nboot=2,nlandmarks=12)
    assert z["status"]=="HOLD_SOURCE_SPECIES_REGION_SPATIAL_COVERAGE"


def test_residual_spatial_Moran_uses_original_photo_sites_and_four_classes():
    n=20
    rng=np.random.default_rng(12)
    positions=pd.DataFrame({
        "latitude":rng.normal(10,.01,n),
        "longitude":rng.normal(100,.01,n)})
    y=np.arange(n)%4
    correct=np.eye(4)[y]
    predicted={
        "mock_baseline":np.full((n,4),.25),
        "mock_perfect":correct,
    }
    result=M.residual_neighbour_autocorrelation(positions,y,predicted)
    assert set(result)=={"100","500"}
    z=result["100"]
    assert z["status"]=="ORIGINAL_NEAREST_PHOTO_RESIDUAL_SPATIAL_AUTOCORRELATION_DIAGNOSTIC"
    assert z["n_nearby_directed_edges"]>0
    assert len(z["residual_methods"]["mock_baseline"]["source_four_colour_class_residual_Moran_I"])==4
    assert z["descriptive_not_a_significance_test"] is True
    bad=positions.copy()
    bad.loc[0,"latitude"]=np.nan
    with pytest.raises(ValueError,match="missing or impossible"):
        M.residual_neighbour_autocorrelation(bad,y,predicted)

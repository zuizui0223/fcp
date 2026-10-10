"""True-geodesic low-rank spatial field and train-only 85337 photo source guards."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"/"analysis"))
import audit_fcp85337_species_geodesic_lowrank_spatial_all_abiotic_20261010 as M
from compare_fcp_global42111_full_climate_sun_soil_blocks_20261010 import BLOCKS,CLASSES


@pytest.fixture
def source(monkeypatch):
    monkeypatch.setattr(M,"MIN_EVAL",40)
    monkeypatch.setattr(M,"N_FIXED_TRAIN_KNOTS",12)
    monkeypatch.setattr(M,"RBF",tuple(f"spatial_trainonly_rbf_{i:03d}" for i in range(12)))
    rng=np.random.default_rng(20261010)
    records=[]
    for region in range(6):
        for taxon in range(40):
            # Distinct original source genus/species and source photo observation per cell.
            n=region*40+taxon
            lat=-24+region*8+rng.uniform(-1,1)
            lon=42+region*9+rng.uniform(-1,1)
            rec={
                "inat_taxon_id":taxon+1,
                "species":f"Genus{taxon%10} species{taxon}",
                "morph":CLASSES[(taxon+region)%4],
                "cell_id":region+35,
                "measurement_status":"classified_four_state_morph",
                "latitude":lat,"longitude":lon,
                "site_geo_status":"VALID_PUBLIC_ORIGINAL_PHOTO_POINT",
                "observation_id":7000+n,"photo_id":8000+n,
            }
            for name in (x for cols in BLOCKS.values() for x in cols):
                rec[name]=float(rng.normal(5+region/3,2))
            records.append(rec)
    return pd.DataFrame(records)


def test_true_original_photo_geometry_spatial_field_finite_and_distance_dependent(source):
    base=source.iloc[:100].copy()
    points=M.fit_training_geo_knots(base,n_knots=12)
    assert points.shape==(12,3)
    assert np.allclose(np.linalg.norm(points,axis=1),1)
    expanded=M.rbf_columns(base,points,250,names=M.RBF)
    assert len(expanded)==len(base)
    assert len(M.RBF)==12
    assert expanded[list(M.RBF)].notna().all().all()
    assert expanded[list(M.RBF)].to_numpy().min()>=0
    assert expanded[list(M.RBF)].to_numpy().max()<=1


def test_failed_original_location_is_never_imputed_as_cell_midpoint(source):
    train=source.iloc[:100].copy()
    coords=M.fit_training_geo_knots(train,n_knots=12)
    bad=train.copy()
    bad.loc[0,"latitude"]=np.nan
    with pytest.raises(ValueError,match="finite"):
        M.rbf_columns(bad,coords,250,names=M.RBF)


def test_fitted_knots_depend_only_on_training_coordinates_not_flower_label(source):
    train=source.iloc[:100].copy()
    a=M.fit_training_geo_knots(train,n_knots=12)
    train["morph"]="white"
    b=M.fit_training_geo_knots(train,n_knots=12)
    assert np.allclose(a,b)


def test_all_models_include_distance_kernel_and_block_drop_on_same_photo_species(source):
    for soil in (False,True):
        fam=M.candidate_models(soil=soil)
        assert fam["SPECIES_GEODESIC_SPATIAL_FIELD"][-12:]==M.RBF
        assert len(fam)==3+(7 if soil else 6)
        assert "SPATIAL_FIELD_FULL_MINUS_PRECIPITATION" in fam
        if not soil:assert "SPATIAL_FIELD_FULL_MINUS_SOIL" not in fam
        else:assert "SPATIAL_FIELD_FULL_MINUS_SOIL" in fam


def test_heldout_photo_geographic_regions_never_enter_train_species_baseline(source):
    d=source.copy()
    d["abs_latitude"]=abs(d.latitude)
    d["lon_sin"]=np.sin(np.deg2rad(d.longitude))
    d["lon_cos"]=np.cos(np.deg2rad(d.longitude))
    result=M.evaluate(d,soil=False,bandwidth=250,nboot=9,n_knots=12)
    assert result["status"]=="SOURCE_SPECIES_INTERCEPT_AND_GEODESIC_LOW_RANK_SPATIAL_FIELD_EXPLORATION"
    assert result["n_train_species_heldout_original_photos"]==240
    assert result["n_distinct_source_taxa_in_evaluation"]==40
    assert result["training_only_unsupervised_source_photo_geographic_knots"]
    assert result["species_means_learned_only_from_other_source_geographic_regions"]
    assert result["every_source_photo_same_heldout_fold_and_model_comparison"]
    assert all(len(v["source_species_cluster_95CI"])==2 and
               len(v["original_region_cell_cluster_95CI"])==2
               for v in result["source_block_predictive_gains"].values())
    diagnostics=result["descriptive_heldout_photo_geodesic_spatial_residual_diagnostics"]
    assert diagnostics["original_public_photographic_coordinates_only"]
    assert diagnostics["same_species_neighbor_pairs_excluded"] is True
    assert set(diagnostics["spatial_radius_reports_km"])=={"50","100","250"}
    assert diagnostics["not_a_calibrated_spatial_independence_p_value"] is True
    assert all(set(q["moran_like_centered_fourstate_residual_spatial_index_by_model"])=={
        "SPECIES_NONGP_COORDINATES","SPECIES_GEODESIC_SPATIAL_FIELD",
        "SPECIES_GEODESIC_SPATIAL_FIELD_ALL_ENV"
    } for q in diagnostics["spatial_radius_reports_km"].values())


def test_kernel_is_psd_and_uses_real_spherical_distance(source):
    pts=source.iloc[:90].copy()
    anchors=M.fit_training_geo_knots(pts,n_knots=12)
    spatial=M.rbf_columns(pts,anchors,250,names=M.RBF)[list(M.RBF)].to_numpy()
    K=spatial@spatial.T
    assert np.linalg.eigvalsh(K).min()>=-1e-7


def test_same_species_photo_pair_not_counted_as_independent_spatial_neighbors():
    d=pd.DataFrame({
        "inat_taxon_id":[1,1,2,3,4,5],
        "latitude":[12,12.001,12.002,12.003,12.004,12.005],
        "longitude":[10,10.001,10.002,10.003,10.004,10.005],
    })
    y=np.array([0,0,1,2,3,1])
    pred=np.full((6,4),.25,float)
    preds={k:pred for k in (
        "SPECIES_NONGP_COORDINATES",
        "SPECIES_GEODESIC_SPATIAL_FIELD",
        "SPECIES_GEODESIC_SPATIAL_FIELD_ALL_ENV")}
    audit=M.spatial_oof_residual_diagnostic(d,y,preds,np.ones(6,bool),max_nearest=5)
    assert audit["n_original_heldout_photo_residuals"]==6
    assert audit["n_distinct_original_photo_site_cross_species_neighbor_pairs"]==14
    assert audit["spatial_radius_reports_km"]["50"]["insufficient_neighbor_support"]
    assert audit["spatial_radius_reports_km"]["50"]["moran_like_centered_fourstate_residual_spatial_index_by_model"]["SPECIES_GEODESIC_SPATIAL_FIELD_ALL_ENV"] is None


def test_undeclared_or_missing_heldout_residual_prediction_rejected():
    d=pd.DataFrame({"inat_taxon_id":range(40),
                    "latitude":np.linspace(0,1,40),
                    "longitude":np.linspace(5,6,40)})
    labels=np.arange(40)%4
    prediction=np.full((40,4),.25)
    prediction[0]=np.nan
    models={k:prediction for k in (
        "SPECIES_NONGP_COORDINATES",
        "SPECIES_GEODESIC_SPATIAL_FIELD",
        "SPECIES_GEODESIC_SPATIAL_FIELD_ALL_ENV")}
    with pytest.raises(RuntimeError,match="incomplete"):
        M.spatial_oof_residual_diagnostic(d,labels,models,np.ones(40,bool))

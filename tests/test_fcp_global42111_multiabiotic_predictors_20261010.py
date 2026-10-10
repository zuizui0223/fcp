"""Mock all requested FCP bioclimate, altitude, radiation and soil blocks."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"/"analysis"))
import compare_fcp_global42111_full_climate_sun_soil_blocks_20261010 as M


@pytest.fixture
def source(monkeypatch):
    monkeypatch.setattr(M,"N_ORIGINAL",240)
    monkeypatch.setattr(M,"N_CLASSIFIED",220)
    monkeypatch.setattr(M,"MIN_COMPLETE",100)
    rng=np.random.default_rng(20261010)
    i=np.arange(240)
    d=pd.DataFrame({
        "inat_taxon_id":i+1,
        "species":[f"Genus{j%50} plant{j}" for j in i],
        "morph":np.where(i<220,np.array(M.CLASSES)[i%4],"UNKNOWN"),
        "measurement_status":np.where(i<220,"classified_four_state_morph","UNCLASSIFIABLE"),
        "site_geo_status":"VALID_PUBLIC_ORIGINAL_PHOTO_POINT",
        "latitude":rng.uniform(-80,80,240),
        "longitude":rng.uniform(-179,179,240),
    })
    for f in M.ALL[3:]:
        d[f]=rng.normal(loc=5,scale=2,size=240)
    d.loc[230:,"wc_srad_annual_kj_m2_day"]=np.nan
    d.loc[200:219,"soil_pH"]=np.nan
    return d


def test_original_source_denominator_and_classification_missingness(source):
    sub,c=M.make_source(source,strict=True)
    assert c["n_original_species_taxa"]==240
    assert c["n_original_classifiable_photo_species"]==220
    assert c["n_original_unclassifiable_photo_species"]==20
    assert c["n_complete_all_environment_and_original_colour"]==200
    assert len(sub)==200
    assert sub.source_cell_162.between(0,161).all()


def test_all_radiation_thermal_precip_soil_elevation_blocks_are_declared():
    families=M.features()
    assert "FULL_ALL_BLOCKS" in families
    assert len(families)==2+len(M.BLOCKS)+sum(map(len,M.BLOCKS.values()))
    assert tuple(M.BLOCKS)==("elevation","temperature","precipitation","solar_radiation","wind","vapor_pressure","soil")
    for v in M.BLOCKS.values():
        for name in v:
            assert "FULL_MINUS_SINGLE_"+name.upper() in families
    assert "wc_srad_annual_kj_m2_day" in M.ALL and "wc_srad_monthly_cv" in M.ALL


def test_identical_species_grouped_folds_under_genus_and_geo(source):
    d,_=M.make_source(source,strict=True)
    for group in ("genus","source_cell_162"):
        out=M.oof(d,group,nboot=29)
        assert out["all_models_same_original_photo_species_and_fixed_folds"]
        assert out["n_original_species_evaluated"]==200
        assert len(out["source_model_scores"])==len(M.features())
        assert all(q["n_original_photo_taxa"]==200 for q in out["source_model_scores"].values())
        assert set(out["full_vs_geography_and_drop_one_predictor_gains"])==(
            {"full_minus_geography_only"}
            | {"conditional_"+x for x in M.BLOCKS}
            | {"conditional_single_"+x for z in M.BLOCKS.values() for x in z}
        )
        assert all(len(g["group_resampled_fixed_oof_95CI"])==2
                   for g in out["full_vs_geography_and_drop_one_predictor_gains"].values())


def test_sun_missingness_cannot_be_replaced_with_original_cell_centroid(source):
    source.loc[0,"wc_srad_monthly_cv"]=np.nan
    d,c=M.make_source(source)
    assert 1 not in set(d.inat_taxon_id)
    assert c["n_original_classifiable_photo_species"]==220
    source.loc[1,"site_geo_status"]="NO_PUBLIC_PHOTO_SITE"
    d,_=M.make_source(source)
    assert 2 not in set(d.inat_taxon_id)


def test_duplicate_source_taxa_fails_and_named_field_missing_fails(source):
    original=source.copy()
    original.loc[1,"inat_taxon_id"]=1
    with pytest.raises(ValueError,match="identities"):
        M.make_source(original)
    bad=source.drop(columns=["wc_bio14"])
    with pytest.raises(ValueError,match="missing"):
        M.make_source(bad)

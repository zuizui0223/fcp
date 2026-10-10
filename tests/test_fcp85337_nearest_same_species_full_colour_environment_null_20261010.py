"""Original photo nearest same-species 4-colour pair sign-null safety and FWER."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"/"analysis"))
import audit_fcp85337_nearest_same_species_full_colour_environment_null_20261010 as M


@pytest.fixture
def original():
    rng=np.random.default_rng(20261010)
    rows=[]
    for taxon in range(1,26):
        for site in range(3):
            # Each species has three old source photo regions. Site 1 is
            # closest to 0; site 2 is much farther. No colour-based selection.
            latitude=taxon*.11
            longitude=(taxon%5)*14.+ (0.04 if site==1 else 4. if site==2 else 0.)
            row={
                "inat_taxon_id":taxon,
                "observation_id":taxon*10+site,
                "photo_id":taxon*100+site,
                "cell_id":50+site,
                "species":f"Genus{taxon%5} species{taxon}",
                "morph":M.CLASSES[(taxon+site)%4],
                "latitude":latitude,
                "longitude":longitude,
                "wc_elevation_m":float(site*100.+taxon)
            }
            for variable in M.FEATURES_SOIL:
                row[variable]=float(taxon+rng.normal()+site)
            rows.append(row)
    return pd.DataFrame(rows)


def test_outcome_independent_closest_photo_pair_and_unique_species(original):
    d=M.deterministic_closest_original_pairs(original)
    assert len(d)==25 and d.inat_taxon_id.nunique()==25
    assert d.photo_id_a.mod(100).eq(0).all()
    assert d.photo_id_b.mod(100).eq(1).all()
    assert d.distance_km.lt(6).all()
    other=original.copy()
    other["morph"]="white"
    z=M.deterministic_closest_original_pairs(other)
    assert d[["inat_taxon_id","photo_id_a","photo_id_b","distance_km"]].equals(
        z[["inat_taxon_id","photo_id_a","photo_id_b","distance_km"]])


def test_original_source_cell_and_photo_reuse_fails_closed(original):
    bad=original.copy()
    bad.loc[1,"cell_id"]=50
    with pytest.raises(ValueError,match="identity"):
        M.deterministic_closest_original_pairs(bad)
    repeat=original.copy()
    repeat.loc[2,"photo_id"]=100
    with pytest.raises(ValueError,match="identity"):
        M.deterministic_closest_original_pairs(repeat)


def test_geo_residualization_does_not_use_original_morph(original):
    a=M.deterministic_closest_original_pairs(original)
    x=M.geography_residualized_gradients(a,M.FEATURES_CLIMATE)
    b=a.copy()
    b["morph_a"]="white"
    b["morph_b"]="red_pink"
    y=M.geography_residualized_gradients(b,M.FEATURES_CLIMATE)
    np.testing.assert_allclose(x,y)
    assert x.shape==(25,len(M.FEATURES_CLIMATE))
    assert np.isfinite(x).all()


def test_all_features_and_maxT_pvalues_not_winner_only(original):
    d=M.deterministic_closest_original_pairs(original)
    z=M.sign_null(d,M.FEATURES_CLIMATE,nperm=19,seed=321)
    assert z["n_fixed_original_source_photo_pairs"]==25
    assert z["one_photo_pair_per_species_fixed_by_nearest_actual_photo_site_distance"]
    assert len(z["environmental_feature_results"])==len(M.FEATURES_CLIMATE)
    assert "temperature" in z["all_named_environmental_block_results"]
    assert "soil" not in z["all_named_environmental_block_results"]
    for v in z["environmental_feature_results"].values():
        assert 0.05<=v["one_sided_label_swap_p_unadjusted"]<=1
        assert v["within_radius_maxT_across_all_individual_features_FWER_p"]>=v["one_sided_label_swap_p_unadjusted"]
        assert v["conservative_8_nested_cohort_radius_tests_Bonferroni_maxT_p"]>=v["within_radius_maxT_across_all_individual_features_FWER_p"]


def test_soil_not_imputed_and_test_uses_one_fixed_pair_per_original_taxon(original):
    d=M.deterministic_closest_original_pairs(original)
    d.loc[0,"d_soil_pH"]=np.nan
    with pytest.raises(ValueError,match="missing"):
        M.geography_residualized_gradients(d,M.FEATURES_SOIL)
    assert d.inat_taxon_id.is_unique


def test_only_all_four_photo_colour_categories_not_pigment_numeric(original):
    d=M.deterministic_closest_original_pairs(original)
    assert set(d.morph_a)==set(M.CLASSES)
    z=M.sign_null(d,("wc_bio12",),nperm=9,seed=1)
    source_components=z["environmental_feature_results"]["wc_bio12"]["signed_photo_colour_components"]
    assert set(source_components)==set(M.CLASSES)
    assert abs(sum(source_components.values()))<1e-7

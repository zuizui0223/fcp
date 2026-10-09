"""Permute original photographed labels only within prequalified genus-cell groups."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"/"analysis"))
import audit_fcp_global42111_local_congener_moisture_null_20261009 as N
import audit_fcp_global42111_local_congeners_20261009 as M


@pytest.fixture
def original():
    rng=np.random.default_rng(20261009)
    rows=[]
    for g in range(24):
        lon=-170+20*(g%10)
        for j in range(15):
            row={"inat_taxon_id":g*15+j+1,
                 "species":f"Genus{g} species{j}",
                 "latitude":12+rng.uniform(-.1,.1),
                 "longitude":lon+rng.uniform(-.1,.1),
                 "site_geo_status":"VALID_PUBLIC_ORIGINAL_PHOTO_POINT",
                 "measurement_status":"classified_four_state_morph",
                 "morph":M.CLASSES[(g+j//3)%4],
                 "environment_climate_complete":True,
                 "environment_soil_complete":True,
                 "wc_elevation_m":rng.uniform(100,400)}
            for v in M.TEMP+M.RAIN+M.SOIL:row[v]=rng.normal(10,3)
            rows.append(row)
    d=pd.DataFrame(rows)
    pop,_=M.photo_populations(d,strict=False)
    return pop["CLIMATE_ALL"]


def test_permutation_sample_uses_true_source_photo_group_distances(original):
    select=N.fixed_source_cohort(original,250)
    assert len(select)==360
    assert select.genus_cell_id.nunique()==24
    assert select.inat_taxon_id.nunique()==360


def test_outofsource_fold_gain_matches_existing_same_local_cell_source(original):
    photo=N.fixed_source_cohort(original,250)
    observed=M.test_local_congeners(photo,include_soil=False)
    expected=observed["incremental_source_photo_prediction"]["moisture_unique_beyond_local_geography_temperature"]["mean_heldout_brier_reduction"]
    found=N.outofsource_photo_moisture_gain(photo,M.fold_plan(photo))
    assert found==pytest.approx(expected,abs=1e-12)


def test_random_shuffling_retests_all_original_species_and_retains_group_composition(original):
    chosen=N.fixed_source_cohort(original,250)
    result=N.null_calibration(chosen,nperm=3,seed=51)
    assert result["n_original_photo_species"]==360
    assert result["n_original_genus_cell_groups"]==24
    assert result["group_colour_composition_preserved_in_every_null"]
    assert result["test_fold_assignments_preserved_in_every_null"]
    assert result["n_null_permutations"]==3
    assert 0.25<=result["one_sided_permutation_p_source_conditional"]<=1
    assert all(np.isfinite(v) for v in result["null_2p5_50_97p5"])


def test_not_enough_source_photo_groups_cannot_be_permutation_rescued(original):
    q=original.iloc[:50].copy()
    with pytest.raises(ValueError,match="Too few"):
        N.null_calibration(q)

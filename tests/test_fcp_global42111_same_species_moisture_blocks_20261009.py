"""Mock-only thermal vs rainfall blocks on fixed photographed same-species pair labels."""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"/"analysis"))
import audit_fcp_global42111_same_species_moisture_blocks_20261009 as M


@pytest.fixture
def source_pair():
    rng=np.random.default_rng(102026)
    n=850
    k=np.arange(n)
    z=pd.DataFrame({
        "inat_taxon_id":k+1,
        "pair_state":np.where(k%4==0,"discordant","same"),
        "genus":[f"Genus{j%100}" for j in k],
        "pair_midpoint_cell_162":k%42,
        "log_geodesic_distance_km":rng.normal(7,1,n),
        "delta_abs_latitude":rng.uniform(0,50,n),
        "delta_elevation_m":rng.uniform(0,2500,n),
    })
    for v in M.THERMAL+M.MOISTURE:
        z[v]=rng.uniform(0,100,n)
    return z


def test_climate_block_specification_preserves_original_four_bio_features():
    assert M.THERMAL==("delta_wc_bio1","delta_wc_bio5")
    assert M.MOISTURE==("delta_wc_bio12","delta_wc_bio15")
    assert set(M.MODELS)=={"GEOGRAPHY","GEOGRAPHY_TEMPERATURE","GEOGRAPHY_MOISTURE","GEOGRAPHY_ALL_CLIMATE"}


def test_source_original_photo_pair_climate_models_are_blocked_and_identical(source_pair):
    for group in ("genus","pair_midpoint_cell_162"):
        out=M.model_group(source_pair,group)
        assert out["status"]=="ORIGINAL_SAME_SPECIES_PHOTOGRAPHED_PAIR_CLIMATE_BLOCK_MODEL"
        assert out["same_species_pairs_and_split_in_all_models"]
        assert out["n_original_source_species_pairs"]==850
        assert set(out["models"])==set(M.MODELS)
        assert {v["source_pairs_evaluated"] for v in out["models"].values()}=={850}
        assert all(len(t["original_genus_or_cell_block_bootstrap_95CI"])==2 for t in out["climate_increments"].values())


def test_source_misclassification_not_inferred_as_null_sample(source_pair):
    with pytest.raises(ValueError,match="One source photographed pair"):
        M.model_group(pd.concat([source_pair,source_pair.iloc[[0]]]),"genus")


def test_environment_missing_not_silently_filled(source_pair):
    source_pair.loc[0,"delta_wc_bio12"]=np.nan
    with pytest.raises(ValueError,match="missing"):
        M.model_group(source_pair,"genus")


def test_unusable_holdout_cells_return_coverage_status(source_pair):
    source_pair["pair_midpoint_cell_162"]=1
    d=M.model_group(source_pair,"pair_midpoint_cell_162")
    assert d["status"]=="HOLD_PAIR_GEOGRAPHIC_OR_TAXON_SUPPORT"

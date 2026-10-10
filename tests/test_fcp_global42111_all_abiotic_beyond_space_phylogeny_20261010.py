"""True phylo+spatial seven-block family assembly and conservative decision matrix."""
from __future__ import annotations
import sys
from pathlib import Path
import pandas as pd
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"/"analysis"))
import run_fcp_global42111_all_abiotic_beyond_space_phylogeny_20261010 as M


def test_full_inventory_is_seven_environment_blocks_not_only_precipitation():
    assert set(M.BLOCKS)=={
        "temperature","precipitation","solar_radiation","wind","vapor_pressure",
        "soil","elevation"
    }
    assert len(M.BLOCKS["temperature"])==11
    assert len(M.BLOCKS["precipitation"])==8
    assert len(M.BLOCKS["soil"])==10
    assert "wc_srad_annual_mean" in M.BLOCKS["solar_radiation"]


def fake_model(blocks, n=300, unexpected=False):
    good=["GEO_SPATIAL_ONLY","GEO_SPATIAL_PHYLOGENY","GEO_SPATIAL_PHYLOGENY_ALL_ENV"]
    for key in blocks:
        good.append("GEO_SPATIAL_PHYLOGENY_MINUS_"+key.upper())
    gains={
        "all_abiotic_beyond_spatial_plus_real_phylogeny":{
            "heldout_brier_reduction":.02,"original_group_cluster_95CI":[.001,.03]
        }
    }
    for key in blocks:
        gains["unique_"+key+"_beyond_space_real_phylogeny_and_other_abiotic"]={
            "heldout_brier_reduction":.0004,
            "original_group_cluster_95CI":[-.0001,.001] if key=="soil" else [.0001,.001],
        }
    return {
        "status":"EXPLORATORY_EXACT_PHYLO_BROWNIAN_KERNEL_PLUS_SPATIAL_KERNEL_ENVIRONMENT",
        "source_environment_blocks":list(blocks)+(["unexpected"] if unexpected else []),
        "modelled_true_phylogeny_BM_shared_branch_covariance":True,
        "modelled_original_photograph_spatial_exponential_covariance":True,
        "source_model_scores":{k:{"n_same_species":n} for k in good},
        "conditional_feature_gains":gains,
    }


def fake_raw():
    output={
        "original_42111_photo_species_denominator":42111,
        "original_one_photo_classifiable":18457,
        "full_872_1761_direct_tip_phylogenetic_coverage_HOLD":True,
        "cohorts":{},
    }
    for cap,n in M.EXPECTED_N.items():
        for mode in M.MODES:
            blocks={k:v for k,v in M.BLOCKS.items() if mode=="soil" or k!="soil"}
            output["cohorts"][cap+"_"+mode]={
                "source":{
                    "n_direct_original_photo_taxa_pre_missingness":n,
                    "n_direct_original_photo_taxa_full_environment_common_case":n,
                },
                "models":{
                    scale+"__"+hold:fake_model(blocks,n=n)
                    for scale in M.SCALES for hold in M.HOLDS
                },
            }
    return output


def test_same_support_all_seven_blocks_across_crossed_heldouts():
    z=M.summary(fake_raw())
    assert "250km_soil__soil" in z
    assert "500km_climate__temperature" in z
    assert "250km_climate__soil" not in z
    assert z["250km_soil__soil"]["all_four_conditional_group_intervals_positive"] is False
    assert z["500km_climate__solar_radiation"]["all_four_conditional_group_intervals_positive"] is True
    assert len(z["500km_soil__precipitation"]["fixed_repeated_spatial_and_heldout_results"])==4


def test_one_bad_original_geographic_or_genus_holdout_blocks_claim():
    raw=fake_raw()
    q=raw["cohorts"]["500km_soil"]["models"]["250km__true_photo_cell"]
    key="unique_precipitation_beyond_space_real_phylogeny_and_other_abiotic"
    q["conditional_feature_gains"][key]["original_group_cluster_95CI"]=[-.01,.001]
    z=M.summary(raw)
    assert not z["500km_soil__precipitation"]["all_four_conditional_group_intervals_positive"]


def test_must_keep_genuine_tree_and_spatial_covariance_and_full_42111_denominator():
    raw=fake_raw()
    raw["cohorts"]["250km_climate"]["models"]["50km__genus"][
        "modelled_true_phylogeny_BM_shared_branch_covariance"]=False
    with pytest.raises(RuntimeError,match="mandatory covariance"):
        M.summary(raw)
    raw=fake_raw()
    raw["original_42111_photo_species_denominator"]=42110
    with pytest.raises(RuntimeError,match="denominator"):
        M.summary(raw)


def test_unknown_extra_environment_not_pass_as_complete_seven_block():
    raw=fake_raw()
    raw["cohorts"]["250km_climate"]["models"]["50km__genus"]["source_environment_blocks"].append("unknown")
    with pytest.raises(RuntimeError,match="omitted"):
        M.summary(raw)

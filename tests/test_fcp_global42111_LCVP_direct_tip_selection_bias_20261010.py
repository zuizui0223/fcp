"""Synthetic exact LCVP tip selection-bias diagnostic, no source photo remeasure."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"/"analysis"))
import audit_fcp_global42111_LCVP_direct_tip_selection_bias_20261010 as M


@pytest.fixture
def cases(monkeypatch):
    monkeypatch.setattr(M,"SAMPLE",{"250":(8,3),"500":(10,4)})
    taxa=np.arange(1,11)
    s=pd.DataFrame({
        "inat_taxon_id":taxa,
        "genus":["G","G","G","G","H","H","H","H","K","K"],
        "photo_cell_162":[9]*8+[10]*2,
        "morph":["white","white","red_pink","blue_purple","yellow_orange",
                 "white","red_pink","blue_purple","white","white"],
    })
    for k in M.FEATURES:
        s[k]=taxa.astype(float)
    direct=[True,True,False,False,True,False,False,False,True,False]
    led=pd.DataFrame({
        "inat_taxon_id":taxa,
        "direct_lcvp_backbone_tip":direct,
        "in_fixed_250km_original_source":[True]*8+[False]*2,
    })
    return s,led


def test_exact_source_missing_tree_tips_do_not_disappear(cases):
    s,led=cases
    a=M.cover(s,led,500,strict=True)
    assert a["n_original_source_photo_species"]==10
    assert a["n_direct_original_LCVP_tip_species"]==4
    assert a["n_not_direct_original_LCVP_tip_species"]==6
    assert a["n_orig_genera_with_both_direct_and_untipped_photo_species"]==3
    assert sum(a["direct_tip_strata"]["DIRECT_LCVP"]["source_photo_four_colour_counts"].values())==4
    assert sum(a["direct_tip_strata"]["NOT_DIRECT_LCVP"]["source_photo_four_colour_counts"].values())==6
    assert a["direct_backbone_tip_selection_did_not_use_original_photo_colour_results"]


def test_original_colour_status_changes_outcome_diagnostics_not_taxon_selection(cases):
    s,led=cases
    a=M.cover(s,led,500,strict=True)
    other=s.copy()
    other["morph"]="white"
    b=M.cover(other,led,500,strict=True)
    assert a["n_direct_original_LCVP_tip_species"]==b["n_direct_original_LCVP_tip_species"]
    assert a["original_photo_site_environment_standardized_difference_direct_minus_untipped"]==b["original_photo_site_environment_standardized_difference_direct_minus_untipped"]
    assert a["max_abs_original_colour_fraction_difference"]>=b["max_abs_original_colour_fraction_difference"]


def test_250_nested_source_only_keeps_original_taxa(cases):
    s,led=cases
    small=M.cover(s.iloc[:8],led,250,strict=True)
    assert small["n_original_source_photo_species"]==8
    assert small["n_direct_original_LCVP_tip_species"]==3
    assert small["original_direct_tree_tip_fraction"]==pytest.approx(3/8)


def test_missing_ledger_photo_species_is_not_mapped_to_its_genus(cases):
    s,led=cases
    bad=led.loc[led.inat_taxon_id!=3]
    with pytest.raises(RuntimeError,match="missing taxonomy"):
        M.cover(s,bad,500,strict=True)


def test_id_reuse_fails_and_climate_missing_fails(cases):
    s,led=cases
    led.loc[1,"inat_taxon_id"]=1
    with pytest.raises(ValueError,match="duplicated"):
        M.cover(s,led,500,strict=True)
    s.loc[1,"wc_bio12"]=np.nan
    with pytest.raises(RuntimeError,match="missing covariates"):
        M.cover(s,cases[1].assign(inat_taxon_id=np.arange(1,11)),500,strict=True)

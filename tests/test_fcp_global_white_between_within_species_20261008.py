"""Global white~50% could mask species-specific all-white/all-coloured sampling."""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT=Path(__file__).resolve().parents[1]/"scripts/analysis"
sys.path.insert(0,str(ROOT))
sp=importlib.util.spec_from_file_location("species_white",ROOT/"audit_fcp_global_white_between_within_species_20261008.py")
m=importlib.util.module_from_spec(sp)
sp.loader.exec_module(m)


def demo_source():
    # Ten distinct cells in one latitude region, two species each 5
    # classified photo-colour states. Global white=5/10, but every
    # source species itself is PURE white or coloured. No 50:50
    # within-species genetic/pigment equilibrium is implied.
    return pd.DataFrame([
        {"inat_taxon_id":1,"cell_id":i,"morph":"white",
         "measurement_status":"classified_four_state_morph"} for i in range(72,77)
    ]+[
        {"inat_taxon_id":2,"cell_id":i,"morph":"blue_purple",
         "measurement_status":"classified_four_state_morph"} for i in range(78,83)
    ])


def test_exact_half_photos_with_zero_mixed_species():
    s=m.compress_species_photo_records(demo_source())
    g=s.loc[s.region=="low_0_30"]
    d=m.decompose(g,123)
    assert d["n_species"]==2
    assert d["n_classifiable_photo_taxon_cell_records"]==10
    assert d["fraction_white_classified_photo_weighted"]==0.5
    assert d["n_species_with_both_white_and_nonwhite_in_sampled_cells"]==0
    assert d["n_species_with_white_only_in_sampled_cells"]==1
    assert d["n_species_with_nonwhite_only_in_sampled_cells"]==1
    assert d["photo_weighted_between_species_share_of_observed_photo_variance"]==1
    assert d["species_equal_between_species_share_of_observed_photo_variance"]==1


def test_all_white_and_mixed_species_differ_but_identity_remains_exact():
    x=demo_source()
    x.loc[(x.inat_taxon_id==1)&(x.cell_id==72),"morph"]="red_pink"
    s=m.compress_species_photo_records(x)
    d=m.decompose(s.loc[s.region=="low_0_30"],124)
    assert d["n_species_with_both_white_and_nonwhite_in_sampled_cells"]==1
    assert 0<d["fraction_species_with_both_observed_states"]<1
    assert 0<d["photo_weighted_between_species_share_of_observed_photo_variance"]<1
    assert 0<d["species_equal_between_species_share_of_observed_photo_variance"]<1


def test_missing_photo_is_in_denominator_not_called_coloured():
    x=demo_source()
    x.loc[0,"morph"]="mixed_uncertain"
    x.loc[0,"measurement_status"]="not_evaluable_roi_or_flip_gate"
    s=m.compress_species_photo_records(x)
    assert int(s.loc[(s.region=="ALL_GLOBAL"),"observed_cells"].sum())==10
    assert int(s.loc[s.region=="ALL_GLOBAL","n_classified"].sum())==9
    assert int(s.loc[s.region=="ALL_GLOBAL","n_white"].sum())==4
    assert int(s.loc[s.region=="low_0_30","n_white"].sum())==4


def test_separate_species_region_and_global_record_grains():
    s=m.compress_species_photo_records(demo_source())
    assert len(s.loc[s.region=="ALL_GLOBAL"])==2
    assert len(s.loc[s.region=="low_0_30"])==2
    assert s.groupby("region").observed_cells.sum()["ALL_GLOBAL"]==10
    assert s.groupby("region").observed_cells.sum()["low_0_30"]==10


def test_zero_differing_observed_labels_not_called_perfect_morph_balance():
    d=pd.DataFrame([{"n_classified":5,"n_white":5},
                    {"n_classified":5,"n_white":5}])
    r=m.decompose(d,125)
    assert r["n_species_with_both_white_and_nonwhite_in_sampled_cells"]==0
    assert r["photo_weighted_between_species_share_of_observed_photo_variance"] is None


def test_duplicate_taxon_cell_raises_before_any_frequency():
    d=pd.concat([demo_source(),demo_source().iloc[[0]]],ignore_index=True)
    with pytest.raises(RuntimeError,match="Repeating same photo"):
        m.compress_species_photo_records(d)


def test_taxon_cell_photo_imbalance_cannot_be_species_equal_by_accident():
    d=pd.DataFrame([{"n_classified":9,"n_white":0},
                    {"n_classified":1,"n_white":1}])
    r=m.decompose(d,128)
    assert r["fraction_white_classified_photo_weighted"]==0.1
    assert r["fraction_white_species_equal"]==0.5

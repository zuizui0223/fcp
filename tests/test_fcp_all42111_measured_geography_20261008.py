"""Check source-frozen global FCP measured palette and independent geography denominators."""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

SRC=Path(__file__).resolve().parents[1]/"scripts/analysis/synthesize_fcp_all42111_measured_geography_20261008.py"
sp=importlib.util.spec_from_file_location("fcp_global_measured",SRC)
m=importlib.util.module_from_spec(sp)
sp.loader.exec_module(m)


def demo_cells():
    return pd.DataFrame([
        {"inat_taxon_id":1,"cell_id":0,"morph":"white","measurement_status":"classified_four_state_morph"},
        {"inat_taxon_id":1,"cell_id":18,"morph":"red_pink","measurement_status":"classified_four_state_morph"},
        {"inat_taxon_id":2,"cell_id":18,"morph":"mixed_uncertain","measurement_status":"not_evaluable_roi_or_flip_gate"},
        {"inat_taxon_id":3,"cell_id":90,"morph":"blue_purple","measurement_status":"classified_four_state_morph"},
        {"inat_taxon_id":4,"cell_id":90,"morph":"yellow_orange","measurement_status":"classified_four_state_morph"},
    ])


def demo_pairs():
    return pd.DataFrame([
        {"inat_taxon_id":1,"cell_id_1":0,"cell_id_2":18,"both_endpoints_classifiable":True,
         "pair_state":"discordant","observer_id_1":"u1","observer_id_2":"u2"},
        {"inat_taxon_id":2,"cell_id_1":18,"cell_id_2":90,"both_endpoints_classifiable":False,
         "pair_state":"unclassifiable","observer_id_1":"u3","observer_id_2":"u4"},
        {"inat_taxon_id":3,"cell_id_1":36,"cell_id_2":126,"both_endpoints_classifiable":True,
         "pair_state":"same","observer_id_1":"u3","observer_id_2":"u8"},
    ])


def test_original_equal_area_cell_midpoints_and_climate_zones():
    assert len(m.row_cell_latitude(np.arange(162)))==162
    lat=m.row_cell_latitude(np.array([0,80,81,161]))
    assert lat[0]<-60 and lat[3]>60
    assert m.classify_region(lat).tolist()==["high_60_90","low_0_30","low_0_30","high_60_90"]
    assert abs(m.row_cell_lon(np.array([0]))[0]+170)<.0001


def test_regional_photo_diversity_does_not_hide_classifier_failures():
    region,cells=m.region_cell_stats(demo_cells())
    assert sum(x["taxon_cell_opportunities"] for x in region)==5
    assert sum(x["classified_taxon_cell"] for x in region)==4
    assert sum(x["unclassified_taxon_cell"] for x in region)==1
    assert cells.n_taxa.sum()==5
    assert all("population_claim_boundary" in x for x in region)


def test_pair_conditional_observed_discordance_and_missingness_bounds(monkeypatch):
    monkeypatch.setattr(m,"BOOTSTRAPS",49)
    zone,check=m.pair_metrics(demo_pairs())
    total=next(x for x in zone if x["region"]=="ALL")
    assert total["fixed_observer_disjoint_species_pairs"]==3
    assert total["both_photos_classified"]==2
    assert total["discordant_photo_colour"]==1
    assert total["one_or_both_unclassifiable"]==1
    assert total["discordance_fraction_among_both_classified"]==0.5
    assert np.allclose(total["no_assumption_full_fixed_pair_discordance_bounds"],[1/3,2/3])
    assert check["n_fixed_pairs"]==3


def test_missing_label_cannot_be_called_same_colour():
    x=demo_pairs()
    x.loc[1,"pair_state"]="same"
    with pytest.raises(RuntimeError,match="Pair state"):
        m.pair_metrics(x)


def test_pair_must_have_distinct_observers_and_cells():
    x=demo_pairs()
    x.loc[0,"observer_id_2"]="u1"
    with pytest.raises(RuntimeError,match="Same-observer"):
        m.pair_metrics(x)
    x=demo_pairs()
    x.loc[0,"cell_id_2"]=x.loc[0,"cell_id_1"]
    with pytest.raises(RuntimeError,match="Same-cell"):
        m.pair_metrics(x)


def test_taxon_cell_duplicate_prohibited():
    x=pd.concat([demo_cells(),demo_cells().iloc[[0]]],ignore_index=True)
    with pytest.raises(RuntimeError,match="more than once"):
        m.region_cell_stats(x)


def test_no_favourable_colour_reuse_and_no_fake_selection():
    source=SRC.read_text()
    assert m.N_WORLD_SPECIES==42111
    assert m.N_WORLD_CELL_ROWS==85337
    assert m.N_CROSSCELL_PAIRS==13416
    assert "photo_url_large" not in source
    assert "image_gen" not in source
    assert "benefit_cost_selection_coefficient" not in source

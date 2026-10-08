"""White ~50% in a world photo atlas need not indicate balancing morph selection."""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT=Path(__file__).resolve().parents[1]/"scripts/analysis"
sys.path.insert(0,str(ROOT))
sp=importlib.util.spec_from_file_location("white_balance",ROOT/"audit_fcp_global_white_balance_20261008.py")
m=importlib.util.module_from_spec(sp)
sp.loader.exec_module(m)


def source():
    # Species 1 is white in the low band and red in middle; species 2 is
    # entirely white across both bands. Species 3 has one unclassified photo.
    return pd.DataFrame([
        {"inat_taxon_id":1,"cell_id":75,"morph":"white","measurement_status":"classified_four_state_morph"},
        {"inat_taxon_id":1,"cell_id":77,"morph":"white","measurement_status":"classified_four_state_morph"},
        {"inat_taxon_id":1,"cell_id":20,"morph":"red_pink","measurement_status":"classified_four_state_morph"},
        {"inat_taxon_id":2,"cell_id":76,"morph":"white","measurement_status":"classified_four_state_morph"},
        {"inat_taxon_id":2,"cell_id":21,"morph":"white","measurement_status":"classified_four_state_morph"},
        {"inat_taxon_id":3,"cell_id":75,"morph":"mixed_uncertain","measurement_status":"not_evaluable_roi_or_flip_gate"},
    ])


def test_unequal_geographical_repeats_do_not_get_equal_species_weight_by_accident():
    s=m.source_species_regions(source())
    low=s.loc[s.region=="low_0_30"]
    assert len(low)==3
    valid=low.loc[low.n_classified>0]
    assert int(valid.n_photo_white.sum())==3
    assert len(valid)==2
    assert float(valid.fraction_white_when_classified.mean())==1


def test_half_observed_is_not_half_all_source_photo_outcomes():
    d=pd.DataFrame([
        {"inat_taxon_id":1,"cell_id":77,"morph":"white","measurement_status":"classified_four_state_morph"},
        {"inat_taxon_id":2,"cell_id":78,"morph":"blue_purple","measurement_status":"classified_four_state_morph"},
        {"inat_taxon_id":3,"cell_id":79,"morph":"mixed_uncertain","measurement_status":"not_evaluable_ambiguous_palette_composition"},
        {"inat_taxon_id":4,"cell_id":80,"morph":"mixed_uncertain","measurement_status":"not_evaluable_roi_or_flip_gate"}
    ])
    s=m.source_species_regions(d)
    x=m.region_coverage(s)
    y=next(v for v in x if v["region"]=="low_0_30")
    assert y["fraction_white_all_classified_cell_records"]==0.5
    assert y["NO_assumption_all_cell_white_frequency_bound"]==[0.25,0.75]
    assert y["n_unclassifiable_cell"]==2
    assert y["species_equal_mean_white_fraction_within_region"]==0.5


def test_paired_species_control_reveals_colour_change_that_cross_species_mix_cannot():
    d=pd.DataFrame([
        {"region":"low_0_30","inat_taxon_id":1,"n_cells":2,
         "n_classified":2,"n_photo_white":2,"fraction_white_when_classified":1.0,
         "white_and_chromatic_images_in_same_region":False},
        {"region":"middle_30_60","inat_taxon_id":1,"n_cells":3,
         "n_classified":3,"n_photo_white":0,"fraction_white_when_classified":0.0,
         "white_and_chromatic_images_in_same_region":False},
        {"region":"low_0_30","inat_taxon_id":2,"n_cells":1,
         "n_classified":1,"n_photo_white":0,"fraction_white_when_classified":0.0,
         "white_and_chromatic_images_in_same_region":False},
        {"region":"middle_30_60","inat_taxon_id":3,"n_cells":1,
         "n_classified":1,"n_photo_white":1,"fraction_white_when_classified":1.0,
         "white_and_chromatic_images_in_same_region":False}
    ])
    result,pair=m.species_pair_region(d,"low_0_30","middle_30_60",1)
    assert result["n_matched_species"]==1
    assert result["species_equal_white_difference_region_b_minus_region_a"]==-1
    assert result["status"]=="HOLD_SPARSE_SAME_SPECIES_REGIONAL_PHOTOS"
    assert int(pair.inat_taxon_id.iloc[0])==1
    x,p=m.species_pair_region(d,"low_0_30","middle_30_60",3)
    assert x["n_matched_species"]==0


def test_mutually_exclusive_region_bins_and_source_date_not_used_as_fitness():
    assert m.original.classify_region(
        m.original.row_cell_latitude(np.arange(162))).shape==(162,)
    src=Path(m.__file__).read_text()
    assert "relative_fitness" not in src
    assert "anthocyanin_loss_genotype" not in src


def test_original_taxon_cell_identifiers_must_be_unique():
    x=source()
    x=pd.concat([x,x.iloc[[0]]],ignore_index=True)
    with pytest.raises(RuntimeError,match="Repeated exact species"):
        m.source_species_regions(x)


def test_region_photo_classification_is_not_imputed_from_missingness():
    s=m.source_species_regions(source())
    x=s.loc[s.inat_taxon_id==3].iloc[0]
    assert x.n_classified==0 and x.n_photo_white==0
    assert np.isnan(x.fraction_white_when_classified)
    low=m.region_coverage(s)[0]
    assert low["n_unclassifiable_cell"]>=1
    assert low["NO_assumption_all_cell_white_frequency_bound"][1]>low["NO_assumption_all_cell_white_frequency_bound"][0]

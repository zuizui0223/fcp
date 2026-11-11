"""Fixed-five-photo white/nonwhite detection, no evolutionary equilibrium assumptions."""
from pathlib import Path
import importlib.util
import sys

import numpy as np
import pandas as pd
import pytest

ROOT=Path(__file__).resolve().parents[1]/"scripts/analysis"
sys.path.insert(0,str(ROOT))
sp=importlib.util.spec_from_file_location("white_five",ROOT/"audit_fcp_global_white_fixed_photo_depth_20261008.py")
m=importlib.util.module_from_spec(sp)
sp.loader.exec_module(m)


def fake_species_region():
    return pd.DataFrame([
       {"region":"ALL_GLOBAL","inat_taxon_id":1,"n_classified":10,"n_white":5,"observed_cells":10},
       {"region":"ALL_GLOBAL","inat_taxon_id":2,"n_classified":5,"n_white":5,"observed_cells":6},
       {"region":"ALL_GLOBAL","inat_taxon_id":3,"n_classified":5,"n_white":0,"observed_cells":5},
       {"region":"low_0_30","inat_taxon_id":1,"n_classified":10,"n_white":5,"observed_cells":10},
       {"region":"low_0_30","inat_taxon_id":2,"n_classified":5,"n_white":5,"observed_cells":6},
       {"region":"low_0_30","inat_taxon_id":3,"n_classified":5,"n_white":0,"observed_cells":5},
    ])


def test_exact_five_photo_subsample_probability_preserves_pure_colour():
    n=np.array([5,5,10,10])
    w=np.array([0,5,0,10])
    assert m.mixed_prob_fixed_sample(n,w,5).tolist()==[0,0,0,0]


def test_five_of_ten_balanced_white_can_almost_always_detect_both():
    p=m.mixed_prob_fixed_sample(np.array([10,10]),np.array([5,1]),5)
    assert np.isclose(p[0],1-2/252)
    assert np.isclose(p[1],.5)


def test_all_observed_five_is_not_50_percent_global_maintenance():
    result=m.analyze(fake_species_region())
    w=next(x for x in result["results"] if x["region"]=="ALL_GLOBAL" and x["minimum_classifiable_cell_photos"]==5)
    assert w["n_species"]==3
    assert w["exact_number_with_both_observed_states"]==1
    assert np.isclose(w["observed_fraction_species_with_both_colours_at_all_available_depth"],1/3)
    assert 0<w["expected_fraction_species_showing_both_colours_in_exactly_5_of_their_existing_photos"]<1/3
    assert w["fraction_species_with_all_observed_labels_white"]==1/3
    assert w["fraction_species_with_all_observed_labels_nonwhite"]==1/3
    assert result["confirmatory_decisions_changed"] is False


def test_invalid_photo_counts_stop_instead_of_infer_color():
    with pytest.raises(ValueError,match="Invalid"):
        m.mixed_prob_fixed_sample(np.array([10]),np.array([11]),5)
    with pytest.raises(ValueError,match="requires >=2"):
        m.finite_photo_partition(np.array([1]),np.array([0]))


def test_independent_photo_correction_reduces_spurious_between_species_share():
    # Even if true underlying white probability is equal, small photo n
    # makes observed per-species colour proportions spuriously different.
    n=np.array([5,5,5,5])
    w=np.array([2,3,1,4])
    c=m.finite_photo_partition(n,w)
    assert c["estimable"]
    assert c["finite_photo_corrected_between_species_share_photo_weighted"]<c["uncorrected_between_species_share_photo_weighted"]
    assert c["finite_photo_corrected_between_species_share_species_equal"]<c["uncorrected_between_species_share_species_equal"]


def test_no_fake_improvement_from_missingness_or_second_species_identity():
    df=fake_species_region()
    df.loc[0,"observed_cells"]=15
    out=m.analyze(df)
    x=next(z for z in out["results"] if z["region"]=="ALL_GLOBAL" and z["minimum_classifiable_cell_photos"]==5)
    assert x["n_species"]==3
    assert x["n_classifiable_original_taxon_cell_photos"]==20
    df=pd.concat([df,df.iloc[[0]]],ignore_index=True)
    with pytest.raises(ValueError,match="counted twice"):
        m.analyze(df)

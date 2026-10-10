"""Finite photo-pair support falsification guards; no mock empirical outcome."""
import importlib.util
from pathlib import Path
import numpy as np
import pandas as pd

path=Path(__file__).resolve().parents[1]/"scripts/analysis/audit_polymorphism_50km_pair_support_ladder_20261010.py"
spec=importlib.util.spec_from_file_location("fcp_pairs",path)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

def species_frame():
    n=60
    return pd.DataFrame({"inat_taxon_id":[13]*n,"species":["X y"]*n,
       "photo_id":np.arange(n),"latitude":[0.]*n,
       "longitude":[0.]*30+[2.]*30,
       "morph":["white"]*30+["blue_purple"]*30})

def test_local_photo_graph_uses_actual_edges_and_no_weight_doubling():
    s=species_frame()
    share=m.maximum_single_photo_edge_share(s)
    # exactly 2 x C(30,2) local edges, 29 adjacent edges at every node
    assert np.isclose(share,29./(2*(30*29//2)),atol=1e-12)

def test_frozen_50km_original_statistic_and_null_agree_on_synthetic():
    s=species_frame()
    row,nul=m.base.analyse_species(s,"discovery")
    assert row["local_pairs_50km"]==2*(30*29//2)
    assert row["depletion_50km"]>0
    assert len(nul[50.])==m.base.PERMUTATIONS

def test_threshold_ladder_preserves_subset_and_reports_hold():
    rows=pd.DataFrame({"local_pairs_50km":[30,49,100,200],
        "depletion_50km":[.10,.10,.10,.10],
        "D_pair":[.5]*4,"D_local_50km":[.4]*4,
        "max_single_photo_edge_share":[.15]*4})
    nulls=[np.zeros(m.base.PERMUTATIONS) for _ in range(4)]
    assert m.summarize_threshold(rows,nulls,30,"discovery")["n_species"]==4
    high=m.summarize_threshold(rows,nulls,100,"discovery")
    assert high["n_species"]==2
    assert high["status"]=="HOLD_INSUFFICIENT_SPECIES"
    assert high["supported_descriptive"] is False
    assert m.summarize_threshold(rows,nulls,500,"discovery")["n_species"]==0

def test_threshold_estimable_gate_and_exact_species_equality():
    n=35
    rows=pd.DataFrame({"local_pairs_50km":[100]*n,
        "depletion_50km":[.10]*n,
        "D_pair":[.5]*n,"D_local_50km":[.4]*n,
        "max_single_photo_edge_share":[.08]*n})
    nulls=[np.zeros(m.base.PERMUTATIONS) for _ in range(n)]
    x=m.summarize_threshold(rows,nulls,100,"validation")
    assert x["n_species"]==35
    assert x["permutation_p_upper"]==1/(m.base.PERMUTATIONS+1)
    assert x["supported_descriptive"] is True
    assert x["species_bootstrap_ci95"]==[.1,.1]

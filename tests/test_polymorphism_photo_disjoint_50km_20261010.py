"""Source-shape synthetic guards for geographically disjoint photographic pairs."""
import importlib.util
from pathlib import Path
import numpy as np
import pandas as pd

script=Path(__file__).resolve().parents[1]/"scripts/analysis/audit_polymorphism_photo_disjoint_50km_20261010.py"
spec=importlib.util.spec_from_file_location("fcp_photo_disjoint",script)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

def two_cluster_photos():
    n=80
    return pd.DataFrame({
      "inat_taxon_id":[47]*n,"species":["X flora"]*n,
      "photo_id":np.arange(n).astype(str),
      "latitude":[0.]*n,
      "longitude":[0.]*40+[2.]*40,
      "observer_id":[str(i//2) for i in range(n)],
      "morph":["white"]*40+["blue_purple"]*40,
    })

def test_photo_disjoint_maximal_pairs_use_each_image_only_once():
    g=two_cluster_photos()
    for method in m.MATCH_SCHEMES:
        a,b,n_edges=m.outcome_blind_matching(g,"discovery",method,"all_photos")
        assert n_edges==2*(40*39//2)
        assert len(a)==40
        assert len(np.unique(np.r_[a,b]))==80
        assert all(a<b)
        assert not np.any(m.base.pairwise_geo_km(g.latitude,g.longitude)[a,b]>50.)

def test_different_observer_policy_never_selfpairs():
    g=two_cluster_photos()
    a,b,_=m.outcome_blind_matching(g,"third","nearest","different_observer")
    assert len(a)>=30
    assert (g.observer_id.to_numpy()[a]!=g.observer_id.to_numpy()[b]).all()
    assert len(np.unique(np.r_[a,b]))==2*len(a)

def test_selection_is_label_blind_and_reproducible():
    g=two_cluster_photos()
    original=g.copy(deep=True)
    a,b,_=m.outcome_blind_matching(g,"discovery","fixed_random","all_photos")
    g["morph"]=g.morph.iloc[::-1].to_numpy()
    aa,bb,_=m.outcome_blind_matching(g,"discovery","fixed_random","all_photos")
    assert np.array_equal(a,aa) and np.array_equal(b,bb)
    assert original.morph.ne(g.morph).any()
    pd.testing.assert_series_equal(g.photo_id,original.photo_id)

def test_segregated_original_colour_has_positive_depletion_under_disjoint_pairs():
    g=two_cluster_photos()
    row, nulls=m.species_test(g,"validation")
    assert row["source_50km_local_pairs"] if "source_50km_local_pairs" in row else True
    assert row["nearest__all_photos__disjoint_local_pairs"]==40
    assert row["nearest__all_photos__depletion"]>.4
    assert len(nulls[("nearest","all_photos")])==m.base.PERMUTATIONS
    assert abs(float(np.mean(nulls[("nearest","all_photos")])))<.10

def test_support_gate_holds_ineligible_species_without_claim():
    rows=pd.DataFrame({
      "nearest__all_photos__disjoint_local_pairs":[5]*10+[10]*5,
      "nearest__all_photos__depletion":[.20]*15,
      "nearest__all_photos__mean_local_pair_discordance":[0.]*15,
      "nearest__all_photos__unique_photos":[10]*10+[20]*5,
      "D_specieswide":[.20]*15,
    })
    nul=[np.zeros(m.base.PERMUTATIONS) for _ in range(15)]
    x=m.summarize(rows,nul,"nearest","all_photos","discovery",10)
    assert x["n_species"]==5
    assert x["status"]=="HOLD_INSUFFICIENT_SPECIES"
    assert x["positive_bounded_evidence"] is False
    assert x["species_bootstrap_ci95"] is None


def test_original_30edge_species_conditioning_preserves_original_population():
    # Source photo-pair baseline eligibility is a DISTINCT population gate
    # from >=10 photo-disjoint pairs. It must not be silently exchanged.
    n=36
    rows=pd.DataFrame({
      "nearest__all_photos__disjoint_local_pairs":[10]*n,
      "nearest__all_photos__depletion":[.12]*n,
      "nearest__all_photos__mean_local_pair_discordance":[.08]*n,
      "nearest__all_photos__unique_photos":[20]*n,
      "D_specieswide":[.20]*n,
      "source_50km_local_pair_edges":[50]*32+[20]*4,
    })
    nulls=[np.zeros(m.base.PERMUTATIONS) for _ in range(n)]
    broad=m.summarize(rows,nulls,"nearest","all_photos","third",10)
    restricted=m.summarize(rows,nulls,"nearest","all_photos","third",10,True)
    assert broad["n_species"]==36
    assert restricted["n_species"]==32
    assert restricted["restricted_to_original_30_edge_species"] is True
    assert broad["restricted_to_original_30_edge_species"] is False
    assert restricted["positive_bounded_evidence"] is True

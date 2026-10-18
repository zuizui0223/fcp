"""Unit and frozen-empirical guardrails for five-cohort within-species FCP commonality."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from scripts.analysis import fcp_five_cohort_polymorphism_commonality_20261018 as f


def test_photo_polymorphism_is_species_equal_and_dominance_is_not_white_preselected():
    rows = []
    for taxon, name, co, states in (
        (1, "Aster alpha", "discovery", ["blue_purple"]*45 + ["white"]*5),
        (2, "Bellis beta", "validation", ["yellow_orange"]*25 + ["red_pink"]*25),
    ):
        for i, morph in enumerate(states):
            rows.append({"inat_taxon_id":taxon,"species":name,"fcp_source_cohort":co,
                         "morph":morph,"global_classifiable":True,"photo_id":taxon*1000+i})
    d = f.all_species_counts(pd.DataFrame(rows))
    assert len(d) == 2
    d = d.set_index("inat_taxon_id")
    assert d.loc[1, "majority_morph"] == "blue_purple"
    assert d.loc[1, "majority_fraction"] == pytest.approx(.9)
    assert d.loc[1, "D"] == pytest.approx(.18)
    assert d.loc[2, "D"] == pytest.approx(.5)
    assert d.loc[2, "majority_morph"] == "yellow_orange"
    assert d.loc[2, "n_observed_states"] == 2
    profile = f.species_profile(d.reset_index())
    assert profile["median_photo_majority_fraction"] == pytest.approx(.70)


def test_original_null_asymmetry_preserved_with_binary_morph_names_ignored():
    # The same vertex assignments must serve both multi-state and majority/other.
    rows = []
    for i in range(50):
        rows.append({"inat_taxon_id":7,"species":"Viola test",
                     "fcp_source_cohort":"prospective","photo_id":900+i,
                     "morph":"white" if i%3 else "red_pink"})
    raw = pd.DataFrame(rows)
    d,y,z=f.fixed_label_perm(raw,"third")
    assert len(y) == 50 and z.shape == (199,50)
    from scripts.analysis.run_polymorphism_distributed_polymorphism_posthoc_20261007 import stable_seed
    rng=np.random.default_rng(stable_seed("third",7))
    assert np.array_equal(z[0],rng.permutation(y))
    assert np.all(z.sum(axis=1) == z[0].sum())
    # Majority does not change inside the permutation because state counts are frozen.
    dominant=np.bincount(y,minlength=4).argmax()
    assert np.all((z==dominant).sum(axis=1)==(y==dominant).sum())


def test_edge_based_statistic_requires_30_local_pairs():
    rng=np.random.default_rng(42)
    y=np.array([0,0,1,1,0,0,1,1,0,1],dtype=np.int8)
    z=np.stack([rng.permutation(y) for _ in range(199)])
    u,v=np.triu_indices(len(y),1)
    global_D=np.mean(y[u]!=y[v])
    assert f.obs_and_null_over_edges(y,z,u,v,global_D,np.arange(len(u))<29) is None
    out=f.obs_and_null_over_edges(y,z,u,v,global_D,np.arange(len(u))<30)
    assert out is not None
    obs,nul=out
    assert np.isfinite(obs) and np.isfinite(nul).all() and nul.shape==(199,)


def test_pooled_null_not_pseudoreplicated_by_photograph_pair():
    d=np.array([.02,.04,-.01])
    z=np.zeros((3,199))
    result=f.mean_null(d,z)
    assert result["n_species"]==3
    assert result["mean_depletion"]==pytest.approx(.05/3)
    assert result["p_upper"]==.005
    assert result["positive_species_fraction"]==pytest.approx(2/3)
    with pytest.raises(ValueError,match="nonfinite"):
        f.mean_null(np.array([np.nan,.01]),np.zeros((2,199)))
    with pytest.raises(ValueError,match="null matrix"):
        f.mean_null(np.array([.01,.02]),np.zeros((2,9)))


def test_release_byte_hashes_are_pinned():
    assert f.PHOTOS_SHA=="8e7f10f0427927208485b2da6fa9c46bb681d6f3611e53b8dfdc606cfc8dabe6"
    assert f.SPATIAL_SHA=="960e8c4bae8a9ecb01fefbef06b6b88cd8fdc63d86af3e4c354cab4eb1aedd35"
    assert f.SPATIAL_RESULT_SHA=="79aeab0555ec2efd865606a9bdc3e06574a5168c38029ee16e22b945e1a9eacf"
    assert f.PERMUTATIONS==199
    assert f.PRIMARY_RADIUS==50.
    assert f.MIN_LOCAL_PAIRS==30

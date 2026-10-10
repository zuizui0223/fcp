"""Synthetic falsification tests for FCP label-swap stress test.

No biological source data or imputed 'true' colour labels in these tests.
"""
import importlib.util
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

path=Path(__file__).resolve().parents[1]/"scripts/analysis/audit_polymorphism_adversarial_label_swaps_20261010.py"
spec=importlib.util.spec_from_file_location("fcp_adversary",path)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def photographs(labels, sites):
    n=len(labels)
    return pd.DataFrame({
        "inat_taxon_id":[1]*n, "species":["A test"]*n,
        "photo_id":np.arange(n)+1000, "latitude":[0.]*n,
        "longitude":sites, "morph":labels})


def spatial_species():
    return m.Species(photographs(["white"]*50+["red_pink"]*50,
                                 [0.]*50+[2.]*50))


def test_exact_swap_gain_formula_matches_recomputed_pair_changes():
    s=spatial_species()
    candidate=s.candidate()
    assert candidate is not None
    gain,delta,i,j=candidate
    local0=s.local
    assert delta>0
    s.swap(i,j,delta)
    assert np.isclose(s.local-local0,gain,atol=1e-12)


def test_swap_preserves_every_colour_count_and_never_reuses_photo():
    s=spatial_species()
    original=np.bincount(s.orig,minlength=4).copy()
    n=0
    for _ in range(4):
        best=s.candidate()
        assert best is not None
        _,delta,i,j=best
        assert i!=j and not s.touched[i] and not s.touched[j]
        s.swap(i,j,delta)
        n+=1
        assert np.array_equal(np.bincount(s.labels,minlength=4),original)
        assert int(s.touched.sum())==2*n


def test_colour_invariant_locality_is_not_artificially_broken():
    s=m.Species(photographs(["white"]*100,[0.]*50+[2.]*50))
    assert s.depletion==0
    assert s.candidate() is None


def test_input_original_labels_are_immutable_and_change_is_monotone():
    g=photographs(["white"]*50+["red_pink"]*50,[0.]*50+[2.]*50)
    original=g.copy(deep=True)
    s=m.Species(g)
    deps=[s.depletion]
    for _ in range(5):
        cand=s.candidate()
        assert cand is not None
        _,delta,i,j=cand
        s.swap(i,j,delta)
        deps.append(s.depletion)
    assert all(a>b for a,b in zip(deps,deps[1:]))
    pd.testing.assert_frame_equal(g,original)


def test_same_colour_neighbours_within_one_site_have_expected_pair_count():
    s=spatial_species()
    assert s.m==2*(50*49//2)
    assert np.isclose(s.local,0)
    assert np.isclose(s.fixed_global,2*50*50/(100*99))
    assert s.depletion>0.5

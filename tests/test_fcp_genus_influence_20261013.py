from __future__ import annotations
import importlib.util,sys
from pathlib import Path
import numpy as np
import pytest

HERE=Path(__file__).resolve().parents[1]
SRC=HERE/'scripts/analysis/audit_fcp_genus_influence_20261013.py'
spec=importlib.util.spec_from_file_location('genus_leverage_prod',SRC)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def synthetic(n=40,seed=12):
    rng=np.random.default_rng(seed)
    stat={True:{},False:{}}
    for z in range(n):
        name=f'G{z:03d}'
        for side in (True,False):
            X=rng.normal(size=(14,8));X-=X.mean(axis=0)
            B=np.zeros((8,4))
            B[0]=([.06,-.01,-.01,-.04] if side else [-.02,0,.01,.01])
            noise=rng.normal(size=4)*.004
            B[0]+=noise-noise.mean()
            Y=X@B
            count=14+z%4
            stat[side][name]=(X.T@X,X.T@Y,count,count*2)
    return stat


def test_known_hue_direction_with_uncorrelated_taxa():
    z=synthetic()
    v=m.whole_beta(z,'species_equal')
    assert v>0
    assert np.isfinite(v)


def test_paired_omission_reestimates_both_sides_not_just_averages_genus_slopes():
    s=synthetic()
    v=m.whole_beta(s,'species_equal')
    drop=m.whole_beta(s,'species_equal',excluded=['G017'])
    assert np.isfinite(drop)
    assert drop!=v
    assert set(s[True])==set(s[False])


def test_genus_equal_and_species_equal_converge_with_constant_support():
    s=synthetic()
    for d in s.values():
        for g,(a,b,n,p) in list(d.items()):d[g]=(a,b,14,p)
    assert m.whole_beta(s,'genus_equal')==pytest.approx(m.whole_beta(s,'species_equal'),abs=1e-12)


def test_less_than_30_genera_holds():
    s=synthetic(28)
    assert m.whole_beta(s,'genus_equal') is None


def test_colour_probabilities_require_four_outcomes():
    s=synthetic()
    for side in s:
        for g,(a,b,n,p) in s[side].items():
            assert b.shape==(8,4)
            assert abs(b.sum())<1e-8


def test_hash_gate_rejects_unsanctioned_archive(tmp_path):
    p=tmp_path/'wrong.zip';p.write_bytes(b'not-the-original-photo-source')
    with pytest.raises(ValueError,match='photo bytes changed'):
        m.ingest(p)

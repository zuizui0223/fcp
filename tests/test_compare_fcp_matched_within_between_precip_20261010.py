"""Synthetic controls for source photo-aligned scale comparisons."""
import importlib.util
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

p=Path(__file__).resolve().parents[1]/"scripts/analysis/compare_fcp_matched_within_between_precip_20261010.py"
s=importlib.util.spec_from_file_location("fcp_aligned",p)
m=importlib.util.module_from_spec(s)
s.loader.exec_module(m)

def fake_train():
    rows=[]
    for genus in ["A","B"]:
        for taxon in range(1,6):
            for point in range(2):
                n=(0 if genus=="A" else 100)+taxon
                rows.append({"inat_taxon_id":n,"genus":genus,"species":f"{genus} species_{taxon}",
                  "morph":m.CLASSES[(taxon+point)%4],
                  "abs_latitude":10.+n*.02+point,
                  "lon_sin":.2+.001*n,
                  "lon_cos":.9-.001*n,
                  "wc_elevation_m":100.+30*taxon,
                  "wc_bio1":10.+point,"wc_bio5":30.+taxon,
                  "wc_bio12":1000.+100*taxon,"wc_bio15":45.+point})
    return pd.DataFrame(rows)


def test_real_four_state_source_dim_and_coarse_environment_schema():
    assert len(m.CLASSES)==4
    assert set(m.MOISTURE)=={"wc_bio12","wc_bio15"}
    assert set(m.TEMPERATURE)=={"wc_bio1","wc_bio5"}
    assert set(m.GEO)=={"abs_latitude","lon_sin","lon_cos","wc_elevation_m"}
    assert not set(m.MOISTURE)&set(m.GEO)


def test_species_fold_is_source_taxon_stable_and_colour_independent():
    ids=[2384,383,140,1123]
    assert [m.phylogeny_free_species_fold(x) for x in ids]==[
        m.phylogeny_free_species_fold(x) for x in ids]
    assert all(0<=m.phylogeny_free_species_fold(x)<5 for x in ids)


def test_focal_species_is_forbidden_from_between_training():
    d=fake_train()
    focal=d[d.inat_taxon_id.eq(1)].iloc[:1]
    with pytest.raises(RuntimeError,match="leaked"):
        m.predict_between_other_species(d,focal,m.FULL_FEATURES)


def test_same_other_genus_prediction_is_four_class_simplex():
    d=fake_train()
    focal=d[d.inat_taxon_id.eq(1)].iloc[:1]
    train=d[d.inat_taxon_id.ne(1)]
    pred=m.predict_between_other_species(train,focal,m.FULL_FEATURES)
    assert pred.shape==(1,4)
    assert np.all(pred>=0)
    assert np.allclose(pred.sum(axis=1),1)


def test_3_other_species_genus_gate_is_fail_closed():
    d=fake_train()
    focal=d[d.inat_taxon_id.eq(1)].iloc[:1]
    train=d[(d.inat_taxon_id!=1)&(d.inat_taxon_id!=2)&(d.inat_taxon_id!=3)]
    with pytest.raises(RuntimeError,match="fewer than 3"):
        m.predict_between_other_species(train,focal,m.FULL_FEATURES)


def test_same_photo_cluster_ci_constant_and_signal_are_finite():
    values=np.full(30,.02)
    groups=np.repeat(np.arange(10),3)
    ci=m.cluster_ci(values,groups,label="synthetic")
    assert np.allclose(ci,[.02,.02])

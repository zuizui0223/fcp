"""Pre-outcome FCP next2000 synthetic pipeline and frozen engine configuration tests."""
import importlib.util
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
PATH=ROOT/"scripts/analysis/qualify_fcp_next2000_synthetic_20261008.py"
spec=importlib.util.spec_from_file_location("future_fcp_synthetic",PATH)
q=importlib.util.module_from_spec(spec)
spec.loader.exec_module(q)


def test_synthetic_geography_not_from_any_selected_real_species():
    s=q.synthetic_species(900000000)
    assert len(s)==80
    assert s.inat_taxon_id.nunique()==1
    assert s.latitude.nunique()==2
    assert s.month.nunique()==2
    assert s.morph.nunique()==2
    assert len(s.observer.unique())>=8
    assert int((s.morph=="white").sum())==40
    assert int((s.morph=="red_pink").sum())==40
    assert 900000000 > 1_000_000


def test_pure_time_effect_has_no_added_spatial_difference():
    d=q.synthetic_species(999999999,seasonal_only=True)
    assert int((d.loc[d.latitude==10,"morph"]=="white").sum())==20
    assert int((d.loc[d.latitude==12,"morph"]=="white").sum())==20
    assert d.groupby("month").morph.nunique().max()==1


def test_single_observer_negative_control():
    x=q.synthetic_species(888888888,one_observer=True)
    assert x.observer.nunique()==1
    assert x.morph.nunique()==2


def test_prospective_contract_source_hashes_and_no_any_live_colour_model():
    assert q.PREOPEN_MANIFEST_SHA256=="e1301095533dfa755c1588cc550522bf31aabaef632083e265a138f0e3af6d2a"
    assert q.HISTORICAL_ENGINE_COMMIT=="969c67d274868b5dc4e8c2c8c9fa5d1f1beb5527"
    assert q.ENGINE_SOURCE_GIT_BLOB_SHA=="c7b277c249fc6121302b35ce2960e2003cfd0930"
    source=PATH.read_text()
    assert "import requests" not in source
    assert "openai" not in source.casefold()
    assert "image_url" not in source
    assert "fcp_next2000_synthetic_end_to_end_qualification_v1" in source


def test_synthetic_photo_ids_are_unique():
    s=[q.synthetic_species(900000000+j) for j in range(4)]
    ids=np.concatenate([x.photo_id.to_numpy() for x in s])
    assert len(np.unique(ids))==len(ids)

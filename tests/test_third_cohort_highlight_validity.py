from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parents[1]

def load(path:str,name:str):
    spec=importlib.util.spec_from_file_location(name,ROOT/path)
    module=importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module

def test_highlight_metrics_exact_fixture():
    m=load("scripts/analysis/measure_third_cohort_highlight_partition_20260922.py","hl_measure")
    rgb=np.array([[255,0,0],[249,249,249],[250,10,10],[0,0,0]],dtype=np.uint8)
    clip,near,q99=m.highlight_metrics(rgb)
    assert clip==0.25
    assert near==0.50
    assert 0.0 <= q99 <= 1.0

def test_highlight_id_and_partition_are_deterministic():
    m=load("scripts/analysis/build_third_cohort_highlight_firewall_20260922.py","hl_firewall")
    a=m.hid(304429869)
    b=m.hid(304429869)
    assert a==b
    assert a.startswith("FCPH2HL-")
    assert m.slot(a)==m.slot(b)
    batch,shard,part=m.slot(a)
    assert 0<=batch<2
    assert 0<=shard<32
    assert 0<=part<4

def test_white_contrast_has_unit_norm_and_zero_sum():
    m=load("scripts/analysis/run_third_cohort_highlight_validity_20260922.py","hl_final")
    q=m.qwhite()
    assert np.isclose(q.sum(),0.0,atol=1e-12)
    assert np.isclose(np.linalg.norm(q),1.0,atol=1e-12)

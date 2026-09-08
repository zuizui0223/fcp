from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
A=ROOT/'scripts/analysis'; sys.path.insert(0,str(A))
spec=importlib.util.spec_from_file_location('wrap',A/'run_rgfca_sharedness_v2_synthetic_qualification_transport.py')
wrap=importlib.util.module_from_spec(spec); assert spec.loader; spec.loader.exec_module(wrap)
v2=wrap.base


def test_transport_geometry_is_exact_150x300():
    _,mapping,frame,xyz,train,test=wrap.load_transport_inputs()
    assert len(frame)==45000 and xyz.shape==(150,300,3)
    assert len(train)==75 and len(test)==75 and not(set(train)&set(test))
    assert mapping['execution']['observed_flower_colour_opened'] is False


def test_fixed_arm_census():
    _,m,*_=wrap.load_transport_inputs()
    assert len(v2.nuisance_arms(m))==21 and len(v2.positive_scenarios(m))==19 and len(v2.evaluation_arms(m))==78


def test_hard_cells_have_two_sd_variants():
    _,m,*_=wrap.load_transport_inputs(); p=v2.positive_scenarios(m)
    for a,f in [(1.,1.),(1.,.5),(2.,.5)]:
        assert sorted(x['threshold_sd'] for x in p if x['amplitude']==a and x['shared_fraction']==f)==[0.,.25]


def test_zero_evidence_score_zero_and_eval_never_changes_posterior():
    _,_,_,_,train,test=wrap.load_transport_inputs(); lr=np.zeros((150,96))
    a=v2.predictive_score(lr,train,test); assert abs(a[0])<1e-12
    rng=np.random.default_rng(4); lr=rng.normal(size=(150,96)); a=v2.predictive_score(lr,train,test)
    lr[test]+=rng.normal(scale=3,size=(75,96)); b=v2.predictive_score(lr,train,test)
    assert np.allclose(a[1:],b[1:],atol=1e-12,rtol=1e-12)


def test_retention_counts_are_exact_and_deterministic():
    _,m,_,xyz,_,_=wrap.load_transport_inputs(); root=int(m['execution']['deterministic_seed_root'])
    for r,k in [(.4,120),(.6,180),(.8,240)]:
        a=v2.retention_indices(root,'calibration',1,r,xyz,False); b=v2.retention_indices(root,'calibration',1,r,xyz,False)
        assert a.shape==(150,k) and np.array_equal(a,b)
        assert all(len(np.unique(x))==k for x in a)


def test_calibration_is_structured_null_only():
    _,m,*_=wrap.load_transport_inputs()
    assert all(x['kind']=='nuisance' and float(x.get('shared_fraction',0))==0 for x in v2.nuisance_arms(m))


def test_transport_wrapper_opens_no_biological_inputs():
    text=(A/'run_rgfca_sharedness_v2_synthetic_qualification_transport.py').read_text()
    assert 'global_monte_carlo_measured_photos' not in text and '34species' not in text and 'six_species' not in text

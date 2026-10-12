import importlib.util
import math
import sys
from pathlib import Path
import pytest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))

def module():
    path=ROOT/'sampling_bound.py'
    assert path.exists(), 'sampling_bound.py is not implemented'
    spec=importlib.util.spec_from_file_location('sampling_bound',path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    return m

def test_threshold_and_closed_form():
    m=module()
    assert m.tail_upper(0,100,3)==1
    assert m.tail_upper(.02,100,3)==1
    assert m.tail_upper(.1,100,3)==pytest.approx(math.exp(-10)*(math.e*5)**2)

@pytest.mark.parametrize('b,n,k',[(-1,100,3),(math.nan,100,3),(math.inf,100,3),(.1,0,3),(.1,100,1),(.1,2.5,3),(.1,True,3)])
def test_invalid_bound_inputs(b,n,k):
    with pytest.raises(ValueError):module().tail_upper(b,n,k)

def test_bound_monotone_and_inversion():
    m=module()
    v=[m.tail_upper(.2,n,3) for n in (1,25,50,75,100)]
    assert v==sorted(v,reverse=True)
    n=m.minimum_n(.2,3,.05/3)
    assert m.tail_upper(.2,n,3)<=.05/3
    assert m.tail_upper(.2,n-1,3)>.05/3

def test_census_not_confused_with_one_observation():
    m=module()
    assert m.validate_design(1,1,100)=='census'
    assert m.validate_design(100,500,100)=='subsample'
    with pytest.raises(ValueError):m.validate_design(1,500,100)
    with pytest.raises(ValueError):m.validate_design(10,5,100)

def test_jensen_kl_identity_and_mixed_census():
    m=module(); q=[[.5,.5,0],[1,0,0]];w=[.25,.75];p=[.5,.5,0]
    b=m.gap(q,w);bar=[sum(a*r[j] for a,r in zip(w,q)) for j in range(3)]
    d=sum(a*m.kl(r,p) for a,r in zip(w,q))
    assert b==pytest.approx(d-m.kl(bar,p))
    assert 0<=b<=d
    with pytest.raises(ValueError):m.gap(q,[1,1])

def test_exact_small_distribution_guards():
    r=module().exact_checks()
    assert r['checks']>=100
    assert r['status']=='PASS'

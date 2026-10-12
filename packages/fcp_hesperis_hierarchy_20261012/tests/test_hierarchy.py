import math
import sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from hierarchy import entropy,vertices,entropy_bounds,decompose,evaluate

def fixture():
    return dict(year=2021,n_all=132,n_poly=118,mono_white=0,mono_pink=2,mono_purple=12,
                all_white=.129,all_pink=.333,all_purple=.538,poly_white=.144,poly_pink=.355,poly_purple=.501,mean_H_all=.71)

def test_entropy_pure():
    assert entropy([1,0,0])==0

def test_entropy_uniform():
    assert entropy([1/3]*3)==pytest.approx(math.log(3))

def test_entropy_invalid_sum():
    with pytest.raises(ValueError): entropy([.1,.2,.3])

def test_entropy_negative():
    with pytest.raises(ValueError): entropy([-.1,.5,.6])

def test_poly_identical_no_withingroup_component():
    d=decompose(4,2,[0,0,2],[.5,.5,0],.5*math.log(2))
    assert d['poly']==pytest.approx(0,abs=1e-14)
    assert d['mono']==pytest.approx(0,abs=1e-14)
    assert d['status']==pytest.approx(math.log(2))

def test_positive_poly_and_chain_identity():
    d=decompose(4,2,[0,0,2],[.5,.5,0],.5*entropy([.9,.1,0]))
    assert d['poly']>0
    assert d['poly']+d['mono']+d['status']==pytest.approx(d['beta'])

def test_all_poly():
    d=decompose(2,2,[0,0,0],[.5,.5,0],entropy([.9,.1,0]))
    assert d['poly']==pytest.approx(d['beta'])
    assert d['poly_share']==pytest.approx(1)
    assert d['mono']==0

def test_pure_zero_diversity_fraction_undefined():
    d=decompose(2,2,[0,0,0],[1,0,0],0)
    assert d['poly_share'] is None

def test_bad_denominator():
    with pytest.raises(ValueError): decompose(4,2,[0,0,1],[.5,.5,0],.2)

def test_noninteger_count():
    with pytest.raises(ValueError): decompose(4,2,[0,0,2.5],[.5,.5,0],.2)

def test_impossible_within_entropy():
    with pytest.raises(ValueError): decompose(4,2,[0,0,2],[.5,.5,0],1)

def test_vertices_known_triangle():
    v=vertices([0,0,0],[1,1,1]); assert len(v)==3
    assert all(sum(p)==pytest.approx(1) for p in v)

def test_infeasible_summary():
    with pytest.raises(ValueError): vertices([.4]*3,[.5]*3)

def test_entropy_bounds_full_simplex():
    b=entropy_bounds([0]*3,[1]*3)
    assert b[0]==pytest.approx(0)
    assert b[1]==pytest.approx(math.log(3))

def test_entropy_bounds_box():
    b=entropy_bounds([.1,.2,.4],[.3,.4,.6]); assert b[0]<=b[1]
    for p in vertices([.1,.2,.4],[.3,.4,.6]): assert b[0]-1e-13<=entropy(p)<=b[1]+1e-13

def test_joint_rounding_representation():
    d=evaluate(fixture()); q=d['representative_poly_mean']; m=d['representative']['mu']
    assert sum(q)==pytest.approx(1)
    for k,c in enumerate(['white','pink','purple']):
        assert abs(q[k]-fixture()['poly_'+c])<=.0005+1e-12
        assert abs(m[k]-fixture()['all_'+c])<=.0005+1e-12
    assert d['interval_type']=='CONSERVATIVE_ROUNDING_OUTER_ENCLOSURE_NOT_CI'

def test_inconsistent_all_and_poly_means_rejected():
    r=fixture();r['all_white']=.7
    with pytest.raises(ValueError): evaluate(r)

def test_poly_shares_sum_0999_not_silently_renormalized():
    r=dict(year=2022,n_all=112,n_poly=102,mono_white=0,mono_pink=1,mono_purple=9,
           all_white=.174,all_pink=.295,all_purple=.531,poly_white=.190,poly_pink=.314,poly_purple=.495,mean_H_all=.77)
    # The published 2022 white frequencies are not jointly rounding-feasible.
    # Rejection, not normalization, is the scientifically required outcome.
    with pytest.raises(ValueError,match='infeasible'): evaluate(r)

def test_enclosures_include_vertices_and_alpha_extremes():
    d=evaluate(fixture()); b=d['enclosures']
    for q in d['feasible_vertices']:
        for a in [.705,.715]:
            x=decompose(132,118,[0,2,12],q,a)
            for key in ['gamma','alpha','beta','poly','mono','status','poly_share']:
                assert b[key][0]-1e-12<=x[key]<=b[key][1]+1e-12


def test_declared_discrepancy_sensitivity_accepts_2022():
    r=dict(year=2022,n_all=112,n_poly=102,mono_white=0,mono_pink=1,mono_purple=9,
           all_white=.174,all_pink=.295,all_purple=.531,poly_white=.190,poly_pink=.314,poly_purple=.495,mean_H_all=.77)
    d=evaluate(r,freq_halfwidth=.001)
    assert d['interval_type']=='PUBLICATION_DISCREPANCY_SENSITIVITY_NOT_CI'
    assert sum(d['representative_poly_mean'])==pytest.approx(1)

def test_rounding_inconsistency_exact_integer_witness():
    # Lower all-stand white mean x112 exceeds upper poly-only mean x102.
    assert 1735*112 > 1905*102

def test_rounded_0999_synthetic_feasible_without_normalizing():
    r=dict(year=2000,n_all=4,n_poly=3,mono_white=0,mono_pink=0,mono_purple=1,
           all_white=.25,all_pink=.25,all_purple=.5,poly_white=.333,poly_pink=.333,poly_purple=.333,mean_H_all=.6)
    d=evaluate(r)
    assert d['literal_poly_sum']==pytest.approx(.999)
    assert sum(d['representative_poly_mean'])==pytest.approx(1)

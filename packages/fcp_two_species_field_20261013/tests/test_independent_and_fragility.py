"""Source-backed independent verifier and *posthoc* site-influence guard tests."""
import sys
from pathlib import Path
import numpy as np
import pytest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from verify_independent import exact_metric,get_original_counts,numerical_verify
from posthoc_site_fragility import audit_one
from field_partition import load_sisyrinchium,load_boechera

def test_exact_example_excess_allocation():
    # Twenty plants meets the frozen eligibility gate. The two-colour blocks
    # are maximally segregated and each quadrat has size ten.
    example=exact_metric([(10,0),(10,10)])
    assert example['B']==1
    assert example['expected_B']==pytest.approx(1/19)
    assert example['delta']==pytest.approx(18/19)
    assert example['expected_mixed']==pytest.approx(1-2/184756)
    # A distinct hand-enumerated N=4 micro-case verifies the null expectation
    # independent of the N>=20 inferential eligibility threshold.
    import itertools
    Bs=[]
    for positions in itertools.combinations(range(4),2):
        x=sum(k<2 for k in positions)
        p=.5
        Bs.append((((x-2*p)**2/2)+(((2-x)-2*p)**2/2))/(4*p*(1-p)))
    assert sum(Bs)/len(Bs)==pytest.approx(1/3)

def test_exact_null_defined_only_for_both_colours():
    assert exact_metric([(30,30),(23,23)]) is None
    assert exact_metric([(30,4),(23,0)]) is None
    assert exact_metric([(0,0),(40,20)]) is None

def test_from_original_source_not_pseudo_replicates():
    src,hashes,qa=get_original_counts(ROOT)
    assert len(src['sisyrinchium'])==16 and len(src['boechera'])==9
    assert sum(len(v) for v in src['sisyrinchium'].values())==512
    assert sum(len(v) for v in src['boechera'].values())==26
    assert sum(sum(x[0] for x in v) for v in src['boechera'].values())==547
    assert qa['sisyrinchium_other_l']+qa['sisyrinchium_other_s']==4492
    assert qa['boechera_id_whitespace_repairs']==3

def test_full_source_to_result_independent_verifier():
    v=numerical_verify(ROOT,9999)
    assert v['status']=='PASS' and v['checks']>=170
    assert v['max_abs_difference']<1e-12

def test_leave_one_judd_falls_is_not_robust():
    _,b,_=load_boechera(ROOT/'source/Boechera_natural_populations_individual.txt',ROOT/'source/Boechera_natural_populations_plot.txt')
    r,ch=audit_one(b,'boechera',R=9999)
    d=next(x for x in r if x['excluded_siteyear'].endswith('Judd_Falls'))
    assert ch['original_p_recomputed']==pytest.approx(.0001)
    assert d['after_leave_one_out_delta_percentage_points']==pytest.approx(.22526370384348796)
    assert d['after_leave_one_out_p_one_sided']==pytest.approx(.3631)

def test_sisyrinchium_leave_any_siteyear_still_positive():
    s=load_sisyrinchium(ROOT/'source/Sisyrinchium_data.zip')
    r,ch=audit_one(s,'sisyrinchium',R=9999)
    assert ch['sites_used']==13
    assert len(r)==13
    assert all(x['after_leave_one_out_delta_percentage_points']>40 for x in r)
    assert all(x['after_leave_one_out_p_one_sided']==pytest.approx(.0001) for x in r)

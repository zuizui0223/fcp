import io
import sys
import zipfile
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from predictive import load_original, score, ridge_objective_gradient, fit_binomial, derive_cohort, evaluate


def test_original_source_and_denominator():
    rows = load_original(ROOT/'source'/'Sisyrinchium_data.zip')
    assert len(rows) == 512
    assert sum(r['v'] for r in rows) == 5081
    assert sum(r['w'] for r in rows) == 8723
    assert sum(r['l'] for r in rows) == 4463
    assert sum(r['s'] for r in rows) == 29
    assert sum(r['v']+r['w']==0 for r in rows) == 126
    assert sum(r['distance']<0 for r in rows) == 3
    assert sum(r['distance']<0 and r['v']+r['w']>0 for r in rows) == 2
    assert len({(r['year'], r['population'], r['line'], r['quad']) for r in rows}) == 512


def test_gradient_matches_numerical_differences():
    n=np.array([5.,12.,7.,3.]);v=np.array([1.,6.,4.,1.]);x=np.array([[0.4,-2.],[-.2,1.],[1.2,0.],[-.9,0.3]])
    theta=np.array([-0.4,.2,.8])
    objective=lambda z:ridge_objective_gradient(z,x,n,v,lam=10.)[0]
    analytical=ridge_objective_gradient(theta,x,n,v,lam=10.)[1]
    step=1e-6
    numerical=np.array([(objective(theta+np.eye(len(theta))[j]*step)-objective(theta-np.eye(len(theta))[j]*step))/(2*step) for j in range(len(theta))])
    assert np.max(np.abs(analytical-numerical)) < 1e-6


def test_fit_baseline_and_no_test_leakage():
    train=[dict(n=5,v=2,distance=0,other_share=0.0),dict(n=5,v=3,distance=10,other_share=.1),dict(n=5,v=4,distance=20,other_share=.2)]
    heldout=[dict(n=1000,v=0,distance=5,other_share=.3)]
    p=fit_binomial(train,heldout,['distance','other_share'])
    assert p.shape==(1,)
    assert 0.01 <=p[0]<=0.99
    q=fit_binomial(train,[dict(heldout[0],v=1000)],['distance','other_share'])
    assert np.array_equal(p,q)  # never use heldout v
    base=fit_binomial(train,heldout,[])
    assert np.isclose(base[0],(9.5/16))


def test_unidentifiable_covariate_gives_baseline():
    train=[dict(n=20,v=5,distance=0,other_share=0),dict(n=20,v=10,distance=0,other_share=0)]
    p0=fit_binomial(train,[dict(n=10,v=1,distance=123,other_share=5)],[])
    p1=fit_binomial(train,[dict(n=10,v=1,distance=123,other_share=5)],['distance','other_share'])
    assert np.allclose(p0,p1,rtol=0,atol=1e-9)


def test_logloss_prefers_correct_fraction():
    rows=[dict(n=10,v=7),dict(n=20,v=3)]
    good=score(rows,[.7,.15]);bad=score(rows,[.2,.8])
    assert good['logloss_per_plant']<bad['logloss_per_plant']
    assert good['fraction_brier_equal_plot']<bad['fraction_brier_equal_plot']


def test_cohort_nested_validity_and_missingness():
    rows=load_original(ROOT/'source'/'Sisyrinchium_data.zip')
    eligible,held=derive_cohort(rows)
    assert len(eligible)>=5
    assert 'Y13|chitei' in held
    assert 'Y11|otto' in held
    assert 'Y13|ootaka2' in held
    assert sum(x['year']=='Y13' and x['population']=='nature' for x in eligible)==1
    assert len([x for x in eligible if x['population']=='nature'])==2
    assert all(len(x['heldout_lines'])>=3 for x in eligible)


def test_full_eval_prediction_rows_and_primary_latest():
    rows=load_original(ROOT/'source'/'Sisyrinchium_data.zip')
    result=evaluate(rows)
    assert len(result['all_siteyears'])==16
    assert len(result['selected_latest'])==len(set(x['population'] for x in result['eligible_siteyears']))
    assert all(p['n']>0 and p['distance']>=0 for p in result['predictions'])
    keys=[(p['year'],p['population'],p['line'],p['quad']) for p in result['predictions']]
    assert len(keys)==len(set(keys))
    for p in result['predictions']:
        for k in ('BASE','DIST','CONGEN','BOTH'):
            assert .01<=p['p_'+k]<=.99


def test_archive_tamper_fails_authentication(tmp_path):
    original=(ROOT/'source'/'Sisyrinchium_data.zip').read_bytes()
    damaged=bytearray(original);damaged[900]^=0x1
    p=tmp_path/'bad.zip';p.write_bytes(damaged)
    with pytest.raises(ValueError,match='HASH_OR_LENGTH_MISMATCH'):
        load_original(p)


def test_heldout_folds_are_mutually_disjoint_and_dates_not_mixed():
    r=evaluate(load_original(ROOT/'source'/'Sisyrinchium_data.zip'))
    known={(x['year'],x['population']):x for x in r['eligible_siteyears']}
    seen=set()
    for x in r['fold_ledger']:
        group=(x['year'],x['population']);assert group in known
        assert x['line'] in known[group]['heldout_lines']
        assert x['test_n']>=10 and x['train_v']>=5 and x['train_w']>=5
    for row in r['predictions']:
        k=(row['year'],row['population'],row['line'],row['quad'])
        assert k not in seen;seen.add(k)
    assert len(r['fold_ledger'])==37
    assert len(r['predictions'])==185
    assert len(seen)==185
    assert r['main']['n_sites']==6
    assert r['main']['n_lines_heldout']==32
    assert r['main']['n_plants_scored']==6692


def test_same_predictions_deterministic_under_rescore():
    # The full analysis is deterministic, and score does not infer new plants.
    rows=load_original(ROOT/'source'/'Sisyrinchium_data.zip')
    a=evaluate(rows)
    b=evaluate(rows)
    import json
    assert json.dumps(a,sort_keys=True)==json.dumps(b,sort_keys=True)


def test_standalone_independent_newton_verifier_reads_original():
    import importlib.util
    p=ROOT/'src'/'verify_independent.py'
    spec=importlib.util.spec_from_file_location('independent_verify',p)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    v=module.verify(ROOT)
    assert v['status']=='PASS' and v['checks']>=1700
    assert v['max_abs_difference'] < 2e-7


def test_exact_hessian_newton_polish_recovers_numerical_optimizer():
    from predictive import hessian_newton_polish
    x=np.array([[-1.0,0.0],[0.2,1.1],[1.4,-0.3],[0.9,1.1]],dtype=float)
    n=np.array([100,300,240,80],float)
    v=np.array([1,45,200,10],float)
    rough=np.array([-0.4,0.2,0.1],float)
    refined=hessian_newton_polish(rough,x,n,v,lam=10)
    assert np.all(np.isfinite(refined))
    assert np.max(np.abs(ridge_objective_gradient(refined,x,n,v,lam=10)[1])) < 1e-6
    assert ridge_objective_gradient(refined,x,n,v,lam=10)[0] <= ridge_objective_gradient(rough,x,n,v,lam=10)[0]

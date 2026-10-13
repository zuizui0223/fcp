"""Focused source and mathematical guards for the source-exposed slope-transfer study."""
import csv,hashlib,json,math,sys
from pathlib import Path
import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
import transport
import posthoc_standardized
import verify_independent


def original():return transport.load_original(ROOT/'source/Sisyrinchium_data.zip')

def results():return json.loads((ROOT/'results/result.json').read_text())


def test_verified_publisher_original_and_colour_denominators():
    x=original()
    assert len(x)==512
    assert [sum(z[k] for z in x) for k in ('v','w','l','s')]==[5081,8723,4463,29]
    assert sum(z['n']==0 for z in x)==126
    assert sum(z['distance']<0 for z in x)==3
    assert sum(z['distance']<0 and z['n']>0 for z in x)==2


def test_wrong_or_tampered_source_bytes_are_held(tmp_path):
    b=bytearray((ROOT/'source/Sisyrinchium_data.zip').read_bytes());b[140]^=1
    q=tmp_path/'bad.zip';q.write_bytes(b)
    with pytest.raises(ValueError,match='SOURCE_BYTES_INVALID'):
        transport.load_original(q)


def test_correct_fixed_site_population_and_32_heldout_lines():
    z=results();s=z['sites'];p=z['predictions'];f=z['folds']
    assert len(s)==6 and len(f)==32 and len(p)==168
    assert sum(x['n'] for x in p)==6692
    assert {x['siteyear'] for x in s}=={'Y13|banpaku','Y13|nature','Y11|ootaka','Y11|ootaka2','Y13|shiga','Y11|toyama'}
    assert all(x['siteyear'] not in x['source_sites'] and len(x['source_sites'])==5 for x in f)
    assert len({(x['siteyear'],x['line'],x['quad']) for x in p})==168


def test_source_slope_orientation_is_heterogeneous():
    r=results()
    assert sum(x['beta_per_distance']<0 for x in r['source_fit'].values())==5
    assert r['source_fit']['Y13|banpaku']['beta_per_distance']>0
    assert r['main']['n_training_fold_sign_agree']==32
    assert all(x['within_site_fold_sign_same']==x['n_target_training_folds'] for x in r['sites'])


def test_transferred_coefficients_exclude_target_and_are_median_or_mean_of_five():
    r=results()
    for s in r['sites']:
        sources=[x['beta_per_distance'] for label,x in r['source_fit'].items() if label!=s['siteyear']]
        assert len(sources)==5
        assert s['source_median_raw_other_sites']==float(np.median(sources))
        assert s['source_mean_raw_other_sites']==float(np.mean(sources))


def test_all_reported_parent_scores_unmodified():
    r=results()
    frozen=json.loads((ROOT/'source/parent_result.json').read_text())
    for s in r['sites']:
        p=next(q for q in frozen['selected_latest'] if q['siteyear']==s['siteyear'])
        assert math.isclose(s['logloss_BASE'],p['logloss_BASE'],abs_tol=1e-11)
        assert math.isclose(s['logloss_LOCAL_DIST'],p['logloss_DIST'],abs_tol=1e-11)
        assert s['n_M']==p['n_scored_plants']


def test_site_equal_scores_not_overweighted_by_big_sites():
    r=results();s=r['sites']
    assert math.isclose(sum(x['logloss_BORROW_MEDIAN'] for x in s)/6,r['main']['mean']['logloss_BORROW_MEDIAN'],abs_tol=1e-13)
    assert r['main']['n_borrow_positive_vs_base']==3
    assert r['main']['n_local_positive_vs_base']==6
    assert math.isclose(r['main']['mean']['gain_borrow_vs_base'],-.12660425209151252,abs_tol=1e-12)
    assert math.isclose(r['main']['mean']['gain_local_vs_base'],.1103369436783896,abs_tol=1e-12)


def test_monomorphic_source_sites_rejected_without_fake_slope():
    d=[{'distance':i*3,'n':8,'v':0,'w':8} for i in range(4)]
    with pytest.raises(ValueError,match='SOURCE_SITE_BOTH_M_REQUIRED'):transport.fit_source_slope(d)


def test_no_distance_spread_source_sites_rejected():
    d=[{'distance':3,'n':8,'v':i*2,'w':8-i*2} for i in range(4)]
    with pytest.raises(ValueError,match='SOURCE_NO_DISTANCE_SPREAD'):transport.fit_source_slope(d)


def test_intercept_calibration_slope_zero_equals_jeffreys_baseline():
    train=[{'distance':1,'v':3,'w':7,'n':10},{'distance':6,'v':8,'w':2,'n':10}]
    a,mu=transport.calibrate_target_intercept(train,0)
    assert math.isclose(math.exp(a)/(1+math.exp(a)),(11+.5)/(20+1),abs_tol=1e-12)
    assert math.isclose(mu,3.5,abs_tol=1e-12)
    assert transport.calibrate_target_intercept(train,-.5)[0]!=a


def test_synthetic_inverse_gradient_slopes():
    d=[]
    for i in range(6):
        n=50
        v=[40,35,30,20,15,10][i]
        d.append({'distance':i*2,'n':n,'v':v,'w':n-v})
    neg=transport.fit_source_slope(d)['beta_per_distance']
    reverse=[{**x,'v':x['w'],'w':x['v']} for x in d]
    pos=transport.fit_source_slope(reverse)['beta_per_distance']
    assert neg<0<pos
    assert math.isclose(neg,-pos,abs_tol=1e-9)


def test_posthoc_normalization_is_separate_and_tracks_opposite_site():
    p=json.loads((ROOT/'results/posthoc_distance_scale.json').read_text())
    assert p['status']=='POSTHOC_DIAGNOSTIC_NOT_PRIMARY'
    assert p['main']['n_gain_median_vs_base']==5
    assert p['main']['n_sites']==6 and p['main']['n_folds']==32 and p['main']['n_quadrats']==168
    assert math.isclose(p['main']['mean_logloss']['BORROW_Z_MEDIAN'],.6581847351121263,abs_tol=1e-12)
    assert next(x for x in p['sites'] if x['siteyear']=='Y13|banpaku')['gain_std_median_vs_base']<0


def test_remove_either_opposite_site_does_not_make_raw_transfer_robust():
    r=results();cuts=r['leave_one_site_out']
    assert len(cuts)==6
    assert all(x['gain_borrow_vs_base']<0 for x in cuts)


def test_analysis_is_deterministic_from_unchanged_source(tmp_path):
    params=(ROOT/'source/Sisyrinchium_data.zip',ROOT/'source/parent_result.json',ROOT/'source/parent_predictions.csv')
    r1=transport.analysis(*params,tmp_path/'a')
    r2=transport.analysis(*params,tmp_path/'b')
    assert r1['main']==r2['main']
    for p in (tmp_path/'a').iterdir():
        assert p.read_bytes()==(tmp_path/'b'/p.name).read_bytes(),p.name


def test_independent_newton_replays_both_primary_and_posthoc():
    z=verify_independent.verify(ROOT/'source/Sisyrinchium_data.zip',ROOT/'results/result.json',ROOT/'results/posthoc_distance_scale.json',ROOT/'results')
    assert z['status']=='PASS' and z['n_independent_arithmetic_checks']>=1000
    assert z['max_absolute_numeric_difference']<3e-9


def test_named_populations_15_day_year_conditional_not_new_site_predictive():
    z=results()
    assert z['scope']=='SOURCE_EXPOSED_EXPLORATORY_NOT_PREREGISTERED'
    assert z['main']['n_sites']==6
    assert not 'p_value' in z['main']
    assert all(z['original_source'][k]>0 for k in ('M_purple','M_white'))

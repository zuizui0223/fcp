import sys
import json
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from transfer import build_strata, select_source_and_test, evaluate, model_specs, source_qa

@pytest.fixture(scope='module')
def original():
    return json.loads((ROOT/'parent/results/result.json').read_text())

@pytest.fixture(scope='module')
def rows(original):
    from transfer import load_original_verified
    return load_original_verified(ROOT/'parent/source/Sisyrinchium_data.zip',ROOT/'parent/src/predictive.py')


def test_expected_source_and_frozen_evaluation_cohort(rows,original):
    qa=source_qa(rows)
    assert qa['n_quadrats']==512
    assert qa['M_v']==5081 and qa['M_w']==8723
    assert qa['L']==4463 and qa['S']==29
    assert qa['zero_M']==126
    assert qa['distance_missing']==3
    strata=build_strata(rows,original)
    assert len(strata)==6
    assert sum(len(s['test']) for s in strata)==168
    assert sum(sum(x['n'] for x in s['test']) for s in strata)==6692
    assert {s['population'] for s in strata}=={'banpaku','nature','ootaka','ootaka2','shiga','toyama'}


def test_never_train_on_heldout_named_site_even_other_years(rows,original):
    strata=build_strata(rows,original)
    for s in strata:
        train,test=select_source_and_test(strata,s['population'])
        assert all(x['population']!=s['population'] for x in train)
        assert all(x['population']==s['population'] for x in test)
        assert len(train)>20 and len(test)>5
        assert all(x['n']>0 and x['distance']>=0 for x in train+test)
        assert not set((x['year'],x['population'],x['line'],x['quad']) for x in train).intersection(
            (x['year'],x['population'],x['line'],x['quad']) for x in test)


def test_four_fixed_models_and_no_extra_hypotheses():
    assert model_specs()=={'BASE':[],'DIST':['distance'],'CONGEN':['other_share'],'BOTH':['distance','other_share']}


def test_same_test_keys_as_prior_heldout_lines(rows,original):
    strata=build_strata(rows,original)
    testkeys={(x['year'],x['population'],x['line'],x['quad']) for s in strata for x in s['test']}
    selected={(s['year'],s['population']) for s in original['selected_latest']}
    parentkeys={(x['year'],x['population'],x['line'],x['quad']) for x in original['predictions'] if (x['year'],x['population']) in selected}
    assert len(testkeys)==len(parentkeys)==168
    assert testkeys==parentkeys


def test_failure_if_parent_support_tampered(rows,original):
    bad=json.loads(json.dumps(original))
    bad['selected_latest'].pop()
    with pytest.raises(ValueError,match='EXACT_SIX_SITES'):
        build_strata(rows,bad)


def test_primary_results_exactly_six_sites_and_168_original_quadrats(rows,original):
    r=evaluate(rows,original,ROOT/'parent/src/predictive.py')
    assert r['main']['n_sites']==6
    assert r['main']['n_quadrats_scored']==168
    assert r['main']['n_plants_scored']==6692
    assert len(r['site_scores'])==6
    assert len(r['predictions'])==168
    assert r['protocol']['new_external_sites'] is False
    assert r['protocol']['whole_site_name_withheld'] is True
    assert set(r['main']['losses'])=={'BASE','DIST','CONGEN','BOTH'}
    assert all(x['test_n']>0 for x in r['site_scores'])


def test_scoring_is_reproducible_and_cannot_overwrite_parent(rows,original):
    from transfer import source_qa
    before=source_qa(rows)
    r1=evaluate(rows,original,ROOT/'parent/src/predictive.py')
    r2=evaluate(rows,original,ROOT/'parent/src/predictive.py')
    assert r1==r2
    assert source_qa(rows)==before
    assert all('p_within_DIST' in x for x in r1['predictions'])
    assert all(x['n']>0 for x in r1['predictions'])


def test_instrumental_same_time_congener_share_uses_M_species_denominator_only(rows,original):
    source=next(x for x in rows if x['n']>0 and x['l']>0)
    expected=(source['l']+source['s'])/(source['n']+source['l']+source['s'])
    assert source['other_share']==expected
    assert source['n']==source['v']+source['w']


def test_training_site_equal_reweighting_is_not_outcome_based(rows,original):
    from transfer import balanced_training
    strata=build_strata(rows,original)
    train,_=select_source_and_test(strata,'shiga')
    revised=balanced_training(train)
    assert len(revised)==len(train)
    assert [r['population'] for r in revised]==[r['population'] for r in train]
    names=sorted(set(x['population'] for x in train))
    counts={name:sum(x['n'] for x in revised if x['population']==name) for name in names}
    assert len(counts)==5
    assert max(counts.values())-min(counts.values())<1e-9
    for old,new in zip(train,revised):
        assert old['other_share']==new['other_share']
        assert old['distance']==new['distance']
        assert abs(new['v']/new['n']-old['v']/old['n'])<1e-14


def test_test_site_color_labels_never_used_to_fit_coefficients(rows,original):
    import numpy as np
    from transfer import parent_module
    strata=build_strata(rows,original)
    train,test=select_source_and_test(strata,'shiga')
    fitter=parent_module(ROOT/'parent/src/predictive.py')
    p=fitter.fit_binomial(train,test,['distance','other_share'],lam=10)
    counterfactual=[dict(q,v=q['w'],w=q['v']) for q in test]
    p_changed_test_labels=fitter.fit_binomial(train,counterfactual,['distance','other_share'],lam=10)
    assert np.array_equal(p,p_changed_test_labels)


def test_source_archive_tamper_rejected(rows,original,tmp_path):
    from transfer import load_original_verified
    data=(ROOT/'parent/source/Sisyrinchium_data.zip').read_bytes()
    broken=tmp_path/'broken.zip';broken.write_bytes(data[:-1]+b'x')
    with pytest.raises(ValueError,match='ORIGINAL_SOURCE_IDENTITY_FAILED'):
        load_original_verified(broken,ROOT/'parent/src/predictive.py')


def test_independent_original_file_and_newton_matches(rows,original):
    from importlib.util import spec_from_file_location,module_from_spec
    spec=spec_from_file_location('independent',ROOT/'src/verify_independent.py')
    m=module_from_spec(spec);spec.loader.exec_module(m)
    record=m.verify(ROOT/'parent/source/Sisyrinchium_data.zip', ROOT/'parent/results/result.json', ROOT/'results/result.json')
    assert record['status']=='PASS'
    assert record['compared_newton_results']==1776
    assert record['max_abs_difference']<1e-9

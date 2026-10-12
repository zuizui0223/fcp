import importlib.util
import json
import math
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('joint_bound', ROOT / 'src' / 'joint_bound.py')
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)


def test_zero_error_recovers_preexisting_sampling_only_pr174():
    parent = json.loads((ROOT/'data'/'parent_entropy_result.json').read_text())
    nums = {x['year']: x['rounding_bounds']['H_between'][0] for x in parent['hesperis']}
    expected = {2021: 5.7051216422e-8, 2022: 4.7620579400e-7, 2023: 9.1683820904e-7}
    for yr, b in nums.items():
        assert M.joint_family_upper(b, 100, 0.0) == pytest.approx(expected[yr], rel=2e-5)


def test_conservative_joint_bound_monotone_and_is_not_independent_p_multiplication():
    b=.22496448574712613
    values = [M.joint_family_upper(b,100,e) for e in (0, .005, .01, .02, .025, .03)]
    assert values == sorted(values)
    assert values[2] > .0001
    assert values[3] > .005


def test_correction_for_failure_of_error_budget_adds_probability_not_multiplies():
    b=.22496448574712613
    base = M.joint_family_upper(b, 100, .01, eta=0)
    risky = M.joint_family_upper(b, 100, .01, eta=.01)
    assert risky == pytest.approx(base+.03)


def test_critical_eps_threshold_and_null_status():
    b=.21776296903477543
    cap = M.max_error_for_family_alpha(b,100,alpha=.05)
    assert .023 < cap < .024
    assert M.joint_family_upper(b,100,cap-.000001) < .05
    assert M.joint_family_upper(b,100,cap+.000001) > .05
    assert M.max_error_for_family_alpha(b,25,alpha=.05) is None


def test_jensen_shannon_recorded_coupling_respects_fano_for_adversarial_stand_specific_errors():
    import random
    rng=random.Random(20261012)
    for rep in range(800):
        true = [[rng.randrange(3) for _ in range(rng.randint(1,30))] for _ in range(3)]
        obs = []
        for g, labels in enumerate(true):
            q=[]
            for z in labels:
                if rng.random() < .35*(g+1)/3:
                    q.append((z+1+rng.randrange(2))%3)
                else:q.append(z)
            obs.append(q)
        e=sum(sum(a!=b for a,b in zip(x,y))/len(x) for x,y in zip(true,obs))/3
        I_x=M.empirical_stand_info(true)
        I_z=M.empirical_stand_info(obs)
        assert I_z <= I_x+M.fano_ceiling(e)+1e-12


def test_complete_census_and_sampling_semantics_are_not_collapsed():
    assert M.sample_null_upper(.22, 100) < M.sample_null_upper(.22,50)
    assert M.sample_null_upper(0,100)==1
    assert M.sample_null_upper(.2,0)==1
    assert M.joint_family_upper(.2,100,.06)==1


def test_input_validation_and_source_lock():
    for eps in [-.01,1.01]:
        with pytest.raises(ValueError):M.fano_ceiling(eps)
    with pytest.raises(ValueError):M.joint_family_upper(.23,100,.01,eta=-.1)
    original=ROOT/'data'/'parent_entropy_result.json'
    assert M.verify_parent(original)['n_annual']==3
    with pytest.raises(ValueError):M.verify_parent(ROOT/'tests'/'test_joint_bound.py')


def test_result_zero_error_regression_and_three_year_support():
    got=M.analyze(ROOT/'data'/'parent_entropy_result.json')
    assert got['null_hypothesis']=='same_true_finite_census_colour_composition_all_stands'
    assert got['sample_size_semantics']=='census_or_simple_random_without_replacement_n_ge_n0'
    assert got['n_years']==3
    assert got['annual'][2021]['critical_error_n100'] > got['annual'][2023]['critical_error_n100']
    for year in (2021,2022,2023):
        assert got['annual'][year]['scenario'][100]['0.0200']['adjusted_bound'] < .05
        assert got['annual'][year]['scenario'][100]['0.0300']['adjusted_bound'] > .05 if year>=2022 else True
        assert got['annual'][year]['scenario'][25]['0.0000']['adjusted_bound'] >= .05


def test_ordinary_fano_budget_is_not_error_validation():
    got=M.analyze(ROOT/'data'/'parent_entropy_result.json')
    assert got['measured_misclassification_rate'] is None
    assert got['raw_stand_counts_recovered'] is False
    assert got['inference_status']=='CONDITIONAL_JOINT_SAMPLING_AND_CLASSIFICATION_BUDGET'


def test_all_scenario_keys_round_trip_exact_prospective_error_grid():
    got=M.analyze(ROOT/'data'/'parent_entropy_result.json')
    for y in (2021,2022,2023):
        for n in M.SAMPLE_MINIMA:
            items=got['annual'][y]['scenario'][n]
            assert len(items)==len(M.ERROR_RATES)
            for key,item in items.items():
                assert float(key)==item['epsilon_upper'], ('misleading epsilon key',key,item)
            assert '0.0025' in items


def test_joint_event_inclusion_under_enumerated_without_replacement_and_adversarial_errors():
    from itertools import combinations
    census=[0,0,1,1,2,2]
    all_samples=[list(s) for s in combinations(range(len(census)),3)]
    examples=0
    for ii in all_samples:
        for jj in all_samples:
            X=[[census[k] for k in ii],[census[k] for k in jj]]
            Z=[x.copy() for x in X]
            # A site-dependent recording artifact in just the second stand.
            Z[1][0]=(Z[1][0]+1)%3
            e=(sum(a!=b for a,b in zip(X[0],Z[0]))/3+
               sum(a!=b for a,b in zip(X[1],Z[1]))/3)/2
            observed=M.empirical_stand_info(Z)
            true_sample=M.empirical_stand_info(X)
            assert observed<=true_sample+M.fano_ceiling(e)+1e-12
            # This is event inclusion, not a bound on the conditional probability
            # P(BZ>=b | realized_error<=epsilon).
            if observed>=.15 and e<=.18:
                assert true_sample >= .15-M.fano_ceiling(.18)-1e-12
                examples+=1
    assert examples>0


def test_eta_union_exposes_probability_budget_failure_without_requiring_independence():
    b=.21776296903477543
    p0=M.joint_year_upper(b,100,.02,eta=0)
    p1=M.joint_year_upper(b,100,.02,eta=.0125)
    assert p1==pytest.approx(p0+.0125)
    assert M.joint_family_upper(b,100,.02,eta=.0125)==pytest.approx(3*(p0+.0125))

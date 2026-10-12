import json, math, pathlib, sys
import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from information_bound import fano_entropy, mutual_info, apply_stand_channels, safe_threshold, analyze_parent, verify_parent


def test_zero_and_saturation():
    assert fano_entropy(0.0,3) == 0
    assert abs(fano_entropy(2/3,3)-math.log(3)) < 1e-13
    assert fano_entropy(0.8,3) == pytest.approx(math.log(3))  # e is an UPPER bound
    with pytest.raises(ValueError): fano_entropy(-0.01,3)
    with pytest.raises(ValueError): fano_entropy(0.5,1)


def test_common_confusion_cannot_increase_information():
    weights=[0.2, 0.8]
    true=[[0.90,0.08,0.02],[0.12,0.20,0.68]]
    channel=[[0.6,0.3,0.1],[0.1,0.6,0.3],[0.25,0.25,0.5]]
    observed=apply_stand_channels(true,[channel,channel])
    assert mutual_info(weights,observed) <= mutual_info(weights,true)+1e-13


def test_stand_specific_confusion_can_generate_information_but_fano_bounds_it():
    true=[[0.4,0.3,0.3],[0.4,0.3,0.3]]
    channels=[[[1,0,0],[0,1,0],[0,0,1]], [[0.9,0.1,0],[0,1,0],[0,0,1]]]
    observed=apply_stand_channels(true,channels)
    p=mutual_info([.5,.5],observed)
    e=.5*.4*.1
    assert p>0 and mutual_info([.5,.5],true)==pytest.approx(0,abs=1e-13)
    assert p <= fano_entropy(e,3)+1e-13


def test_asymmetric_but_low_error_recovers_positive_lower_bound():
    # Raw publication error rate is NOT observed; this only checks a scenario.
    p=0.25386012428895355
    assert p - fano_entropy(.02,3) > 0
    assert p - fano_entropy(.1,3) < 0


def test_threshold_solves_one_entropy_bound_not_two():
    for b in [.24814489869438183,.22496448574712613,.21776296903477543]:
        t=safe_threshold(b,3)
        assert .04<t<.07
        assert abs(fano_entropy(t,3)-b)<1e-12
        assert fano_entropy(t-1e-6,3)<b<fano_entropy(t+1e-6,3)


def test_summary_verification_and_scope():
    source=ROOT/'data'/'parent_entropy_result.json'
    verified=verify_parent(source)
    assert verified['status']=='PASS'
    assert verified['n_years']==3
    assert verified['source_not_raw_stand_counts']
    output=analyze_parent(source,[0,.005,.01,.02,.03,.04,.05])
    assert output['scope']=='DETERMINISTIC_LABEL_ERROR_SENSITIVITY_NOT_FIELD_DATA'
    assert output['verified_empirical_error_rate'] is None
    assert len(output['annual'])==3
    for year in output['annual']:
        assert year['B_original_rounding_min']>0
        assert year['threshold_error_rate']>0
        for scenario in year['scenarios']:
            assert 0 <= scenario['error_rate_upper'] <= .05
            assert 0 <= scenario['guaranteed_true_sample_B_min'] <= year['B_original_rounding_min']
        assert year['scenarios'][0]['guaranteed_true_sample_B_min']==pytest.approx(year['B_original_rounding_min'])


def test_refuse_non_parent_provenance(tmp_path):
    q=tmp_path/'fake.json';q.write_text(json.dumps({'schema':'other','hesperis':[]}))
    with pytest.raises(ValueError,match='provenance'):verify_parent(q)


def test_general_information_bound_under_400_variable_stand_channels():
    import random
    rng = random.Random(20261012)
    for _ in range(400):
        n = rng.randint(2, 5)
        weights = [rng.random()+.01 for _ in range(n)]
        weights = [v/sum(weights) for v in weights]
        true = []
        channels = []
        for _ in range(n):
            z=[rng.random() for _ in range(3)];true.append([v/sum(z) for v in z])
            C=[]
            for i in range(3):
                c=[rng.random()+(.01 if j==i else 0) for j in range(3)]
                C.append([v/sum(c) for v in c])
            channels.append(C)
        observed=apply_stand_channels(true,channels)
        error=sum(weights[g]*sum(true[g][i]*(1-channels[g][i][i]) for i in range(3)) for g in range(n))
        assert mutual_info(weights,observed)-mutual_info(weights,true) <= fano_entropy(error,3)+1e-12
        assert mutual_info(weights,true)-mutual_info(weights,observed) <= fano_entropy(error,3)+1e-12

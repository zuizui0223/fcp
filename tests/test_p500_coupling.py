import numpy as np
import pytest

from fcp_pipeline.p500_coupling import fit_coupling, prepare_coupling_inputs


def fixture():
    mapping = {1: 'A', 2: 'A', 3: 'B', 4: 'B', 5: 'C', 6: 'C', 7: 'D'}
    values = [.1, .3, .5, .9, .2, .2, .3]
    states = [False, True, False, True, False, True, True]
    rows = [dict(photo_id=i, species=s, white=states[i-1], near_clip_fraction=values[i-1]) for i,s in mapping.items()]
    return rows, mapping


def test_species_standardization_and_complete_exclusion_audit():
    rows, mapping = fixture()
    result = prepare_coupling_inputs(rows, mapping)
    assert [r['photo_id'] for r in result['rows']] == [1, 2, 3, 4]
    assert np.allclose([r['near_clip_z'] for r in result['rows']], [-2**-.5, 2**-.5]*2)
    assert result['species_audit'][2]['reasons'] == ['zero_predictor_variance']
    assert result['species_audit'][3]['reasons'] == ['single_response_state', 'zero_predictor_variance']
    assert sum(a['n'] for a in result['species_audit']) == len(mapping)
    assert result['model_fitted'] is result['opening_authorized'] is False
    assert result == prepare_coupling_inputs(rows[::-1], mapping)


@pytest.mark.parametrize('mutation', ['missing', 'duplicate', 'species', 'response', 'nan', 'extra'])
def test_invalid_input_fails_instead_of_silent_dropping(mutation):
    rows, mapping = fixture()
    if mutation == 'missing': rows.pop()
    elif mutation == 'duplicate': rows.append(rows[0])
    elif mutation == 'species': rows[0]['species'] = 'B'
    elif mutation == 'response': rows[0]['white'] = None
    elif mutation == 'nan': rows[0]['near_clip_fraction'] = float('nan')
    elif mutation == 'extra': rows[0]['latitude'] = 1
    with pytest.raises(ValueError):
        prepare_coupling_inputs(rows, mapping)


def matched_pairs(positive, negative):
    rows, mapping = [], {}
    for pair in range(positive + negative):
        for side in (0, 1):
            identity = 2 * pair + side + 1
            species = f'synthetic-{pair:03}'
            mapping[identity] = species
            rows.append(dict(photo_id=identity, species=species,
                             near_clip_fraction=.1 + .2 * side,
                             white=bool(side if pair < positive else 1 - side)))
    return rows, mapping


@pytest.mark.parametrize('positive,negative', [(75, 25), (50, 50), (25, 75)])
def test_fit_matches_analytic_matched_pairs_and_order_invariance(positive, negative):
    rows, mapping = matched_pairs(positive, negative)
    result = fit_coupling(rows, mapping)
    assert result['status'] == 'QUALIFIED_NUMERIC_FIT', result
    expected_beta = np.log(positive / negative) / np.sqrt(2)
    probability = positive / (positive + negative)
    expected_information = 2 * (positive + negative) * probability * (1 - probability)
    expected_se = expected_information ** -.5
    assert result['beta'] == pytest.approx(expected_beta, abs=1e-9)
    assert result['information'] == pytest.approx(expected_information, rel=1e-5)
    assert result['ci95'] == pytest.approx(np.exp([expected_beta - 1.959963984540054 * expected_se,
                                                 expected_beta + 1.959963984540054 * expected_se]), rel=1e-5)
    assert result['opening_authorized'] is result['execution_chronology_verified'] is False
    assert result == fit_coupling(rows[::-1], dict(reversed(list(mapping.items()))))


@pytest.mark.parametrize('positive,negative', [(10, 0), (0, 10)])
def test_separation_has_no_estimate(positive, negative):
    result = fit_coupling(*matched_pairs(positive, negative))
    assert result['status'] == 'INDETERMINATE'
    assert result['reason'] == 'conditional_support_boundary_or_near_boundary'
    assert result['model_fitted'] is False
    assert 'odds_ratio' not in result


def test_no_informative_species_has_no_model():
    rows, mapping = matched_pairs(5, 5)
    for row in rows:
        row['white'] = True
    result = fit_coupling(rows, mapping)
    assert result['reason'] == 'no_informative_species'
    assert len(result['species_audit']) == 10


def test_runtime_drift_fails_closed(monkeypatch):
    import scipy
    monkeypatch.setattr(scipy, '__version__', 'unqualified')
    assert fit_coupling(*matched_pairs(5, 5))['reason'] == 'unqualified_runtime'


def test_numeric_warning_has_no_estimate(monkeypatch):
    import warnings
    from statsmodels.discrete.conditional_models import ConditionalLogit

    def warn(*args, **kwargs):
        warnings.warn('Synthetic numerical failure', RuntimeWarning)

    monkeypatch.setattr(ConditionalLogit, 'score', warn)
    result = fit_coupling(*matched_pairs(5, 5))
    assert result['reason'] == 'numeric_failure:RuntimeWarning'
    assert result['model_fitted'] is False


def test_invalid_full_cohort_cannot_be_fitted():
    rows, mapping = matched_pairs(5, 5)
    with pytest.raises(ValueError, match='full expected cohort'):
        fit_coupling(rows[:-1], mapping)


def test_multiphoto_strata_match_independent_enumerated_likelihood():
    from itertools import combinations
    from scipy.optimize import brentq
    from scipy.special import logsumexp, softmax

    rows, mapping = [], {}
    cases = [([.1, .2, .4, .9], [0, 0, 1, 1]),
             ([.0, .3, .5, .8, 1.], [1, 0, 1, 0, 0]),
             ([.1, .2, .5, .9], [0, 1, 0, 0])]
    supports, observed = [], 0.0
    for group, (values, outcomes) in enumerate(cases):
        z = (np.array(values) - np.mean(values)) / np.std(values, ddof=1)
        supports.append(np.array([sum(c) for c in combinations(z, sum(outcomes))]))
        observed += float(z @ outcomes)
        for value, outcome in zip(values, outcomes):
            identity = len(rows) + 1
            mapping[identity] = str(group)
            rows.append(dict(photo_id=identity, species=str(group), white=bool(outcome),
                             near_clip_fraction=value))

    def independent_score(beta):
        return observed - sum(float(softmax(beta * support) @ support) for support in supports)

    beta = brentq(independent_score, -5, 5)
    information = sum(float(softmax(beta * s) @ (s ** 2) - (softmax(beta * s) @ s) ** 2)
                      for s in supports)
    loglike = beta * observed - sum(logsumexp(beta * s) for s in supports)
    result = fit_coupling(rows, mapping)
    assert result['status'] == 'QUALIFIED_NUMERIC_FIT', result
    assert result['beta'] == pytest.approx(beta, abs=1e-9)
    assert result['information'] == pytest.approx(information, rel=1e-5)
    assert result['loglike'] == pytest.approx(loglike, abs=1e-9)


def test_tied_support_boundary_is_not_rescued():
    values, outcomes = [.1, .1, .3, .3], [0, 0, 1, 1]
    mapping = {i: 'tied' for i in range(1, 5)}
    rows = [dict(photo_id=i, species='tied', white=bool(y), near_clip_fraction=x)
            for i, (x, y) in enumerate(zip(values, outcomes), 1)]
    assert fit_coupling(rows, mapping)['reason'] == 'conditional_support_boundary_or_near_boundary'


def test_bad_curvature_is_not_rescued(monkeypatch):
    from statsmodels.discrete.conditional_models import ConditionalLogit
    monkeypatch.setattr(ConditionalLogit, 'hessian', lambda *args: np.array([[1.0]]))
    result = fit_coupling(*matched_pairs(5, 5))
    assert result['reason'] == 'curvature_or_likelihood_not_qualified'
    assert 'odds_ratio' not in result

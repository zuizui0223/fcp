import numpy as np
import pytest

from fcp_pipeline.p500_coupling import prepare_coupling_inputs


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

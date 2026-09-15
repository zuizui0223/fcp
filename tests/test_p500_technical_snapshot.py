import hashlib
import json

import pytest

from fcp_pipeline.p500_technical_snapshot import seal_high_clip_snapshot, verify_high_clip_snapshot


def rows():
    return [{'photo_id': i + 1, 'near_clip_fraction': i / 100} for i in range(20)]


def test_round_trip_and_order_invariance():
    data = rows()
    encoded, digest = seal_high_clip_snapshot(data, range(1, 21))
    assert seal_high_clip_snapshot(data[::-1], range(20, 0, -1)) == (encoded, digest)
    assert verify_high_clip_snapshot(encoded, digest, range(1, 21)) == (20,)
    payload = json.loads(encoded)
    assert payload['opening_authorized'] is False
    assert payload['execution_chronology_verified'] is False


@pytest.mark.parametrize('mutation', ['extra_response', 'missing', 'duplicate', 'nan', 'bool', 'outside'])
def test_reject_invalid_or_response_bearing_rows(mutation):
    data = rows()
    if mutation == 'extra_response': data[0]['white'] = True
    elif mutation == 'missing': data.pop()
    elif mutation == 'duplicate': data[-1]['photo_id'] = 1
    elif mutation == 'nan': data[0]['near_clip_fraction'] = float('nan')
    elif mutation == 'bool': data[0]['near_clip_fraction'] = True
    elif mutation == 'outside': data[0]['near_clip_fraction'] = 1.1
    with pytest.raises(ValueError):
        seal_high_clip_snapshot(data, range(1, 21))


def test_tampering_and_reissued_digest_do_not_override_membership():
    encoded, digest = seal_high_clip_snapshot(rows(), range(1, 21))
    payload = json.loads(encoded)
    payload['high_clip_photo_ids'] = []
    changed = json.dumps(payload, sort_keys=True, separators=(',', ':')).encode()
    with pytest.raises(ValueError):
        verify_high_clip_snapshot(changed, digest, range(1, 21))
    with pytest.raises(ValueError):
        verify_high_clip_snapshot(changed, hashlib.sha256(changed).hexdigest(), range(1, 21))


def test_exact_threshold_is_not_excluded_and_identity_change_rejected():
    encoded, digest = seal_high_clip_snapshot([{'photo_id': 1, 'near_clip_fraction': .01}], [1])
    assert verify_high_clip_snapshot(encoded, digest, [1]) == ()
    with pytest.raises(ValueError):
        verify_high_clip_snapshot(encoded, digest, [2])

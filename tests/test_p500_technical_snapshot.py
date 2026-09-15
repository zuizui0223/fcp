import hashlib
import json

import pytest

from fcp_pipeline.p500_technical_snapshot import seal_high_clip_snapshot, verify_high_clip_snapshot
from fcp_pipeline.p500_technical_snapshot import write_high_clip_snapshot, join_verified_technical_snapshot


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


def test_exclusive_writer_roundtrip_and_no_overwrite(tmp_path):
    encoded, digest = seal_high_clip_snapshot(rows(), range(1, 21))
    target = tmp_path / 'technical.json'
    write_high_clip_snapshot(target, encoded, digest, range(1, 21))
    assert target.read_bytes() == encoded
    with pytest.raises(FileExistsError):
        write_high_clip_snapshot(target, encoded, digest, range(1, 21))
    assert target.read_bytes() == encoded
    with pytest.raises(ValueError):
        write_high_clip_snapshot(tmp_path / 'invalid.json', encoded, 'invalid', range(1, 21))
    assert not (tmp_path / 'invalid.json').exists()


def test_join_retains_all_rows_and_fixed_membership():
    encoded, digest = seal_high_clip_snapshot(rows(), range(1, 21))
    response = [{'photo_id': i, 'white': i % 2 == 0} for i in range(20, 0, -1)]
    joined = join_verified_technical_snapshot(encoded, digest, range(1, 21), response)
    assert len(joined) == 20
    assert [r['photo_id'] for r in joined if r['high_clip']] == [20]
    assert all(r['white'] == (r['photo_id'] % 2 == 0) for r in joined)


@pytest.mark.parametrize('mutation', ['missing', 'duplicate', 'unknown', 'score', 'extra'])
def test_response_join_refuses_partial_or_ambiguous_data(mutation):
    encoded, digest = seal_high_clip_snapshot(rows(), range(1, 21))
    response = [{'photo_id': i, 'white': False} for i in range(1, 21)]
    if mutation == 'missing': response.pop()
    elif mutation == 'duplicate': response[-1]['photo_id'] = 1
    elif mutation == 'unknown': response[-1]['photo_id'] = 30
    elif mutation == 'score': response[0]['white'] = 0
    elif mutation == 'extra': response[0]['species'] = 'unfrozen'
    with pytest.raises(ValueError):
        join_verified_technical_snapshot(encoded, digest, range(1, 21), response)


def test_bad_snapshot_is_rejected_before_response_iteration():
    encoded, digest = seal_high_clip_snapshot(rows(), range(1, 21))
    def forbidden_response_read():
        raise AssertionError('Response iterator must not be opened')
        yield
    with pytest.raises(ValueError):
        join_verified_technical_snapshot(encoded, 'invalid', range(1, 21), forbidden_response_read())

"""Technical-only snapshots; integrity checks are not opening authorization.

Pure in-memory operations. No files, images, network, colour responses or
historical timestamps are read or written. Trusted digest custody is external.
"""
import hashlib
import json
import math
from collections.abc import Mapping
from numbers import Real

from .p500_white_measurement_control import response_blind_high_clip_mask


def _canonical(records, expected_photo_ids):
    expected = list(expected_photo_ids)
    if not expected or any(type(i) is not int or i <= 0 for i in expected):
        raise ValueError('Expected photo IDs must be positive integers')
    if len(expected) != len(set(expected)):
        raise ValueError('Duplicate expected photo IDs')
    rows = []
    for record in records:
        if not isinstance(record, Mapping) or set(record) != {'photo_id', 'near_clip_fraction'}:
            raise ValueError('Only photo_id and near_clip_fraction are allowed')
        identity, value = record['photo_id'], record['near_clip_fraction']
        if type(identity) is not int or identity <= 0:
            raise ValueError('Invalid photo ID')
        if isinstance(value, bool) or not isinstance(value, Real):
            raise ValueError('Invalid technical fraction type')
        value = float(value)
        if not math.isfinite(value) or not 0 <= value <= 1:
            raise ValueError('Invalid technical fraction range')
        rows.append({'photo_id': identity, 'near_clip_fraction': value})
    ids = [r['photo_id'] for r in rows]
    if len(ids) != len(set(ids)) or set(ids) != set(expected):
        raise ValueError('Technical rows must exactly match expected identities')
    return sorted(rows, key=lambda r: r['photo_id'])


def seal_high_clip_snapshot(records, expected_photo_ids):
    """Return canonical bytes and SHA256 for all supplied evaluable identities.

    The caller must establish and record the evaluable identity set separately;
    this function neither chooses the cohort nor permits dropping failed rows.
    """
    rows = _canonical(records, expected_photo_ids)
    mask, threshold = response_blind_high_clip_mask([r['near_clip_fraction'] for r in rows])
    payload = {
        'schema': 'p500-technical-snapshot-v1',
        'rule': 'near_clip_fraction > max(0.01, q95_linear)',
        'rows': rows,
        'threshold': threshold,
        'high_clip_photo_ids': [r['photo_id'] for r, high in zip(rows, mask) if high],
        'opening_authorized': False,
        'execution_chronology_verified': False,
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()
    return encoded, hashlib.sha256(encoded).hexdigest()


def verify_high_clip_snapshot(encoded, trusted_sha256, expected_photo_ids):
    """Validate exact bytes and recompute membership; never authorize opening.

    trusted_sha256 must be supplied from a separately preserved pre-response
    record, not accepted from the same untrusted payload. This does not establish
    when that record was created or whether responses were already inspected.
    """
    if type(encoded) is not bytes or hashlib.sha256(encoded).hexdigest() != trusted_sha256:
        raise ValueError('Snapshot digest mismatch')
    try:
        payload = json.loads(encoded)
        rebuilt, _ = seal_high_clip_snapshot(payload['rows'], expected_photo_ids)
    except (TypeError, KeyError, json.JSONDecodeError) as exc:
        raise ValueError('Malformed snapshot') from exc
    if rebuilt != encoded:
        raise ValueError('Snapshot content does not match frozen technical rule')
    return tuple(payload['high_clip_photo_ids'])

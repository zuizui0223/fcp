"""Technical-only snapshots; integrity checks are not opening authorization.

No images, network or historical timestamps are read. The optional writer uses
exclusive creation. Joining accepts already supplied response rows only after
snapshot integrity checks. Trusted digest custody and authorization are external.
"""
import hashlib
import json
import math
import os
from pathlib import Path
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


def write_high_clip_snapshot(path, encoded, trusted_sha256, expected_photo_ids):
    """Exclusively create a validated snapshot; never overwrite an existing path.

    Failure during writing may leave a partial file, which must remain a failure
    record and will not pass digest verification. This is not WORM storage and
    cannot prevent other programs from rewriting files. No chronology is inferred.
    """
    verify_high_clip_snapshot(encoded, trusted_sha256, expected_photo_ids)
    with Path(path).open('xb') as stream:
        stream.write(encoded)
        stream.flush()
        os.fsync(stream.fileno())


def join_verified_technical_snapshot(encoded, trusted_sha256, expected_photo_ids, responses):
    """One-to-one complete join, using fixed exclusion membership only.

    This is an integrity operation, not permission to read a biological response.
    No loader, image path, download or real-data entry point is provided here.
    Missing response measurements must be handled by a separately frozen
    attrition contract; silently joining only complete cases is forbidden.
    """
    expected = tuple(expected_photo_ids)
    excluded = set(verify_high_clip_snapshot(encoded, trusted_sha256, expected))
    supplied = {}
    for response in responses:
        if not isinstance(response, Mapping) or set(response) != {'photo_id', 'white'}:
            raise ValueError('Only photo_id and boolean white response are allowed')
        identity = response['photo_id']
        if type(identity) is not int or identity <= 0 or identity in supplied:
            raise ValueError('Invalid or duplicate response identity')
        if type(response['white']) is not bool:
            raise ValueError('Response must be a boolean, not a missing value or score')
        supplied[identity] = response['white']
    if set(supplied) != set(expected):
        raise ValueError('Response identities must exactly match the technical snapshot')
    return [dict(row, white=supplied[row['photo_id']],
                 high_clip=row['photo_id'] in excluded)
            for row in json.loads(encoded)['rows']]

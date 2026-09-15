"""Pre-response coupling implementation: strict inputs and audited strata only.

No real-data entry point or opening authorization. Input construction is not
model fitting, calibration, or independent biological validation.
"""
from collections.abc import Mapping
from numbers import Real

import numpy as np


def prepare_coupling_inputs(records, expected_species_by_photo):
    """Validate a complete cohort; sample-SD standardize within eligible species.

    Expected identities and species mapping must come from an independently
    frozen source. Exclusion reasons are returned for every non-informative
    species. Missing individual rows are errors, not silent complete-case drops.
    """
    if not isinstance(expected_species_by_photo, Mapping) or not expected_species_by_photo:
        raise ValueError('Expected cohort mapping is required')
    for identity, species in expected_species_by_photo.items():
        if type(identity) is not int or identity <= 0 or type(species) is not str or not species.strip():
            raise ValueError('Invalid expected identity or species')
    groups, seen = {}, set()
    for row in records:
        if not isinstance(row, Mapping) or set(row) != {'photo_id', 'species', 'white', 'near_clip_fraction'}:
            raise ValueError('Unexpected coupling fields')
        identity = row['photo_id']
        if type(identity) is not int or identity not in expected_species_by_photo or identity in seen:
            raise ValueError('Missing, duplicate or unexpected identity')
        if row['species'] != expected_species_by_photo[identity]:
            raise ValueError('Species differs from frozen mapping')
        if type(row['white']) is not bool:
            raise ValueError('White response must be boolean')
        value = row['near_clip_fraction']
        if isinstance(value, bool) or not isinstance(value, Real) or not np.isfinite(value) or not 0 <= value <= 1:
            raise ValueError('Invalid near-clip fraction')
        seen.add(identity)
        groups.setdefault(row['species'], []).append((identity, float(value), int(row['white'])))
    if seen != set(expected_species_by_photo):
        raise ValueError('Coupling input does not cover the full expected cohort')
    prepared, audit = [], []
    for species, items in sorted(groups.items()):
        items.sort()
        x = np.array([item[1] for item in items])
        y = np.array([item[2] for item in items])
        reasons = []
        if len(set(y)) != 2:
            reasons.append('single_response_state')
        sd = float(np.std(x, ddof=1)) if len(x) > 1 else 0.0
        if sd == 0:
            reasons.append('zero_predictor_variance')
        if not np.isfinite(sd):
            raise ValueError('Nonfinite within-species SD')
        audit.append({'species': species, 'n': len(items), 'included': not reasons, 'reasons': reasons})
        if not reasons:
            z = (x - x.mean()) / sd
            for item, standardized in zip(items, z):
                prepared.append({'photo_id': item[0], 'species': species,
                                 'white': item[2], 'near_clip_z': float(standardized)})
    return {'rows': prepared, 'species_audit': audit, 'sd_ddof': 1,
            'model_fitted': False, 'opening_authorized': False,
            'execution_chronology_verified': False}

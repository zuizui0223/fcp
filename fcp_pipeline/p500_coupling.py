"""Pre-response coupling implementation: strict inputs and audited strata only.

No real-data entry point or opening authorization. Input construction is not
model fitting, calibration, or independent biological validation.
"""
from collections.abc import Mapping
from numbers import Real
import warnings

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


def fit_coupling(records, expected_species_by_photo):
    """Fit one conditional slope, or return INDETERMINATE without an estimate.

    Artificial-data qualified only. This pure function does not grant access to
    responses. Invalid cohort inputs raise ValueError before model construction.
    Numeric conventions are fixed in P500_COUPLING_IMPLEMENTATION_SPEC.md.
    """
    import scipy
    import statsmodels
    from scipy.optimize import brentq
    from statsmodels.discrete.conditional_models import ConditionalLogit

    prepared = prepare_coupling_inputs(records, expected_species_by_photo)
    rows = prepared['rows']
    result = {key: value for key, value in prepared.items() if key != 'rows'}
    result.update(status='INDETERMINATE', n_included=len(rows),
                  runtime={'scipy': scipy.__version__, 'statsmodels': statsmodels.__version__})

    def fail(reason):
        return dict(result, reason=reason)

    if scipy.__version__ != '1.15.3' or statsmodels.__version__ != '0.14.6':
        return fail('unqualified_runtime')
    if not rows:
        return fail('no_informative_species')
    x = np.array([r['near_clip_z'] for r in rows])
    y = np.array([r['white'] for r in rows])
    groups = np.array([r['species'] for r in rows])
    lower, upper, observed = 0.0, 0.0, float(y @ x)
    for group in sorted(set(groups)):
        selected = groups == group
        values = np.sort(x[selected])
        count = int(y[selected].sum())
        lower += float(values[:count].sum())
        upper += float(values[-count:].sum())
    tolerance = 1e-10 * max(1.0, abs(lower), abs(upper))
    result['support'] = dict(lower=lower, observed=observed, upper=upper, tolerance=tolerance)
    if observed - lower <= tolerance or upper - observed <= tolerance:
        return fail('conditional_support_boundary_or_near_boundary')
    try:
        with warnings.catch_warnings(), np.errstate(all='raise'):
            warnings.simplefilter('error')
            model = ConditionalLogit(y, x[:, None], groups=groups, missing='raise')

            def score(beta):
                value = float(model.score(np.array([beta]))[0])
                if not np.isfinite(value):
                    raise FloatingPointError('Nonfinite score')
                return value

            bracket = None
            for bound in (1.0, 2.0, 4.0, 8.0, 16.0, 20.0):
                if score(-bound) > 0 and score(bound) < 0:
                    bracket = (-bound, bound)
                    break
            if bracket is None:
                return fail('no_strict_score_bracket_within_fixed_limits')
            beta, convergence = brentq(score, *bracket, xtol=1e-10, rtol=1e-12,
                                       maxiter=100, full_output=True, disp=False)
            residual = score(beta)
            result['optimizer'] = dict(converged=bool(convergence.converged),
                                       function_calls=int(convergence.function_calls),
                                       bracket=list(bracket), score=residual)
            if not convergence.converged or abs(residual) > 1e-7:
                return fail('root_or_score_not_qualified')
            information = float(-model.hessian(np.array([beta]))[0, 0])
            h = 1e-4
            check_information = -(score(beta + h) - score(beta - h)) / (2 * h)
            loglike = float(model.loglike(np.array([beta])))
            if (not np.isfinite([beta, information, check_information, loglike]).all()
                    or information <= 1e-8 or check_information <= 1e-8
                    or not np.isclose(information, check_information, rtol=1e-3, atol=1e-8)):
                return fail('curvature_or_likelihood_not_qualified')
            se = float(1 / np.sqrt(information))
            interval = np.exp([beta - 1.959963984540054 * se,
                               beta + 1.959963984540054 * se])
            odds_ratio = float(np.exp(beta))
            if not np.isfinite(interval).all() or (interval <= 0).any():
                return fail('nonfinite_or_zero_interval')
    except (Warning, FloatingPointError, OverflowError, ValueError, RuntimeError) as exc:
        # No penalized fallback, hidden exclusion, or second optimizer.
        return fail('numeric_failure:' + type(exc).__name__)
    return dict(result, status='QUALIFIED_NUMERIC_FIT', model_fitted=True,
                beta=float(beta), standard_error=se, odds_ratio=odds_ratio,
                ci95=interval.tolist(), information=information,
                check_information=check_information, loglike=loglike,
                covariance='model_based_not_observer_cluster_robust')

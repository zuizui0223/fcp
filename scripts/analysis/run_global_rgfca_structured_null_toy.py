#!/usr/bin/env python3
"""Frozen ring proof-of-concept; never reads or reclassifies observed colours.

Specification: docs/supporting/global_rgfca_structured_null_toy_contract_v1.json
Frozen in commit b86a5967fd53d998e0915070e107e5f870ce0ec9 before these outputs.
Run: python scripts/analysis/run_global_rgfca_structured_null_toy.py --output-dir OUT
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import numpy as np

PROTOCOL = 'global-rgfca-structured-null-specificity-toy-v1'
SPEC_COMMIT = 'b86a5967fd53d998e0915070e107e5f870ce0ec9'
SEED, SPECIES, CELLS, WIDTH = 2026090615, 40, 64, 3.0
CALIBRATION, EVALUATION, BATCH = 2000, 1000, 100
AMPLITUDES = (0.0, 0.5, 1.0, 2.0)
FRACTIONS = (0.0, 0.1, 0.25, 0.5, 1.0)
NOISE_SD = 1.0


def wilson(k: int, n: int) -> tuple[float, float]:
    """95% interval for MC rejection probability, conditional on calibration."""
    if n <= 0 or not 0 <= k <= n:
        raise ValueError('require 0 <= k <= n and n > 0')
    z = 1.959963984540054
    p = k / n
    den = 1 + z*z/n
    center = (p + z*z/(2*n)) / den
    half = z * math.sqrt(p*(1-p)/n + z*z/(4*n*n)) / den
    return max(0.0, center-half), min(1.0, center+half)


def shape() -> np.ndarray:
    d = np.minimum(np.arange(CELLS), CELLS-np.arange(CELLS))
    x = np.exp(-0.5*(d/WIDTH)**2)
    return (x-x.mean()) / x.std(ddof=0)


def generate(rng: np.random.Generator, n: int, amplitude: float, q: int) -> np.ndarray:
    """Species share a phase, not the independent observational noise."""
    if n <= 0 or not 0 <= q <= SPECIES or amplitude < 0:
        raise ValueError('invalid simulation dimensions/amplitude')
    phases = rng.integers(0, CELLS, size=(n, SPECIES))
    common = rng.integers(0, CELLS, size=n)
    phases[:, :q] = common[:, None]
    indices = (np.arange(CELLS)[None, None, :]-phases[:, :, None]) % CELLS
    x = amplitude*shape()[indices] + NOISE_SD*rng.normal(size=(n, SPECIES, CELLS))
    return x - x.mean(axis=-1, keepdims=True)


def decomposition(x: np.ndarray) -> np.ndarray:
    """Return total, self, and cross; uniform cells and fixed equal species weights."""
    x = np.asarray(x, dtype=np.float64)
    if x.ndim != 3 or x.shape[1] < 2 or x.shape[2] < 2 or not np.isfinite(x).all():
        raise ValueError('expected finite [replicate, species>=2, cells>=2] array')
    centered = x-x.mean(axis=-1, keepdims=True)
    s = centered.shape[1]
    total = np.mean(centered.mean(axis=1)**2, axis=1)
    self_term = np.sum(np.mean(centered**2, axis=2), axis=1)/(s*s)
    return np.column_stack((total, self_term, total-self_term))


def expected(amplitude: float, q: int) -> np.ndarray:
    self_term = (amplitude**2 + NOISE_SD**2*(CELLS-1)/CELLS)/SPECIES
    cross = q*(q-1)*amplitude**2/(SPECIES*SPECIES)
    return np.array([self_term+cross, self_term, cross])


def stream(amplitude_index: int, stage: int, fraction_index: int = 0):
    return np.random.default_rng(np.random.SeedSequence([SEED, amplitude_index, stage, fraction_index]))


def run(output: Path) -> dict:
    result_path = output/'global_rgfca_structured_null_toy_result_v1.json'
    rates_path = output/'global_rgfca_structured_null_toy_rates_v1.csv'
    if result_path.exists() or rates_path.exists():
        raise FileExistsError('refusing to overwrite completed toy outputs')
    output.mkdir(parents=True, exist_ok=True)
    rates, moments, thresholds = [], [], []
    maximum_identity_error = 0.0
    for ai, amplitude in enumerate(AMPLITUDES):
        rng, scramble_rng = stream(ai, 0), stream(ai, 1)
        structured, scrambled = [], []
        for _ in range(CALIBRATION//BATCH):
            x = generate(rng, BATCH, amplitude, 0)
            structured.append(decomposition(x))
            # Independent permutations along each species' last (cell) axis.
            scrambled.append(decomposition(scramble_rng.permuted(x, axis=2)))
            gram = x[0] @ x[0].T / CELLS
            direct_cross = 2*np.triu(gram, k=1).sum()/(SPECIES*SPECIES)
            err = abs(direct_cross-structured[-1][0, 2])
            maximum_identity_error = max(maximum_identity_error, float(err))
        models = {'structured_independent_phases': np.vstack(structured),
                  'within_species_cell_scramble': np.vstack(scrambled)}
        cuts = {key: np.quantile(val, 0.95, axis=0, method='higher') for key, val in models.items()}
        for name, cut in cuts.items():
            thresholds.append({'amplitude': amplitude, 'calibration': name,
                               'total_threshold': float(cut[0]), 'cross_threshold': float(cut[2])})
        fractions = FRACTIONS if amplitude else (0.0,)
        for fi, fraction in enumerate(fractions):
            q = round(fraction*SPECIES)
            rng = stream(ai, 2, fi)
            values = np.vstack([decomposition(generate(rng, BATCH, amplitude, q))
                                for _ in range(EVALUATION//BATCH)])
            exp = expected(amplitude, q)
            for j, stat in enumerate(('total', 'self', 'cross')):
                moments.append({'amplitude': amplitude, 'shared_fraction': fraction, 'statistic': stat,
                                'empirical_mean': float(values[:, j].mean()), 'analytic_mean': float(exp[j]),
                                'mean_mcse': float(values[:, j].std(ddof=1)/math.sqrt(EVALUATION))})
            for model, cut in cuts.items():
                for j, stat in ((0, 'total'), (2, 'cross')):
                    count = int(np.count_nonzero(values[:, j] >= cut[j]))
                    lo, hi = wilson(count, EVALUATION)
                    rates.append({'amplitude': amplitude, 'shared_fraction': fraction,
                                  'shared_species': q, 'calibration': model, 'statistic': stat,
                                  'rejections': count, 'replicates': EVALUATION,
                                  'rejection_fraction': count/EVALUATION, 'mc_wilson95_low': lo,
                                  'mc_wilson95_high': hi})
    if maximum_identity_error > 1e-12:
        raise ArithmeticError('covariance decomposition identity failed')
    with rates_path.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(rates[0])); writer.writeheader(); writer.writerows(rates)
    result = {'protocol': PROTOCOL, 'status': 'complete_synthetic_ring_proof_of_concept',
              'specification_commit': SPEC_COMMIT, 'seed': SEED, 'n_species': SPECIES, 'n_cells': CELLS,
              'calibration_replicates_per_amplitude': CALIBRATION, 'evaluation_replicates_per_scenario': EVALUATION,
              'maximum_decomposition_identity_error': maximum_identity_error,
              'thresholds': thresholds, 'moments': moments, 'rates': rates,
              'runner_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'rates_sha256': hashlib.sha256(rates_path.read_bytes()).hexdigest(),
              'numpy_version': np.__version__, 'observed_biological_data_used': False,
              'observed_decisions_reclassified': False,
              'claim_ceiling': 'Uniform periodic-ring methodological example only; not calibrated real-world RGFCA commonness, observed sharing fraction, geographic zone, or biological evidence.',
              'uncertainty': 'Wilson intervals quantify evaluation Monte Carlo error conditional on fixed empirical thresholds; they omit threshold calibration uncertainty.'}
    result_path.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    r = run(args.output_dir)
    print(json.dumps({'status': r['status'], 'rows': len(r['rates']),
                      'maximum_decomposition_identity_error': r['maximum_decomposition_identity_error']}))

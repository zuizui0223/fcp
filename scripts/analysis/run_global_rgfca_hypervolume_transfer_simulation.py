#!/usr/bin/env python3
"""Frozen, synthetic Gaussian-hypervolume / species-disjoint transfer benchmark.

Not real-climate inference, a full RGFCA calibration, or causal identification.
Uses numpy only. Definitions are fixed in the companion simulation contract.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import platform

import numpy as np

SEED = 2026090621
SPECIFICATION_COMMIT = "d309aa12529cdfd229f73d8edfc7877437f02373"
SPECIFICATION_PATH = "docs/supporting/global_rgfca_hypervolume_transfer_simulation_contract_v1.json"
DESIGN = dict(species_per_world=60, photos_per_species=20, training_species=30,
              heldout_species=30, calibration_worlds_per_nuisance=1000,
              evaluation_worlds_per_scenario=1000, amplitudes=[0.5, 1.0, 2.0],
              noise_sd=1.0, transition_width=0.5, gaussian_shrinkage=0.2,
              covariance_ridge=0.000001, minimum_photos_per_colour_class=4,
              alpha_two_domain_family=0.05)
NUISANCES = ("unstructured", "species_specific_geography", "species_specific_environment")
WORLDS = NUISANCES + ("geography_shared", "environment_shared")
DOMAINS = ("geography", "environment")


def unit_vectors(rng: np.random.Generator, shape: tuple[int, ...]) -> np.ndarray:
    x = rng.normal(size=shape)
    lengths = np.linalg.norm(x, axis=-1, keepdims=True)
    if np.any(lengths == 0):
        raise RuntimeError("Zero normal vector; refuse silent resampling.")
    return x / lengths


def world_inputs(stage: int, replicate: int, species: int = 60, photos: int = 20) -> dict:
    if stage not in (0, 1) or replicate < 0 or species < 4 or photos < 8:
        raise ValueError("Invalid stage/replicate or insufficient design size")
    rng = np.random.default_rng(np.random.SeedSequence([SEED, stage, replicate]))
    xyz = unit_vectors(rng, (species, photos, 3))
    x, y = xyz[..., 0], xyz[..., 1]
    geo = np.sqrt(3.0) * xyz
    env = np.stack((x**4 - 6*x*x*y*y + y**4, 4*x*y*(x*x-y*y)), axis=-1)
    env *= np.sqrt(315.0)/8.0  # Unit variance under the uniform-sphere design.
    gnormal = unit_vectors(rng, (3,))
    enormal = unit_vectors(rng, (2,))
    gspecific = unit_vectors(rng, (species, 3))
    especific = unit_vectors(rng, (species, 2))
    return dict(geography=geo, environment=env,
                geography_shared=geo @ gnormal, environment_shared=env @ enormal,
                species_specific_geography=np.einsum('snd,sd->sn', geo, gspecific),
                species_specific_environment=np.einsum('snd,sd->sn', env, especific),
                sign=rng.choice([-1.0, 1.0], size=(species, 1)),
                noise=rng.normal(size=(species, photos)))


def synthetic_labels(inputs: dict, world: str, amplitude: float) -> np.ndarray:
    if not np.isfinite(amplitude) or amplitude < 0:
        raise ValueError("Invalid amplitude")
    if world == 'unstructured':
        signal = 0.0
    else:
        projection = inputs['geography'][..., 0] if world == 'confounded_shared' else inputs[world]
        signal = amplitude * inputs['sign'] * np.tanh(projection / DESIGN['transition_width'])
    return (signal + DESIGN['noise_sd'] * inputs['noise'] > 0).astype(np.int8)


def gaussian_classes(features: np.ndarray, labels: np.ndarray) -> dict:
    """Fit each species separately. Arrays: batch x species x photo x dimension."""
    x = np.asarray(features, dtype=float)
    y = np.asarray(labels)
    if x.ndim != 4 or x.shape[:-1] != y.shape or not np.isfinite(x).all():
        raise ValueError("Features must be complete and match the label shape")
    if not np.isin(y, [0, 1]).all():
        raise ValueError("Binary synthetic colours required")
    d = x.shape[-1]
    means, covs, counts = [], [], []
    for k in (0, 1):
        w = (y == k).astype(float)
        count = w.sum(axis=-1)
        mu = np.einsum('bsn,bsnd->bsd', w, x) / np.maximum(count[..., None], 1)
        diff = x - mu[..., None, :]
        cov = np.einsum('bsn,bsni,bsnj->bsij', w, diff, diff, optimize=True)
        cov /= np.maximum(count[..., None, None]-1, 1)
        target = np.trace(cov, axis1=-2, axis2=-1)/d
        lam = DESIGN['gaussian_shrinkage']
        cov = (1-lam)*cov + (lam*target + DESIGN['covariance_ridge'])[..., None, None]*np.eye(d)
        means.append(mu); covs.append(cov); counts.append(count)
    mu = np.stack(means, axis=2)
    cov = np.stack(covs, axis=2)
    counts = np.stack(counts, axis=-1)
    sign, logdet = np.linalg.slogdet(cov)
    if not (sign > 0).all():
        raise RuntimeError("Nonpositive shrinkage covariance")
    valid = (counts >= DESIGN['minimum_photos_per_colour_class']).all(axis=-1)
    return dict(mean=mu, covariance=cov, precision=np.linalg.inv(cov),
                logdet=logdet, valid=valid)


def bhattacharyya_affinity(fit: dict) -> np.ndarray:
    """Gaussian distribution affinity, not ellipsoid intersection volume."""
    cov = fit['covariance']
    pooled = (cov[:, :, 0] + cov[:, :, 1])/2
    delta = fit['mean'][:, :, 0] - fit['mean'][:, :, 1]
    logdet = np.linalg.slogdet(pooled)[1]
    distance = np.einsum('bsi,bsij,bsj->bs', delta, np.linalg.inv(pooled), delta)
    logbc = .25*fit['logdet'].sum(axis=-1) - .5*logdet - .125*distance
    return np.exp(np.minimum(logbc, 0))


def ellipsoid_log_volume(features: np.ndarray) -> np.ndarray:
    """95% Gaussian ellipsoid volume, dimensions 2 or 3, colour-independent."""
    x = np.asarray(features, dtype=float)
    d = x.shape[-1]
    if d not in (2, 3) or x.shape[-2] < 2 or not np.isfinite(x).all():
        raise ValueError("Complete 2D/3D features and >=2 observations required")
    diff = x - x.mean(axis=-2, keepdims=True)
    cov = np.einsum('...ni,...nj->...ij', diff, diff)/(x.shape[-2]-1)
    target = np.trace(cov, axis1=-2, axis2=-1)/d
    lam = DESIGN['gaussian_shrinkage']
    cov = (1-lam)*cov + (lam*target + DESIGN['covariance_ridge'])[..., None, None]*np.eye(d)
    chi2_95 = {2: 5.991464547107979, 3: 7.814727903251179}[d]
    return (d/2*math.log(math.pi) - math.lgamma(d/2+1)
            + d/2*math.log(chi2_95) + .5*np.linalg.slogdet(cov)[1])


def binary_ari(predicted: np.ndarray, truth: np.ndarray) -> np.ndarray:
    """ARI for last-axis binary partitions, with no-information pairs set to zero."""
    pred, true = np.broadcast_arrays(predicted, truth)
    n = pred.shape[-1]
    a = np.sum((pred == 1) & (true == 1), axis=-1).astype(float)
    p = np.sum(pred == 1, axis=-1).astype(float)
    t = np.sum(true == 1, axis=-1).astype(float)
    b, c, d = p-a, t-a, n-p-t+a
    choose = lambda v: v*(v-1)/2
    obs = choose(a)+choose(b)+choose(c)+choose(d)
    r, q = choose(p)+choose(n-p), choose(t)+choose(n-t)
    expected = r*q/choose(float(n))
    denominator = (r+q)/2-expected
    out = np.zeros_like(expected)
    np.divide(obs-expected, denominator, out=out, where=denominator > 0)
    minimum = DESIGN['minimum_photos_per_colour_class']
    return np.where((t >= minimum) & (n-t >= minimum), out, 0)


def transfer_score(features: np.ndarray, labels: np.ndarray, train_species: int | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Fit ONLY first-half species; evaluate all other species with equal weights.

    Training/target colour identities can be swapped independently without
    changing ARI. Invalid training models contribute zero, not exclusions.
    """
    if train_species is None:
        train_species = DESIGN['training_species']
    if not 0 < train_species < features.shape[1]:
        raise ValueError("Training and heldout species must both be nonempty")
    fit = gaussian_classes(features[:, :train_species], labels[:, :train_species])
    target = features[:, train_species:]
    scores = []
    for k in (0, 1):
        diff = target[:, None] - fit['mean'][:, :, k, None, None, :]
        quadratic = np.sum((diff @ fit['precision'][:, :, k, None]) * diff, axis=-1)
        scores.append(-.5*(quadratic + fit['logdet'][:, :, k, None, None]))
    predicted = scores[1] > scores[0]
    ari = binary_ari(predicted, labels[:, None, train_species:])
    ari = np.where(fit['valid'][:, :, None], ari, 0)
    score = ari.mean(axis=(1, 2))
    separation = np.where(fit['valid'], 1-bhattacharyya_affinity(fit), 0).mean(axis=1)
    return score, separation, fit['valid'].mean(axis=1)


def batch_metrics(stage: int, indices: range, worlds: tuple[str, ...], amplitudes: list[float]) -> list[dict]:
    inputs = [world_inputs(stage, r) for r in indices]
    g = np.stack([x['geography'] for x in inputs])
    e = np.stack([x['environment'] for x in inputs])
    rows = []
    for amp in amplitudes:
        for world in worlds:
            labels = np.stack([synthetic_labels(x, world, amp) for x in inputs])
            domains = (g, g[..., :2] if world == 'confounded_shared' else e)
            out = [transfer_score(x, labels) for x in domains]
            volumes = [ellipsoid_log_volume(x).mean(axis=1) for x in domains]
            true_n = labels[:, DESIGN['training_species']:].sum(axis=-1)
            m = DESIGN['minimum_photos_per_colour_class']
            valid_target = ((true_n >= m) & (true_n <= DESIGN['photos_per_species']-m)).mean(axis=1)
            for j, replicate in enumerate(indices):
                row = dict(stage=('calibration' if stage == 0 else 'evaluation'),
                           replicate=replicate, amplitude=amp, world=world,
                           heldout_class_support_fraction=float(valid_target[j]))
                for domain, result, v in zip(DOMAINS, out, volumes):
                    row[domain+'_transfer'] = float(result[0][j])
                    row[domain+'_separation'] = float(result[1][j])
                    row[domain+'_train_support_fraction'] = float(result[2][j])
                    row[domain+'_log_total_volume'] = float(v[j])
                rows.append(row)
    return rows


def wilson(successes: int, n: int) -> list[float]:
    if n <= 0 or not 0 <= successes <= n:
        raise ValueError("Invalid binomial counts")
    z = 1.959963984540054
    p = successes/n
    den = 1+z*z/n
    center = (p+z*z/(2*n))/den
    radius = z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return [max(0., center-radius), min(1., center+radius)]


def verify_contract(path: Path) -> None:
    c = json.loads(path.read_text())
    if c.get('status') != 'fixed_before_first_benchmark_outcome' or c.get('seed') != SEED or c.get('design') != DESIGN:
        raise RuntimeError("Frozen contract values do not match implementation")


def summarize(rows: list[dict]) -> dict:
    thresholds, rates, details = {}, [], []
    for amp in DESIGN['amplitudes']:
        thresholds[str(amp)] = {}
        for domain in DOMAINS:
            q = {}
            for world in NUISANCES:
                values = [r[domain+'_transfer'] for r in rows if r['stage']=='calibration' and r['amplitude']==amp and r['world']==world]
                q[world] = float(np.quantile(values, .975, method='higher'))
            thresholds[str(amp)][domain] = dict(value=max(q.values()), nuisance_quantiles=q)
        for world in WORLDS + (('confounded_shared',) if amp == 2.0 else ()):
            subset = [r for r in rows if r['stage']=='evaluation' and r['amplitude']==amp and r['world']==world]
            gs = np.array([r['geography_transfer'] for r in subset])
            es = np.array([r['environment_transfer'] for r in subset])
            passed_g = gs > thresholds[str(amp)]['geography']['value']
            passed_e = es > thresholds[str(amp)]['environment']['value']
            groups = {'geography_only': passed_g & ~passed_e, 'environment_only': passed_e & ~passed_g,
                      'both': passed_g & passed_e, 'neither': ~passed_g & ~passed_e,
                      'any': passed_g | passed_e, 'geography': passed_g, 'environment': passed_e}
            for decision, values in groups.items():
                n, k = len(values), int(values.sum())
                ci = wilson(k, n)
                rates.append(dict(amplitude=amp, world=world, decision=decision, successes=k, replicates=n,
                                  rate=k/n, wilson95_low=ci[0], wilson95_high=ci[1]))
            details.append(dict(amplitude=amp, world=world,
                                geography_transfer_mean=float(gs.mean()), environment_transfer_mean=float(es.mean()),
                                geography_separation_mean=float(np.mean([r['geography_separation'] for r in subset])),
                                environment_separation_mean=float(np.mean([r['environment_separation'] for r in subset])),
                                heldout_class_support_fraction=float(np.mean([r['heldout_class_support_fraction'] for r in subset]))))
    # Matched geometry must give exactly the same total volume under all ordinary worlds.
    groups = {}
    for row in rows:
        if row['world'] == 'confounded_shared':
            continue
        key = (row['stage'], row['replicate'])
        vector = tuple(row[d+'_log_total_volume'] for d in DOMAINS)
        if key in groups and vector != groups[key]:
            raise RuntimeError("Colour-independent total-volume invariant failed")
        groups[key] = vector
    return dict(protocol='global-rgfca-hypervolume-transfer-simulation-v1',
                status='complete_synthetic_proof_of_concept', specification_commit=SPECIFICATION_COMMIT,
                seed=SEED, design=DESIGN, thresholds=thresholds, rates=rates, diagnostics=details,
                total_volume_label_invariant=True, matched_geometry_worlds=len(groups),
                world_rows=len(rows), observed_biological_data_used=False,
                real_climate_data_used=False, actual_geometry_power_calibrated=False,
                geographic_block_holdout_completed=False, parent_decisions_modified=False,
                uncertainty='Wilson intervals condition on the estimated thresholds; calibration uncertainty is not included.',
                claim_ceiling='Synthetic representation-specific partition transfer only. Both-domain and neither-domain results are retained; not causal identification or observed flower-colour support.',
                versions=dict(python=platform.python_version(), numpy=np.__version__))


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open('w', newline='') as fout:
        w = csv.DictWriter(fout, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir', type=Path, required=True)
    p.add_argument('--contract', type=Path)
    p.add_argument('--batch-size', type=int, default=10)
    args = p.parse_args()
    if args.batch_size <= 0:
        p.error('batch-size must be positive')
    if args.contract is not None:
        verify_contract(args.contract)
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    result_path = out/'global_rgfca_hypervolume_transfer_result_v1.json'
    if result_path.exists():
        raise SystemExit('Refusing to overwrite existing result; use a fresh output directory for exact reproducibility checks.')
    rows = []
    for stage, n, worlds in ((0, DESIGN['calibration_worlds_per_nuisance'], NUISANCES),
                             (1, DESIGN['evaluation_worlds_per_scenario'], WORLDS)):
        for start in range(0, n, args.batch_size):
            idx = range(start, min(start+args.batch_size, n))
            rows.extend(batch_metrics(stage, idx, worlds, DESIGN['amplitudes']))
            if stage == 1:
                rows.extend(batch_metrics(stage, idx, ('confounded_shared',), [2.0]))
            if start % 100 == 0:
                print(json.dumps(dict(stage=stage, completed=min(start+args.batch_size, n), total=n)), flush=True)
    summary = summarize(rows)
    draws_path = out/'global_rgfca_hypervolume_transfer_draws_v1.csv'
    rates_path = out/'global_rgfca_hypervolume_transfer_rates_v1.csv'
    write_csv(draws_path, rows); write_csv(rates_path, summary['rates'])
    summary['hashes'] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__), draws_path, rates_path)}
    result_path.write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(status=summary['status'], world_rows=len(rows), result=str(result_path))), flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

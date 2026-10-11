#!/usr/bin/env python3
"""Fixed shared nonlinear environmental transfer, conditional on original folds.

Retrospective photographic prediction, with no causal mechanism or fitness claim.
The frozen protocol fixes the learner; original aggregation, ridge, deterministic
folds and whole-species raw-photo footprint exclusions are imported unchanged.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
import sklearn
from sklearn.ensemble import ExtraTreesRegressor
import fcp_audit_36_preflight_20261011 as src
from fcp_all36_global_region_cv_20261011 import aggregate_species
from fcp_joint_environment_taxongeo_stable_20261011 import StableGenusKFold, ALGORITHM
from fcp_joint_environment_taxongeo_20261011 import (
    sector_ids, spatial_exclusions, select_train, environmental_support, SCHEMES)
from fcp_joint_environment_unseen_genera_20261011 import (
    predict, loss_components, verify_file, SHA_PHOTO, SHA_ENV, FULL, GEO, CLASSES, ALPHA)

SCHEMA = 'fcp_fixed_nonlinear_environment_transport_20261011_v1'
COMPONENTS = ('total', 'white_axis', 'remaining_hue')
MODELS = ('ridge_geo', 'ridge_full', 'nonlinear_geo', 'nonlinear_full')
TREE_PARAMS = dict(n_estimators=128, max_depth=16, min_samples_leaf=20,
    max_features=1.0, bootstrap=False, criterion='squared_error',
    random_state=2026101161, n_jobs=2)
CONTRASTS = {
    'nonlinear_environment_gain': (0., 0., 1., -1.),
    'ridge_environment_gain': (1., -1., 0., 0.),
    'nonlinear_minus_ridge_environment_gain': (-1., 1., 1., -1.),
    'full_model_improvement': (0., 1., 0., -1.),
    'geography_model_improvement': (1., 0., -1., 0.),
}
TARGET_KEY = ['inat_taxon_id', 'sector', 'scheme']
TARGET_ID = TARGET_KEY + ['fold', 'genus']
SEED = 2026101151
PARENT = 'ae1136e7d43bc6c5711778a5a7c727b0b7ed8d7e'
ARCHIVES = {
    'fcp_original_photo_site_climate_soil_85337_20261010.zip': SHA_PHOTO,
    'fcp_100543_all36_expanded_photo_site_environment_20261011.zip': SHA_ENV,
}


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):h.update(chunk)
    return h.hexdigest()


def composition(array):
    a = np.asarray(array, float)
    if (a.ndim != 2 or a.shape[1] != 4 or not np.isfinite(a).all()
            or (a < -1e-12).any() or (a > 1 + 1e-12).any()
            or not np.allclose(a.sum(1), 1., atol=1e-10, rtol=0)):
        raise ValueError('finite four-colour simplex composition required')
    return a


def check_split(train, target):
    if not len(train):raise ValueError('empty training support')
    if not train.inat_taxon_id.is_unique:raise ValueError('training species must be unique')
    if not target.inat_taxon_id.is_unique:raise ValueError('target species must be unique within sector')
    if set(train.genus) & set(target.genus):raise ValueError('genus leakage into training')
    if set(train.inat_taxon_id) & set(target.inat_taxon_id):raise ValueError('species leakage into training')


def predict_nonlinear(train, target, features):
    """Joint four-output tree; target labels are never read or used for selection."""
    check_split(train, target)
    x = train[list(features)].to_numpy(float)
    xt = target[list(features)].to_numpy(float)
    y = composition(train[['hue_' + h for h in CLASSES]].to_numpy(float))
    if not np.isfinite(x).all() or not np.isfinite(xt).all():
        raise ValueError('nonfinite training or target features')
    model = ExtraTreesRegressor(**TREE_PARAMS)
    model.fit(x, y, sample_weight=np.ones(len(train)))
    # Raw tree averages already lie on the simplex; do not project or recalibrate.
    return composition(model.predict(xt))


def predict_models(train, target, environment=FULL):
    check_split(train, target)
    result = {
        'ridge_geo': predict(train, target, GEO, 'absolute'),
        'ridge_full': predict(train, target, GEO + tuple(environment), 'absolute'),
        'nonlinear_geo': predict_nonlinear(train, target, GEO),
        'nonlinear_full': predict_nonlinear(train, target, GEO + tuple(environment)),
    }
    for value in result.values():composition(value)
    return result


def validate_records(records, expected=None):
    """Fail on coverage loss, duplicate keys, label misalignment or damaged losses."""
    if records.empty:raise ValueError('no target records')
    if records[TARGET_ID].isna().any().any():raise ValueError('missing target identifiers')
    if records.duplicated(TARGET_KEY).any():raise ValueError('duplicate species-sector target rows')
    if records.groupby('inat_taxon_id').genus.nunique().max() != 1:
        raise ValueError('target species assigned inconsistent genera')
    if records.groupby('genus').fold.nunique().max() != 1:
        raise ValueError('target genus crosses folds')
    if expected is not None:
        if expected.duplicated(TARGET_KEY).any():raise ValueError('duplicate expected target rows')
        got = records[TARGET_ID].sort_values(TARGET_KEY).reset_index(drop=True)
        want = expected[TARGET_ID].sort_values(TARGET_KEY).reset_index(drop=True)
        if not got.equals(want):raise ValueError('mismatched or missing comparison target rows')
    truth = composition(records[['truth__' + h for h in CLASSES]].to_numpy(float))
    for model in MODELS:
        prediction = composition(records[[model + '__p__' + h for h in CLASSES]].to_numpy(float))
        losses = records[[model + '__' + c for c in COMPONENTS]].to_numpy(float)
        if not np.isfinite(losses).all():raise ValueError('nonfinite model loss')
        if not np.allclose(losses, loss_components(truth - prediction), atol=1e-12, rtol=0):
            raise ValueError('saved loss disagrees with probability prediction')


def interval(point, draws):
    return {'estimate': float(point), 'CI95': np.quantile(draws, [.025, .975]).tolist(),
            'CI97_5': np.quantile(draws, [.0125, .9875]).tolist()}


def summarize_records(records, reps=1999, seed=SEED):
    """Average sectors per species first, then pair and bootstrap whole genera."""
    if reps < 2:raise ValueError('at least two bootstrap draws required')
    if records.empty:return {'status': 'HOLD_NO_TARGET_RECORDS', 'bootstrap_seed': seed}
    if records.scheme.nunique() != 1:raise ValueError('summarize one training regime at a time')
    validate_records(records)
    cols = [model + '__' + c for model in MODELS for c in COMPONENTS]
    sp = records.groupby('inat_taxon_id', sort=True).agg({**dict.fromkeys(cols, 'mean'), 'genus': 'first'})
    genera = sp.groupby('genus', sort=True)
    sums = genera[cols].sum().to_numpy(float)
    counts = genera.size().to_numpy(float)
    ng = len(counts)
    out = {'status': 'OK', 'n_species': len(sp), 'n_genera': ng,
           'n_species_sector_records': len(records), 'n_sectors': int(records.sector.nunique()),
           'bootstrap_replicates': reps, 'bootstrap_seed': seed,
           'uncertainty': 'paired whole-genus resampling conditional on fixed fits and folds'}
    if ng < 5:
        out['status'] = 'HOLD_FEWER_THAN_5_GENERA'
        return out
    means = sums / counts[:, None]
    rng = np.random.default_rng(seed)
    bootstrap = {'genus_equal': [], 'species_equal': []}
    # Bounded batches: fits are never rerun inside bootstrap.
    for start in range(0, reps, 80):
        ix = rng.integers(ng, size=(min(80, reps - start), ng))
        bootstrap['genus_equal'].append(means[ix].mean(1))
        bootstrap['species_equal'].append(sums[ix].sum(1) / counts[ix].sum(1)[:, None])
    contrast_matrix = np.asarray(list(CONTRASTS.values()))
    for weight, batches in bootstrap.items():
        point = means.mean(0) if weight == 'genus_equal' else sp[cols].mean().to_numpy(float)
        point = point.reshape(4, 3)
        draws = np.concatenate(batches).reshape(reps, 4, 3)
        contrasts = contrast_matrix @ point
        contrast_draws = np.einsum('km,bmc->bkc', contrast_matrix, draws)
        out[weight] = {}
        for j, c in enumerate(COMPONENTS):
            out[weight][c] = {
                'losses': {model: interval(point[i, j], draws[:, i, j]) for i, model in enumerate(MODELS)},
                'contrasts': {name: interval(contrasts[i, j], contrast_draws[:, i, j])
                              for i, name in enumerate(CONTRASTS)},
            }
    return out


def compare_archived_ridge(records, archive):
    """Read-only audit for independent parent comparison; never discovers inputs."""
    old = pd.read_csv(archive)
    for d in (records, old):
        if d.duplicated(TARGET_KEY).any():raise ValueError('duplicate archived target')
    a = records.sort_values(TARGET_KEY).reset_index(drop=True)
    b = old.sort_values(TARGET_KEY).reset_index(drop=True)
    if not a[TARGET_ID].equals(b[TARGET_ID]):raise ValueError('archived target rows mismatch')
    differences = {}
    for new, previous in [('ridge_geo', 'geography'), ('ridge_full', 'environment36')]:
        for c in COMPONENTS:
            aa = a[new + '__' + c].to_numpy(float)
            bb = b[previous + '__' + c].to_numpy(float)
            if not np.isfinite(bb).all():raise ValueError('nonfinite archived ridge loss')
            difference = float(np.max(np.abs(aa - bb)))
            differences[new + '__' + c] = difference
            if not np.allclose(aa, bb, atol=1e-11, rtol=0):
                raise ValueError('recomputed ridge differs from archived loss')
    return {'status': 'MATCH', 'n_target_rows': len(a), 'max_absolute_loss_difference': differences,
            'archive_sha256': sha256(archive)}


def fold_mapping(species):
    mapping = {}
    for fold, (_, te) in enumerate(StableGenusKFold(5).split(species, groups=species.genus)):
        for genus in species.iloc[te].genus.unique():mapping[genus] = fold
    raw = json.dumps(mapping, sort_keys=True, ensure_ascii=True, separators=(',', ':')).encode()
    return mapping, hashlib.sha256(raw).hexdigest()


def identity_hash(values):
    return hashlib.sha256(json.dumps(sorted(set(values)), separators=(',', ':')).encode()).hexdigest()


def evaluate_case(raw, quality, exclusions, reps, outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    d = raw.copy()
    if quality:
        d = d.loc[(d.flower_effective_pixels >= 23826) &
                  (d.nuisance_pixel_fraction <= 0.1375717742195135)].copy()
    d['sector'] = sector_ids(d.longitude)
    species = aggregate_species(d).sort_values('inat_taxon_id', kind='stable').reset_index(drop=True)
    mapping, fold_hash = fold_mapping(species)
    species['fold'] = species.genus.map(mapping)
    case = 'soil_quality' if quality else 'soil_all'
    pd.DataFrame({'genus': sorted(mapping), 'fold': [mapping[g] for g in sorted(mapping)]}).to_csv(
        outdir / (case + '_genus_fold_map.csv'), index=False)
    rows, audits, expected = [], [], []
    for sector in range(4):
        target_sp = aggregate_species(d.loc[d.sector.eq(sector)]).sort_values('inat_taxon_id').reset_index(drop=True)
        target_sp['fold'] = target_sp.genus.map(mapping)
        for scheme in SCHEMES:
            e = target_sp[['inat_taxon_id', 'genus', 'fold']].assign(sector=sector, scheme=scheme)
            expected.append(e)
        for fold in range(5):
            target = target_sp.loc[target_sp.fold.eq(fold)].copy()
            if not len(target):continue
            truth = composition(target[['hue_' + h for h in CLASSES]].to_numpy(float))
            for scheme in SCHEMES:
                forbidden = set() if scheme == 'genus_only' else exclusions[sector][scheme]
                train = select_train(species, fold, forbidden)
                guard = {'sector': sector, 'fold': fold, 'scheme': scheme,
                    'n_train_species': len(train), 'n_train_genera': int(train.genus.nunique()),
                    'n_target_species': len(target), 'n_target_photos': int(target.n_photos.sum()),
                    'shared_train_test_genus': len(set(train.genus) & set(target.genus)),
                    'shared_train_test_species': len(set(train.inat_taxon_id) & set(target.inat_taxon_id)),
                    'training_species_intersect_spatial_blacklist': len(set(train.inat_taxon_id) & forbidden),
                    'train_species_sha256': identity_hash(train.inat_taxon_id.astype(int).tolist()),
                    'target_species_sha256': identity_hash(target.inat_taxon_id.astype(int).tolist()),
                    'n_models_identical_split': len(MODELS)}
                if any(guard[k] for k in ['shared_train_test_genus', 'shared_train_test_species',
                                         'training_species_intersect_spatial_blacklist']):
                    raise ValueError('source isolation failure: ' + str(guard))
                if len(train) < 300 or train.genus.nunique() < 30:
                    guard['status'] = 'HOLD_TRAINING_SUPPORT'
                    audits.append(guard)
                    pd.DataFrame(audits).to_csv(outdir / (case + '_fold_audit.csv'), index=False)
                    raise ValueError('insufficient training support; target rows cannot be dropped: ' + str(guard))
                prediction = predict_models(train, target)
                support, distance, cutoff = environmental_support(train[list(FULL)], target[list(FULL)])
                guard.update(status='OK', n_target_within_training_environment_envelope=int(support.sum()),
                             train_environment_d2_99_cutoff=cutoff)
                audits.append(guard)
                r = target[['inat_taxon_id', 'genus', 'n_photos', 'cell']].copy()
                r['sector'] = sector;r['fold'] = fold;r['scheme'] = scheme
                r['environment_support'] = support;r['environment_squared_distance'] = distance
                for j, h in enumerate(CLASSES):r['truth__' + h] = truth[:, j]
                for model, values in prediction.items():
                    losses = loss_components(truth - values)
                    for j, h in enumerate(CLASSES):r[model + '__p__' + h] = values[:, j]
                    for j, c in enumerate(COMPONENTS):r[model + '__' + c] = losses[:, j]
                rows.append(r)
        print(json.dumps({'case': case, 'sector_completed': sector, 'n_species_sector_targets': len(target_sp)}), flush=True)
    records = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    validate_records(records, pd.concat(expected, ignore_index=True))
    records.to_csv(outdir / (case + '_oof_predictions.csv.gz'), index=False, compression='gzip')
    pd.DataFrame(audits).to_csv(outdir / (case + '_fold_audit.csv'), index=False)
    result = {'n_input_photos': len(d), 'n_input_species': len(species),
              'n_input_genera': int(species.genus.nunique()), 'genus_fold_map_sha256': fold_hash,
              'all_eligible_source_photos_accounted': True, 'identical_complete_target_rows_all_models_and_regimes': True,
              'fold_status_counts': pd.Series([a['status'] for a in audits]).value_counts().to_dict(), 'schemes': {}}
    for scheme in SCHEMES:
        rr = records.loc[records.scheme.eq(scheme)]
        if (rr.inat_taxon_id.nunique() != len(species) or rr.n_photos.sum() != len(d)):
            raise ValueError('incomplete source target species/photo coverage')
        seed = SEED + sum(map(ord, case + scheme))
        result['schemes'][scheme] = {
            'n_evaluated_species': int(rr.inat_taxon_id.nunique()),
            'n_evaluated_source_photos': int(rr.n_photos.sum()), 'n_species_sector_records': len(rr),
            'overall': summarize_records(rr, reps, seed),
            'regions': {str(s): summarize_records(rr.loc[rr.sector.eq(s)], reps, seed + 11 + s)
                        for s in range(4)},
        }
    return result


def run(source_root, outdir, bootstrap=1999, quality='both'):
    if bootstrap < 99:raise ValueError('at least 99 bootstrap replicates required')
    if quality not in ('both', 'all', 'strict'):raise ValueError('invalid quality cohort')
    root, outdir = Path(source_root), Path(outdir)
    for name, expected in ARCHIVES.items():verify_file(root / name, expected)
    src.A, src.B = [root / name for name in ARCHIVES]
    data, counts, masks = src.run()
    if (len(data) != 100543 or data.photo_id.nunique() != 100543 or
            data.inat_taxon_id.nunique() != 42111 or counts['classified'] != 45953):
        raise ValueError('source opportunity drift')
    if counts['soil_complete'] != 31925 or counts['soil_species'] != 20786:
        raise ValueError('all36 complete source cohort drift')
    exclusions, geography_audit = spatial_exclusions(data, buffer_km=500.)
    outdir.mkdir(parents=True, exist_ok=True)
    directory = Path(__file__).resolve().parent
    hashes = {p.name: sha256(p) for p in [Path(__file__).resolve(), *[
        directory / name for name in ('fcp_audit_36_preflight_20261011.py',
        'fcp_all36_global_region_cv_20261011.py', 'fcp_joint_environment_unseen_genera_20261011.py',
        'fcp_joint_environment_taxongeo_20261011.py', 'fcp_joint_environment_taxongeo_stable_20261011.py')]]}
    protocol = directory.parents[1] / 'docs/FCP_NONLINEAR_ENVIRONMENT_TRANSPORT_PROTOCOL_20261011.md'
    hashes[protocol.name] = sha256(protocol)
    result = {'schema': SCHEMA, 'source_commit_parent': PARENT, 'source_hashes': ARCHIVES,
        'script_and_helper_hashes': hashes, 'source_counts': counts, 'full36_predictors': list(FULL),
        'geography_predictors': list(GEO), 'classes': list(CLASSES), 'components': list(COMPONENTS),
        'models': list(MODELS), 'contrast_coefficients_in_model_order': CONTRASTS,
        'model_settings': {'nonlinear': TREE_PARAMS, 'ridge_alpha': ALPHA,
            'training_weights': 'equal species', 'ridge_preprocessing': 'training mean/std only; original simplex projection',
            'tree_preprocessing': 'none; raw training features; joint four-output averages'},
        'genus_fold_algorithm': ALGORITHM, 'n_genus_folds': 5, 'bootstrap': bootstrap,
        'primary': {'case': 'soil_all', 'scheme': 'buffer500', 'weight': 'genus_equal',
                    'component': 'total', 'contrast': 'nonlinear_environment_gain'},
        'buffer_km': 500., 'longitude_sectors': [[-180, -90], [-90, 0], [0, 90], [90, 180]],
        'strict_quality_thresholds': {'flower_effective_pixels_min': 23826,
                                    'nuisance_pixel_fraction_max': 0.1375717742195135},
        'source_raw_footprint_exclusion_audit': geography_audit,
        'ridge_archive_comparison': 'pending independent comparison against archived PR169 OOF losses',
        'main_H1H2_changed': False, 'new_colour_measurements': False,
        'execution_status': 'LOCAL_SOURCE_REPLAY_IN_PROGRESS',
        'runtime': {'python': platform.python_version(), 'numpy': np.__version__, 'pandas': pd.__version__,
                    'scipy': scipy.__version__, 'scikit_learn': sklearn.__version__},
        'limitations': [
            'Retrospective source-exposed fixed model-form test, not preregistration or independent biological replication.',
            'Four longitude sectors are artificial domains, not independent evolutionary replicates.',
            'GEO is a specified coordinate baseline, not comprehensive spatial adjustment.',
            'Bootstrap conditions on existing fits and folds and does not refit learners.',
            'Photographic colour proportions do not identify pigments, selection, pollination, herbivory, fitness or genotypes.',
            'Missing classification, measurement error and unobserved real populations remain.',
            'Predictive failure cannot establish universal absence of environmental responses.'], 'cases': {}}
    for strict in ((False, True) if quality == 'both' else (quality == 'strict',)):
        name = 'soil_quality' if strict else 'soil_all'
        case = evaluate_case(data.loc[masks['soil']], strict, exclusions, bootstrap, outdir)
        if case['n_input_species'] != (6457 if strict else 20786):raise ValueError('fixed cohort species count drift')
        result['cases'][name] = case
        (outdir / 'result.json').write_text(json.dumps(result, indent=2, ensure_ascii=False) + '\n')
    result['execution_status'] = 'SOURCE_REPLAY_COMPLETE'
    (outdir / 'result.json').write_text(json.dumps(result, indent=2, ensure_ascii=False) + '\n')
    print('SOURCE_REPLAY_COMPLETE', flush=True)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-root', type=Path, required=True)
    p.add_argument('--outdir', type=Path, required=True)
    p.add_argument('--bootstrap', type=int, default=1999)
    p.add_argument('--quality', choices=('both', 'all', 'strict'), default='both')
    a = p.parse_args()
    run(a.source_root, a.outdir, a.bootstrap, a.quality)


if __name__ == '__main__':main()

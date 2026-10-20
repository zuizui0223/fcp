#!/usr/bin/env python3
"""Execute the frozen FCP held-out-genus and held-out-spatial-cell analyses.

This is post hoc internal transport evaluation, not external validation or a
causal analysis. See docs/FCP_ENVIRONMENT_TRANSPORT_PROTOCOL_20261019.md.
All predictor choices, folds, buffers, and comparison families are fixed.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import sys
import time

import numpy as np
import pandas as pd
import scipy
from scipy.optimize import brentq
from scipy.special import expit

from scripts.analysis.run_fcp_447500_environment_model_20261018 import (
    PREDICTORS, fit_firth, load_panel, reml_meta,
)
from scripts.analysis.summarize_fcp_447500_within_species_environment_20261018 import (
    PHOTO_SHA, ZIP_SHA, demand, json_safe, sha,
)

GENUS_SEED = 'FCP-20261019-genus10-v1'
CELL_SEED = 'FCP-20261019-global1deg5-v1'
INCREMENTS = np.array([1., 10., 1000.])
EARTH_RADIUS_KM = 6371.0088
MODELS = ('intercept', 'climate', 'geography', 'geoclimate')
COMPARISONS = (('climate_vs_intercept', 0, 1), ('geoclimate_vs_geography', 2, 3))
BUFFERS_KM = (0., 50.)
BOOTSTRAP_REPETITIONS = 2000
GENUS_BOOTSTRAP_SEED = 2026101901
SPATIAL_BOOTSTRAP_SEED = 2026101902


def stable_fold(key, seed, n_folds):
    """SHA256 of UTF-8 seed:exact-key; whole digest interpreted big-endian."""
    digest = hashlib.sha256((seed + ':' + str(key)).encode('utf-8')).digest()
    return int.from_bytes(digest, 'big') % n_folds


def convert_units(beta, covariance, species_sd, increments=INCREMENTS):
    factor = np.asarray(increments, float) / np.asarray(species_sd, float)
    demand(np.isfinite(factor).all() and (factor > 0).all(), 'invalid physical-unit scale')
    return np.asarray(beta) * factor, np.asarray(covariance) * np.outer(factor, factor)


def learn_fold_coefficients(slopes, heldout_fold):
    """Only training-genus joint slopes enter the three marginal REML fits."""
    genus = slopes.species.str.split().str[0]
    folds = genus.map(lambda g: stable_fold(g, GENUS_SEED, 10))
    train = slopes.loc[folds.ne(heldout_fold)].copy()
    held = slopes.loc[folds.eq(heldout_fold)]
    training_genera = set(train.species.str.split().str[0])
    heldout_genera = set(held.species.str.split().str[0])
    demand(training_genera.isdisjoint(heldout_genera), 'genus fold leakage')
    records = []
    for variable, increment in zip(PREDICTORS, INCREMENTS):
        g = train.loc[train.variable.eq(variable)]
        factor = increment / g.primary_within_species_sd_native.to_numpy(float)
        y = g.slope.to_numpy(float) * factor
        variance = g.variance_approx.to_numpy(float) * factor ** 2
        result, _, _, _ = reml_meta(y, variance)
        records.append({'heldout_genus_fold':heldout_fold, 'variable':variable,
                        'physical_increment':increment, 'coefficient':result['reml_mean'],
                        'n_training_species':int(g.inat_taxon_id.nunique()),
                        'n_training_genera':len(training_genera),
                        'n_heldout_species':int(held.inat_taxon_id.nunique()),
                        'n_heldout_genera':len(heldout_genera), **result})
    return pd.DataFrame(records)


def pairwise_metrics(white_score, chromatic_score):
    """Every white/chromatic pair, with bounded temporary memory and half ties."""
    white = np.asarray(white_score, float)
    chromatic = np.asarray(chromatic_score, float)
    demand(len(white) > 0 and len(chromatic) > 0, 'pair metric needs both states')
    demand(np.isfinite(white).all() and np.isfinite(chromatic).all(), 'nonfinite pair score')
    n_pairs = len(white) * len(chromatic)
    loss = 0.
    wins = 0.
    block = max(1, 1_000_000 // len(chromatic))
    for start in range(0, len(white), block):
        delta = white[start:start + block, None] - chromatic[None, :]
        loss += float(np.logaddexp(0., -delta).sum())
        wins += float((delta > 0).sum()) + .5 * float((delta == 0).sum())
    loss /= n_pairs
    return {'n_pairs':n_pairs, 'pair_logloss':loss,
            'pair_logloss_gain':float(np.log(2.) - loss), 'pair_auc':wins / n_pairs}


def cell_folds(latitude, longitude):
    """Canonical global one-degree cells; +180 and -180 share their cell."""
    lat = np.asarray(latitude, float)
    lon = (np.asarray(longitude, float) + 180.) % 360. - 180.
    demand(np.isfinite(lat).all() and np.isfinite(lon).all(), 'nonfinite coordinates')
    demand(((lat >= -90.) & (lat <= 90.)).all(), 'latitude outside globe')
    latcell = np.minimum(np.floor(lat).astype(int), 89)
    loncell = np.floor(lon).astype(int)
    cells = np.array([f'{a},{b}' for a, b in zip(latcell, loncell)])
    mapping = {c:stable_fold(c, CELL_SEED, 5) for c in np.unique(cells)}
    return cells, np.array([mapping[c] for c in cells], dtype=np.int8)


def great_circle_km(a_latlon, b_latlon):
    """Haversine distance, including antimeridian and polar geometry."""
    a = np.deg2rad(np.asarray(a_latlon, float))
    b = np.deg2rad(np.asarray(b_latlon, float))
    dlat = a[:, None, 0] - b[None, :, 0]
    dlon = a[:, None, 1] - b[None, :, 1]
    hav = (np.sin(dlat / 2.) ** 2
           + np.cos(a[:, None, 0]) * np.cos(b[None, :, 0]) * np.sin(dlon / 2.) ** 2)
    return 2. * EARTH_RADIUS_KM * np.arcsin(np.sqrt(np.clip(hav, 0., 1.)))


def spatial_split_indices(coords, folds, test_fold, buffer_km):
    """Counts, cells and coordinates alone select rows; no response is accepted."""
    train = np.flatnonzero(np.asarray(folds) != test_fold)
    test = np.flatnonzero(np.asarray(folds) == test_fold)
    if not len(test):
        return train, test, np.full(len(train), np.nan)
    if not len(train):
        return train, test, np.array([], float)
    nearest = great_circle_km(np.asarray(coords)[train], np.asarray(coords)[test]).min(axis=1)
    if buffer_km > 0.:
        keep = nearest > buffer_km
        train, nearest = train[keep], nearest[keep]
    demand(not set(train).intersection(test), 'spatial row leakage')
    return train, test, nearest


def fit_feature_basis(train_climate, train_coords, model):
    """Fit centering, SD and circular longitude origin on training features only."""
    demand(model in MODELS, 'unknown spatial model')
    climate = np.asarray(train_climate, float)
    coords = np.asarray(train_coords, float)
    names = []
    columns = []
    basis = {'status':'ok', 'model':model}
    if model in ('climate', 'geoclimate'):
        sd = climate.std(axis=0)
        threshold = 32. * np.finfo(float).eps * np.maximum(1., np.abs(climate).max(axis=0))
        basis['training_climate_sd'] = sd
        if np.any(sd <= threshold):
            basis['status'] = 'zero_or_near_zero_training_predictor_sd'
            return basis
        names.extend(PREDICTORS)
        columns.append(climate)
    if model in ('geography', 'geoclimate'):
        circular = np.exp(1j * np.deg2rad(coords[:, 1])).mean()
        origin = float(np.angle(circular))
        delta = np.rad2deg(np.angle(np.exp(1j * (np.deg2rad(coords[:, 1]) - origin))))
        geographic = np.column_stack([coords[:, 0], delta])
        basis.update({'longitude_origin_radians':origin,
                      'longitude_circular_resultant':float(np.abs(circular)),
                      'unwrapped_longitude_range_degrees':float(np.ptp(delta)),
                      'training_geographic_sd':geographic.std(axis=0)})
        if np.abs(circular) < .1 or np.ptp(delta) > 180. or np.any(geographic.std(axis=0) < 1e-8):
            basis['status'] = 'geographic_basis_not_identifiable'
            return basis
        names.extend(('latitude', 'unwrapped_longitude'))
        columns.append(geographic)
    raw = np.column_stack(columns) if columns else np.empty((len(coords), 0))
    basis.update({'names':names, 'mean':raw.mean(axis=0), 'scale':raw.std(axis=0),
                  'min':raw.min(axis=0), 'max':raw.max(axis=0)})
    return basis


def apply_feature_basis(basis, climate, coords):
    columns = []
    if basis['model'] in ('climate', 'geoclimate'):
        columns.append(np.asarray(climate, float))
    if basis['model'] in ('geography', 'geoclimate'):
        coords = np.asarray(coords, float)
        delta = np.rad2deg(np.angle(np.exp(1j * (
            np.deg2rad(coords[:, 1]) - basis['longitude_origin_radians']))))
        columns.append(np.column_stack([coords[:, 0], delta]))
    raw = np.column_stack(columns) if columns else np.empty((len(coords), 0))
    return np.column_stack([np.ones(len(coords)), (raw - basis['mean']) / basis['scale']])


def recalibrate_intercept(train_eta, train_y):
    """Fixed smoothed train-only target, keeping all Firth slopes unchanged.

    This is the pre-result protocol clarification, not standard FLIC. A common
    Jeffreys-smoothed target avoids different Firth model dimensions changing
    mean training calibration and remains finite for constant responses.
    """
    eta = np.asarray(train_eta, float)
    y = np.asarray(train_y, float)
    demand(len(y) > 0 and len(y) == len(eta), 'recalibration training shape')
    demand(np.isfinite(eta).all() and np.isin(y, [0., 1.]).all(), 'invalid recalibration training data')
    target = float((y.sum() + .5) / (len(y) + 1.))
    objective = lambda shift: float(expit(eta + shift).mean() - target)
    shift = float(brentq(objective, -float(eta.max()) - 40.,
                         -float(eta.min()) + 40., xtol=1e-12, rtol=1e-14))
    achieved = float(expit(eta + shift).mean())
    demand(abs(achieved - target) < 1e-10, 'recalibration root did not attain training target')
    return {'shift':shift, 'target':target, 'achieved':achieved,
            'raw_training_mean':float(expit(eta).mean())}


def fit_one_model(train_climate, train_coords, train_y, test_climate, test_coords, model):
    """Train one fixed Firth model; the function never accepts held-out labels."""
    y = np.asarray(train_y, float)
    out = {'model':model, 'n_train':len(y), 'n_train_white':int(y.sum()),
           'constant_training_response':bool(len(y) > 0 and (np.all(y == 0) or np.all(y == 1)))}
    if len(y) < 20:
        return {**out, 'status':'training_n_below_20'}
    if len(test_coords) == 0:
        return {**out, 'status':'no_test_photos'}
    basis = fit_feature_basis(train_climate, train_coords, model)
    out['basis'] = basis
    if basis['status'] != 'ok':
        return {**out, 'status':basis['status']}
    xtrain = apply_feature_basis(basis, train_climate, train_coords)
    xtest = apply_feature_basis(basis, test_climate, test_coords)
    fitted = fit_firth(xtrain, y)
    out.update(fitted)
    if fitted['status'] != 'ok':
        return out
    recalibration = recalibrate_intercept(xtrain @ fitted['beta'], y)
    raw_eta = xtest @ fitted['beta']
    eta = raw_eta + recalibration['shift']
    out.update({'raw_intercept':float(fitted['beta'][0]),
                'adjusted_intercept':float(fitted['beta'][0] + recalibration['shift']),
                'intercept_recalibration_shift':recalibration['shift'],
                'training_probability_target':recalibration['target'],
                'recalibrated_training_mean_probability':recalibration['achieved'],
                'raw_firth_training_mean_probability':recalibration['raw_training_mean']})
    if not np.isfinite(eta).all():
        out['status'] = 'nonfinite_test_prediction'
        return out
    if xtest.shape[1] > 1:
        rawtest = xtest[:, 1:] * basis['scale'] + basis['mean']
        outside = np.any((rawtest < basis['min']) | (rawtest > basis['max']), axis=1)
        zmax = np.max(np.abs(xtest[:, 1:]), axis=1)
    else:
        outside = np.zeros(len(xtest), bool)
        zmax = np.zeros(len(xtest), float)
    out.update({'test_eta':eta, 'test_raw_eta':raw_eta, 'test_max_abs_train_z':zmax,
                'test_outside_train_feature_range':outside,
                'n_test_outside_train_feature_range':int(outside.sum()),
                'max_test_abs_train_z':float(zmax.max()),
                'n_test_abs_train_z_gt_3':int((zmax > 3.).sum()),
                'max_test_linear_predictor_variance':float(np.max(np.einsum(
                    'ij,jk,ik->i', xtest, fitted['cov'], xtest)))})
    return out


def paired_gains(y, baseline_eta, augmented_eta):
    """Stable binary scores on the exact common successful prediction set."""
    y = np.asarray(y, float)
    baseline = np.asarray(baseline_eta, float)
    augmented = np.asarray(augmented_eta, float)
    keep = np.isfinite(y) & np.isfinite(baseline) & np.isfinite(augmented)
    y, baseline, augmented = y[keep], baseline[keep], augmented[keep]
    base_loss = np.logaddexp(0., baseline) - y * baseline
    augmented_loss = np.logaddexp(0., augmented) - y * augmented
    return keep, base_loss - augmented_loss, (expit(baseline) - y) ** 2 - (expit(augmented) - y) ** 2


def cluster_weighted_means(values, genus_codes, weights):
    """Resample genera, retain every member species and equal-species means."""
    values = np.asarray(values, float)
    codes = np.asarray(genus_codes, int)
    weights = np.asarray(weights)
    sums = np.zeros((weights.shape[1], values.shape[1]))
    counts = np.zeros_like(sums)
    np.add.at(sums, codes, np.nan_to_num(values, nan=0.))
    np.add.at(counts, codes, np.isfinite(values))
    numerator = weights @ sums
    denominator = weights @ counts
    return np.divide(numerator, denominator, out=np.full_like(numerator, np.nan), where=denominator > 0)


def summarize_species_metrics(metrics, panel_species, key_columns, metric_names, family_size, seed):
    """Shared genus draws across all contrasts, including differing coverage."""
    species = panel_species.sort_values('inat_taxon_id').reset_index(drop=True)
    genera, codes = np.unique(species.genus.to_numpy(str), return_inverse=True)
    rng = np.random.default_rng(seed)
    weights = rng.multinomial(len(genera), np.full(len(genera), 1. / len(genera)), size=BOOTSTRAP_REPETITIONS)
    records = []
    draw_columns = []
    draw_names = []
    for group, rows in metrics.groupby(key_columns, sort=True):
        group = group if isinstance(group, tuple) else (group,)
        key = dict(zip(key_columns, group))
        joined = species[['inat_taxon_id', 'genus']].merge(
            rows[['inat_taxon_id', *metric_names]], on='inat_taxon_id', how='left', validate='one_to_one')
        draws = cluster_weighted_means(joined[metric_names].to_numpy(float), codes, weights)
        record = {**key, 'n_species':len(rows), 'n_genera':rows.genus.nunique(),
                  'n_photos':int(rows.n_photos.sum()), 'bootstrap_repetitions':BOOTSTRAP_REPETITIONS,
                  'bootstrap_master_genera':len(genera), 'bootstrap_seed':seed,
                  'primary_comparison_family_size':family_size}
        for j, metric in enumerate(metric_names):
            estimate = float(rows[metric].mean())
            values = draws[:, j]
            nominal = np.nanquantile(values, [.025, .975])
            family = np.nanquantile(values, [.025 / family_size, 1. - .025 / family_size])
            record.update({metric:estimate, metric + '_ci95_low':float(nominal[0]),
                           metric + '_ci95_high':float(nominal[1]),
                           metric + f'_family{family_size}_low':float(family[0]),
                           metric + f'_family{family_size}_high':float(family[1]),
                           metric + ('_above_half_species_fraction' if metric == 'pair_auc' else '_positive_species_fraction'):float((rows[metric] > (.5 if metric == 'pair_auc' else 0.)).mean()),
                           metric + '_bootstrap_finite_draws':int(np.isfinite(values).sum())})
            draw_columns.append(values)
            draw_names.append('|'.join(map(str, group)) + '|' + metric)
        records.append(record)
    return pd.DataFrame(records), np.column_stack(draw_columns), np.array(draw_names)


def run_genus_transfer(panel, slopes, outdir):
    coefficients = pd.concat([learn_fold_coefficients(slopes, f) for f in range(10)], ignore_index=True)
    coefficients.to_csv(outdir / 'genus_fold_coefficients.csv', index=False)
    scores = np.full((len(panel), 2), np.nan)
    records = []
    for tid, g in panel.groupby('inat_taxon_id', sort=True):
        fold = int(g.genus_fold.iloc[0])
        beta = coefficients.loc[coefficients.heldout_genus_fold.eq(fold)].set_index('variable').loc[list(PREDICTORS), 'coefficient'].to_numpy(float)
        units = g[list(PREDICTORS)].to_numpy(float) / INCREMENTS
        # Score location is arbitrary for conditional pairs. This fixed origin is
        # physical zero, and neither test labels nor test prevalence enter it.
        current = np.column_stack([units[:, 0] * beta[0], units @ beta])
        scores[g.index] = current
        white = g.white.to_numpy(bool)
        for j, model in enumerate(('BIO5_component', 'joint_three')):
            records.append({'inat_taxon_id':int(tid), 'species':g.species.iloc[0],
                            'genus':g.genus.iloc[0], 'source_cohort':g.fcp_source_cohort.iloc[0],
                            'heldout_genus_fold':fold, 'model':model, 'n_photos':len(g),
                            'n_white':int(white.sum()), 'n_chromatic':int((~white).sum()),
                            **pairwise_metrics(current[white, j], current[~white, j])})
    demand(np.isfinite(scores).all(), 'incomplete genus OOF scores')
    metrics = pd.DataFrame(records)
    metrics.to_csv(outdir / 'genus_species_metrics.csv', index=False)
    np.savez_compressed(outdir / 'genus_oof_scores.npz', photo_index=panel.index.to_numpy(np.int32),
                        models=np.array(['BIO5_component', 'joint_three']), score=scores)
    species = panel[['inat_taxon_id', 'genus']].drop_duplicates()
    summary, draws, names = summarize_species_metrics(
        metrics, species, ['model'], ['pair_logloss_gain', 'pair_auc'], 2, GENUS_BOOTSTRAP_SEED)
    summary.to_csv(outdir / 'genus_summary.csv', index=False)
    np.savez_compressed(outdir / 'genus_bootstrap_draws.npz', names=names, means=draws)
    return summary


def flatten_fit_diagnostics(result):
    arrays = {'beta', 'cov', 'basis', 'test_eta', 'test_raw_eta', 'test_max_abs_train_z', 'test_outside_train_feature_range'}
    record = {k:v for k, v in result.items() if k not in arrays}
    basis = result.get('basis', {})
    if 'training_climate_sd' in basis:
        for variable, sd in zip(PREDICTORS, basis['training_climate_sd']):
            record['train_sd_' + variable] = float(sd)
    if 'training_geographic_sd' in basis:
        for variable, sd in zip(('latitude', 'unwrapped_longitude'), basis['training_geographic_sd']):
            record['train_sd_' + variable] = float(sd)
    for key in ('longitude_origin_radians', 'longitude_circular_resultant', 'unwrapped_longitude_range_degrees'):
        if key in basis:
            record[key] = basis[key]
    if 'names' in basis:
        for j, name in enumerate(basis['names']):
            for label, attr in (('mean', 'mean'), ('sd', 'scale'), ('min', 'min'), ('max', 'max')):
                record[f'train_{label}_{name}'] = float(basis[attr][j])
    if 'beta' in result:
        for name, beta in zip(['intercept', *basis.get('names', [])], result['beta']):
            record['coefficient_' + name] = float(beta)
    return record


def run_spatial_transfer(panel, outdir):
    n = len(panel)
    etas = np.full((len(BUFFERS_KM), n, len(MODELS)), np.nan)
    raw_etas = np.full_like(etas, np.nan)
    zmax = np.full_like(etas, np.nan)
    outside = np.full(etas.shape, -1, dtype=np.int8)
    nearest = np.full((len(BUFFERS_KM), n), np.nan)
    train_n = np.zeros((len(BUFFERS_KM), n), dtype=np.int16)
    diagnostics = []
    fold_metrics = []
    fold_selection = []
    start = time.monotonic()
    for species_number, (tid, g) in enumerate(panel.groupby('inat_taxon_id', sort=True), 1):
        climate = g[list(PREDICTORS)].to_numpy(float)
        coords = g[['latitude', 'longitude']].to_numpy(float)
        y = g.white.to_numpy(float)
        folds = g.cell_fold.to_numpy(int)
        global_indices = g.index.to_numpy(int)
        common = {'inat_taxon_id':int(tid), 'species':g.species.iloc[0], 'genus':g.genus.iloc[0],
                  'source_cohort':g.fcp_source_cohort.iloc[0], 'n_species_panel':len(g)}
        for bi, buffer_km in enumerate(BUFFERS_KM):
            for fold in range(5):
                train, test, train_dist = spatial_split_indices(coords, folds, fold, buffer_km)
                base = {**common, 'buffer_km':buffer_km, 'heldout_cell_fold':fold,
                        'n_test':len(test), 'n_train_before_buffer':int((folds != fold).sum()),
                        'n_train_retained':len(train), 'n_train_removed_by_buffer':int((folds != fold).sum()) - len(train),
                        'n_train_white':int(y[train].sum()), 'n_test_white':int(y[test].sum()),
                        'constant_training_response':bool(len(train) > 0 and (np.all(y[train] == 0) or np.all(y[train] == 1))),
                        'min_train_to_test_distance_km':float(train_dist.min()) if len(train_dist) and len(test) else np.nan}
                if len(test) and len(train):
                    test_distance = great_circle_km(coords[test], coords[train]).min(axis=1)
                    nearest[bi, global_indices[test]] = test_distance
                    base.update({'median_test_to_train_distance_km':float(np.median(test_distance)),
                                 'max_test_to_train_distance_km':float(test_distance.max())})
                fold_selection.append({**base, 'selection_status':'no_test_photos' if not len(test) else
                                       ('training_n_below_20' if len(train) < 20 else 'training_count_eligible')})
                train_n[bi, global_indices[test]] = len(train)
                if len(test):
                    demand(set(g.cell_id.iloc[train]).isdisjoint(g.cell_id.iloc[test]), 'shared heldout cell')
                for mi, model in enumerate(MODELS):
                    if not len(test):
                        diagnostics.append({**base, 'model':model, 'status':'no_test_photos'})
                        continue
                    result = fit_one_model(climate[train], coords[train], y[train], climate[test], coords[test], model)
                    diagnostics.append({**base, **flatten_fit_diagnostics(result)})
                    if result['status'] != 'ok':
                        continue
                    etas[bi, global_indices[test], mi] = result['test_eta']
                    raw_etas[bi, global_indices[test], mi] = result['test_raw_eta']
                    zmax[bi, global_indices[test], mi] = result['test_max_abs_train_z']
                    outside[bi, global_indices[test], mi] = result['test_outside_train_feature_range'].astype(np.int8)
                for comparison, baseline, augmented in COMPARISONS:
                    keep, loss, brier = paired_gains(y[test], etas[bi, global_indices[test], baseline], etas[bi, global_indices[test], augmented])
                    if keep.any():
                        fold_metrics.append({**base, 'comparison':comparison, 'n_matched_test':int(keep.sum()),
                                             'logloss_gain_sum':float(loss.sum()), 'brier_gain_sum':float(brier.sum()),
                                             'logloss_gain':float(loss.mean()), 'brier_gain':float(brier.mean())})
        if species_number % 100 == 0:
            success = sum(d['status'] == 'ok' for d in diagnostics)
            print('SPATIAL_SPECIES', species_number, 'successful_fits', success,
                  'elapsed_seconds', round(time.monotonic() - start, 1), flush=True)
    pd.DataFrame(diagnostics).to_csv(outdir / 'spatial_fit_diagnostics.csv.gz', index=False,
                                    compression={'method':'gzip', 'mtime':0})
    pd.DataFrame(fold_selection).to_csv(outdir / 'spatial_fold_selection.csv.gz', index=False,
                                       compression={'method':'gzip', 'mtime':0})
    pd.DataFrame(fold_metrics).to_csv(outdir / 'spatial_fold_metrics.csv.gz', index=False,
                                     compression={'method':'gzip', 'mtime':0})
    np.savez_compressed(outdir / 'spatial_oof_predictions.npz',
                        photo_index=panel.index.to_numpy(np.int32), buffers_km=np.array(BUFFERS_KM),
                        models=np.array(MODELS), eta=etas, probability=expit(etas),
                        raw_firth_eta=raw_etas,
                        max_abs_train_z=zmax, outside_train_feature_range=outside,
                        nearest_retained_training_km=nearest, n_training_photos=train_n)
    records = []
    for tid, g in panel.groupby('inat_taxon_id', sort=True):
        indices = g.index.to_numpy(int)
        y = g.white.to_numpy(float)
        for bi, buffer_km in enumerate(BUFFERS_KM):
            for comparison, baseline, augmented in COMPARISONS:
                keep, loss, brier = paired_gains(y, etas[bi, indices, baseline], etas[bi, indices, augmented])
                if not keep.any():
                    continue
                eta0, eta1 = etas[bi, indices[keep], baseline], etas[bi, indices[keep], augmented]
                target = y[keep]
                records.append({'inat_taxon_id':int(tid), 'species':g.species.iloc[0], 'genus':g.genus.iloc[0],
                                'source_cohort':g.fcp_source_cohort.iloc[0], 'buffer_km':buffer_km,
                                'comparison':comparison, 'n_photos':int(keep.sum()), 'n_white':int(target.sum()),
                                'n_evaluable_folds':int(g.cell_fold.iloc[np.flatnonzero(keep)].nunique()),
                                'logloss_gain':float(loss.mean()), 'brier_gain':float(brier.mean()),
                                'baseline_logloss':float(np.mean(np.logaddexp(0., eta0) - target * eta0)),
                                'augmented_logloss':float(np.mean(np.logaddexp(0., eta1) - target * eta1)),
                                'baseline_brier':float(np.mean((expit(eta0) - target) ** 2)),
                                'augmented_brier':float(np.mean((expit(eta1) - target) ** 2)),
                                'augmented_outside_train_feature_range_photos':int((outside[bi, indices[keep], augmented] == 1).sum()),
                                'augmented_abs_train_z_gt_3_photos':int((zmax[bi, indices[keep], augmented] > 3.).sum()),
                                'max_augmented_abs_train_z':float(np.max(zmax[bi, indices[keep], augmented])),
                                'min_test_to_train_km':float(np.min(nearest[bi, indices[keep]])),
                                'median_test_to_train_km':float(np.median(nearest[bi, indices[keep]]))})
    metrics = pd.DataFrame(records)
    metrics.to_csv(outdir / 'spatial_species_metrics.csv', index=False)
    species = panel[['inat_taxon_id', 'genus']].drop_duplicates()
    summary, draws, names = summarize_species_metrics(
        metrics, species, ['buffer_km', 'comparison'], ['logloss_gain', 'brier_gain'], 4, SPATIAL_BOOTSTRAP_SEED)
    for i, row in summary.iterrows():
        subset = metrics.loc[metrics.buffer_km.eq(row.buffer_km) & metrics.comparison.eq(row.comparison)]
        summary.loc[i, 'outside_training_range_fraction_photos'] = subset.augmented_outside_train_feature_range_photos.sum() / subset.n_photos.sum()
        summary.loc[i, 'abs_train_z_gt_3_fraction_photos'] = subset.augmented_abs_train_z_gt_3_photos.sum() / subset.n_photos.sum()
        summary.loc[i, 'baseline_logloss_equal_species'] = subset.baseline_logloss.mean()
        summary.loc[i, 'augmented_logloss_equal_species'] = subset.augmented_logloss.mean()
    summary.to_csv(outdir / 'spatial_summary.csv', index=False)
    np.savez_compressed(outdir / 'spatial_bootstrap_draws.npz', names=names, means=draws)
    return summary, pd.DataFrame(diagnostics)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--photos', type=Path, required=True)
    parser.add_argument('--worldclim-root', type=Path, required=True)
    parser.add_argument('--slopes', type=Path, default=Path('results/fcp_real_447500_environment_model_20261018/species_slopes.csv'))
    parser.add_argument('--protocol', type=Path, default=Path('docs/FCP_ENVIRONMENT_TRANSPORT_PROTOCOL_20261019.md'))
    parser.add_argument('--outdir', type=Path, required=True)
    args = parser.parse_args()
    started = datetime.now(timezone.utc).isoformat()
    start_clock = time.monotonic()
    demand(args.protocol.is_file(), 'frozen protocol missing')
    protocol_sha256 = sha(args.protocol)
    args.outdir.mkdir(parents=True, exist_ok=True)
    panel, raster_meta, missing = load_panel(args.photos, args.worldclim_root)
    panel = panel.sort_values(['inat_taxon_id', 'photo_id'], kind='stable').reset_index(drop=True)
    panel['genus'] = panel.species.str.split().str[0]
    panel['genus_fold'] = panel.genus.map(lambda g:stable_fold(g, GENUS_SEED, 10))
    panel['cell_id'], panel['cell_fold'] = cell_folds(panel.latitude, panel.longitude)
    panel.index.name = 'photo_index'
    slopes = pd.read_csv(args.slopes)
    slopes = slopes.loc[slopes.model.eq('primary')].copy()
    demand(slopes.groupby('inat_taxon_id').size().eq(3).all(), 'primary slope completeness')
    demand(set(slopes.inat_taxon_id) == set(panel.inat_taxon_id), 'panel and slopes identities differ')
    demand(slopes.set_index('inat_taxon_id').species.to_dict() == panel.set_index('inat_taxon_id').species.to_dict(), 'species name mismatch')
    photo_columns = ['photo_id', 'observation_id', 'inat_taxon_id', 'species', 'genus', 'fcp_source_cohort',
                     'white', 'latitude', 'longitude', *PREDICTORS, 'genus_fold', 'cell_id', 'cell_fold']
    panel[photo_columns].to_csv(args.outdir / 'photo_index.csv.gz', index=True,
                                compression={'method':'gzip', 'mtime':0})
    print('PANEL_READY', panel.inat_taxon_id.nunique(), len(panel), 'genera', panel.genus.nunique(), flush=True)
    genus_start = time.monotonic()
    genus_summary = run_genus_transfer(panel, slopes, args.outdir)
    genus_seconds = time.monotonic() - genus_start
    print('GENUS_TRANSFER_COMPLETE', flush=True)
    print(genus_summary[['model', 'n_species', 'pair_logloss_gain', 'pair_auc']].to_string(index=False), flush=True)
    spatial_start = time.monotonic()
    spatial_summary, diagnostics = run_spatial_transfer(panel, args.outdir)
    spatial_seconds = time.monotonic() - spatial_start
    profile = {
        'schema':'fcp_environment_transport_20261019_v1',
        'status':'complete_post_hoc_internal_transport_evaluation',
        'n_panel_species':panel.inat_taxon_id.nunique(), 'n_panel_photos':len(panel),
        'n_panel_genera':panel.genus.nunique(), 'n_global_1deg_cells':panel.cell_id.nunique(),
        'missing_any_raster_classifiable_photos':missing,
        'genus_seed':GENUS_SEED, 'genus_folds':10, 'cell_seed':CELL_SEED, 'cell_folds':5,
        'fold_hash':'SHA256(UTF8(seed + colon + key)); whole digest big-endian integer modulo fold count',
        'cell_key':'floor(latitude),floor(canonical longitude in [-180,180)); north pole belongs to latitude cell 89',
        'genus_key':'exact first whitespace-separated species-name word',
        'predictors':list(PREDICTORS), 'physical_increments':INCREMENTS.tolist(),
        'buffers_km':list(BUFFERS_KM), 'buffer_exclusion':'exclude training photos at distance <= buffer around any conspecific test photo; buffer 0 has no additional distance filter',
        'distance':'haversine, mean Earth radius 6371.0088 km',
        'spatial_models':list(MODELS), 'minimum_training_photos':20,
        'constant_training_response_is_eligible':True, 'adaptive_predictor_dropping':False,
        'spatial_recalibration':{'method':'pre-result common smoothed train-only intercept recalibration, not standard FLIC',
                                 'target':'(training_white + 0.5)/(training_n + 1)',
                                 'slopes_unchanged':True, 'test_information_used':False,
                                 'raw_and_adjusted_intercepts_saved':True,
                                 'all_spatial_scores_use_recalibrated_predictions':True},
        'genus_coefficient_training':'marginal REML of training-genus primary joint species slopes in common physical units',
        'genus_primary_metric':'equal-species mean of log(2) minus conditional loss over every white/chromatic pair',
        'genus_secondary_metric':'equal-species mean within-species pair AUC with half ties',
        'spatial_primary_metric':'binary log loss baseline minus augmented; exact matched successful predictions; pooled test photos within species then equal species weight',
        'spatial_secondary_metric':'same paired Brier-score gain',
        'bootstrap':{'repetitions':BOOTSTRAP_REPETITIONS, 'genus_seed':GENUS_BOOTSTRAP_SEED,
                     'spatial_seed':SPATIAL_BOOTSTRAP_SEED, 'sampling':'727 panel genus clusters sampled with replacement; retain all member species, equal species weights',
                     'shared_draws_across_contrasts':True, 'nominal_quantiles':[.025, .975],
                     'family2_quantiles':[.0125, .9875], 'family4_quantiles':[.00625, .99375],
                     'conditional_on_fixed_OOF_fits':True, 'refits_inside_bootstrap':False},
        'prediction_storage':{'join_key':'photo_index shared across NPZ and photo_index.csv.gz',
                              'spatial_array_axes':['buffer', 'photo_index', 'model'],
                              'genus_array_axes':['photo_index', 'model'],
                              'failed_model_prediction':'NaN', 'outside_range_missing':-1,
                              'photo_ids_are_actual_source_ids':True},
        'spatial_fit_status_counts':diagnostics.groupby(['buffer_km', 'model', 'status']).size().rename('n').reset_index().to_dict(orient='records'),
        'constant_response_fit_status_counts':diagnostics.loc[diagnostics.constant_training_response.fillna(False)].groupby(['buffer_km', 'model', 'status']).size().rename('n').reset_index().to_dict(orient='records'),
        'genus_summary':genus_summary.to_dict(orient='records'),
        'spatial_summary':spatial_summary.to_dict(orient='records'),
        'limitations':[
            'The panel and predictor choice were previously examined; this is post hoc internal validation.',
            'Held-out-genera scores discriminate photos within selected species and do not calibrate absolute white frequency.',
            'Genus-cluster percentile intervals condition on fitted cross-validation; cross-fold dependence and training uncertainty remain.',
            'Genus blocking does not constitute complete phylogenetic control.',
            'Buffered spatial comparisons have changing coverage and environmental extrapolation.',
            'Photos, classification, spatial sampling, cultivation, and observer biases remain; no adaptation, pigment-cost, or fitness inference.',
        ],
    }
    (args.outdir / 'profile.json').write_text(json.dumps(json_safe(profile), indent=2, sort_keys=True, allow_nan=False) + '\n')
    outputs = sorted(p for p in args.outdir.iterdir() if p.is_file() and p.name != 'execution_receipt.json')
    demand(sha(args.protocol) == protocol_sha256, 'protocol changed during execution')
    receipt = {'started_utc':started, 'completed_utc':datetime.now(timezone.utc).isoformat(),
               'command':sys.argv, 'python':sys.version, 'platform':platform.platform(),
               'packages':{'numpy':np.__version__, 'pandas':pd.__version__, 'scipy':scipy.__version__},
               'source_photo_sha256':PHOTO_SHA, 'worldclim_zip_sha256':ZIP_SHA, 'raster_metadata':raster_meta,
               'slopes_sha256':sha(args.slopes), 'protocol_sha256':sha(args.protocol), 'code_sha256':sha(Path(__file__)),
               'imported_environment_model_code_sha256':sha(Path('scripts/analysis/run_fcp_447500_environment_model_20261018.py')),
               'imported_panel_sampler_code_sha256':sha(Path('scripts/analysis/summarize_fcp_447500_within_species_environment_20261018.py')),
               'seconds':{'genus_analysis':genus_seconds, 'spatial_analysis':spatial_seconds,
                          'total':time.monotonic() - start_clock},
               'output_sha256':{p.name:sha(p) for p in outputs}}
    (args.outdir / 'execution_receipt.json').write_text(json.dumps(json_safe(receipt), indent=2, sort_keys=True, allow_nan=False) + '\n')
    print(spatial_summary[['buffer_km', 'comparison', 'n_species', 'n_photos', 'logloss_gain', 'brier_gain']].to_string(index=False), flush=True)
    print('FCP_ENVIRONMENT_TRANSPORT_COMPLETE', flush=True)


if __name__ == '__main__':
    main()

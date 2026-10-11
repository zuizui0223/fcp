#!/usr/bin/env python3
"""Independent saved-prediction verification; imports no analysis implementation.

Reconstruct losses, species/genus estimands and the primary paired genus
bootstrap using only NumPy/Pandas and the archived parent CSVs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

CLASSES = ('white', 'yellow_orange', 'red_pink', 'blue_purple')
MODELS = ('ridge_geo', 'ridge_full', 'nonlinear_geo', 'nonlinear_full')
COMPONENTS = ('total', 'white_axis', 'remaining_hue')
SCHEMES = ('genus_only', 'sector_purged', 'buffer500')
CASES = ('soil_all', 'soil_quality')
KEYS = ['inat_taxon_id', 'sector', 'fold', 'scheme']
LOSS_COLS = [m + '__' + c for m in MODELS for c in COMPONENTS]
CONTRASTS = {
    'nonlinear_environment_gain': (0, 0, 1, -1),
    'ridge_environment_gain': (1, -1, 0, 0),
    'nonlinear_minus_ridge_environment_gain': (-1, 1, 1, -1),
    'full_model_improvement': (0, 1, 0, -1),
    'geography_model_improvement': (1, 0, -1, 0),
}
TOL = 1e-11


class Verification:
    def __init__(self):
        self.receipt = {'status': 'RUNNING', 'tolerance': TOL, 'cases': {},
                        'checks': 0, 'max_abs_errors': {}}

    def require(self, condition, label):
        self.receipt['checks'] += 1
        if not condition:
            raise ValueError(label)

    def close(self, actual, expected, label, category):
        a, b = np.asarray(actual, float), np.asarray(expected, float)
        self.require(a.shape == b.shape, label + ': shape mismatch')
        self.require(np.isfinite(a).all() and np.isfinite(b).all(), label + ': nonfinite')
        error = float(np.max(np.abs(a - b))) if a.size else 0.0
        old = self.receipt['max_abs_errors'].get(category, 0.0)
        self.receipt['max_abs_errors'][category] = max(old, error)
        self.require(error <= TOL, f'{label}: max absolute error {error:.17g}')


def read_unique(path, keys, v):
    data = pd.read_csv(path)
    v.require(not data.empty, str(path) + ': empty')
    v.require(not data[keys].isna().any().any(), str(path) + ': missing identity')
    v.require(not data.duplicated(keys).any(), str(path) + ': duplicate identity')
    return data.set_index(keys).sort_index()


def compare_parent(rows, parent, case, v):
    old = read_unique(parent / (case + '_oof_losses.csv.gz'), KEYS, v)
    new = rows.set_index(KEYS).sort_index()
    v.require(new.index.equals(old.index), case + ': parent target/fold/scheme mismatch')
    for col in ('genus', 'n_photos', 'cell'):
        v.require(new[col].equals(old[col]), case + ': parent ' + col + ' mismatch')
    for model, previous in [('ridge_geo', 'geography'), ('ridge_full', 'environment36')]:
        for component in COMPONENTS:
            v.close(new[model + '__' + component], old[previous + '__' + component],
                    case + ': archived ' + model + ' ' + component, 'archived_ridge')
    return len(old)


def validate_predictions(rows, case, v):
    v.require(set(rows.scheme.unique()) == set(SCHEMES), case + ': incorrect schemes')
    v.require(set(rows.sector.unique()) == set(range(4)), case + ': incorrect sectors')
    v.require(set(rows.fold.unique()) == set(range(5)), case + ': incorrect folds')
    v.require(not rows.duplicated(['inat_taxon_id', 'sector', 'scheme']).any(),
              case + ': duplicate species-sector-scheme')
    v.require(not rows[['inat_taxon_id', 'genus', 'sector', 'fold', 'scheme', 'n_photos', 'cell']].isna().any().any(),
              case + ': missing target identity')
    v.require((rows.n_photos > 0).all() and (rows.n_photos == rows.n_photos.astype(int)).all(),
              case + ': invalid photo counts')
    v.require(rows.groupby('inat_taxon_id').genus.nunique().max() == 1,
              case + ': one species has multiple genera')
    v.require(rows.groupby('genus').fold.nunique().max() == 1,
              case + ': one genus has multiple folds')
    truth_cols = ['truth__' + c for c in CLASSES]
    fixed = ['genus', 'n_photos', 'cell', 'fold', *truth_cols]
    base = rows[rows.scheme.eq(SCHEMES[0])].set_index(['inat_taxon_id', 'sector']).sort_index()
    for scheme in SCHEMES[1:]:
        other = rows[rows.scheme.eq(scheme)].set_index(['inat_taxon_id', 'sector']).sort_index()
        v.require(base.index.equals(other.index), case + ': different targets across schemes')
        v.require(base[fixed].equals(other[fixed]), case + ': changed truth/photos/fold across schemes')
    truth = rows[truth_cols].to_numpy(float)
    for name, prob in [('truth', truth)] + [(m, rows[[m + '__p__' + c for c in CLASSES]].to_numpy(float)) for m in MODELS]:
        v.require(np.isfinite(prob).all(), case + ': nonfinite ' + name)
        v.require((prob >= -TOL).all() and (prob <= 1 + TOL).all(), case + ': out-of-simplex ' + name)
        v.close(prob.sum(axis=1), np.ones(len(rows)), case + ': ' + name + ' normalization', 'simplex')
        if name == 'truth':
            continue
        residual = truth - prob
        total = np.einsum('ij,ij->i', residual, residual)
        white = (4.0 / 3.0) * residual[:, 0] ** 2
        for component, loss in zip(COMPONENTS, (total, white, total - white)):
            column = name + '__' + component
            v.close(rows[column].to_numpy(float), loss, case + ': ' + column, 'saved_loss')
            rows[column] = loss  # All downstream estimates use the recomputed losses.
    return base


def group_values(rows):
    species = rows.groupby('inat_taxon_id', sort=True).agg({**dict.fromkeys(LOSS_COLS, 'mean'), 'genus': 'first'})
    genera = species.groupby('genus', sort=True)
    sums = genera[LOSS_COLS].sum().to_numpy(float)
    counts = genera.size().to_numpy(float)
    means = sums / counts[:, None]
    return species, sums, counts, means


def metric_values(losses):
    # Supports a point vector or a bootstrap matrix.
    shaped = losses.reshape(*losses.shape[:-1], len(MODELS), len(COMPONENTS))
    values = {('losses', model): shaped[..., i, :] for i, model in enumerate(MODELS)}
    for name, coefficients in CONTRASTS.items():
        values[('contrasts', name)] = np.einsum('...mc,m->...c', shaped, coefficients)
    return values


def verify_summary(rows, summary, label, v, bootstrap=False):
    species, sums, counts, means = group_values(rows)
    v.require(summary.get('status', 'OK') == 'OK', label + ': summary not OK')
    for field, expected in [('n_species', len(species)), ('n_genera', len(counts)),
                            ('n_species_sector_records', len(rows)), ('n_sectors', rows.sector.nunique())]:
        if field in summary:
            v.require(summary[field] == expected, label + ': ' + field + ' mismatch')
    for weight, point in [('genus_equal', means.mean(axis=0)),
                           ('species_equal', species[LOSS_COLS].mean().to_numpy(float))]:
        values = metric_values(point)
        for i, component in enumerate(COMPONENTS):
            reported = summary[weight][component]
            v.require(set(reported['losses']) == set(MODELS), label + ': model summary mismatch')
            v.require(set(reported['contrasts']) == set(CONTRASTS), label + ': contrast summary mismatch')
            for (kind, name), value in values.items():
                v.close(reported[kind][name]['estimate'], value[i],
                        label + ': ' + weight + '/' + component + '/' + name, 'point_estimates')
    if not bootstrap:
        return
    v.require('bootstrap_seed' in summary, label + ': bootstrap seed absent')
    v.require(summary.get('bootstrap_replicates') == 1999, label + ': bootstrap must be 1999')
    rng = np.random.default_rng(summary['bootstrap_seed'])
    draws = {'genus_equal': [], 'species_equal': []}
    for start in range(0, 1999, 32):
        ix = rng.integers(len(counts), size=(min(32, 1999 - start), len(counts)))
        draws['genus_equal'].append(means[ix].mean(axis=1))
        draws['species_equal'].append(sums[ix].sum(axis=1) / counts[ix].sum(axis=1)[:, None])
    for weight, batches in draws.items():
        values = metric_values(np.concatenate(batches))
        for i, component in enumerate(COMPONENTS):
            for (kind, name), value in values.items():
                for ci, quantiles in [('CI95', [.025, .975]), ('CI97_5', [.0125, .9875])]:
                    v.close(summary[weight][component][kind][name][ci], np.quantile(value[:, i], quantiles),
                            label + ': ' + weight + '/' + component + '/' + name + '/' + ci, 'primary_bootstrap')


def verify_audits(rows, folder, parent, case, metadata, v):
    audit_keys = ['sector', 'fold', 'scheme']
    audit = read_unique(folder / (case + '_fold_audit.csv'), audit_keys, v)
    old_audit = read_unique(parent / (case + '_fold_audit.csv'), audit_keys, v)
    expected_keys = rows[audit_keys].drop_duplicates().set_index(audit_keys).sort_index().index
    v.require(audit.index.equals(expected_keys) and audit.index.equals(old_audit.index),
              case + ': audit/fold identity mismatch')
    for column in ['shared_train_test_genus', 'shared_train_test_species', 'training_species_intersect_spatial_blacklist']:
        v.require(audit[column].eq(0).all(), case + ': nonzero leakage audit ' + column)
    v.require(audit.status.eq('OK').all(), case + ': unsuccessful fold audit')
    for column in ['n_train_species', 'n_train_genera', 'n_target_species', 'n_target_photos']:
        v.require(audit[column].equals(old_audit[column]), case + ': parent audit ' + column + ' mismatch')
    target_counts = rows.groupby(audit_keys, sort=True).agg(n_target_species=('inat_taxon_id', 'size'), n_target_photos=('n_photos', 'sum'))
    for column in target_counts:
        v.require(audit[column].equals(target_counts[column]), case + ': target audit ' + column + ' mismatch')
    fold_map = read_unique(folder / (case + '_genus_fold_map.csv'), ['genus'], v)
    parent_map = read_unique(parent / (case + '_genus_fold_map.csv'), ['genus'], v)
    actual = rows.groupby('genus', sort=True).fold.first().astype(int)
    v.require(fold_map.equals(parent_map), case + ': parent fold map mismatch')
    v.require(fold_map.fold.equals(actual), case + ': prediction fold map mismatch')
    encoded = json.dumps(actual.to_dict(), sort_keys=True, ensure_ascii=True, separators=(',', ':')).encode()
    digest = hashlib.sha256(encoded).hexdigest()
    v.require(metadata['genus_fold_map_sha256'] == digest, case + ': fold map hash mismatch')
    return len(audit), digest


def verify(result_dir, parent_results, v):
    result = json.loads((result_dir / 'result.json').read_text())
    parent_result = json.loads((parent_results / 'result.json').read_text())
    v.require(set(result['cases']) == set(CASES), 'unexpected or missing cases')
    v.require(result.get('bootstrap') == 1999, 'result bootstrap must be 1999')
    for case in CASES:
        rows = pd.read_csv(result_dir / (case + '_oof_predictions.csv.gz'))
        base = validate_predictions(rows, case, v)
        archived_matches = compare_parent(rows, parent_results, case, v)
        metadata = result['cases'][case]
        for field, actual in [('n_input_photos', int(base.n_photos.sum())),
                              ('n_input_species', int(base.index.get_level_values('inat_taxon_id').nunique())),
                              ('n_input_genera', int(base.genus.nunique()))]:
            v.require(metadata[field] == actual == parent_result['cases'][case][field],
                      case + ': cohort ' + field + ' mismatch')
        v.require(metadata['genus_fold_map_sha256'] == parent_result['cases'][case]['genus_fold_map_sha256'],
                  case + ': changed parent fold-map hash')
        audits, digest = verify_audits(rows, result_dir, parent_results, case, metadata, v)
        v.require(set(metadata['schemes']) == set(SCHEMES), case + ': summary schemes mismatch')
        summaries = 0
        for scheme in SCHEMES:
            subset = rows.loc[rows.scheme.eq(scheme)]
            schema = metadata['schemes'][scheme]
            label = case + '/' + scheme
            verify_summary(subset, schema['overall'], label + '/overall', v,
                           bootstrap=(case == 'soil_all' and scheme == 'buffer500'))
            v.require(set(schema['regions']) == {'0', '1', '2', '3'}, label + ': summary sectors mismatch')
            summaries += 1
            for sector in range(4):
                verify_summary(subset.loc[subset.sector.eq(sector)], schema['regions'][str(sector)],
                               label + '/sector' + str(sector), v)
                summaries += 1
        v.receipt['cases'][case] = {'prediction_rows': len(rows), 'species_sector_targets': len(base),
            'species': int(base.index.get_level_values('inat_taxon_id').nunique()),
            'genera': int(rows.genus.nunique()), 'archived_ridge_rows_matched': archived_matches,
            'archived_ridge_loss_values_matched': archived_matches * 6,
            'fold_audits': audits, 'fold_map_sha256': digest, 'point_summaries_verified': summaries}
    v.receipt['primary_bootstrap'] = {'case': 'soil_all', 'scheme': 'buffer500', 'replicates': 1999,
        'weights': ['genus_equal', 'species_equal'], 'components': list(COMPONENTS),
        'loss_models': list(MODELS), 'contrasts': list(CONTRASTS), 'intervals': ['CI95', 'CI97_5']}
    v.receipt['status'] = 'PASS'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--result-dir', required=True, type=Path)
    parser.add_argument('--parent-results', required=True, type=Path)
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    v = Verification()
    try:
        verify(args.result_dir, args.parent_results, v)
    except Exception as exc:
        v.receipt.update(status='FAIL', error=str(exc), error_type=type(exc).__name__)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(v.receipt, indent=2, sort_keys=True) + '\n')
    print(json.dumps(v.receipt, sort_keys=True))
    return 0 if v.receipt['status'] == 'PASS' else 1


if __name__ == '__main__':
    sys.exit(main())

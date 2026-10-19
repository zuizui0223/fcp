#!/usr/bin/env python3
"""Correct only numerical range-boundary diagnostics in completed FCP outputs.

Run from the repository root. Fits, selection, eta, probabilities, primary and
secondary scores, and bootstrap draws are retained exactly. Original scientific
execution bytes are archived in script_at_execution.py and authenticated by the
original receipt. The isolated photo-index repair is separately recorded there.
"""
from collections import defaultdict
from datetime import datetime, timezone
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path

import numpy as np
import pandas as pd

from scripts.analysis.run_fcp_environment_transport_20261019 import (
    MODELS, PREDICTORS, COMPARISONS, outside_training_range,
)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def atomic_blob(path, data):
    staged = path.with_name(path.name + '.correction_tmp')
    staged.write_bytes(data)
    assert staged.read_bytes() == data
    staged.replace(path)
    assert path.read_bytes() == data


def rewrite_columns(path, columns, update):
    original = gzip.decompress(path.read_bytes()) if path.suffix == '.gz' else path.read_bytes()
    reader = csv.DictReader(io.StringIO(original.decode()))
    stream = io.StringIO(newline='')
    writer = csv.DictWriter(stream, fieldnames=reader.fieldnames, lineterminator='\n')
    writer.writeheader()
    n = 0
    for row in reader:
        saved = dict(row)
        row.update(update(row))
        assert all(row[k] == value for k, value in saved.items() if k not in columns)
        writer.writerow(row)
        n += 1
    data = stream.getvalue().encode()
    if path.suffix == '.gz':
        data = gzip.compress(data, mtime=0)
        assert gzip.decompress(data).count(b'\n') == n + 1
    atomic_blob(path, data)
    return n


def main():
    root = Path('results/fcp_environment_transport_20261019')
    receipt_path = root / 'execution_receipt.json'
    receipt = json.loads(receipt_path.read_text())
    assert not receipt.get('range_diagnostic_correction'), 'correction already applied'
    assert digest((root / 'script_at_execution.py').read_bytes()) == receipt['code_sha256']
    for name, expected in receipt['output_sha256'].items():
        assert digest((root / name).read_bytes()) == expected, name
    old_artifact_hashes = dict(receipt['output_sha256'])
    photo = pd.read_csv(root / 'photo_index.csv.gz')
    diagnostics = pd.read_csv(root / 'spatial_fit_diagnostics.csv.gz', low_memory=False)
    npz_path = root / 'spatial_oof_predictions.npz'
    with np.load(npz_path) as source:
        arrays = {key:source[key] for key in source.files}
    old_component_hashes = {key:digest(value.tobytes()) for key, value in arrays.items()}
    indices = photo.groupby(['inat_taxon_id', 'cell_fold'], sort=False).indices
    climate = photo[list(PREDICTORS)].to_numpy(float)
    coords = photo[['latitude', 'longitude']].to_numpy(float)
    counts = defaultdict(lambda: {'old_outside_photos':0, 'new_outside_photos':0,
                                  'changed_flags':0, 'successful_test_photos':0})
    fold_counts = {}
    for row in diagnostics.loc[diagnostics.status.eq('ok')].itertuples(index=False):
        bi = 0 if row.buffer_km == 0. else 1
        mi = MODELS.index(row.model)
        idx = indices[(row.inat_taxon_id, row.heldout_cell_fold)]
        names = []
        if row.model in ('climate', 'geoclimate'):
            names.extend(PREDICTORS)
        if row.model in ('geography', 'geoclimate'):
            names.extend(('latitude', 'unwrapped_longitude'))
        basis = {'model':row.model,
                 'min':np.array([getattr(row, 'train_min_' + name) for name in names]),
                 'max':np.array([getattr(row, 'train_max_' + name) for name in names]),
                 'longitude_origin_radians':row.longitude_origin_radians}
        corrected = outside_training_range(basis, climate[idx], coords[idx])
        old = arrays['outside_train_feature_range'][bi, idx, mi]
        assert np.all(old >= 0)
        group = counts[(row.buffer_km, row.model)]
        group['old_outside_photos'] += int((old == 1).sum())
        group['new_outside_photos'] += int(corrected.sum())
        group['changed_flags'] += int((corrected != (old == 1)).sum())
        group['successful_test_photos'] += len(idx)
        arrays['outside_train_feature_range'][bi, idx, mi] = corrected.astype(np.int8)
        fold_counts[(row.buffer_km, row.inat_taxon_id, row.heldout_cell_fold, row.model)] = int(corrected.sum())
    stream = io.BytesIO()
    np.savez_compressed(stream, **arrays)
    atomic_blob(npz_path, stream.getvalue())
    with np.load(npz_path) as saved:
        component_hashes = {key:digest(saved[key].tobytes()) for key in saved.files}
    for key in old_component_hashes:
        if key != 'outside_train_feature_range':
            assert old_component_hashes[key] == component_hashes[key], key
    def update_fit(row):
        if row['status'] != 'ok':
            return {}
        key = (float(row['buffer_km']), int(row['inat_taxon_id']), int(row['heldout_cell_fold']), row['model'])
        return {'n_test_outside_train_feature_range':str(fold_counts[key])}
    rewrite_columns(root / 'spatial_fit_diagnostics.csv.gz', {'n_test_outside_train_feature_range'}, update_fit)
    by_species = photo.groupby('inat_taxon_id', sort=False).indices
    comparison_map = {name:(baseline, augmented) for name, baseline, augmented in COMPARISONS}
    summary_counts = defaultdict(lambda: [0, 0])
    def update_species(row):
        buffer = float(row['buffer_km'])
        bi = 0 if buffer == 0. else 1
        baseline, augmented = comparison_map[row['comparison']]
        idx = by_species[int(row['inat_taxon_id'])]
        keep = np.isfinite(arrays['eta'][bi, idx, baseline]) & np.isfinite(arrays['eta'][bi, idx, augmented])
        assert int(keep.sum()) == int(row['n_photos'])
        n = int((arrays['outside_train_feature_range'][bi, idx[keep], augmented] == 1).sum())
        key = (buffer, row['comparison'])
        summary_counts[key][0] += n
        summary_counts[key][1] += int(keep.sum())
        return {'augmented_outside_train_feature_range_photos':str(n)}
    rewrite_columns(root / 'spatial_species_metrics.csv', {'augmented_outside_train_feature_range_photos'}, update_species)
    fractions = {key:count / n for key, (count, n) in summary_counts.items()}
    def update_summary(row):
        key = (float(row['buffer_km']), row['comparison'])
        return {'outside_training_range_fraction_photos':repr(fractions[key])}
    rewrite_columns(root / 'spatial_summary.csv', {'outside_training_range_fraction_photos'}, update_summary)
    profile_path = root / 'profile.json'
    profile = json.loads(profile_path.read_text())
    for row in profile['spatial_summary']:
        row['outside_training_range_fraction_photos'] = fractions[(row['buffer_km'], row['comparison'])]
    profile['extrapolation_range_boundary_tolerance'] = '64 * float64 epsilon * max(1, abs(training feature minimum), abs(training feature maximum)); direct physical features, not inverse standardization'
    profile['raw_firth_variance_diagnostic'] = 'max_test_linear_predictor_variance uses raw Firth covariance; it is not recalibrated prediction variance or an interval'
    atomic_blob(profile_path, (json.dumps(profile, indent=2, sort_keys=True, allow_nan=False) + '\n').encode())
    modified = {'spatial_oof_predictions.npz', 'spatial_fit_diagnostics.csv.gz',
                'spatial_species_metrics.csv', 'spatial_summary.csv', 'profile.json'}
    for name, expected in old_artifact_hashes.items():
        if name not in modified:
            assert digest((root / name).read_bytes()) == expected, name
    correction = {
        'completed_utc':datetime.now(timezone.utc).isoformat(),
        'scope':'Only floating-point outside-training-range flags and linked descriptive counts/fractions; no fitting, selection, scores, or bootstrap changes.',
        'reason':'Inverse standardization could classify an exact training-boundary value as extrapolation at machine precision.',
        'rule':profile['extrapolation_range_boundary_tolerance'],
        'scientific_execution_code_sha256':receipt['code_sha256'],
        'corrected_analysis_code_sha256':digest(Path('scripts/analysis/run_fcp_environment_transport_20261019.py').read_bytes()),
        'correction_script_sha256':digest(Path(__file__).read_bytes()),
        'regression_test_sha256':digest(Path('tests/test_fcp_environment_transport_20261019.py').read_bytes()),
        'counts':[{'buffer_km':buffer, 'model':model, **value} for (buffer, model), value in sorted(counts.items())],
        'paired_comparison_fractions':[{'buffer_km':buffer, 'comparison':comparison,
                                       'outside_training_range_fraction_photos':fraction}
                                      for (buffer, comparison), fraction in sorted(fractions.items())],
        'original_spatial_npz_component_sha256':old_component_hashes,
        'corrected_spatial_npz_component_sha256':component_hashes,
        'all_nonrange_npz_components_bytewise_unchanged':True,
        'all_nondiagnostic_csv_columns_stringwise_unchanged':True,
        'all_bootstrap_artifacts_bytewise_unchanged':True,
        'original_modified_artifact_sha256':{name:old_artifact_hashes[name] for name in sorted(modified)},
    }
    receipt['range_diagnostic_correction'] = correction
    for name in (*modified, 'script_at_execution.py', Path(__file__).name):
        receipt['output_sha256'][name] = digest((root / name).read_bytes())
    atomic_blob(receipt_path, (json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + '\n').encode())
    print(json.dumps({'counts':correction['counts'], 'fractions':correction['paired_comparison_fractions'],
                      'nonrange_components_unchanged':True}, indent=2))


if __name__ == '__main__':
    main()

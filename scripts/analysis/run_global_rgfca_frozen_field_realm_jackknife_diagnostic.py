#!/usr/bin/env python3
"""Audit geographic leverage in frozen G1 fields without refitting or new nulls.

This is an output-domain deletion diagnostic, NOT a full leave-one-realm-out
refit and NOT a way to change the frozen G1/G2 decisions.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

REALMS = ('Afrotropic','Antarctic','Australasia','Indomalaya','Nearctic',
          'Neotropic','Oceania','Palearctic','UNASSIGNED')


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def weighted_variance(values: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Last-axis weighted population variance; never silently drop data."""
    x = np.asarray(values, dtype=np.float64)
    w = np.asarray(weights, dtype=np.float64)
    if w.ndim != 1 or x.ndim not in (1, 2) or x.shape[-1] != len(w):
        raise ValueError('incompatible field and opportunity shapes')
    if len(w) < 3 or not np.isfinite(x).all() or not np.isfinite(w).all():
        raise ValueError('need at least three completely finite cells')
    if np.any(w <= 0):
        raise ValueError('opportunity must be strictly positive')
    mu = np.sum(x * w, axis=-1, keepdims=True) / w.sum()
    return np.sum((x - mu)**2 * w, axis=-1) / w.sum()


def deletion_rows(observed, null_fields, opportunity, evaluable, realms):
    """All groups are fixed in advance; an empty group is not a passed check."""
    y = np.asarray(observed, dtype=float)
    null = np.asarray(null_fields, dtype=float)
    w = np.asarray(opportunity, dtype=float)
    mask = np.asarray(evaluable, dtype=bool)
    labels = np.asarray(realms, dtype=str)
    if y.ndim != 1 or null.shape != (999, len(y)):
        raise ValueError('expected all 999 frozen aggregate null fields')
    if w.shape != y.shape or mask.shape != y.shape or labels.shape != y.shape:
        raise ValueError('unaligned cells')
    if not set(labels).issubset(REALMS):
        raise ValueError('unregistered realm labels')
    base_obs = float(weighted_variance(y[mask], w[mask]))
    base_null = weighted_variance(null[:, mask], w[mask])
    base_excess = base_obs - float(base_null.mean())
    rows = []
    for label in ('FULL_DOMAIN',) + REALMS:
        removed = mask & (labels == label) if label != 'FULL_DOMAIN' else np.zeros_like(mask)
        keep = mask & ~removed
        row = {'deleted_group': label, 'removed_cells': int(removed.sum()),
               'retained_cells': int(keep.sum()),
               'removed_opportunity_fraction': float(w[removed].sum()/w[mask].sum())}
        if label != 'FULL_DOMAIN' and not removed.any():
            row['status'] = 'not_evaluable_no_evaluable_cells'
        elif keep.sum() < 3:
            row['status'] = 'not_evaluable_fewer_than_three_retained_cells'
        else:
            obs = float(weighted_variance(y[keep], w[keep]))
            nstat = weighted_variance(null[:, keep], w[keep])
            nmean = float(nstat.mean())
            excess = obs - nmean
            row.update(status='complete_descriptive_only', observed_concentration=obs,
                       mean_null_concentration=nmean, concentration_excess=excess,
                       observed_to_null_mean_ratio=obs/nmean if nmean > 0 else None,
                       excess_change_from_full_domain=excess-base_excess,
                       positive_excess=bool(excess > 0))
        rows.append(row)
    baseline = {'reconstructed_observed_concentration':base_obs,
                'reconstructed_null_mean_concentration':float(base_null.mean()),
                'reconstructed_observed_to_null_mean_ratio':base_obs/float(base_null.mean()),
                'reconstructed_existing_primary_p_upper':float((1 + np.count_nonzero(base_null >= base_obs))/1000),
                'original_evaluable_cells':int(mask.sum())}
    return rows, baseline, base_null


def main() -> int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--primary-zip',type=Path,required=True)
    ap.add_argument('--resolve-zip',type=Path,required=True)
    ap.add_argument('--primary-dir',type=Path,required=True)
    ap.add_argument('--resolve-dir',type=Path,required=True)
    ap.add_argument('--contract',type=Path,required=True)
    ap.add_argument('--output-dir',type=Path,required=True)
    args=ap.parse_args()
    contract=json.loads(args.contract.read_text())
    frozen=contract['inputs']
    for p,key in ((args.primary_zip,'g1_artifact_zip_sha256'),(args.resolve_zip,'resolve_artifact_zip_sha256')):
        if sha256(p)!=frozen[key]:
            raise RuntimeError(f'frozen artifact mismatch: {p}')
    raster=args.resolve_dir/'resolve_ecoregion2017_zones_v1.tif'
    gpkg=args.resolve_dir/'Ecoregions2017.gpkg'
    for p,key in ((raster,'resolve_raster_sha256'),(gpkg,'resolve_gpkg_sha256')):
        if sha256(p)!=frozen[key]:raise RuntimeError(f'frozen source mismatch: {p}')
    import geopandas as gpd
    import rasterio
    source=gpd.read_file(gpkg)
    if source.crs is None or source.crs.to_epsg()!=4326:
        raise RuntimeError('unexpected GPKG CRS')
    arrays=np.load(args.primary_dir/'global_rgfca_g1_primary_v1.npz',allow_pickle=False)
    grid_lat=arrays['grid_latitude'];grid_lon=arrays['grid_longitude']
    with rasterio.open(raster) as ds:
        if ds.crs.to_epsg()!=4326:raise RuntimeError('unexpected raster CRS')
        zone=np.array([int(x[0]) for x in ds.sample(zip(grid_lon,grid_lat))])
        nodata=int(ds.nodata)
    if not set(zone).issubset(set(range(1,len(source)+1))|{nodata}):
        raise RuntimeError('invalid frozen raster row IDs')
    labels=np.full(len(zone),'UNASSIGNED',dtype=object)
    valid=zone!=nodata
    labels[valid]=source.iloc[zone[valid]-1]['REALM'].astype(str).to_numpy()
    rows,baseline,reconstructed_null=deletion_rows(
        arrays['aggregate_field'], arrays['null_aggregate_fields'],
        arrays['aggregate_opportunity'], arrays['evaluable'], labels)
    original=json.loads((args.primary_dir/'global_rgfca_g1_result_v1.json').read_text())
    if original['global_recurrence_supported'] is not False:
        raise RuntimeError('unexpected original support state')
    if not np.isclose(baseline['reconstructed_observed_concentration'],original['observed_concentration'],rtol=1e-12,atol=1e-15):
        raise RuntimeError('original observed statistic not reproduced')
    if not np.allclose(reconstructed_null,arrays['null_concentration'],rtol=1e-12,atol=1e-15):
        raise RuntimeError('original null statistics not reproduced')
    if baseline['reconstructed_existing_primary_p_upper']!=original['p_upper'] or original['p_upper']!=0.07:
        raise RuntimeError('original p-value not reproduced')
    out=args.output_dir;out.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(rows).to_csv(out/'global_rgfca_frozen_field_realm_jackknife_diagnostic_v1.csv',index=False)
    pd.DataFrame({'cell_index':np.arange(len(zone)), 'latitude':grid_lat,'longitude':grid_lon,
                  'frozen_resolve_row_id':zone,'realm':labels,
                  'original_evaluable':arrays['evaluable'],
                  'original_opportunity':arrays['aggregate_opportunity']}).to_csv(
                      out/'global_rgfca_frozen_field_realm_assignments_v1.csv',index=False)
    result={'protocol':contract['protocol'],
            'status':'complete_descriptive_output_domain_realm_deletion_no_refit',
            'contract_sha256':sha256(args.contract),
            'source_lineage':frozen,
            'baseline_verification':baseline,'deletions':rows,
            'interpretation_guard':contract['guard'],
            'limitations':contract['limitations']}
    (out/'global_rgfca_frozen_field_realm_jackknife_diagnostic_result_v1.json').write_text(
        json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'baseline':baseline,'deletions':rows},indent=2))
    return 0

if __name__=='__main__':
    raise SystemExit(main())

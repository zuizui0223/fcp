"""Source-photo matched environment from pinned rasters, never nearest mois photos.

The manifest identifies the same 36 covariate *names* as mois. Units, depth,
transforms and raster checksums must be stated explicitly for each supplied
field: equal names alone do not certify identical soil measurements. Missing
values are retained. This script measures covariates only, not ecological effects.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio.warp import transform as warp_coordinates

MOIS_NONSOIL = (
    "wc_elevation_m", *(f"wc_bio{i}" for i in range(1, 20)),
    "wc_srad_annual_mean", "wc_srad_monthly_cv",
    "wc_wind_annual_mean", "wc_wind_monthly_cv",
    "wc_vapr_annual_mean", "wc_vapr_monthly_cv",
)
MOIS_SOIL = (
    "soil_pH", "soil_SOC", "soil_N", "soil_clay", "soil_available_water_proxy",
    "soil_cec_0_30cm_source_raw", "soil_sand_0_30cm_source_raw",
    "soil_silt_0_30cm_source_raw", "soil_bdod_0_30cm_source_raw",
    "soil_cfvo_0_30cm_source_raw",
)
MOIS_36 = MOIS_NONSOIL + MOIS_SOIL
assert len(MOIS_NONSOIL) == 26 and len(MOIS_SOIL) == 10 and len(MOIS_36) == 36


def file_hash(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def read_manifest(path: Path, *, require_all36: bool = False) -> dict:
    obj = json.loads(path.read_text(encoding='utf-8'))
    if obj.get('schema') != 'fcp_mois36_direct_photo_raster_manifest_v1':
        raise ValueError('Unknown provenance manifest schema')
    specs = obj.get('features')
    if not isinstance(specs, dict) or not specs or set(specs) - set(MOIS_36):
        raise ValueError('Unknown or missing source environment feature names')
    if require_all36 and set(specs) != set(MOIS_36):
        raise ValueError('All 36 mois-aligned raster features are required')
    base = path.parent.resolve()
    for feature, rec in specs.items():
        if not isinstance(rec, dict) or not rec.get('unit') or not rec.get('source'):
            raise ValueError(f'{feature}: unit and original source are mandatory')
        method = rec.get('reducer')
        raw = rec.get('rasters')
        if method not in ('identity', 'monthly_mean', 'monthly_cv') or not isinstance(raw, list):
            raise ValueError(f'{feature}: invalid raster reducer')
        expected = 1 if method == 'identity' else 12
        if len(raw) != expected:
            raise ValueError(f'{feature}: expected {expected} rasters, got {len(raw)}')
        if feature.startswith('soil_') and method != 'identity':
            raise ValueError('Soil properties require explicit processed rasters')
        for raster in raw:
            if not isinstance(raster, dict) or not raster.get('path') or not raster.get('sha256'):
                raise ValueError(f'{feature}: source raster path/hash required')
            f = (base / raster['path']).resolve()
            if not f.is_file() or file_hash(f) != raster['sha256']:
                raise ValueError(f'{feature}: raster SHA256 mismatch or unavailable')
            if not f.suffix.lower() in ('.tif', '.tiff'):
                raise ValueError('Only georeferenced TIFF raster sources supported')
            raster['_verified_path'] = str(f)
    return obj


def values_at_points(spec: dict, lon: np.ndarray, lat: np.ndarray) -> np.ndarray:
    monthly=[]
    for record in spec['rasters']:
        with rasterio.open(record['_verified_path']) as ds:
            if ds.count != 1 or ds.crs is None:
                raise ValueError('Raster must have one band and a declared CRS')
            if ds.crs.to_epsg() == 4326:
                xs, ys = lon.tolist(), lat.tolist()
            else:
                xs, ys = warp_coordinates('EPSG:4326', ds.crs, lon.tolist(), lat.tolist())
            from rasterio.transform import rowcol
            rows, cols = rowcol(ds.transform, xs, ys)
            inside = (np.asarray(rows) >= 0) & (np.asarray(rows) < ds.height) & (np.asarray(cols) >= 0) & (np.asarray(cols) < ds.width)
            sampled = np.asarray([
                np.nan if np.ma.is_masked(v[0]) else float(v[0])
                for v in ds.sample(zip(xs,ys), masked=True)
            ], dtype=float)
            sampled[~inside] = np.nan
            # Rasterio's masked sampling does not mask every raster-specific
            # ocean sentinel; explicit nodata and raster bounds are mandatory.
            if ds.nodata is not None and np.isfinite(ds.nodata):
                sampled[np.isclose(sampled, ds.nodata, rtol=0, atol=1e-10)] = np.nan
            sampled[~np.isfinite(sampled)] = np.nan
            monthly.append(sampled)
    arr=np.vstack(monthly)
    if spec['reducer'] == 'identity':
        return arr[0]
    valid = np.isfinite(arr).all(axis=0)
    result = np.full(arr.shape[1], np.nan, dtype=float)
    if np.any(valid):
        x=arr[:,valid]
        if spec['reducer'] == 'monthly_mean':
            result[valid] = x.mean(axis=0)
        else:
            means=x.mean(axis=0)
            good=np.abs(means)>1e-12
            z=np.full(len(means),np.nan)
            z[good]=x[:,good].std(axis=0,ddof=0)/np.abs(means[good])
            result[valid]=z
    return result


def sample_rows(data: pd.DataFrame, manifest: dict) -> pd.DataFrame:
    required={'inat_taxon_id','photo_id','latitude','longitude'}
    if required - set(data):
        raise ValueError(f'Missing source-photo keys/locations {required-set(data)}')
    if data['photo_id'].isna().any() or data['photo_id'].astype(str).duplicated().any():
        raise ValueError('Every photo_id must be present and unique in this partition')
    out = data[['inat_taxon_id','photo_id','latitude','longitude'] +
               (['positional_accuracy_m'] if 'positional_accuracy_m' in data else [])].copy()
    out['photo_id']=out.photo_id.astype(str)
    lat=pd.to_numeric(out.latitude,errors='coerce').to_numpy(float)
    lon=pd.to_numeric(out.longitude,errors='coerce').to_numpy(float)
    valid=np.isfinite(lat)&np.isfinite(lon)&(np.abs(lat)<=90)&(np.abs(lon)<=180)
    out['georeferenced']=valid
    if 'positional_accuracy_m' in out:
        acc=pd.to_numeric(out.positional_accuracy_m,errors='coerce')
        out['accuracy_le1000m']=valid & acc.between(0,1000).to_numpy(bool)
        out['accuracy_le5000m']=valid & acc.between(0,5000).to_numpy(bool)
    else:
        out['accuracy_le1000m']=False
        out['accuracy_le5000m']=False
    for feat,spec in manifest['features'].items():
        value=np.full(len(data),np.nan,dtype=float)
        if valid.any():
            value[valid]=values_at_points(spec,lon[valid],lat[valid])
        out[feat]=value
    out['location_source']='direct_source_photo_coordinate'
    return out


def stream_all_photos(source_paths: list[Path], manifest_path: Path, out: Path,
                      *, require_all36: bool = False, chunksize: int = 5000) -> dict:
    manifest=read_manifest(manifest_path,require_all36=require_all36)
    out.parent.mkdir(parents=True,exist_ok=True)
    db=sqlite3.connect(out.with_suffix('.seen.sqlite'))
    db.execute('CREATE TABLE IF NOT EXISTS ids(photo_id TEXT PRIMARY KEY)')
    n=0
    n_geo=0
    complete={k:0 for k in manifest['features']}
    if out.exists():
        raise ValueError('Refuse to overwrite existing direct-photo environmental table')
    first=True
    try:
        for source in source_paths:
            for chunk in pd.read_csv(source,chunksize=chunksize,dtype={'photo_id':str},low_memory=False):
                if chunk.photo_id.isna().any():
                    raise ValueError('Missing photo IDs')
                try:
                    db.executemany('INSERT INTO ids(photo_id) VALUES (?)',((x,) for x in chunk.photo_id))
                except sqlite3.IntegrityError as ex:
                    raise ValueError('Repeated photo ID: new-source overlap must be resolved before environmental inference') from ex
                sampled=sample_rows(chunk,manifest)
                sampled.to_csv(out,mode='w' if first else 'a',index=False,header=first)
                first=False
                n+=len(sampled)
                n_geo+=int(sampled.georeferenced.sum())
                for k in complete:
                    complete[k]+=int(sampled[k].notna().sum())
        if first:  # no photographs
            pd.DataFrame(columns=['inat_taxon_id','photo_id','latitude','longitude',
                    'georeferenced','accuracy_le1000m','accuracy_le5000m',
                    *manifest['features'],'location_source']).to_csv(out,index=False)
    finally:
        db.close()
        db_path=out.with_suffix('.seen.sqlite')
        if db_path.exists():db_path.unlink()
    receipt={'schema':'fcp_direct_photo_environment_receipt_v1',
            'n_photo_ids':n,'n_georeferenced':n_geo,
            'per_feature_available_photo_count':complete,
            'raster_manifest_sha256':file_hash(manifest_path),
            'source_photo_csv_sha256':{str(p):file_hash(p) for p in source_paths},
            'features_requested':list(manifest['features']),
            'is_full_mois36':set(manifest['features'])==set(MOIS_36),
            'all_raster_values_sampled_at_source_photo_coordinates':True,
            'proxy_nearest_mois_photo_values_used':False,
            'missing_not_imputed':True,'landscape_or_genetic_selection_identified':False}
    out.with_suffix('.receipt.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    return receipt


def main(argv: list[str] | None=None) -> int:
    p=argparse.ArgumentParser()
    p.add_argument('--measured',action='append',required=True,type=Path)
    p.add_argument('--manifest',required=True,type=Path)
    p.add_argument('--out',required=True,type=Path)
    p.add_argument('--require-all36',action='store_true')
    a=p.parse_args(argv)
    report=stream_all_photos(a.measured,a.manifest,a.out,require_all36=a.require_all36)
    print(json.dumps({'photos':report['n_photo_ids'],
                      'features':len(report['features_requested']),
                      'all36':report['is_full_mois36']}))
    return 0


if __name__=='__main__':
    raise SystemExit(main())

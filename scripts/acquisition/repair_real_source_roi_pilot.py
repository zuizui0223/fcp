#!/usr/bin/env python3
"""One-shot repair of EXACTLY two case-sensitive extension failures in a REAL FCP pilot.

Original 94 real-image measurements from the immutable October 15 artifact are
preserved byte/provenance checked. Only image IDs 840402 and 3201197, whose
original source URL uses `large.JPG`, are newly downloaded and measured with the
same frozen location-blind historical ROI-v4 implementation. This is STILL NOT
an all-42,111-species result, not a sample extension, and not FCP inference.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from fcp_pipeline.fcp_allphotos import sha256_file
from fcp_pipeline.fcp_fresh_study import BIO
from scripts.acquisition.run_fcp_fresh_image_partition import verify_frozen_method
from scripts.acquisition.run_real_source_roi_pilot import (
    SALT, SOURCE_ARTIFACT, SOURCE_RUN, selected_metadata, original_photo_bytes,
)

FIRST_REAL_CSV_SHA = '772dda9c51cf24fa4fdabf6408d16b2ba96fc965b532af82533c150cf0144c09'
FIX_IDS = frozenset({'840402', '3201197'})
EXPECTED_TOTAL = 96


def frozen_previous(prior_csv: Path, prior_receipt: Path, source_meta: Path):
    if sha256_file(prior_csv) != FIRST_REAL_CSV_SHA:
        raise ValueError('Original 96-real-photo measurement artifact has changed')
    prior = json.loads(prior_receipt.read_text(encoding='utf-8'))
    if (prior.get('real_fresh_photo_terminal_csv_sha256') != FIRST_REAL_CSV_SHA
            or prior.get('n_real_source_original_photo_ids') != EXPECTED_TOTAL
            or prior.get('n_real_original_images_processed_by_frozen_roi_v4') != 94
            or prior.get('source', {}).get('source_artifact_id') != SOURCE_ARTIFACT
            or prior.get('source', {}).get('source_run_id') != SOURCE_RUN
            or prior.get('eligible_for_main_fcp_W_null_or_climate_inference') is not False):
        raise ValueError('Prior real-photo pilot provenance or quarantine was altered')
    selected, source = selected_metadata(source_meta, 24)
    if prior.get('source', {}).get('selected_photo_id_pairs_sha256') != source['selected_photo_id_pairs_sha256']:
        raise ValueError('Exact outcome-blind source selection changed since original image pilot')
    data = pd.read_csv(prior_csv, dtype=str, keep_default_na=False)
    if len(data) != EXPECTED_TOTAL or data.photo_id.duplicated().any() or data.measurement_id.duplicated().any():
        raise ValueError('Prior original-photo terminal source IDs are missing or duplicated')
    original_ids = {str(r['photo_id']): r for r in selected}
    if set(data.photo_id) != set(original_ids):
        raise ValueError('Previous pilot photo IDs differ from the frozen provider manifest')
    if set(data.loc[data.measurement_status.eq('image_acquisition_failed'), 'photo_id']) != FIX_IDS:
        raise ValueError('Original source rejected a different set of photo IDs')
    for row in data.itertuples(index=False):
        src = original_ids[str(row.photo_id)]
        expected = 'FCPR-PILOT-' + hashlib.sha256(f'{SALT}:{row.photo_id}'.encode()).hexdigest().upper()[:24]
        if (str(row.measurement_id) != expected
                or str(row.source_url) != str(src['photo_url_large'])
                or str(row.species) != str(src['species'])
                or str(row.observation_id) != str(src['observation_id'])
                or str(row.photo_license) != str(src['photo_license'])):
            raise ValueError('Previous pilot row does not match the original provider photo ID')
    for id_ in FIX_IDS:
        src = original_ids[id_]
        if not str(src['photo_url_large']).endswith('/large.JPG'):
            raise ValueError('Only the two documented uppercase-extension photos may be repaired')
    return data, prior, source, original_ids


def repair(prior_csv: Path, prior_receipt: Path, metadata_root: Path, output_root: Path,
           source: Path, detector: Path, sam: Path, *, delay_seconds: float = 1.1):
    if delay_seconds < 1.0:
        raise ValueError('Provider image requests require >=1 second intervals')
    if output_root.exists() and any(output_root.iterdir()):
        raise ValueError('Repair output directory must be empty')
    data, prior, selected_info, id_map = frozen_previous(prior_csv, prior_receipt, metadata_root)
    method = verify_frozen_method(source, detector, sam)
    if method != prior['original_roi_v4_method_hashes']:
        raise ValueError('Frozen ROI-v4 runtime differs from first 94 real images')
    output_root.mkdir(parents=True, exist_ok=True)
    new_hashes = {}
    failures = {}
    actual_worker = []
    with tempfile.TemporaryDirectory(prefix='fcp_two_original_photo_extension_repair_') as tmp:
        scratch = Path(tmp)
        images = scratch/'images'
        images.mkdir()
        for i, pid in enumerate(sorted(FIX_IDS, key=int)):
            if i:
                time.sleep(delay_seconds)
            src = id_map[pid]
            mid = data.loc[data.photo_id.eq(pid), 'measurement_id'].iloc[0]
            try:
                blob = original_photo_bytes(src['photo_url_large'], pid)
                name = mid+'.jpg'
                (images/name).write_bytes(blob)
                new_hashes[mid] = hashlib.sha256(blob).hexdigest()
                actual_worker.append({'measurement_id':mid, 'image_filename':name,
                                      'photo_license':src['photo_license']})
                print(f'UPPERCASE_REAL_PHOTO_RECOVERED {pid} size={len(blob)}', flush=True)
            except Exception as e:
                failures[mid] = f'{type(e).__name__}:{str(e)[:180]}'
                print(f'UPPERCASE_REAL_PHOTO_STILL_UNAVAILABLE {pid} {type(e).__name__}', flush=True)
        measured = {}
        if actual_worker:
            worker = scratch/'blind_worker.csv'
            with worker.open('w', encoding='utf-8', newline='') as f:
                w=csv.DictWriter(f, fieldnames=('measurement_id','image_filename','photo_license'))
                w.writeheader()
                w.writerows(actual_worker)
            terminal = scratch/'original_frozen_roi_v4_two_photos.csv'
            subprocess.run([sys.executable, str(source), '--worker-manifest', str(worker),
                            '--images-dir', str(images), '--detector-weight', str(detector),
                            '--efficient-sam-dir', str(sam), '--output-csv', str(terminal)],
                           check=True, timeout=1000)
            df = pd.read_csv(terminal, dtype=str, keep_default_na=False)
            if len(df) != len(actual_worker) or df.measurement_id.duplicated().any() or set(df.measurement_id) != {x['measurement_id'] for x in actual_worker}:
                raise ValueError('Original exact ROI-v4 implementation changed its terminal IDs')
            measured = df.set_index('measurement_id').to_dict(orient='index')
        for pid in FIX_IDS:
            mask = data.photo_id.eq(pid)
            mid = data.loc[mask,'measurement_id'].iloc[0]
            if mid in measured:
                val = measured[mid]
                if (val['measurement_status']=='classified_four_state_morph') != (val['morph'] in ('white','yellow_orange','red_pink','blue_purple')):
                    raise ValueError('Original ROI colour decision inconsistent')
                for col, value in val.items():
                    if col in data:
                        data.loc[mask, col] = value
                data.loc[mask,'original_image_bytes_sha256'] = new_hashes[mid]
                data.loc[mask,'global_classifiable'] = str(val['measurement_status']=='classified_four_state_morph')
            else:
                data.loc[mask,'failure_reasons'] = failures[mid]
        counts = data.measurement_status.value_counts(dropna=False).to_dict()
        classes = int(data.measurement_status.eq('classified_four_state_morph').sum())
        failed = int(data.measurement_status.eq('image_acquisition_failed').sum())
        n_roi = int(EXPECTED_TOTAL - classes - failed)
        out = output_root/'real_original_roi_v4_new_photo_pilot_repaired.csv'
        data.to_csv(out, index=False, lineterminator='\n')
        report = {
            'schema':'fcp_real_original_roi_v4_case_sensitive_extension_repair_v1',
            'status':'REPAIRED_REAL_PHOTO_TECHNICAL_PILOT_NOT_42111_CENSUS',
            'utc':datetime.now(timezone.utc).isoformat(),
            'prior_real_photo_pilot_run':38671341143,
            'prior_real_photo_pilot_artifact':11944359701,
            'prior_real_photo_pilot_csv_sha256': FIRST_REAL_CSV_SHA,
            'source_artifact':SOURCE_ARTIFACT,
            'identical_outcome_blind_selected_photo_id_pairs_sha256':selected_info['selected_photo_id_pairs_sha256'],
            'source_photo_id_count':EXPECTED_TOTAL,
            'id_repair_candidates':sorted(FIX_IDS, key=int),
            'n_original_images_retried':len(FIX_IDS),
            'n_original_images_recovered':len(actual_worker),
            'n_original_image_acquisition_failed':failed,
            'n_actual_original_image_measurements_preserved_from_prior':94,
            'n_actual_original_image_measurements_new':len(actual_worker),
            'n_actual_original_image_measurements_total':94+len(actual_worker),
            'n_classified_four_state_photo_rows': classes,
            'n_unclassifiable_roi_palette_photo_rows':n_roi,
            'morph_counts': data.loc[data.measurement_status.eq('classified_four_state_morph'),'morph'].value_counts().to_dict(),
            'status_counts':counts,
            'per_species': {name:{'source_photo_ids':len(group),
                                  'n_classified':int(group.measurement_status.eq('classified_four_state_morph').sum()),
                                  'n_acquisition_failed':int(group.measurement_status.eq('image_acquisition_failed').sum())}
                            for name,group in data.groupby('species')},
            'terminal_csv_sha256':sha256_file(out),
            'original_method':method,
            'previous_94_new_pixel_measurements_reused_verbatim':True,
            'older_149900_photo_measurements_read':False,
            'fully_original_source_42111_taxa_collected':False,
            'eligible_for_main_fcp_W_null_or_climate_inference':False,
            'original_images_persisted':False
        }
        (output_root/'real_original_roi_v4_new_photo_pilot_repaired_receipt.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
        print('REAL_TWO_UPPERCASE_PHOTO_REPAIR_RESULT '+json.dumps({
            'successes':len(actual_worker),'new_classified_total':classes,'failures':failed,
            'status_counts':counts,'sha256':report['terminal_csv_sha256']}, sort_keys=True),flush=True)
        return report


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--original-pilot-csv',type=Path,required=True)
    p.add_argument('--original-pilot-receipt',type=Path,required=True)
    p.add_argument('--source-metadata-root',type=Path,required=True)
    p.add_argument('--original-measure-script',type=Path,required=True)
    p.add_argument('--original-detector',type=Path,required=True)
    p.add_argument('--original-efficient-sam',type=Path,required=True)
    p.add_argument('--output-root',type=Path,required=True)
    p.add_argument('--delay-seconds',type=float,default=1.1)
    a=p.parse_args(argv)
    return repair(a.original_pilot_csv,a.original_pilot_receipt,a.source_metadata_root,a.output_root,
                  a.original_measure_script,a.original_detector,a.original_efficient_sam,
                  delay_seconds=a.delay_seconds)

if __name__=='__main__':
    main()

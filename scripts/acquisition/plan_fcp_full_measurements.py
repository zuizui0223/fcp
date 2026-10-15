#!/usr/bin/env python3
"""Enumerate EVERY original-photo ROI-v4 work unit from one COMPLETE FCP freeze.

No species/sample/photo cap; batch_size caps TEMPORARY work-unit memory only.
All 128 source-blind partitions are inspected without opening source images,
old colour labels, or environmental values.  Existing successful terminal
batches are validated before being marked complete for resumable execution.
The output is a task ledger, NOT flower-colour or ecological evidence.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import pandas as pd

from fcp_pipeline.fcp_allphotos import sha256_file
from fcp_pipeline.fcp_fresh_study import WORKER, validate_method_attestation


def plan_jobs(firewall: Path, completed: Path, out: Path,
              *, batch_size: int = 64, _test_only_small_source: bool = False) -> dict:
    if batch_size < 1 or batch_size > 1024:
        raise ValueError('Compute work unit must be between 1 and 1024, NEVER an inclusion cap')
    if out.exists() and any(out.iterdir()):
        raise ValueError('Fresh distributed processing queue output must be empty')
    frozen_path=firewall/'fresh_study_freeze.json'
    doc=json.loads(frozen_path.read_text(encoding='utf-8'))
    if (doc.get('status') != 'FRESH_BLIND_METADATA_FROZEN'
        or doc.get('new_measurements_completed') != 0
        or doc.get('historical_measured_CSV_read') is not False):
        raise ValueError('Main source census must be fully frozen before planning any original photo pixels')
    if _test_only_small_source:
        if (not 1 <= int(doc.get('source_species',0)) <= 20 or
                not 1 <= int(doc.get('worker_packets',0)) <= 8):
            raise ValueError('Synthetic short-source bypass is only for tiny regression fixtures')
    elif doc.get('source_species') != 42111 or doc.get('worker_packets') != 128:
        raise ValueError('Main 42111 species / 128-partition source must be fully frozen')
    receipt=firewall/'firewall_receipt.json'
    if doc.get('blind_firewall_sha256') != sha256_file(receipt):
        raise ValueError('Frozen source firewall bytes changed')
    output_counts={'complete':0,'pending':0}
    n_photos=0;n_jobs=0
    out.mkdir(parents=True,exist_ok=True)
    existing_files=set()
    if completed.exists():
        for e in completed.iterdir():
            if e.name.startswith('terminal_') and e.suffix in ('.json','.csv'):
                existing_files.add(e.name)
            elif e.name == 'frozen_method_attestation.json':
                validate_method_attestation(e)
            else:
                raise ValueError('Unexpected output in sealed terminal batch directory')
    queue=out/'all_photo_ROI_v4_work_units.ndjson'
    with queue.open('w',encoding='utf-8') as h:
        for part in range(int(doc['worker_packets'])):
            worker_name=f'worker_partition_{part:03}.csv'
            key_name=f'acquisition_partition_{part:03}.csv'
            for name in (worker_name,key_name):
                source=firewall/name
                if sha256_file(source)!=doc['frozen_file_sha256'][name]:
                    raise ValueError(f'Worker/source URI partition changed after full census: {name}')
            worker=firewall/worker_name
            seen=set()
            buffer=[]
            def flush(sub):
                nonlocal n_jobs, n_photos
                if not sub:
                    return
                batch=n_jobs_by_part[part]
                n_jobs_by_part[part]+=1
                current=f'terminal_p{part:03}_b{batch:07}'
                csv_file=completed/(current+'.csv')
                rec_file=completed/(current+'.json')
                worker_hash=hashlib.sha256(pd.DataFrame(sub,columns=list(WORKER)).to_csv(index=False).encode()).hexdigest()
                exists=(csv_file.exists(),rec_file.exists())
                if exists[0]!=exists[1]:
                    raise ValueError('A terminal CSV and its source-bound receipt must be present together')
                state='pending'
                if all(exists):
                    rec=json.loads(rec_file.read_text(encoding='utf-8'))
                    if (rec.get('status')!='COMPLETE_WORK_UNIT'
                        or rec.get('partition')!=part or rec.get('batch')!=batch
                        or rec.get('compute_batch_size_not_photo_cap')!=batch_size
                        or rec.get('selected_source_worker_sha256')!=worker_hash
                        or rec.get('n_terminal_photo_ids')!=len(sub)
                        or rec.get('n_source_photo_ids')!=len(sub)
                        or rec.get('terminal_csv_sha256')!=sha256_file(csv_file)
                        or rec.get('reused_historical_measurement_csv') is not False
                        or rec.get('image_bytes_persisted') is not False
                        or rec.get('full_study_complete') is not False):
                        raise ValueError('Cannot mark a changed or synthetic terminal batch as completed')
                    method=rec.get('method')
                    if not isinstance(method,dict):
                        raise ValueError('Missing original exact ROI-v4 method attestation')
                    method_path=completed/'frozen_method_attestation.json'
                    if not method_path.is_file() or method!=validate_method_attestation(method_path):
                        raise ValueError('ROI source or model method changed between original photo batches')
                    table=pd.read_csv(csv_file,dtype={'measurement_id':str})
                    if (len(table)!=len(sub) or table.measurement_id.duplicated().any()
                        or set(table.measurement_id)!= {r['measurement_id'] for r in sub}):
                        raise ValueError('Sealed terminal result contains wrong original photos')
                    state='complete'
                    existing_files.discard(csv_file.name)
                    existing_files.discard(rec_file.name)
                output_counts[state]+=1
                n_jobs+=1
                n_photos+=len(sub)
                h.write(json.dumps({'schema':'fcp_all_photos_immutable_work_unit_v1',
                                    'partition':part,'batch':batch,
                                    'compute_batch_size_not_sampling_cap':batch_size,
                                    'n_source_photo_ids':len(sub),
                                    'source_worker_sha256':worker_hash,
                                    'status':state},sort_keys=True)+'\n')
            if part==0:
                n_jobs_by_part={i:0 for i in range(int(doc['worker_packets']))}
            with worker.open('r',encoding='utf-8',newline='') as source:
                reader=csv.DictReader(source)
                if tuple(reader.fieldnames or ())!=WORKER:
                    raise ValueError('Worker packet is not original exactly three-column blind photo input')
                for row in reader:
                    mid=row['measurement_id']
                    if not mid or mid in seen:
                        raise ValueError('Duplicate or missing frozen original photo measurement ID')
                    seen.add(mid)
                    buffer.append(row)
                    if len(buffer)==batch_size:
                        flush(buffer)
                        buffer=[]
                flush(buffer)
    if existing_files:
        raise ValueError(f'Unknown sealed real ROI-v4 terminal files: {sorted(existing_files)[:6]}')
    if n_photos!=doc['all_available_licensed_photo_ids']:
        raise ValueError('All-photo work units do not cover the complete source census')
    result={'schema':'fcp_all_42111_original_photo_work_plan_v1',
            'status':'COMPLETE_UNCAPPED_SOURCE_PHOTO_WORK_PLAN_NOT_ALL_IMAGES_MEASURED',
            'source_species':doc['source_species'],'source_photo_ids':n_photos,
            'work_units':n_jobs,'n_completed_work_units':output_counts['complete'],
            'n_pending_work_units':output_counts['pending'],
            'source_firewall_sha256':doc['blind_firewall_sha256'],
            'work_unit_ledger_sha256':sha256_file(queue),
            'all_photo_ids_included_once':True,
            'target_classifiable_min_per_species':40,
            'this_queue_not_genetic_or_environmental_inference':True,
            'full_biological_analysis_complete':False,
            'synthetic_test_only':_test_only_small_source}
    (out/'work_plan_receipt.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    return result


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--firewall',type=Path,required=True)
    p.add_argument('--completed',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--batch-size',type=int,default=64)
    a=p.parse_args(argv)
    print(json.dumps(plan_jobs(a.firewall,a.completed,a.out,batch_size=a.batch_size),sort_keys=True))
    return 0

if __name__=='__main__':
    raise SystemExit(main())

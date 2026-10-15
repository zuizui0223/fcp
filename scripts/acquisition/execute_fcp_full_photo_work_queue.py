#!/usr/bin/env python3
"""Run EVERY frozen FCP original photograph through ROI-v4, resumably.

Only call after the 42,111-taxon complete source PHOTO census has passed.
The budget below limits a single compute invocation, never taxon/photo inclusion.
A batch is published ONLY after its original model attestation, exact worker
subset fingerprint and terminal result bytes are checked. An interrupted batch
is not a result and cannot satisfy reassembly. No previous colour CSV is read.

On normal completion a sealed terminal set is READY FOR REASSEMBLY, not a
statistical finding. The final complete original-photo CSV, H2, 999 null worlds,
and direct environment joins all remain separate hard gates.
"""
from __future__ import annotations
import argparse
import json
import os
import shutil
import tempfile
from datetime import datetime,timezone
from pathlib import Path

from fcp_pipeline.fcp_allphotos import sha256_file
from fcp_pipeline.fcp_fresh_study import validate_method_attestation
from scripts.acquisition.plan_fcp_full_measurements import plan_jobs
from scripts.acquisition import run_fcp_fresh_image_partition as roi


def _json(path:Path)->dict:
    return json.loads(path.read_text(encoding='utf-8'))


def _utc()->str:
    return datetime.now(timezone.utc).isoformat()


def _write_atomic(path:Path,obj:dict)->None:
    path.parent.mkdir(parents=True,exist_ok=True)
    candidate=path.with_name(path.name+'.tmp')
    if candidate.exists():
        raise ValueError('Previous interrupted progress receipt must be reconciled')
    candidate.write_text(json.dumps(obj,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    os.replace(candidate,path)


def _sealed_status(firewall:Path,completed:Path,batch_size:int,
                   *,partition:int|None=None,
                   _test_only_small_source:bool=False)->tuple[dict,list[dict]]:
    with tempfile.TemporaryDirectory(prefix='fcp_all_photo_queue_check_') as tmp:
        checked=plan_jobs(firewall,completed,Path(tmp)/'plan',
                          batch_size=batch_size,partition=partition,
                          _test_only_small_source=_test_only_small_source)
        path=Path(tmp)/'plan/all_photo_ROI_v4_work_units.ndjson'
        queue=[json.loads(line) for line in path.read_text().splitlines()]
        if (sum(j['n_source_photo_ids'] for j in queue)!=checked['source_photo_ids'] or
                len(queue)!=checked['work_units']):
            raise ValueError('Full frozen original photo work queue inconsistent')
        return checked,queue


def execute_pending(firewall:Path,output:Path,*,source:Path,detector:Path,sam:Path,
                    batch_size:int=64,max_work_units:int=32,partition:int|None=None,
                    _test_only_small_source:bool=False)->dict:
    if max_work_units<0:
        raise ValueError('Work-unit runtime budget must be >=0, with 0 meaning unlimited')
    if partition is not None and partition<0:
        raise ValueError('Selected compute partition must be non-negative')
    # Original method bytes must match BEFORE opening any new source pixels.
    fingerprint=roi.verify_frozen_method(source,detector,sam)
    completed=output/'sealed_terminal_batches'
    stages=output/'unsealed_inflight_batches'
    report_path=output/'original_photo_work_queue_progress.json'
    checked,queue=_sealed_status(firewall,completed,batch_size,partition=partition,
                                 _test_only_small_source=_test_only_small_source)
    if not _test_only_small_source and checked['source_species']!=42111:
        raise ValueError('Only the actual 42,111-taxon complete frozen source is admissible')
    if partition is not None and partition>=_json(firewall/'fresh_study_freeze.json')['worker_packets']:
        raise ValueError('Compute partition lies outside frozen worker count')
    if completed.exists() and (completed/'frozen_method_attestation.json').exists():
        if validate_method_attestation(completed/'frozen_method_attestation.json')!=fingerprint:
            raise ValueError('Frozen production model drifted relative to sealed results')
    pending=[job for job in queue if job['status']=='pending' and
             (partition is None or job['partition']==partition)]
    selected=pending if max_work_units==0 else pending[:max_work_units]
    output.mkdir(parents=True,exist_ok=True)
    completed.mkdir(exist_ok=True)
    stages.mkdir(exist_ok=True)
    state={'schema':'fcp_full_original_photo_roi_v4_execution_progress_v1',
           'status':'INCOMPLETE_FRESH_ORIGINAL_IMAGE_TERMINALS_NOT_H2',
           'started_utc':_utc(), 'source_species':checked['source_species'],
           'all_frozen_source_photo_ids':checked['all_frozen_source_photo_ids'],
           'selected_partition_source_photo_ids':(checked['source_photo_ids']
                                                 if partition is not None else None),
           'completion_audit_scope':('global' if partition is None else 'one_partition'),
           'compute_batch_size_not_photo_cap':batch_size,
           'invocation_work_unit_budget_not_sample_cap':max_work_units,
           'invocation_partition':partition,
           'source_firewall_sha256':checked['source_firewall_sha256'],
           'frozen_original_method':fingerprint,
           'n_prior_completed_work_units':checked['n_completed_work_units'],
           'n_processed_in_this_invocation':0,
           'full_study_original_photo_measurement_complete':False,
           'reused_historical_flower_colours':False,
           'W_or_environment_analyses_done':False}
    _write_atomic(report_path,state)
    for job in selected:
        part,batch=int(job['partition']),int(job['batch'])
        stem=f'terminal_p{part:03}_b{batch:07}'
        stage=stages/stem
        if stage.exists():
            # The previous attempt never produced a complete sealed unit.
            # Preserve its untrusted files rather than silently overwriting it.
            raise ValueError(f'Unsealed interrupted work unit requires review before retry: {stage}')
        stage.mkdir()
        try:
            result=roi.run_batch(firewall,part,batch,batch_size,stage,
                                 source=source,detector=detector,sam=sam)
            terminal=stage/(stem+'.csv')
            receipt=stage/(stem+'.json')
            if (not terminal.is_file() or not receipt.is_file() or
                    result.get('status')!='COMPLETE_WORK_UNIT' or
                    result.get('selected_source_worker_sha256')!=job['source_worker_sha256'] or
                    result.get('n_terminal_photo_ids')!=job['n_source_photo_ids'] or
                    result.get('terminal_csv_sha256')!=sha256_file(terminal) or
                    result.get('method')!=fingerprint):
                raise ValueError('Original image terminal does not match its frozen exact source work unit')
            att=stage/'frozen_method_attestation.json'
            if not att.is_file() or validate_method_attestation(att)!=fingerprint:
                raise ValueError('Frozen original model not actually attested in ROI output')
            target_att=completed/'frozen_method_attestation.json'
            if not target_att.exists():
                shutil.copy2(att,target_att)
            elif _json(target_att)!=fingerprint:
                raise ValueError('Mixed original ROI model versions among sealed work units')
            target_csv=completed/terminal.name
            target_receipt=completed/receipt.name
            if target_csv.exists() or target_receipt.exists():
                raise ValueError('Existing source terminal cannot be overwritten')
            # Terminal receipt is published last. Any crash during the two-file
            # commit yields an orphan, rejected by plan_jobs on the next call.
            os.replace(terminal,target_csv)
            os.replace(receipt,target_receipt)
            state['n_processed_in_this_invocation']+=1
            state['last_complete_unit']={'partition':part,'batch':batch,
                                         'terminal_csv_sha256':result['terminal_csv_sha256'],
                                         'n_source_photo_ids':result['n_source_photo_ids'],
                                         'n_classifiable':result['n_classifiable']}
            _write_atomic(report_path,state)
            print('FCP_ORIGINAL_FULL_PHOTO_TERMINAL_COMPLETE',json.dumps(state['last_complete_unit']),flush=True)
        finally:
            # Remove stage ONLY when its exact checked receipt has reached the
            # sealed set; preserve all incomplete files for explicit diagnosis.
            if (completed/(stem+'.csv')).is_file() and (completed/(stem+'.json')).is_file():
                shutil.rmtree(stage)
    after,_=_sealed_status(firewall,completed,batch_size,partition=partition,
                           _test_only_small_source=_test_only_small_source)
    state.update(n_total_completed_work_units=after['n_completed_work_units'],
                 n_pending_work_units=after['n_pending_work_units'],
                 n_completed_original_photo_terminals=(
                     after['source_photo_ids'] if partition is None and
                     not after['n_pending_work_units'] else None),
                 checked_at_utc=_utc())
    if not after['n_pending_work_units']:
        if partition is None:
            state['status']='ALL_FROZEN_SOURCE_IMAGES_TERMINAL_READY_FOR_REASSEMBLY_NOT_H2'
            state['full_study_original_photo_measurement_complete']=True
        else:
            # Finishing one compute partition does NOT finish 42,111 taxa.
            # Only the independent all-128 collector can verify full coverage.
            state['status']='SELECTED_PARTITION_ORIGINAL_PHOTO_TERMINALS_COMPLETE_NOT_FULL_STUDY'
    _write_atomic(report_path,state)
    return state


def main(argv=None)->int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--firewall',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path)
    p.add_argument('--original-measure-script',required=True,type=Path)
    p.add_argument('--original-detector',required=True,type=Path)
    p.add_argument('--original-efficient-sam',required=True,type=Path)
    p.add_argument('--batch-size',type=int,default=64)
    p.add_argument('--max-work-units',type=int,default=32,
                   help='Compute per invocation; 0 runs all remaining work units (NO photo or species cap)')
    p.add_argument('--partition',type=int,default=None)
    a=p.parse_args(argv)
    r=execute_pending(a.firewall,a.output,
                      source=a.original_measure_script,detector=a.original_detector,
                      sam=a.original_efficient_sam,batch_size=a.batch_size,
                      max_work_units=a.max_work_units,partition=a.partition)
    print('FCP_ORIGINAL_FULL_PHOTO_EXECUTION',json.dumps({
        k:r[k] for k in ('status','source_species','all_frozen_source_photo_ids',
                        'n_total_completed_work_units','n_pending_work_units',
                        'full_study_original_photo_measurement_complete')},sort_keys=True))
    return 0


if __name__=='__main__':
    raise SystemExit(main())

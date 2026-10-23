"""Synthetic-only checks for exact original-photo work-unit execution/resumption."""
import csv
import json

import pandas as pd
import pytest

from fcp_pipeline.fcp_fresh_study import BIO,freeze_fresh,reassemble_fresh
from test_fcp_fresh_study import setup_source,attestation
from scripts.acquisition.execute_fcp_full_photo_work_queue import execute_pending
from scripts.acquisition import execute_fcp_full_photo_work_queue as executor


def sample(tmp_path):
    p,root=setup_source(tmp_path)
    firewall=tmp_path/'firewall'
    freeze_fresh(p,root,firewall,shards=2,workers=2,expected=4)
    return firewall


def patch_full_method(monkeypatch,tmp_path):
    method=json.loads(attestation(tmp_path).read_text())
    monkeypatch.setattr(executor.roi,'verify_frozen_method',lambda *args:method)
    monkeypatch.setattr(executor.roi,'image_loader',lambda url:b'real-test-bytes')
    monkeypatch.setattr(executor.roi.time,'sleep',lambda _:None)
    def evaluate(worker,images,source,detector,sam,out):
        source_rows=list(csv.DictReader(worker.open()))
        assert source_rows and set(source_rows[0])=={'measurement_id','image_filename','photo_license'}
        assert all((images/r['image_filename']).exists() for r in source_rows)
        return pd.DataFrame([{'measurement_id':r['measurement_id'],
                              'morph':'white','measurement_status':'classified_four_state_morph',
                              **{f'flower_fraction_{c}':float(c=='white') for c in BIO}}
                             for r in source_rows])
    monkeypatch.setattr(executor.roi,'evaluate_with_frozen_source',evaluate)
    return method


def call(firewall,out,*,budget,partition=None):
    return execute_pending(firewall,out,
        source=out/'original_measure_source.py',detector=out/'original_detector.pt',
        sam=out/'original_sam',batch_size=2,max_work_units=budget,
        partition=partition,_test_only_small_source=True)


def test_all_photos_are_measured_across_restarts_and_none_used_as_H2(tmp_path,monkeypatch):
    firewall=sample(tmp_path)
    method=patch_full_method(monkeypatch,tmp_path)
    out=tmp_path/'queue'
    first=call(firewall,out,budget=2)
    assert first['n_processed_in_this_invocation']==2
    assert first['status']=='INCOMPLETE_FRESH_ORIGINAL_IMAGE_TERMINALS_NOT_H2'
    assert not first['full_study_original_photo_measurement_complete']
    assert len(list((out/'sealed_terminal_batches').glob('terminal_*.csv')))==2
    last=call(firewall,out,budget=0)
    assert last['status']=='ALL_FROZEN_SOURCE_IMAGES_TERMINAL_READY_FOR_REASSEMBLY_NOT_H2'
    assert last['n_pending_work_units']==0
    assert last['n_completed_original_photo_terminals']==8
    assert last['all_frozen_source_photo_ids']==8
    assert last['W_or_environment_analyses_done'] is False
    result=reassemble_fresh(firewall,out/'sealed_terminal_batches',attestation(tmp_path),
                            tmp_path/'fresh_outcomes')
    assert result['status']=='COMPLETE_FRESH_MEASURED_PHOTO_CENSUS'
    assert result['n_terminal_measurements']==8
    assert result['n_classifiable']==8
    assert result['frozen_H2_null_computed'] is False
    assert json.loads((out/'original_photo_work_queue_progress.json').read_text())['n_pending_work_units']==0
    assert not list((out/'unsealed_inflight_batches').iterdir())


def test_worker_subset_forgery_never_relabels_completion(tmp_path,monkeypatch):
    firewall=sample(tmp_path)
    patch_full_method(monkeypatch,tmp_path)
    out=tmp_path/'queue'
    call(firewall,out,budget=1)
    receipt=next((out/'sealed_terminal_batches').glob('terminal_*.json'))
    data=json.loads(receipt.read_text())
    data['selected_source_worker_sha256']='0'*64
    receipt.write_text(json.dumps(data))
    with pytest.raises(ValueError,match='changed or synthetic'):
        call(firewall,out,budget=0)
    assert len(list((out/'sealed_terminal_batches').glob('terminal_*.csv')))==1


def test_interrupted_inflight_never_silently_overwritten(tmp_path,monkeypatch):
    firewall=sample(tmp_path)
    patch_full_method(monkeypatch,tmp_path)
    out=tmp_path/'queue'
    first=call(firewall,out,budget=1)
    assert first['n_pending_work_units']>0
    from scripts.acquisition.plan_fcp_full_measurements import plan_jobs
    plan_jobs(firewall,out/'sealed_terminal_batches',tmp_path/'check',batch_size=2,
              _test_only_small_source=True)
    queue=[json.loads(s) for s in (tmp_path/'check/all_photo_ROI_v4_work_units.ndjson').read_text().splitlines()]
    pending=next(j for j in queue if j['status']=='pending')
    stem=f"terminal_p{pending['partition']:03}_b{pending['batch']:07}"
    quarantined=out/'unsealed_inflight_batches'/stem
    quarantined.mkdir()
    (quarantined/'crashed_note.txt').write_text('unsealed interrupted compute marker')
    with pytest.raises(ValueError,match='requires review'):
        call(firewall,out,budget=0)
    assert (quarantined/'crashed_note.txt').is_file()


def test_systemic_provider_failure_keeps_work_unit_unsealed_for_explicit_recovery(tmp_path,monkeypatch):
    """Not even a fully enumerated photo corpus may count 429 as completion."""
    firewall=sample(tmp_path)
    patch_full_method(monkeypatch,tmp_path)
    output=tmp_path/'queue'
    def outage(_url):
        raise executor.roi.ProviderTransientAccessError('synthetic HTTP 429 exhaustion')
    monkeypatch.setattr(executor.roi,'image_loader',outage)
    with pytest.raises(executor.roi.ProviderTransientAccessError,match='HTTP 429'):
        call(firewall,output,budget=1)
    assert not list((output/'sealed_terminal_batches').glob('terminal_*.csv'))
    assert not list((output/'sealed_terminal_batches').glob('terminal_*.json'))
    pending=list((output/'unsealed_inflight_batches').glob('terminal_*'))
    assert len(pending)==1
    progress=json.loads((output/'original_photo_work_queue_progress.json').read_text())
    assert progress['full_study_original_photo_measurement_complete'] is False
    # A retry must first explicitly reconcile the unsealed interrupted stage.
    with pytest.raises(ValueError,match='requires review'):
        call(firewall,output,budget=1)



def test_partition_executor_cannot_claim_all_photos_after_one_partition(tmp_path,monkeypatch):
    firewall=sample(tmp_path)
    patch_full_method(monkeypatch,tmp_path)
    only=call(firewall,tmp_path/'partition_queue',budget=0,partition=0)
    assert only['completion_audit_scope']=='one_partition'
    assert only['n_pending_work_units']==0
    assert only['status']=='SELECTED_PARTITION_ORIGINAL_PHOTO_TERMINALS_COMPLETE_NOT_FULL_STUDY'
    assert only['full_study_original_photo_measurement_complete'] is False
    assert only['n_completed_original_photo_terminals'] is None
    assert only['all_frozen_source_photo_ids']==8
    assert 0<only['selected_partition_source_photo_ids']<8
    whole=call(firewall,tmp_path/'global_queue',budget=0)
    assert whole['completion_audit_scope']=='global'
    assert whole['full_study_original_photo_measurement_complete'] is True
    assert whole['n_completed_original_photo_terminals']==8

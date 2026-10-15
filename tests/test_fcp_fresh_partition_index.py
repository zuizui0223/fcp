"""No-network source-index tests; all photographs and taxa are synthetic."""
from __future__ import annotations

import csv
import json
import math

import pytest

from fcp_pipeline.fcp_fresh_study import freeze_fresh
from scripts.acquisition.fcp_fresh_partition_index import FrozenPartitionSourceIndex
from scripts.acquisition import run_fcp_fresh_image_partition as roi
from scripts.acquisition import execute_fcp_full_photo_work_queue as executor
from test_fcp_fresh_study import setup_source
from test_execute_fcp_full_photo_work_queue import patch_full_method, call


def sample(tmp_path):
    frame,source=setup_source(tmp_path)
    firewall=tmp_path/'frozen'
    freeze_fresh(frame,source,firewall,shards=2,workers=2,expected=4)
    return firewall


def test_indexed_photo_batches_match_traditional_frozen_metadata_preflight(tmp_path):
    firewall=sample(tmp_path)
    for part in range(2):
        index=FrozenPartitionSourceIndex(firewall,part)
        for batch_size in (1,2,3):
            for batch in range(math.ceil(index.total_photo_ids/batch_size)):
                worker,acq,count=index.get_batch(part,batch,batch_size)
                old_worker,old_acq,old_count=roi.preflight(
                    firewall,part,batch,batch_size)
                assert worker.equals(old_worker)
                assert acq.equals(old_acq)
                assert count==old_count
                assert set(worker.columns)=={'measurement_id','image_filename','photo_license'}
                assert 'photo_url_large' not in worker.columns
                assert 'photo_url_large' in acq.columns
        if index.total_photo_ids:
            with pytest.raises(ValueError,match='past last source photograph'):
                index.get_batch(part,index.total_photo_ids,1)
    with pytest.raises(ValueError,match='outside frozen partition'):
        FrozenPartitionSourceIndex(firewall,2)


def test_index_rejects_modified_url_and_crosswired_source_identity(tmp_path):
    firewall=sample(tmp_path)
    acquired=firewall/'acquisition_partition_000.csv'
    old=acquired.read_bytes()
    acquired.write_bytes(old.replace(b'example.org',b'other.example.org'))
    with pytest.raises(ValueError,match='SHA256 changed'):
        FrozenPartitionSourceIndex(firewall,0)
    acquired.write_bytes(old)
    assert FrozenPartitionSourceIndex(firewall,0).total_photo_ids>=2
    with acquired.open(newline='') as h:
        rows=list(csv.DictReader(h))
    rows[0]['measurement_id']='SYNTHETIC_TAMPERED_PHOTO_ID'
    with acquired.open('w',newline='') as h:
        out=csv.DictWriter(h,fieldnames=list(rows[0]))
        out.writeheader();out.writerows(rows)
    with pytest.raises(ValueError,match='identity not source-bound'):
        FrozenPartitionSourceIndex(firewall,0)


def test_execution_indexes_one_partition_once_for_multiple_genuine_work_units(tmp_path,monkeypatch):
    firewall=sample(tmp_path)
    patch_full_method(monkeypatch,tmp_path)
    invocations=[]
    original_class=executor.roi.FrozenPartitionSourceIndex
    class CountedIndex(original_class):
        def __init__(self,source,partition):
            invocations.append(partition)
            super().__init__(source,partition)
    monkeypatch.setattr(executor.roi,'FrozenPartitionSourceIndex',CountedIndex)
    def forbid_reparse(*a,**k):
        raise AssertionError('Never reload entire worker/acquisition pandas CSV per batch')
    monkeypatch.setattr(executor.roi,'preflight',forbid_reparse)
    out=tmp_path/'partition_outputs'
    result=call(firewall,out,budget=0,partition=0)
    assert result['n_processed_in_this_invocation']>=2
    assert result['completion_audit_scope']=='one_partition'
    assert invocations==[0]
    assert result['full_study_original_photo_measurement_complete'] is False
    assert result['n_pending_work_units']==0
    terminals=list((out/'sealed_terminal_batches').glob('terminal_*.csv'))
    assert len(terminals)==result['n_processed_in_this_invocation']
    assert sum(sum(1 for _ in csv.DictReader(p.open()))
               for p in terminals)==result['selected_partition_source_photo_ids']


def test_index_does_not_trust_empty_source_or_invalid_work_unit(tmp_path):
    firewall=sample(tmp_path)
    index=FrozenPartitionSourceIndex(firewall,0)
    with pytest.raises(ValueError,match='outside frozen photo source'):
        index.get_batch(1,0,1)
    with pytest.raises(ValueError,match='outside frozen photo source'):
        index.get_batch(0,-1,1)
    with pytest.raises(ValueError,match='outside frozen photo source'):
        index.get_batch(0,0,1025)

"""Resumable all-photo work-unit planner, with synthetic fixture only."""
from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

from test_fcp_fresh_study import setup_source, terminal, attestation
from fcp_pipeline.fcp_fresh_study import freeze_fresh
from scripts.acquisition.plan_fcp_full_measurements import plan_jobs


def make_firewall(tmp_path):
    p,root=setup_source(tmp_path)
    f=tmp_path/'frozen'
    freeze_fresh(p,root,f,shards=2,expected=4,workers=2)
    return f


def test_all_source_images_in_work_units_without_species_or_photo_cap(tmp_path):
    f=make_firewall(tmp_path)
    q=plan_jobs(f,tmp_path/'nothing_completed',tmp_path/'plan',batch_size=3,
                _test_only_small_source=True)
    assert q['source_species']==4
    assert q['source_photo_ids']==8
    assert q['work_units']>=3
    assert q['n_completed_work_units']==0
    assert q['n_pending_work_units']==q['work_units']
    jobs=[json.loads(x) for x in (tmp_path/'plan/all_photo_ROI_v4_work_units.ndjson').read_text().splitlines()]
    assert sum(x['n_source_photo_ids'] for x in jobs)==8
    assert all('photo_id' not in j and 'latitude' not in j and 'species' not in j for j in jobs)
    with pytest.raises(ValueError,match='Main 42111 species'):
        plan_jobs(f,tmp_path/'nothing_completed',tmp_path/'invalid',batch_size=3)


def test_resume_rejects_forged_or_mixed_terminal_outcomes(tmp_path):
    f=make_firewall(tmp_path)
    t=terminal(tmp_path,f)
    shutil.copy2(attestation(tmp_path),t/'frozen_method_attestation.json')
    q=plan_jobs(f,t,tmp_path/'complete_plan',batch_size=100,_test_only_small_source=True)
    assert q['n_completed_work_units']==q['work_units']
    assert q['n_pending_work_units']==0
    # A bogus completed receipt never enters a resumable queue as completed.
    p=next(t.glob('terminal_*.json'))
    d=json.loads(p.read_text())
    d['selected_source_worker_sha256']='0'*64
    p.write_text(json.dumps(d))
    with pytest.raises(ValueError,match='synthetic terminal batch|changed or synthetic'):
        plan_jobs(f,t,tmp_path/'tampered_plan',batch_size=100,_test_only_small_source=True)



def test_partition_only_plan_is_fast_scope_not_global_completion(tmp_path):
    f=make_firewall(tmp_path)
    part=plan_jobs(f,tmp_path/'sealed',tmp_path/'part_only',batch_size=2,
                   partition=1,_test_only_small_source=True)
    assert part['partition_scope']==1
    assert part['is_global_all_photo_plan'] is False
    assert part['all_photo_ids_included_once'] is False
    assert part['status']=='SELECTED_ORIGINAL_PHOTO_PARTITION_PLAN_NOT_GLOBAL_COMPLETION'
    assert part['all_frozen_source_photo_ids']==8
    assert 0<part['source_photo_ids']<8
    assert part['n_pending_work_units']==part['work_units']
    work=[json.loads(x) for x in
          (tmp_path/'part_only/all_photo_ROI_v4_work_units.ndjson').read_text().splitlines()]
    assert work and all(x['partition']==1 for x in work)
    assert sum(x['n_source_photo_ids'] for x in work)==part['source_photo_ids']
    full=plan_jobs(f,tmp_path/'sealed',tmp_path/'whole_source',batch_size=2,
                   _test_only_small_source=True)
    assert full['is_global_all_photo_plan'] is True
    assert full['all_photo_ids_included_once'] is True
    assert full['source_photo_ids']==8
    assert full['n_pending_work_units']>part['n_pending_work_units']
    with pytest.raises(ValueError,match='outside frozen partition'):
        plan_jobs(f,tmp_path/'sealed',tmp_path/'bad_partition',batch_size=2,
                  partition=2,_test_only_small_source=True)

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

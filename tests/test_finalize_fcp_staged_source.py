"""Real-shaped synthetic four-table stage -> independent photo shards -> admission."""
from __future__ import annotations

import json
import tarfile

import pytest

from scripts.acquisition.execute_fcp_42111_census import execute
from scripts.acquisition.emit_fcp_staged_source_shards import emit_shards
from scripts.acquisition.finalize_fcp_staged_source import finalize_staged_work
from scripts.acquisition.admit_fcp_complete_photo_source import admit_source
from test_execute_fcp_42111_census import pack_archive
from test_fcp_open_data_bulk_import import fixture


def source_work(tmp_path):
    frame,taxa,obs,photos,observers,_=fixture(tmp_path)
    provider_archive=pack_archive(tmp_path/'original_20260927.tar.gz',{
        'taxa':taxa,'observations':obs,'photos':photos,'observers':observers})
    root=tmp_path/'durable_source'
    staged=execute(frame,root,source_archive=provider_archive,
                   species=4,shards=2,workers=2,
                   stage_only=True,skip_package=True)
    assert staged['status']=='SOURCE_SQLITE_STAGED_FOR_INDEPENDENT_SHARDS_NOT_FULL_CENSUS'
    assert staged['complete_metadata_census'] is False
    assert staged['staged_original_source_photo_ids']==3
    assert staged['staged_sqlite_bytes']>0
    return frame,root


def test_complete_split_source_ready_for_real_photo_ingest_not_old_colours(tmp_path):
    frame,root=source_work(tmp_path)
    metadata=root/'complete_42111_metadata_source'
    with pytest.raises(FileNotFoundError):
        finalize_staged_work(root,frame,species=4,shards=2,workers=2,
                             skip_package=True)
    emit_shards(frame,metadata,start=0,stop=1,expected=4,shards=2)
    with pytest.raises(FileNotFoundError):
        finalize_staged_work(root,frame,species=4,shards=2,workers=2,
                             skip_package=True)
    emit_shards(frame,metadata,start=1,stop=2,expected=4,shards=2)
    receipt=finalize_staged_work(root,frame,species=4,shards=2,workers=2,
                                 skip_package=False)
    assert receipt['complete_metadata_census'] is True
    assert receipt['complete_original_images_measured'] is False
    assert receipt['n_all_available_licensed_original_photo_ids']==3
    assert receipt['species_count_checkpoints_audited']==4
    assert receipt['frozen_location_blind_worker_partitions']==2
    assert receipt['source_shards_sealed_from_one_original_sqlite']==2
    package=root/'fcp_complete_42111_source_and_blind_firewall.tar.gz'
    with tarfile.open(package,'r:gz') as t:
        names=t.getnames()
        assert not any('_working_fcp_snapshot.sqlite' in s for s in names)
        assert any(s.endswith('frozen_location_blind_all_photo_firewall/fresh_study_freeze.json') for s in names)
        assert any(s.endswith('complete_42111_metadata_source/bulk_source_receipt.json') for s in names)
    assert (metadata/'_working_fcp_snapshot.sqlite').exists()
    admission=admit_source(
        root,root/'full_source_execution_receipt.json',
        root/'frozen_4_metadata_tables/fcp_source_archive_pin.json',
        root/'all_42111_species_source_photo_counts.csv',
        tmp_path/'admitted',batch_size=2,_test_only_small_source=True)
    assert admission['status']=='FULL_SOURCE_PHOTO_WORK_PLAN_READY_NO_PIXELS_YET'
    assert admission['source_photo_ids']==3
    assert admission['new_original_photo_colours_measured']==0
    assert admission['n_species_checkpoints_verified']==4


def test_staged_source_cannot_claim_complete_before_all_shards(tmp_path):
    frame,root=source_work(tmp_path)
    meta=root/'complete_42111_metadata_source'
    emit_shards(frame,meta,start=0,stop=1,expected=4,shards=2)
    with pytest.raises(FileNotFoundError):
        finalize_staged_work(root,frame,species=4,shards=2,workers=2)
    record=json.loads((root/'full_source_execution_receipt.json').read_text())
    assert record['complete_metadata_census'] is False
    assert not (root/'frozen_location_blind_all_photo_firewall').exists()
    assert not (root/'all_42111_species_source_photo_counts.csv').exists()
    assert not (root/'fcp_complete_42111_source_and_blind_firewall.tar.gz').exists()

"""Synthetic split-source census parity; no provider downloads and no flower pixels."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from fcp_pipeline.fcp_allphotos import sha256_file
from fcp_pipeline.fcp_fresh_study import required_source_dirs, freeze_fresh
from scripts.acquisition.import_inat_open_data_snapshot import ingest
from scripts.acquisition.emit_fcp_staged_source_shards import (
    emit_shards,finalize_shards,
)
from test_fcp_open_data_bulk_import import fixture


def stage_and_full(tmp_path):
    frame,taxa,observations,photos,observers,hashes=fixture(tmp_path)
    pin=tmp_path/'original_dated_archive_pin.json'
    pin.write_text(json.dumps({
        'schema':'fcp_dated_inaturalist_metadata_bundle_pin_v1',
        'status':'FOUR_ORIGINAL_METADATA_TABLES_EXTRACTED_AND_PINNED',
        'snapshot_date':'2026-09-27',
        's3_object_key':'metadata/inaturalist-open-data-20260927.tar.gz',
        'metadata_files_sha256':hashes,
        'full_bundle_sha256':'a'*64,
        'source_images_opened':False,
        'historical_photo_colours_read':False,
    }))
    original=tmp_path/'monolithic_expected'
    mono=ingest(frame,taxa,observations,photos,original,'2026-09-27',
                observers=observers,source_archive_pin=pin,
                shards=2,expected_species=4,input_sha256=hashes)
    staged=tmp_path/'durable_staged_sqlite'
    stage=ingest(frame,taxa,observations,photos,staged,'2026-09-27',
                 observers=observers,source_archive_pin=pin,
                 shards=2,expected_species=4,input_sha256=hashes,
                 stage_sqlite_only=True)
    return frame,original,mono,staged,stage


def test_disjoint_shards_match_monolithic_4table_original_source_bytes(tmp_path):
    frame,original,mono,staged,stage=stage_and_full(tmp_path)
    assert stage['status']=='SOURCE_SQLITE_STAGED_SHARDS_NOT_YET_EMITTED'
    assert stage['full_source_census_complete'] is False
    assert stage['source_image_pixels_opened'] is False
    assert (staged/'_working_fcp_snapshot.sqlite').is_file()
    assert not list(staged.glob('shard_*/species/*/checkpoint.json'))
    left=emit_shards(frame,staged,start=0,stop=1,expected=4,shards=2)
    assert left['completed_here']==1
    with pytest.raises(FileNotFoundError):
        finalize_shards(frame,staged,expected=4,shards=2)
    right=emit_shards(frame,staged,start=1,stop=2,expected=4,shards=2)
    assert right['completed_here']==1
    assert left['newly_or_previously_sealed_photo_ids_in_range']+right[
        'newly_or_previously_sealed_photo_ids_in_range']==mono['n_target_taxon_licensed_photos']
    # Completed shards can be rechecked without re-extracting 36GB or
    # writing a new species/photo manifest.
    repeated=emit_shards(frame,staged,start=0,stop=1,expected=4,shards=2)
    assert repeated['newly_or_previously_sealed_photo_ids_in_range']==left[
        'newly_or_previously_sealed_photo_ids_in_range']
    sealed=finalize_shards(frame,staged,expected=4,shards=2)
    assert sealed['status']=='COMPLETE_SNAPSHOT_METADATA_ALL_NOMINATED_TAXA'
    assert sealed['n_target_taxon_licensed_photos']==mono['n_target_taxon_licensed_photos']==3
    assert sealed['n_conflicting_photo_ids_quarantined']==mono[
        'n_conflicting_photo_ids_quarantined']
    assert sealed['number_independently_sealed_source_shards']==2
    assert sealed['old_photo_measurements_read'] is False
    for src in original.glob('shard_*_of_*/species/*/photo_manifest.ndjson'):
        rel=src.relative_to(original)
        assert sha256_file(src)==sha256_file(staged/rel)
        assert (staged/rel.parent/'checkpoint.json').is_file()
    dirs,audit=required_source_dirs(frame,staged,shards=2,expected=4)
    assert audit['status']=='SOURCE_CENSUS_COMPLETE'
    assert len(dirs)==4
    fw=tmp_path/'frozen'
    frozen=freeze_fresh(frame,staged,fw,shards=2,expected=4,workers=2)
    assert frozen['all_available_licensed_photo_ids']==3


def test_staged_source_rejects_corrupted_db_and_incomplete_shards(tmp_path):
    frame,original,mono,staged,stage=stage_and_full(tmp_path)
    db=staged/'_working_fcp_snapshot.sqlite'
    with db.open('ab') as handle:
        handle.write(b'BULK_SOURCE_TAMPERING')
    with pytest.raises(ValueError,match='changed'):
        emit_shards(frame,staged,start=0,stop=1,expected=4,shards=2)
    # Restaging is required before any shard can be accepted, never
    # resolving a damaged database by substituting old colour results.


def test_partial_shard_is_not_silently_resumed_or_global_complete(tmp_path):
    frame,original,mono,staged,stage=stage_and_full(tmp_path)
    orphan=staged/'shard_000_of_002'
    orphan.mkdir()
    (orphan/'crashed_metadata.txt').write_text('incomplete source shard')
    with pytest.raises(ValueError,match='quarantine'):
        emit_shards(frame,staged,start=0,stop=1,expected=4,shards=2)
    assert not (staged/'bulk_source_receipt.json').exists()
    assert not any(staged.glob('shard_*/shard_source_receipt.json'))


def test_source_shard_receipt_tampering_fails_final_audit(tmp_path):
    frame,original,mono,staged,stage=stage_and_full(tmp_path)
    emit_shards(frame,staged,start=0,stop=2,expected=4,shards=2)
    shard=staged/'shard_001_of_002/shard_source_receipt.json'
    record=json.loads(shard.read_text())
    record['n_original_source_photo_ids']+=1
    shard.write_text(json.dumps(record))
    with pytest.raises(ValueError,match='counts do not match'):
        finalize_shards(frame,staged,expected=4,shards=2)


def test_out_of_range_partition_emission_fails_before_any_source_change(tmp_path):
    frame,original,mono,staged,stage=stage_and_full(tmp_path)
    with pytest.raises(ValueError,match='outside full immutable'):
        emit_shards(frame,staged,start=2,stop=3,expected=4,shards=2)
    assert not list(staged.glob('shard_*/shard_source_receipt.json'))



def test_full_shard_range_hashes_and_integrity_checks_global_db_once(tmp_path,monkeypatch):
    """A 128-shard emission must not rehash the huge SQLite 128/8 times."""
    import scripts.acquisition.emit_fcp_staged_source_shards as shard_module
    frame,original,mono,staged,stage=stage_and_full(tmp_path)
    checks=[]
    verify=shard_module._validate_stage

    def count_source_integrity(*args,**kwargs):
        checks.append(True)
        return verify(*args,**kwargs)

    monkeypatch.setattr(shard_module,'_validate_stage',count_source_integrity)
    result=shard_module.emit_shards(
        frame,staged,start=0,stop=2,expected=4,shards=2)
    assert len(checks)==1
    assert result['completed_here']==2
    assert result['newly_or_previously_sealed_photo_ids_in_range']==3
    assert len(list(staged.glob('shard_*_of_*/shard_source_receipt.json')))==2
    assert shard_module.finalize_shards(
        frame,staged,expected=4,shards=2
    )['n_target_taxon_licensed_photos']==3
    assert len(checks)==2   # Independent full-source finalization also audits.


def test_reusing_all_sealed_shards_never_rewrites_original_photo_manifests(tmp_path):
    frame,original,mono,staged,stage=stage_and_full(tmp_path)
    emit_shards(frame,staged,start=0,stop=2,expected=4,shards=2)
    frozen={str(p.relative_to(staged)):sha256_file(p)
            for p in staged.glob('shard_*_of_*/species/*/photo_manifest.ndjson')}
    assert frozen
    result=emit_shards(frame,staged,start=0,stop=2,expected=4,shards=2)
    assert result['completed_here']==2
    assert frozen=={
        str(p.relative_to(staged)):sha256_file(p)
        for p in staged.glob('shard_*_of_*/species/*/photo_manifest.ndjson')}



def test_real_source_progress_receipt_is_diagnostic_only_and_preserves_stage(tmp_path):
    frame,original,mono,staged,stage=stage_and_full(tmp_path)
    trace=staged/'fcp_source_progress.ndjson'
    events=[json.loads(line) for line in trace.read_text().splitlines()]
    phases=[e['phase'] for e in events]
    assert phases[0]=='begin_verified_provider_scan'
    assert 'all_original_observations_scanned' in phases
    assert 'all_original_photos_scanned' in phases
    assert 'original_photo_ids_globally_deduplicated' in phases
    assert phases[-1]=='source_sqlite_staged_ready_for_independent_shards'
    assert events[-1]['n_candidate_rows_retained']==3
    assert all(e['source_census_complete'] is False for e in events)
    assert all(e['original_pixels_measured'] is False for e in events)
    assert all(e['source_frame_sha256']==sha256_file(frame) for e in events)
    assert all(e['scratch_free_bytes']>=0 for e in events)
    assert all(e['sqlite_file_bytes']>=0 for e in events)
    assert stage['full_source_census_complete'] is False
    assert not (staged/'bulk_source_receipt.json').exists()

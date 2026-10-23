"""Synthetic full four-table tar end-to-end; no biological images or previous outcomes."""
from __future__ import annotations

import gzip
import io
import json
import tarfile
from pathlib import Path

import pytest

from test_fcp_open_data_bulk_import import fixture
from scripts.acquisition.execute_fcp_42111_census import execute


def pack_archive(destination: Path, tables: dict[str, Path]) -> Path:
    with tarfile.open(destination, 'w:gz') as tar:
        folder = tarfile.TarInfo('inaturalist-open-data-20260927/')
        folder.type = tarfile.DIRTYPE
        tar.addfile(folder)
        for name, source in tables.items():
            with gzip.open(source, 'rb') as fh:
                raw=fh.read()
            member=tarfile.TarInfo(f'inaturalist-open-data-20260927/{name}.csv')
            member.size=len(raw)
            tar.addfile(member, io.BytesIO(raw))
    return destination


def test_source_census_end_to_end_from_4_table_archive_without_any_photo_pixels(tmp_path):
    frame, taxa, observations, photos, observers, _ = fixture(tmp_path)
    original = pack_archive(tmp_path/'20260927_source.tar.gz', {
        'observations':observations,'photos':photos,'taxa':taxa,'observers':observers})
    result=execute(frame, tmp_path/'complete_source', source_archive=original,
                   species=4, shards=2, workers=2, skip_package=True)
    assert result['complete_metadata_census'] is True
    assert result['complete_original_images_measured'] is False
    assert result['n_all_available_licensed_original_photo_ids']==3
    assert result['species_count_checkpoints_audited']==4
    assert result['frozen_location_blind_worker_partitions']==2
    assert result['status']=='COMPLETE_42111_SOURCE_METADATA_FROZEN_PIXEL_MEASUREMENT_STILL_NOT_DONE'
    assert result['old_measured_flower_colours_used'] is False
    frozen=tmp_path/'complete_source/frozen_location_blind_all_photo_firewall'
    assert (frozen/'fresh_study_freeze.json').exists()
    for path in frozen.glob('worker_partition_*.csv'):
        assert path.read_text().splitlines()[0]=='measurement_id,image_filename,photo_license'
    ledger=(tmp_path/'complete_source/all_42111_species_source_photo_counts.csv').read_text().splitlines()
    assert len(ledger)==5
    assert original.exists()  # user-supplied archive never deleted


def test_census_failure_leaves_explicit_incomplete_receipt_without_fake_results(tmp_path):
    frame, taxa, observations, photos, observers, _ = fixture(tmp_path)
    original=pack_archive(tmp_path/'20260927_source.tar.gz', {
        'observations':observations,'photos':photos,'taxa':taxa,'observers':observers})
    # A wrong nominated frame count must never be converted into a partial success.
    with pytest.raises(ValueError):
        execute(frame,tmp_path/'not_ready',source_archive=original,species=5,
                shards=2,workers=2,skip_package=True)
    receipt=json.loads((tmp_path/'not_ready/full_source_execution_receipt.json').read_text())
    assert receipt['complete_metadata_census'] is False
    assert receipt['complete_original_images_measured'] is False
    assert receipt['status']=='FAILED_OR_UNFINISHED_SOURCE_METADATA_NOT_SCIENTIFIC_RESULTS'


def test_verified_four_table_scratch_is_evicted_only_after_full_import(tmp_path):
    frame,taxa,observations,photos,observers,_=fixture(tmp_path)
    original=pack_archive(tmp_path/'dated_original.tar.gz',{
        'taxa':taxa,'observations':observations,'photos':photos,'observers':observers})
    out=tmp_path/'complete_pruned_source'
    got=execute(frame,out,source_archive=original,species=4,shards=2,workers=2,
                skip_package=True,release_verified_temporary_tables=True)
    assert got['complete_metadata_census'] is True
    assert got['pruned_verified_scratch_tables_after_scan'] is True
    assert got['pruned_verified_scratch_table_bytes']>0
    extracted=out/'frozen_4_metadata_tables'
    assert (extracted/'fcp_source_archive_pin.json').is_file()
    assert not list(extracted.glob('*.csv.gz'))
    assert original.is_file()  # Never delete caller-owned source archive.
    assert got['complete_original_images_measured'] is False
    from fcp_pipeline.fcp_fresh_study import required_source_dirs
    dirs,audit=required_source_dirs(
        frame,out/'complete_42111_metadata_source',shards=2,expected=4)
    assert audit['status']=='SOURCE_CENSUS_COMPLETE'
    assert len(dirs)==4
    assert json.loads((out/'complete_42111_metadata_source/bulk_source_receipt.json').read_text())[
        'pruned_pinned_temporary_metadata_after_full_scan'] is True


def test_aborted_long_stage_persists_true_in_progress_phase(tmp_path, monkeypatch):
    # Reproduces the 2026-10-15 355-minute GitHub-hosted timeout. The saved
    # source receipt must NOT falsely claim that only tar extraction happened.
    from scripts.acquisition.execute_fcp_42111_census import run
    import scripts.acquisition.execute_fcp_42111_census as module
    receipt=tmp_path/'full_source_execution_receipt.json'
    status={'status':'FAILED_OR_UNFINISHED_SOURCE_METADATA_NOT_SCIENTIFIC_RESULTS',
            'last_completed_stage':'extract_all_four_and_pin_exact_SHA256',
            'complete_metadata_census':False,
            'complete_original_images_measured':False}
    def persist():
        temp=receipt.with_suffix('.json.pending')
        temp.write_text(json.dumps(status))
        temp.replace(receipt)
    def interrupted(*args, **kwargs):
        observed=json.loads(receipt.read_text())
        assert observed['stage']=='scan_entire_dated_source_and_emit_all_photo_manifests'
        assert observed['stage_in_progress'] is True
        assert observed['last_completed_stage']=='extract_all_four_and_pin_exact_SHA256'
        raise KeyboardInterrupt('simulation: hosted runner SIGTERM/timeout')
    monkeypatch.setattr(module.subprocess,'run',interrupted)
    with pytest.raises(KeyboardInterrupt):
        run(['never-executed'],cwd=tmp_path,
            stage='scan_entire_dated_source_and_emit_all_photo_manifests',
            status=status,persist_stage=persist)
    observed=json.loads(receipt.read_text())
    assert observed['complete_metadata_census'] is False
    assert observed['complete_original_images_measured'] is False
    assert observed['stage_in_progress'] is True
    assert not receipt.with_suffix('.json.pending').exists()


def test_successful_long_stage_persists_completed_stage(tmp_path,monkeypatch):
    from scripts.acquisition.execute_fcp_42111_census import run
    import scripts.acquisition.execute_fcp_42111_census as module
    receipt=tmp_path/'full_source_execution_receipt.json'
    status={'last_completed_stage':'extract_all_four_and_pin_exact_SHA256'}
    def persist():
        receipt.write_text(json.dumps(status))
    monkeypatch.setattr(module.subprocess,'run',lambda *args,**kwargs:None)
    run(['a-noop'],cwd=tmp_path,stage='scan_entire_dated_source_and_emit_all_photo_manifests',
        status=status,persist_stage=persist)
    observed=json.loads(receipt.read_text())
    assert observed['stage_in_progress'] is False
    assert observed['last_completed_stage']==observed['stage']


def test_photo_opportunities_cover_all_nominated_taxa_without_opening_original_pixels(tmp_path):
    # Real production's first bounded output: exact all-species source counts.
    # A complete COUNT output is not the frozen photo-ID firewall or W.
    import csv
    frame,taxa,observations,photos,observers,_=fixture(tmp_path)
    original=pack_archive(tmp_path/'source_counts_archive.tar.gz',{
        'taxa':taxa,'observations':observations,'photos':photos,'observers':observers})
    out=tmp_path/'counts_study'
    got=execute(frame,out,source_archive=original,species=4,shards=2,workers=2,
                photo_counts_only=True)
    assert got['status']=='COMPLETE_REAL_42111_PHOTO_OPPORTUNITY_COUNTS_ONLY_NOT_FULL_SOURCE'
    assert got['source_photo_opportunity_counts_completed'] is True
    assert got['complete_metadata_census'] is False
    assert got['complete_original_images_measured'] is False
    assert got['photo_id_manifest_universe_frozen'] is False
    assert got['n_frozen_source_photo_opportunities']==3
    with (out/'complete_42111_metadata_source/all_42111_species_photo_opportunity_counts.csv').open() as hdl:
        rows=list(csv.DictReader(hdl))
    assert len(rows)==4
    assert {row['inat_taxon_id']:int(row['n_unique_source_photo_ids']) for row in rows}=={
        '11':2,'12':1,'13':0,'14':0}
    assert not (out/'frozen_location_blind_all_photo_firewall').exists()
    assert original.is_file()
    assert not (out/'complete_42111_metadata_source/_working_fcp_snapshot.sqlite').exists()
    receipt=json.loads((out/'full_source_execution_receipt.json').read_text())
    assert receipt['status']==got['status']
    assert receipt['new_W_or_environment_results_computed'] is False
    assert receipt['last_completed_stage']=='scan_complete_source_and_count_all_species_photo_opportunities'

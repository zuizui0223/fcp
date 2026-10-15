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

"""Synthetic admission of complete fresh photo universe; no original images."""
from __future__ import annotations

import csv
import json

import pytest

from test_fcp_open_data_bulk_import import fixture
from test_execute_fcp_42111_census import pack_archive
from scripts.acquisition.execute_fcp_42111_census import execute
from scripts.acquisition.admit_fcp_complete_photo_source import admit_source


def completed(tmp_path):
    frame,taxa,obs,photos,observers,_=fixture(tmp_path)
    source=pack_archive(tmp_path/'dated_provider.tar.gz',{
        'taxa':taxa,'observations':obs,'photos':photos,'observers':observers})
    package=tmp_path/'source_work'
    result=execute(frame,package,source_archive=source,
                   species=4,shards=2,workers=2,skip_package=True)
    assert result['status']=='COMPLETE_42111_SOURCE_METADATA_FROZEN_PIXEL_MEASUREMENT_STILL_NOT_DONE'
    return (package,package/'full_source_execution_receipt.json',
            package/'frozen_4_metadata_tables/fcp_source_archive_pin.json',
            package/'all_42111_species_source_photo_counts.csv')


def call(parts,out):
    package,status,pin,counts=parts
    return admit_source(package,status,pin,counts,out,batch_size=2,
                        _test_only_small_source=True)


def test_full_source_admitted_only_as_future_photo_work(tmp_path):
    parts=completed(tmp_path)
    result=call(parts,tmp_path/'admitted')
    assert result['status']=='FULL_SOURCE_PHOTO_WORK_PLAN_READY_NO_PIXELS_YET'
    assert result['source_species']==4
    assert result['source_photo_ids']==3
    assert result['n_species_checkpoints_verified']==4
    assert result['worker_partitions']==2
    assert result['pending_work_units']==result['work_units']
    assert result['new_original_photo_colours_measured']==0
    assert not result['w_h2_environment_inference_complete']
    assert list((tmp_path/'admitted'/'plan').glob('all_photo_ROI_v4_work_units.ndjson'))
    with pytest.raises(ValueError,match='empty'):
        call(parts,tmp_path/'admitted')
    with pytest.raises(ValueError,match='42,111'):
        admit_source(*parts,tmp_path/'production_claim',batch_size=2)


def test_partial_source_receipt_forbids_photographic_opening(tmp_path):
    parts=completed(tmp_path)
    _,status,_,_=parts
    doc=json.loads(status.read_text())
    doc['complete_metadata_census']=False
    status.write_text(json.dumps(doc))
    with pytest.raises(ValueError,match='incomplete'):
        call(parts,tmp_path/'unauthorized')


def test_wrong_table_hash_or_archive_identity_forbids_handoff(tmp_path):
    parts=completed(tmp_path)
    _,_,pin,_=parts
    doc=json.loads(pin.read_text())
    doc['metadata_files_sha256']['photos']='0'*64
    pin.write_text(json.dumps(doc))
    with pytest.raises(ValueError,match='source mismatch'):
        call(parts,tmp_path/'mixed_source')


def test_missing_checkpoint_refuses_optimistic_species_completion(tmp_path):
    parts=completed(tmp_path)
    package,_,_,_=parts
    checkpoint=next((package/'complete_42111_metadata_source').glob(
        'shard_*_of_*/species/*/checkpoint.json'))
    checkpoint.unlink()
    with pytest.raises(ValueError,match='Missing'):
        call(parts,tmp_path/'missing_source')


def test_photo_count_ledger_cannot_hide_extra_or_lost_source_photo(tmp_path):
    parts=completed(tmp_path)
    _,_,_,counts=parts
    with counts.open('r',newline='') as h:
        rows=list(csv.DictReader(h))
    one=next(r for r in rows if int(r['n_unique_target_photo_ids'])>0)
    one['n_unique_target_photo_ids']=str(int(one['n_unique_target_photo_ids'])+1)
    with counts.open('w',newline='') as h:
        w=csv.DictWriter(h,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)
    with pytest.raises(ValueError,match='photo denominator'):
        call(parts,tmp_path/'altered_counts')

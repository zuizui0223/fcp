#!/usr/bin/env python3
"""Fail-closed handoff of the complete dated FCP source to original-photo work.

This reads ONLY source metadata/checksums and blind worker IDs. It never opens
image pixels, historic photo-colour outcomes, species W or environmental data.
The source-status artifact and complete source-package artifact are separate
GitHub Actions artifacts; both are REQUIRED to admit an actual photo census.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

from fcp_pipeline.fcp_allphotos import sha256_file
from fcp_pipeline.fcp_fresh_study import STUDY
from scripts.acquisition.execute_fcp_42111_census import ARCHIVE_KEY, FROZEN_DATE, FRAME_SHA256
from scripts.acquisition.plan_fcp_full_measurements import plan_jobs

COMPLETE_SOURCE='COMPLETE_42111_SOURCE_METADATA_FROZEN_PIXEL_MEASUREMENT_STILL_NOT_DONE'
COMPLETE_BULK='COMPLETE_SNAPSHOT_METADATA_ALL_NOMINATED_TAXA'
COMPLETE_FREEZE='FRESH_BLIND_METADATA_FROZEN'
TABLES={'taxa','observations','observers','photos'}
HEX=re.compile(r'^[a-f0-9]{64}$')


def read_doc(path:Path)->dict:
    data=json.loads(path.read_text(encoding='utf-8'))
    if not isinstance(data,dict):
        raise ValueError('FCP source certificate must be a JSON object')
    return data


def assert_equal(label:str,*values):
    if len(set(str(x) for x in values))!=1:
        raise ValueError(f'FCP handoff source mismatch: {label}')


def admit_source(package_root:Path, source_status:Path, original_pin:Path,
                 source_counts:Path, output:Path, *,batch_size:int=64,
                 _test_only_small_source:bool=False)->dict:
    """Admit only a fully documented, complete, source-blind photo universe."""
    status=read_doc(source_status)
    pin=read_doc(original_pin)
    bulk=read_doc(package_root/'complete_42111_metadata_source'/'bulk_source_receipt.json')
    firewall=package_root/'frozen_location_blind_all_photo_firewall'
    freeze=read_doc(firewall/'fresh_study_freeze.json')
    count_expected=(4 if _test_only_small_source else 42111)
    packet_expected=(2 if _test_only_small_source else 128)
    if (status.get('schema')!='fcp_full_42111_dated_source_execution_v1'
            or status.get('status')!=COMPLETE_SOURCE
            or status.get('complete_metadata_census') is not True
            or status.get('complete_original_images_measured') is not False
            or status.get('old_measured_flower_colours_used') is not False
            or status.get('target_frame_species')!=count_expected
            or status.get('species_count_checkpoints_audited')!=count_expected
            or status.get('worker_packets')!=packet_expected):
        raise ValueError('42,111 species source execution receipt is incomplete or contradicts full census')
    if (bulk.get('status')!=COMPLETE_BULK
            or bulk.get('taxa_in_initial_frame')!=count_expected
            or bulk.get('old_photo_measurements_read') is not False
            or bulk.get('image_pixels_opened') is not False
            or bulk.get('photo_cap','INVALID') is not None):
        raise ValueError('Four-table source importer never completed an uncapped census')
    if (freeze.get('status')!=COMPLETE_FREEZE
            or freeze.get('study_id')!=STUDY
            or freeze.get('source_species')!=count_expected
            or freeze.get('worker_packets')!=packet_expected
            or freeze.get('historical_measured_CSV_read') is not False
            or freeze.get('new_measurements_completed')!=0
            or freeze.get('reused_old_colour_labels') is not False):
        raise ValueError('Fresh location-blind photo firewall is not fully frozen')
    if (pin.get('schema')!='fcp_dated_inaturalist_metadata_bundle_pin_v1'
            or pin.get('status')!='FOUR_ORIGINAL_METADATA_TABLES_EXTRACTED_AND_PINNED'
            or pin.get('source_images_opened') is not False
            or pin.get('historical_photo_colours_read') is not False):
        raise ValueError('Dated original archive and four source-table hashes are not pinned')
    expected_hashes=pin.get('metadata_files_sha256')
    if (not isinstance(expected_hashes,dict) or set(expected_hashes)!=TABLES
            or any(not isinstance(v,str) or not HEX.fullmatch(v)
                   for v in expected_hashes.values())):
        raise ValueError('The four original table SHA256 values are incomplete')
    if not HEX.fullmatch(str(pin.get('full_bundle_sha256',''))):
        raise ValueError('Original archive digest missing')
    assert_equal('original archive date',status.get('immutable_source_date'),
                 pin.get('snapshot_date'),bulk.get('source_snapshot_date'),
                 freeze.get('frozen_created_d2'),FROZEN_DATE)
    assert_equal('original archive path',status.get('immutable_source_key'),
                 pin.get('s3_object_key'),ARCHIVE_KEY)
    archive_claim=freeze.get('source_archive_claim')
    if not _test_only_small_source and not isinstance(archive_claim,dict):
        raise ValueError('42,111 species frozen source lacks original archive claim')
    assert_equal('original archive SHA256',status.get('original_archive_sha256'),
                 pin.get('full_bundle_sha256'),
                 bulk.get('verified_original_dated_archive',{}).get('source_archive_sha256'),
                 ((archive_claim or {}).get('source_archive_sha256')
                  if not _test_only_small_source else pin.get('full_bundle_sha256')))
    assert_equal('source tables SHA256',json.dumps(status.get('all_four_extracted_table_sha256'),sort_keys=True),
                 json.dumps(expected_hashes,sort_keys=True),
                 json.dumps(bulk.get('metadata_snapshots_sha256'),sort_keys=True),
                 json.dumps(freeze.get('source_input_sha256'),sort_keys=True))
    assert_equal('original 42,111 taxon opportunity frame SHA256',
                 bulk.get('original_species_frame_sha256'),
                 freeze.get('source_frame_sha256'),
                 FRAME_SHA256 if not _test_only_small_source else bulk.get('original_species_frame_sha256'))
    assert_equal('all licensed photo IDs',status.get('n_all_available_licensed_original_photo_ids'),
                 bulk.get('n_target_taxon_licensed_photos'),
                 freeze.get('all_available_licensed_photo_ids'))
    if source_counts.is_file() is False:
        raise ValueError('The complete 42,111-species photo count ledger was not uploaded')
    with source_counts.open('r',encoding='utf-8',newline='') as f:
        reader=csv.DictReader(f)
        expected_columns={'inat_taxon_id','species','n_unique_target_photo_ids','source_status','source_manifest_sha256'}
        if set(reader.fieldnames or ())!=expected_columns:
            raise ValueError('Source taxon count ledger has unexpected columns')
        ledger={}
        n_ledger_photos=0
        for row in reader:
            ident=int(row['inat_taxon_id'])
            if ident<1 or ident in ledger or row['source_status']!='complete':
                raise ValueError('Duplicate/incomplete taxon in source count ledger')
            n=int(row['n_unique_target_photo_ids'])
            digest=row['source_manifest_sha256']
            if n<0 or not HEX.fullmatch(str(digest)):
                raise ValueError('Invalid source photo count/hash')
            ledger[ident]=(row['species'],n,digest)
            n_ledger_photos+=n
    if len(ledger)!=count_expected:
        raise ValueError('Not all nominated taxa are in the photo count ledger')
    assert_equal('count ledger photo denominator',n_ledger_photos,
                 freeze.get('all_available_licensed_photo_ids'))
    # Validate all independently materialized checkpoint/manifests and detect
    # omissions/duplicates, without ever looking at historic colour labels.
    source=package_root/'complete_42111_metadata_source'
    recovered=set()
    for checkpoint in source.glob('shard_*_of_*/species/*/checkpoint.json'):
        ck=read_doc(checkpoint)
        tax=ck.get('taxon_id')
        if tax not in ledger or tax in recovered:
            raise ValueError('Unexpected/duplicate species checkpoint in source package')
        species,n,digest=ledger[tax]
        manifest=checkpoint.parent/'photo_manifest.ndjson'
        if (ck.get('status')!='complete' or ck.get('photo_cap','INVALID') is not None
                or ck.get('colour_pixels_opened') is not False
                or ck.get('species')!=species or ck.get('source_photo_rows')!=n
                or ck.get('metadata_sha256')!=digest
                or sha256_file(manifest)!=digest):
            raise ValueError('Complete taxon metadata source did not match photo-count ledger')
        recovered.add(tax)
    if recovered!=set(ledger):
        raise ValueError(f'Missing {len(set(ledger)-recovered)} nominated source checkpoints')
    if n_ledger_photos<1:
        raise ValueError('No photographic opportunity: cannot start original-photo measurement')
    if output.exists() and any(output.iterdir()):
        raise ValueError('New source admission output must be empty: no silent overwrite')
    output.mkdir(parents=True,exist_ok=True)
    plan=plan_jobs(firewall,output/'sealed_terminal_batches',output/'plan',
                   batch_size=batch_size,_test_only_small_source=_test_only_small_source)
    assert_equal('photo work units cover source without cap',plan['source_photo_ids'],n_ledger_photos)
    if plan['n_completed_work_units'] or plan['n_pending_work_units']!=plan['work_units']:
        raise ValueError('Source-only admission must not reuse photo measurement outcomes')
    result={
        'schema':'fcp_complete_photo_source_admission_v1',
        'status':'FULL_SOURCE_PHOTO_WORK_PLAN_READY_NO_PIXELS_YET',
        'source_species':count_expected,
        'source_photo_ids':plan['source_photo_ids'],
        'n_species_checkpoints_verified':len(recovered),
        'worker_partitions':packet_expected,
        'work_units':plan['work_units'],
        'pending_work_units':plan['n_pending_work_units'],
        'compute_batch_size_not_biological_cap':batch_size,
        'source_firewall_sha256':freeze['blind_firewall_sha256'],
        'original_archive_sha256':pin['full_bundle_sha256'],
        'source_status_receipt_sha256':sha256_file(source_status),
        'four_original_table_sha256':expected_hashes,
        'single_coordinated_request_cadence_seconds':1.1,
        'one_request_per_photo_minimum_hours_at_assumed_single_global_cadence':
            round(max(0,n_ledger_photos-1)*1.1/3600,3),
        'provider_bulk_access_and_durable_compute_still_required':True,
        'new_original_photo_colours_measured':0,
        'w_h2_environment_inference_complete':False,
        'synthetic_test_only':_test_only_small_source,
    }
    (output/'source_admission_receipt.json').write_text(
        json.dumps(result,indent=2,ensure_ascii=False,sort_keys=True)+'\n',encoding='utf-8')
    return result


def main()->int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-package-root',type=Path,required=True,
                   help='Root of extracted complete source + 128-partition firewall archive')
    p.add_argument('--source-status-receipt',type=Path,required=True,
                   help='full_source_execution_receipt.json from separate status artifact')
    p.add_argument('--original-pin',type=Path,required=True,
                   help='fcp_source_archive_pin.json from separate source status artifact')
    p.add_argument('--source-counts',type=Path,required=True,
                   help='all_42111_species_source_photo_counts.csv from status artifact')
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--batch-size',type=int,default=64)
    args=p.parse_args()
    result=admit_source(args.source_package_root,args.source_status_receipt,
                        args.original_pin,args.source_counts,args.out,
                        batch_size=args.batch_size)
    print(json.dumps(result,ensure_ascii=False,sort_keys=True))
    return 0

if __name__=='__main__':
    raise SystemExit(main())

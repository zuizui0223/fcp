#!/usr/bin/env python3
"""Seal source-only FCP census AFTER all 128 immutable SQLite source shards.

All relevant original provider metadata are fetched/joined once beforehand.
This script never downloads a photo, opens a pixel, calculates flower colour,
or silently takes a subset of the original 42,111 nominated taxa.

For production use durable scratch with enough capacity for the staged source
database, all 128 incremental metadata shards and the final sealed artifact.
The packaged source intentionally EXCLUDES the large temporary SQLite DB so
downstream workers can consume the smaller complete manifests and firewall.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import tarfile
from pathlib import Path

from fcp_pipeline.fcp_allphotos import sha256_file
from fcp_pipeline.fcp_fresh_study import freeze_fresh
from scripts.acquisition.emit_fcp_staged_source_shards import finalize_shards
from scripts.acquisition.execute_fcp_42111_census import (
    FRAME_SHA256,FROZEN_DATE,ARCHIVE_KEY,sha256_stream,
)

STAGED='SOURCE_SQLITE_STAGED_FOR_INDEPENDENT_SHARDS_NOT_FULL_CENSUS'
FINAL='COMPLETE_42111_SOURCE_METADATA_FROZEN_PIXEL_MEASUREMENT_STILL_NOT_DONE'


def _read(path:Path)->dict:
    return json.loads(path.read_text(encoding='utf-8'))


def finalize_staged_work(work:Path,frame:Path,*,species:int=42111,shards:int=128,
                         workers:int=128,skip_package:bool=False)->dict:
    """Prove all photo IDs were emitted exactly once before sealing source."""
    if species<1 or shards<1 or workers<1:
        raise ValueError('Invalid source taxon/shard/worker counts')
    if species==42111 and sha256_file(frame)!=FRAME_SHA256:
        raise ValueError('Original 42111 species frame changed before source finalization')
    status_path=work/'full_source_execution_receipt.json'
    record=_read(status_path)
    if (record.get('status')!=STAGED
            or record.get('complete_metadata_census') is not False
            or record.get('complete_original_images_measured') is not False
            or record.get('old_measured_flower_colours_used') is not False
            or record.get('target_frame_species')!=species
            or record.get('shards')!=shards
            or record.get('worker_packets')!=workers
            or record.get('immutable_source_date')!=FROZEN_DATE
            or record.get('immutable_source_key')!=ARCHIVE_KEY):
        raise ValueError('Cannot finalize incomplete, changed, or previous-colour source stage')
    source=work/'complete_42111_metadata_source'
    stage=_read(source/'fcp_staged_sqlite_source_receipt.json')
    if (stage.get('sqlite_sha256')!=record.get('staged_source_sqlite_sha256')
            or stage.get('source_input_sha256')!=record.get('all_four_extracted_table_sha256')
            or stage.get('verified_original_dated_archive',{}).get('source_archive_sha256')!=
                record.get('original_archive_sha256')
            or stage.get('source_species_frame_sha256')!=sha256_file(frame)):
        raise ValueError('Staged SQLite source does not match original four-table archive')
    pin=_read(work/'frozen_4_metadata_tables/fcp_source_archive_pin.json')
    if (pin.get('full_bundle_sha256')!=record.get('original_archive_sha256')
            or pin.get('metadata_files_sha256')!=stage['source_input_sha256']):
        raise ValueError('Original full archive/metadata SHA pins changed')
    census=finalize_shards(frame,source,expected=species,shards=shards)
    if census['n_target_taxon_licensed_photos']!=stage['n_target_taxon_licensed_photos']:
        raise ValueError('Disjoint source-photo partitions lost a photo ID')

    counts=work/'all_42111_species_source_photo_counts.csv'
    if counts.exists():
        raise ValueError('Refuse overwrite of already existing source-count ledger')
    temporary=counts.with_suffix('.csv.tmp')
    actual_taxa=0
    with temporary.open('w',encoding='utf-8',newline='') as out:
        w=csv.writer(out)
        w.writerow(['inat_taxon_id','species','n_unique_target_photo_ids',
                    'source_status','source_manifest_sha256'])
        for checkpoint in source.glob('shard_*_of_*/species/*/checkpoint.json'):
            ck=_read(checkpoint)
            if ck.get('status')!='complete' or ck.get('photo_cap','INVALID') is not None:
                raise ValueError('Incomplete source checkpoint while writing photo counts')
            actual_taxa+=1
            w.writerow([ck['taxon_id'],ck['species'],ck['source_photo_rows'],
                        ck['status'],ck['metadata_sha256']])
    if actual_taxa!=species:
        temporary.unlink(missing_ok=True)
        raise ValueError('Source counts do not cover every nominated taxon')
    os.replace(temporary,counts)

    firewall=work/'frozen_location_blind_all_photo_firewall'
    if firewall.exists():
        raise ValueError('A previous frozen firewall requires explicit reconciliation')
    frozen=freeze_fresh(frame,source,firewall,
                        shards=shards,expected=species,workers=workers)
    if (frozen.get('status')!='FRESH_BLIND_METADATA_FROZEN'
            or frozen.get('source_species')!=species
            or frozen.get('all_available_licensed_photo_ids')!=
                census['n_target_taxon_licensed_photos']):
        raise ValueError('Original-photo source-blind firewall did not freeze completely')

    package_file=work/'fcp_complete_42111_source_and_blind_firewall.tar.gz'
    extras={}
    if not skip_package:
        if package_file.exists():
            raise ValueError('Cannot overwrite previously sealed source package')
        temp=package_file.with_suffix('.tar.gz.tmp')
        # SQLite is kept separately on durable disk for possible source audit.
        # Do not put it into the portable manifests/worker artifact.
        with tarfile.open(temp,'w:gz',compresslevel=3) as t:
            for f in sorted(source.glob('shard_*_of_*')):
                if not f.is_dir(): continue
                t.add(f,arcname=f'{source.name}/{f.name}',recursive=True)
            t.add(source/'bulk_source_receipt.json',
                  arcname=f'{source.name}/bulk_source_receipt.json')
            t.add(firewall,arcname=firewall.name,recursive=True)
        os.replace(temp,package_file)
        extras={
            'packed_all_source_manifests_sha256':sha256_stream(package_file),
            'packed_all_source_manifests_bytes':package_file.stat().st_size,
            'packed_all_source_manifests_path':str(package_file)
        }
    record.update(
        status=FINAL,
        complete_metadata_census=True,
        no_per_species_photo_cap=True,
        full_source_census_receipt=census,
        species_count_checkpoints_audited=actual_taxa,
        n_all_available_licensed_original_photo_ids=
            frozen['all_available_licensed_photo_ids'],
        frozen_location_blind_worker_partitions=frozen['worker_packets'],
        source_shards_sealed_from_one_original_sqlite=shards,
        original_staged_sqlite_excluded_from_portable_source_archive=True,
        full_pixel_measurement_complete=False,
        biological_h2_and_environment_inference_ready=False,
        **extras
    )
    temporary_status=status_path.with_suffix('.json.tmp')
    temporary_status.write_text(json.dumps(
        record,ensure_ascii=False,sort_keys=True,indent=2)+'\n',encoding='utf-8')
    os.replace(temporary_status,status_path)
    return record


def main(argv=None)->int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--frame',required=True,type=Path)
    p.add_argument('--work',required=True,type=Path)
    p.add_argument('--expected-species',type=int,default=42111)
    p.add_argument('--shards',type=int,default=128)
    p.add_argument('--worker-partitions',type=int,default=128)
    p.add_argument('--skip-package',action='store_true')
    a=p.parse_args(argv)
    d=finalize_staged_work(a.work,a.frame,species=a.expected_species,
                          shards=a.shards,workers=a.worker_partitions,
                          skip_package=a.skip_package)
    print('FCP_STAGED_ALL_SOURCE_FINALIZED',json.dumps({
        k:d[k] for k in ('status','complete_metadata_census',
                         'species_count_checkpoints_audited',
                         'n_all_available_licensed_original_photo_ids')},
        sort_keys=True))
    return 0


if __name__=='__main__':
    raise SystemExit(main())

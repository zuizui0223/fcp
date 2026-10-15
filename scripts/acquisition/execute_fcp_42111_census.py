#!/usr/bin/env python3
"""One-shot, source-frozen FCP complete iNaturalist metadata census (NO pixels).

Runs a user-selected, exact dated public Open Data tarball through four-table
SHA verification, full 42,111-target taxon census and fresh, location-blind
photo-manifest freeze. Never imports historical colour labels. This entry point
may require >80 GB of scratch disk and hours of computing time. It does not
produce W, infer polymorphism, or download original photograph pixels.

On failure, a fail-closed JSON receipt records the stage and complete-source
status; no partial metadata result is published as a frozen source census.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ARCHIVE_KEY='metadata/inaturalist-open-data-20260927.tar.gz'
ARCHIVE_URL='https://inaturalist-open-data.s3.amazonaws.com/'+ARCHIVE_KEY
ARCHIVE_SIZE=35_966_549_880
FRAME_SHA256='5d871bd1f1190c68fd513bb527e6c93de1606dc35ff8cfff5f6187850dc382cc'
REPRODUCED_ARCHIVE_SHA256='bfc1a97c436bf5989fccd43f12392953128489b76b92cea24da335ad73da4d8e'
FROZEN_DATE='2026-09-27'


def sha256_stream(path:Path)->str:
    d=hashlib.sha256()
    with path.open('rb') as h:
        for part in iter(lambda:h.read(8*1024*1024),b''):
            d.update(part)
    return d.hexdigest()


def disk_free(path:Path)->int:
    return shutil.disk_usage(path).free


def run(argv:list[str],*,cwd:Path,stage:str,status:dict):
    status['stage']=stage
    status['stage_started_utc']=datetime.now(timezone.utc).isoformat()
    print('FCP_42111_FULL_SOURCE_STAGE',stage,flush=True)
    env=os.environ.copy()
    env['PYTHONPATH']=str(cwd)+(os.pathsep+env['PYTHONPATH'] if env.get('PYTHONPATH') else '')
    subprocess.run(argv,cwd=cwd,check=True,env=env)
    status['last_completed_stage']=stage


def download_exact_archive(path:Path,status:dict,min_space:int):
    path.parent.mkdir(parents=True,exist_ok=True)
    free=disk_free(path.parent)
    status['initial_scratch_disk_free_bytes']=free
    if free<min_space:
        raise ValueError(f'Not enough persistent scratch disk for exact archive and four extracts: free={free} required={min_space}')
    req=urllib.request.Request(ARCHIVE_URL,method='HEAD')
    with urllib.request.urlopen(req,timeout=35) as h:
        remote=int(h.headers.get('Content-Length','-1'))
    if remote!=ARCHIVE_SIZE:
        raise ValueError('Provider dated archive object length changed; original fixed-date source is not admitted')
    # Stream to a separate partially named file. curl resumes only a bounded,
    # immutable date-selected object; end-to-end SHA/length is always verified.
    partial=path.with_suffix('.tar.gz.part')
    status['stage']='download_exact_dated_archive'
    subprocess.run(['curl','--fail','--silent','--show-error','--location','--retry','4','--retry-delay','10',
                    '--retry-all-errors','--connect-timeout','60',
                    '--continue-at','-','--output',str(partial),ARCHIVE_URL],check=True)
    if partial.stat().st_size!=ARCHIVE_SIZE:
        raise ValueError('Incomplete iNaturalist dated archive transfer; refuse all photo inference')
    partial.rename(path)
    return path


def execute(frame:Path,work:Path,*,source_archive:Path|None=None,
            species:int=42111,shards:int=128,workers:int=128,
            skip_package:bool=False,required_scratch_gb:int=80,
            release_verified_temporary_tables:bool=False,
            stage_only:bool=False)->dict:
    if workers<1 or shards<1 or species<1:
        raise ValueError('Shards, workers and species must be positive')
    if work.exists() and any(work.iterdir()):
        raise ValueError('Fresh FCP work output must be empty to prevent mixing science inputs')
    work.mkdir(parents=True,exist_ok=True)
    repo=Path(__file__).resolve().parents[2]
    record={'schema':'fcp_full_42111_dated_source_execution_v1',
            'status':'FAILED_OR_UNFINISHED_SOURCE_METADATA_NOT_SCIENTIFIC_RESULTS',
            'complete_metadata_census':False,'complete_original_images_measured':False,
            'new_W_or_environment_results_computed':False,
            'old_measured_flower_colours_used':False,
            'immutable_source_date':FROZEN_DATE,'immutable_source_key':ARCHIVE_KEY,
            'target_frame_species':species,'shards':shards,'worker_packets':workers,
            'started_utc':datetime.now(timezone.utc).isoformat()}
    receipt=work/'full_source_execution_receipt.json'
    def store():
        receipt.write_text(json.dumps(record,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    store()
    try:
        if species==42111 and sha256_stream(frame)!=FRAME_SHA256:
            raise ValueError('Full census species opportunity frame SHA256 mismatched')
        with frame.open('r',encoding='utf-8',newline='') as handle:
            reader=csv.DictReader(handle)
            if not {'inat_taxon_id','species'}.issubset(reader.fieldnames or ()):
                raise ValueError('Complete frozen species frame lacks required taxon fields')
            actual_species=sum(1 for _ in reader)
        if actual_species!=species:
            raise ValueError(f'Frozen opportunity frame count {actual_species} differs from declared {species}')
        original=source_archive
        if original is None:
            original=download_exact_archive(work/'source'/'inaturalist-open-data-20260927.tar.gz',record,
                                            required_scratch_gb*1_000_000_000)
            record['archive_downloaded_from_verified_dated_public_s3']=True
        else:
            if not original.is_file():
                raise FileNotFoundError(original)
            record['archive_downloaded_from_verified_dated_public_s3']=False
        record['stage']='hash_complete_original_archive'
        observed=sha256_stream(original)
        if species==42111 and observed!=REPRODUCED_ARCHIVE_SHA256:
            raise ValueError('2026-09-27 source archive differs from previously observed verified 36GB SHA')
        record['original_archive_sha256']=observed
        record['original_archive_bytes']=original.stat().st_size
        if source_archive is None and record['original_archive_bytes']!=ARCHIVE_SIZE:
            raise ValueError('Dated Open Data archive changed length after transfer')
        # SHA is recorded before any four-table data are opened. Date and
        # source key are fixed by the predeclared scientific protocol.
        (work/'source_archive_sha256_before_source_data.txt').write_text(observed+'\n',encoding='utf-8')
        store()
        frozen=work/'frozen_4_metadata_tables'
        run([sys.executable,'scripts/acquisition/pin_inat_open_data_snapshot.py',
             '--archive',str(original),'--s3-key',ARCHIVE_KEY,
             '--archive-sha256',observed,'--outdir',str(frozen)],
             cwd=repo,stage='extract_all_four_and_pin_exact_SHA256',status=record)
        if source_archive is None:
            original.unlink()  # Verified original digest retained; frees ~36GB.
            record['source_archive_deleted_after_full_SHA_verification']=True
        pin=frozen/'fcp_source_archive_pin.json'
        metadata=json.loads(pin.read_text(encoding='utf-8'))
        h=metadata['metadata_files_sha256']
        assert set(h)=={'taxa','observations','observers','photos'}
        record['all_four_extracted_table_sha256']=h
        store()
        dest=work/'complete_42111_metadata_source'
        cmd=[sys.executable,'scripts/acquisition/import_inat_open_data_snapshot.py',
             '--species-frame',str(frame),'--snapshot-date',FROZEN_DATE,
             '--taxa',str(frozen/'taxa.csv.gz'),
             '--observations',str(frozen/'observations.csv.gz'),
             '--observers',str(frozen/'observers.csv.gz'),
             '--photos',str(frozen/'photos.csv.gz'),
             '--snapshot-pin',str(pin),
             '--taxa-sha256',h['taxa'],'--observations-sha256',h['observations'],
             '--observers-sha256',h['observers'],'--photos-sha256',h['photos'],
             '--expected-species',str(species),'--shards',str(shards),
             '--output-root',str(dest)]
        if release_verified_temporary_tables:
            cmd.append('--prune-verified-temporary-tables')
        if stage_only:
            cmd.append('--stage-sqlite-only')
        run(cmd,cwd=repo,stage=('scan_source_once_and_stage_sqlite_for_sharded_emission'
                               if stage_only else
                               'scan_entire_dated_source_and_emit_all_photo_manifests'),status=record)
        if stage_only:
            staged=json.loads((dest/'fcp_staged_sqlite_source_receipt.json').read_text(encoding='utf-8'))
            if (staged.get('status')!='SOURCE_SQLITE_STAGED_SHARDS_NOT_YET_EMITTED'
                    or staged.get('source_species')!=species
                    or staged.get('source_input_sha256')!=h
                    or staged.get('n_target_taxon_licensed_photos',-1)<0):
                raise ValueError('Provider snapshot did not yield one complete immutable SQLite stage')
            record.update(
                status='SOURCE_SQLITE_STAGED_FOR_INDEPENDENT_SHARDS_NOT_FULL_CENSUS',
                complete_metadata_census=False,
                staged_original_source_photo_ids=staged['n_target_taxon_licensed_photos'],
                staged_source_sqlite_sha256=staged['sqlite_sha256'],
                staged_sqlite_bytes=staged['sqlite_bytes'],
                n_shards_still_to_emit=shards,
                full_pixel_measurement_complete=False,
                biological_h2_and_environment_inference_ready=False)
            store()
            return record
        census=json.loads((dest/'bulk_source_receipt.json').read_text(encoding='utf-8'))
        if census['status']!='COMPLETE_SNAPSHOT_METADATA_ALL_NOMINATED_TAXA' or census['taxa_in_initial_frame']!=species:
            raise ValueError('Not every target species was source-censused')
        record['full_source_census_receipt']=census
        record['pruned_verified_scratch_tables_after_scan']=census.get(
            'pruned_pinned_temporary_metadata_after_full_scan',False)
        record['pruned_verified_scratch_table_bytes']=census.get('pruned_temporary_table_bytes',0)
        store()
        from csv import writer
        counts=work/'all_42111_species_source_photo_counts.csv'
        with counts.open('w',encoding='utf-8',newline='') as hdl:
            csvwriter=writer(hdl)
            csvwriter.writerow(['inat_taxon_id','species','n_unique_target_photo_ids','source_status','source_manifest_sha256'])
            count=0
            for ck in dest.glob('shard_*_of_*/species/*/checkpoint.json'):
                row=json.loads(ck.read_text(encoding='utf-8'))
                assert row['status']=='complete' and row['photo_cap'] is None
                count+=1
                csvwriter.writerow([row['taxon_id'],row['species'],row['source_photo_rows'],
                                    row['status'],row['metadata_sha256']])
        if count!=species:
            raise ValueError(f'Only {count} complete per-species photo source manifests')
        record['species_count_checkpoints_audited']=count
        store()
        firewall=work/'frozen_location_blind_all_photo_firewall'
        run([sys.executable,'-m','fcp_pipeline.fcp_fresh_study','freeze',
             '--source-species',str(frame),'--metadata-root',str(dest),
             '--shards',str(shards),'--worker-partitions',str(workers),
             '--expected-species',str(species),'--outdir',str(firewall)],
             cwd=repo,stage='seal_every_real_photo_ID_without_opening_a_pixel',status=record)
        fz=json.loads((firewall/'fresh_study_freeze.json').read_text())
        assert fz['status']=='FRESH_BLIND_METADATA_FROZEN'
        assert fz['source_species']==species
        assert fz['all_available_licensed_photo_ids']==census['n_target_taxon_licensed_photos']
        record.update(status='COMPLETE_42111_SOURCE_METADATA_FROZEN_PIXEL_MEASUREMENT_STILL_NOT_DONE',
                      complete_metadata_census=True,
                      no_per_species_photo_cap=True,
                      n_all_available_licensed_original_photo_ids=fz['all_available_licensed_photo_ids'],
                      frozen_location_blind_worker_partitions=fz['worker_packets'],
                      full_pixel_measurement_complete=False,
                      biological_h2_and_environment_inference_ready=False)
        store()
        if not skip_package:
            # IMPORTANT: GitHub ephemeral disk is not a durable source store;
            # package the exact manifests + blind packets for artifact retention.
            archive=work/'fcp_complete_42111_source_and_blind_firewall.tar.gz'
            run(['tar','-C',str(work),'-czf',str(archive),
                 dest.name,firewall.name],cwd=repo,
                stage='archive_verbatim_complete_source_manifests_and_firewall',status=record)
            record['packed_all_source_manifests_sha256']=sha256_stream(archive)
            record['packed_all_source_manifests_bytes']=archive.stat().st_size
            record['packed_all_source_manifests_path']=str(archive)
            store()
        return record
    except BaseException as exc:
        record['error_type']=type(exc).__name__
        record['error_message']=str(exc)[:750]
        record['status']='FAILED_OR_UNFINISHED_SOURCE_METADATA_NOT_SCIENTIFIC_RESULTS'
        store()
        raise


def main(argv=None)->int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--frame',required=True,type=Path)
    ap.add_argument('--work',required=True,type=Path)
    ap.add_argument('--archive',type=Path,default=None)
    ap.add_argument('--expected-species',type=int,default=42111)
    ap.add_argument('--shards',type=int,default=128)
    ap.add_argument('--worker-partitions',type=int,default=128)
    ap.add_argument('--skip-package',action='store_true')
    ap.add_argument('--stage-only',action='store_true',
                    help='Stage full source once on durable storage, then emit photo manifests by shard')
    ap.add_argument('--minimum-free-disk-gb',type=int,default=80)
    ap.add_argument('--release-verified-temporary-tables',action='store_true',
                    help='Reclaim only SHA-verified 4-table scratch copies after full source loading')
    args=ap.parse_args(argv)
    out=execute(args.frame.resolve(),args.work.resolve(),source_archive=args.archive.resolve() if args.archive else None,
                species=args.expected_species,shards=args.shards,workers=args.worker_partitions,
                skip_package=args.skip_package,required_scratch_gb=args.minimum_free_disk_gb,
                release_verified_temporary_tables=args.release_verified_temporary_tables,
                stage_only=args.stage_only)
    print('FCP_ALL_42111_SOURCE_EXECUTION_RESULT',json.dumps(out,ensure_ascii=False,sort_keys=True))
    return 0

if __name__=='__main__':
    raise SystemExit(main())

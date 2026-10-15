#!/usr/bin/env python3
"""Incrementally emit all FCP source manifests from ONE durable, verified SQLite import.

Usage: python scripts/acquisition/emit_fcp_staged_source_shards.py \
 --frame /persistent/fcp/42111_species.csv \
 --staged-metadata /persistent/fcp/complete_42111_metadata_source \
 --shard-start 0 --shard-stop 8 --expected-species 42111 --shards 128

Each disjoint source shard has independent, immutable complete checkpoints,
including zero-photo taxa. The SQL photo universe was globally de-duplicated
once when imported; emitters NEVER re-download the 36GB archive, rescan raw
monthly TSVs, open image pixels, or read historical flower colour outcomes.

A failed shard is NOT resumable without explicit inspection/quarantine.
Only finalize_source_shards() after all expected shard receipts independently
verify makes a SOURCE METADATA census, not a photo measurement or W result.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import shutil
import sqlite3
from pathlib import Path

from fcp_pipeline.fcp_allphotos import opportunity_frame, opportunity_schedule, sha256_file
from scripts.acquisition.import_inat_open_data_snapshot import (
    FIELDS, SNAPSHOT_SOURCE, AUTHORITATIVE_42111_SHA256, integer,
)

STAGE_SCHEMA='fcp_staged_sqlite_photo_universe_v1'
STAGE_STATUS='SOURCE_SQLITE_STAGED_SHARDS_NOT_YET_EMITTED'
SHARD_STATUS='SOURCE_SHARD_COMPLETE_FROM_FROZEN_SQLITE_NO_IMAGES'
TABLE_KEYS={'taxa','observations','observers','photos'}


def _read(path:Path)->dict:
    obj=json.loads(path.read_text(encoding='utf-8'))
    if not isinstance(obj,dict):
        raise ValueError('FCP staged source receipt must be a JSON object')
    return obj


def _validate_stage(stage_root:Path,frame:Path,*,expected:int,shards:int,
                    verify_sqlite:bool=True):
    receipt=_read(stage_root/'fcp_staged_sqlite_source_receipt.json')
    if (receipt.get('schema')!=STAGE_SCHEMA
            or receipt.get('status')!=STAGE_STATUS
            or receipt.get('source_kind')!=SNAPSHOT_SOURCE
            or receipt.get('full_source_census_complete') is not False
            or receipt.get('full_study_original_images_measured') is not False
            or receipt.get('historical_flower_colour_csv_used') is not False
            or receipt.get('source_image_pixels_opened') is not False
            or receipt.get('original_photo_ids_deduplicated_globally') is not True
            or receipt.get('source_species')!=expected
            or receipt.get('shards')!=shards
            or receipt.get('source_species_frame_sha256')!=sha256_file(frame)
            or set(receipt.get('source_input_sha256',{}))!=TABLE_KEYS
            or (expected==42111 and
                receipt.get('source_species_frame_sha256')!=AUTHORITATIVE_42111_SHA256)):
        raise ValueError('Incomplete or changed immutable global FCP source stage')
    if not receipt.get('verified_original_dated_archive',{}).get('source_archive_sha256'):
        raise ValueError('Staged monthly source lacks verified original dated archive')
    name=receipt.get('sqlite_basename')
    if name!='_working_fcp_snapshot.sqlite':
        raise ValueError('Unexpected staged source SQLite basename')
    db_file=stage_root/name
    if not db_file.is_file() or db_file.stat().st_size!=receipt.get('sqlite_bytes'):
        raise ValueError('Staged source SQLite absent/changed size')
    if verify_sqlite and sha256_file(db_file)!=receipt.get('sqlite_sha256'):
        raise ValueError('Staged source SQLite changed after original provider ingestion')
    frame_data=opportunity_frame(frame,expected)
    schedule=opportunity_schedule(frame_data,shards)
    source=sqlite3.connect(db_file.as_uri()+'?mode=ro&immutable=1',uri=True)
    source.execute('PRAGMA query_only=ON')
    integrity=source.execute('PRAGMA quick_check').fetchone()
    if integrity is None or integrity[0]!='ok':
        source.close()
        raise ValueError('Staged source SQLite is not structurally intact')
    n=source.execute('SELECT COUNT(*) FROM photos WHERE ambiguous=0').fetchone()[0]
    if n!=receipt.get('n_target_taxon_licensed_photos'):
        source.close()
        raise ValueError('Original photo ID count changed within staged SQLite')
    return receipt,schedule,source


def _photo_query()->str:
    return """SELECT observations.uuid,observations.observer,
          observations.observed_on,observations.lat,observations.lon,observations.accuracy,
          photos.photo_id,photos.license,photos.extension,photos.photographer,
          observer_labels.login,observer_labels.name
       FROM photos INDEXED BY photos_taxon_numeric_order
       JOIN observations ON photos.uuid=observations.uuid
       LEFT JOIN observer_labels ON photos.photographer=observer_labels.observer_id
       WHERE photos.ambiguous=0 AND photos.taxon=?
       ORDER BY CAST(photos.photo_id AS INTEGER)"""


def _write_shard(stage_root:Path,frame_sha:str,doc:dict,db:sqlite3.Connection,
                 names:dict[int,str],taxa:list[int],part:int,shards:int)->dict:
    folder=stage_root/f'shard_{part:03}_of_{shards:03}'
    if folder.exists():
        receipt=folder/'shard_source_receipt.json'
        if not receipt.is_file():
            raise ValueError(f'Unsealed interrupted source shard {part} requires explicit quarantine')
        completed=_read(receipt)
        if (completed.get('status')!=SHARD_STATUS
                or completed.get('stage_sqlite_sha256')!=doc['sqlite_sha256']
                or completed.get('source_species_frame_sha256')!=frame_sha
                or completed.get('partition')!=part
                or completed.get('n_source_species')!=len(taxa)
                or completed.get('source_species_ids')!=taxa):
            raise ValueError(f'Previously emitted source shard {part} belongs to another source')
        count=0
        for tid in taxa:
            checkpoint=folder/'species'/str(tid)/'checkpoint.json'
            record=_read(checkpoint)
            manifest=checkpoint.parent/'photo_manifest.ndjson'
            if (record.get('status')!='complete' or record.get('taxon_id')!=tid
                    or record.get('metadata_sha256')!=sha256_file(manifest)
                    or record.get('source_photo_rows',-1)<0):
                raise ValueError(f'Existing source shard {part} has a corrupted checkpoint')
            count+=record['source_photo_rows']
        if count!=completed['n_original_source_photo_ids']:
            raise ValueError(f'Existing source shard {part} no longer matches source photo counts')
        return completed

    folder.mkdir(parents=True,exist_ok=False)
    src_claim=doc['verified_original_dated_archive']
    contract={
        'schema':'fcp_uncapped_shard_source_contract_v1',
        'source_species_frame_sha256':frame_sha,
        'expected_source_species':doc['source_species'],
        'n_source_species_in_shard':len(taxa),
        'shard':part,'nshards':shards,
        'created_d2':doc['source_snapshot_date'],
        'source_kind':SNAPSHOT_SOURCE,
        'created_d2_field_semantics':'immutable_source_snapshot_date_not_API_created_d2_query',
        'input_metadata_sha256':doc['source_input_sha256'],
        **src_claim,
        'photo_cap':None,'outcome_blind':True
    }
    (folder/'source_contract.json').write_text(json.dumps(contract,indent=2)+'\n',encoding='utf-8')
    query=_photo_query()
    plan=[str(t[3]) for t in db.execute('EXPLAIN QUERY PLAN '+query,(taxa[0],))] if taxa else []
    if any('USE TEMP B-TREE FOR ORDER BY' in p.upper() for p in plan):
        raise RuntimeError('Partitioned photo census would create an unbounded SQLite sort')
    if taxa and not any('photos_taxon_numeric_order' in p for p in plan):
        raise RuntimeError('Frozen photo source covering index not used for species')
    n_photo=0
    attributable=unresolved=0
    taxa_with_photos=0
    for tid in taxa:
        species_folder=folder/'species'/str(tid)
        species_folder.mkdir(parents=True,exist_ok=False)
        manifest=species_folder/'photo_manifest.ndjson'
        local=0
        with manifest.open('x',encoding='utf-8') as out:
            for (uuid,observer,observed_on,lat,lon,accuracy,pid,lic,ext,photographer,
                 login,person_name) in db.execute(query,(tid,)):
                author=integer(photographer)
                person=str(person_name or '').strip() or str(login or '').strip()
                if person:
                    attributable+=1
                else:
                    unresolved+=1
                    person=(f'iNaturalist photo creator {author}' if author
                            else 'iNaturalist photographer not resolved')
                profile=(f'https://www.inaturalist.org/users/{author}' if author
                         else 'https://www.inaturalist.org/')
                attribution=(f'{person} ({profile}), CC0 public domain' if lic=='cc0'
                             else f'© {person} ({profile}), {lic.upper()}')
                obj=dict(zip(FIELDS,(tid,names[tid],str(uuid),str(pid),
                    f'https://inaturalist-open-data.s3.amazonaws.com/photos/{pid}/large.{ext}',
                    lic,attribution,observer,str(observed_on or ''),lat,lon,accuracy,
                    doc['source_snapshot_date'],SNAPSHOT_SOURCE)))
                out.write(json.dumps(obj,sort_keys=True,ensure_ascii=False)+'\n')
                local+=1
                n_photo+=1
                if n_photo % 25000==0 and shutil.disk_usage(folder).free<6_000_000_000:
                    raise RuntimeError(
                        f'FCP source shard disk-reserve exhausted at shard={part}, photo={n_photo}; '
                        'preserve all earlier sealed shards and move to durable storage')
        taxa_with_photos+=bool(local)
        checkpoint={
            'schema':'fcp_uncapped_metadata_census_v1',
            'taxon_id':tid,'species':names[tid],
            'created_d2':doc['source_snapshot_date'],
            'source_kind':SNAPSHOT_SOURCE,
            'snapshot_source_sha256':doc['source_input_sha256'],
            'cursor':0,'batches':0,'status':'complete','page_size':None,
            'photo_cap':None,'colour_pixels_opened':False,
            'source_photo_rows':local,'metadata_sha256':sha256_file(manifest)
        }
        (species_folder/'checkpoint.json').write_text(
            json.dumps(checkpoint,indent=2)+'\n',encoding='utf-8')
    complete={
        'schema':'fcp_frozen_source_shard_emission_v1',
        'status':SHARD_STATUS,'stage_sqlite_sha256':doc['sqlite_sha256'],
        'source_species_frame_sha256':frame_sha,
        'partition':part,'n_shards':shards,
        'source_species_ids':taxa,
        'n_source_species':len(taxa),
        'n_original_source_photo_ids':n_photo,
        'n_taxa_with_source_photos':taxa_with_photos,
        'n_taxa_without_source_photos':len(taxa)-taxa_with_photos,
        'n_photos_with_named_creator':attributable,
        'n_photos_creator_unresolved':unresolved,
        'historical_colour_outcomes_used':False,
        'photo_pixels_opened':False,
        'per_species_photo_cap':None
    }
    candidate=folder/'shard_source_receipt.json.tmp'
    candidate.write_text(json.dumps(complete,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    os.replace(candidate,folder/'shard_source_receipt.json')
    print('FCP_SOURCE_PARTITION_SEALED',json.dumps({
        'partition':part,'taxa':len(taxa),'photos':n_photo}),flush=True)
    return complete


def emit_shards(frame:Path,stage_root:Path,*,start:int,stop:int,
                expected:int=42111,shards:int=128)->dict:
    """Run disjoint source emitters one partition block at a time on durable disk."""
    if not 0<=start<stop<=shards:
        raise ValueError('Source shard start/stop outside full immutable source universe')
    stage,schedule,db=_validate_stage(stage_root,frame,expected=expected,shards=shards)
    try:
        names=dict(zip(schedule.inat_taxon_id.astype(int),schedule.species.astype(str)))
        frame_sha=stage['source_species_frame_sha256']
        receipts=[]
        for i in range(start,stop):
            taxa=schedule.loc[schedule.partition.eq(i),'inat_taxon_id'].astype(int).tolist()
            receipts.append(_write_shard(stage_root,frame_sha,stage,db,names,taxa,i,shards))
        return {
            'status':'INCREMENTAL_SOURCE_SHARDS_SEALED_NOT_GLOBAL_CENSUS',
            'source_species':expected,
            'all_shards_expected':shards,
            'completed_here':len(receipts),
            'partition_start':start,
            'partition_stop_exclusive':stop,
            'newly_or_previously_sealed_photo_ids_in_range':
                sum(x['n_original_source_photo_ids'] for x in receipts),
            'source_sha256':stage['sqlite_sha256'],
            'all_photos_measured':False
        }
    finally:
        db.close()


def finalize_shards(frame:Path,stage_root:Path,*,expected:int=42111,
                    shards:int=128)->dict:
    """Only all complete taxon/photo receipts can become a bulk source receipt."""
    stage,schedule,db=_validate_stage(stage_root,frame,expected=expected,shards=shards)
    try:
        source_frame_sha=stage['source_species_frame_sha256']
        total=attributed=unresolved=with_photos=without_photos=0
        sealed_taxa=set()
        for part in range(shards):
            folder=stage_root/f'shard_{part:03}_of_{shards:03}'
            complete=_read(folder/'shard_source_receipt.json')
            expected_taxa=schedule.loc[schedule.partition.eq(part),'inat_taxon_id'].astype(int).tolist()
            if (complete.get('status')!=SHARD_STATUS
                    or complete.get('stage_sqlite_sha256')!=stage['sqlite_sha256']
                    or complete.get('source_species_frame_sha256')!=source_frame_sha
                    or complete.get('partition')!=part
                    or complete.get('source_species_ids')!=expected_taxa
                    or complete.get('n_source_species')!=len(expected_taxa)
                    or complete.get('historical_colour_outcomes_used') is not False
                    or complete.get('photo_pixels_opened') is not False
                    or complete.get('per_species_photo_cap','INVALID') is not None):
                raise ValueError(f'Source partition {part} not source-complete')
            count=0
            n_with=0
            for tid in expected_taxa:
                if tid in sealed_taxa:
                    raise ValueError('Same frozen taxon occurs in more than one source partition')
                sealed_taxa.add(tid)
                check=_read(folder/'species'/str(tid)/'checkpoint.json')
                manifest=folder/'species'/str(tid)/'photo_manifest.ndjson'
                if (check.get('taxon_id')!=tid
                        or check.get('status')!='complete'
                        or check.get('photo_cap','INVALID') is not None
                        or check.get('colour_pixels_opened') is not False
                        or check.get('snapshot_source_sha256')!=stage['source_input_sha256']
                        or check.get('metadata_sha256')!=sha256_file(manifest)):
                    raise ValueError(f'Source photo manifest altered/incomplete for taxon {tid}')
                n=check.get('source_photo_rows')
                if not isinstance(n,int) or n<0:
                    raise ValueError('Source checkpoint photo count invalid')
                count+=n
                n_with+=int(n>0)
            if (count!=complete['n_original_source_photo_ids']
                    or n_with!=complete['n_taxa_with_source_photos']
                    or len(expected_taxa)-n_with!=complete['n_taxa_without_source_photos']):
                raise ValueError(f'Source partition {part} photo/taxon counts do not match checkpoint ledger')
            total+=count
            with_photos+=n_with
            without_photos+=len(expected_taxa)-n_with
            attributed+=complete['n_photos_with_named_creator']
            unresolved+=complete['n_photos_creator_unresolved']
        if (len(sealed_taxa)!=expected or total!=stage['n_target_taxon_licensed_photos']
                or attributed+unresolved!=total or with_photos+without_photos!=expected):
            raise ValueError('All original photo/taxon denominators are not preserved across source shards')
        receipt={
            'schema':'fcp_open_data_bulk_complete_source_v1',
            'status':'COMPLETE_SNAPSHOT_METADATA_ALL_NOMINATED_TAXA',
            'source_kind':SNAPSHOT_SOURCE,
            'source_snapshot_date':stage['source_snapshot_date'],
            'outcome_blind':True,'old_photo_measurements_read':False,
            'image_pixels_opened':False,
            'original_species_frame_sha256':source_frame_sha,
            'taxa_in_initial_frame':expected,
            'n_source_observations_scanned':stage['n_source_observations_scanned'],
            'n_target_taxon_research_grade_observations':
                stage['n_target_taxon_research_grade_observations'],
            'n_source_observers_scanned':stage['n_source_observers_scanned'],
            'n_source_photos_scanned':stage['n_source_photos_scanned'],
            'n_target_taxon_licensed_photo_metadata_rows':
                stage['n_target_taxon_licensed_photo_metadata_rows'],
            'n_target_taxon_licensed_photos':total,
            'n_photo_manifest_rows_written':total,
            'manifest_indexed_stream_without_sqlite_order_temp_btree':True,
            'n_conflicting_photo_ids_quarantined':
                stage['n_conflicting_photo_ids_quarantined'],
            'n_target_photos_with_named_creator':attributed,
            'n_target_photos_creator_unresolved':unresolved,
            'n_taxa_with_source_photos':with_photos,
            'n_taxa_without_source_photos':without_photos,
            'metadata_snapshots_sha256':stage['source_input_sha256'],
            'verified_original_dated_archive':
                stage['verified_original_dated_archive'],
            'photo_cap':None,
            'all_available_snapshot_photo_ids_included':True,
            'source_has_immutable_monthly_cutoff':True,
            'not_literally_live_API_created_d2_window':True,
            'unavailable_or_unlicensed_photo_not_called_polymorphic':True,
            'original_S3_extension_case_preserved':True,
            'conflicting_reused_photo_IDs_quarantined_without_morph_selection':True,
            'pruned_pinned_temporary_metadata_after_full_scan':
                stage['pruned_pinned_temporary_metadata_after_full_scan'],
            'pruned_temporary_table_bytes':stage['pruned_temporary_table_bytes'],
            'source_stage_sqlite_sha256':stage['sqlite_sha256'],
            'number_independently_sealed_source_shards':shards
        }
        output=stage_root/'bulk_source_receipt.json'
        if output.exists():
            if _read(output)!=receipt:
                raise ValueError('Finalized immutable source receipt changed')
        else:
            tmp=output.with_suffix('.json.tmp')
            tmp.write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n',encoding='utf-8')
            os.replace(tmp,output)
        return receipt
    finally:
        db.close()


def main(argv=None)->int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--frame',required=True,type=Path)
    p.add_argument('--staged-metadata',required=True,type=Path)
    p.add_argument('--expected-species',type=int,default=42111)
    p.add_argument('--shards',type=int,default=128)
    g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--shard-start',type=int)
    g.add_argument('--finalize',action='store_true')
    p.add_argument('--shard-stop',type=int)
    args=p.parse_args(argv)
    if args.finalize:
        if args.shard_stop is not None:
            p.error('--shard-stop cannot accompany --finalize')
        result=finalize_shards(args.frame,args.staged_metadata,
                               expected=args.expected_species,shards=args.shards)
    else:
        if args.shard_stop is None:
            p.error('--shard-stop is required with --shard-start')
        result=emit_shards(args.frame,args.staged_metadata,
                           start=args.shard_start,stop=args.shard_stop,
                           expected=args.expected_species,shards=args.shards)
    print(json.dumps(result,ensure_ascii=False,sort_keys=True))
    return 0

if __name__=='__main__':
    raise SystemExit(main())

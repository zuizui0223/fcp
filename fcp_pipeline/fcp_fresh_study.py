"""New FCP study from ZERO colour measurements, with a frozen source pipeline.

Historical photo-measurement CSVs are not inputs. The historical algorithm / ROI
weights are method specifications ONLY, never biological observations here.

Two independent gates:
1. finish the entire frozen outcome-blind taxon/photo metadata census, then
   produce blind worker packets containing *all* photo IDs, including IDs that
   may have appeared in historical studies (they must be measured anew).
2. reassemble sealed fresh ROI-v4 outputs and refuse missing/extra/duplicate
   measurements or algorithm fingerprints. An incomplete run is not a dataset.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd

from .fcp_allphotos import opportunity_frame, opportunity_schedule, sha256_file

STUDY = 'fcp_fresh_uncapped_20261014_v1'
SOURCE_METHOD_COMMIT = '9fae6ccdf684a46026f72ba12e98de2c5c54bf2a'
SOURCE_MEASUREMENT_BLOB = 'b1d7091db44682559bbe4b48a9f3b8db648501a6'
SOURCE_DETECTOR_SHA256 = 'f1aaeec4664fe2c178e5cf2bc1f508977bef3e4aa7b40613026cb8ae3de789d5'
SOURCE_ENCODER_SHA256 = '84ed466ffcc5c1f8d08409bc34a23bb364ab2c15e402cb12d4335a42be0e0951'
SOURCE_DECODER_SHA256 = 'a62f8fa5ea080447c0689418d69e58f1e83e0b7adf9c142e2bd9bcc8045c0b11'
MORPHS=('white','yellow_orange','red_pink','blue_purple')
BIO=('white','yellow','orange','red','pink','magenta','purple','blue','bronze')
WORKER=('measurement_id','image_filename','photo_license')
ANNOTATION=('inat_taxon_id','species','observation_id','photo_id','photo_license',
            'observer_id','observed_on','latitude','longitude','positional_accuracy_m','attribution')


def _read_json(path:Path)->dict:
    return json.loads(path.read_text(encoding='utf-8'))


def _builder():
    path=Path(__file__).resolve().parents[1]/'scripts/acquisition/build_fcp_uncapped_blind_batches.py'
    if not path.is_file():
        raise ValueError('Original-interface blind builder missing')
    spec=importlib.util.spec_from_file_location('_fresh_blind_builder',path)
    if spec is None or spec.loader is None:
        raise ValueError('Cannot load fresh blind builder')
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def now_compact_date(value) -> str:
    """Convert the source-dated ISO cutoff, never an invented observation filter."""
    from datetime import date
    d=date.fromisoformat(str(value))
    return d.strftime('%Y%m%d')


def required_source_dirs(frame_path:Path, metadata_root:Path,*,shards:int,expected:int)->tuple[list[Path],dict]:
    """A complete source census means EVERY canonical species, not a chosen subset."""
    if shards<1:
        raise ValueError('shards must be >=1')
    frame=opportunity_frame(frame_path,expected)
    plan=opportunity_schedule(frame,shards)
    frame_sha=sha256_file(frame_path)
    date=None
    source_kind=None
    source_metadata_hashes=None
    source_archive_claim=None
    dirs=[]
    pending=[]
    for i in range(shards):
        folder=metadata_root/f'shard_{i:03}_of_{shards:03}'
        contract_file=folder/'source_contract.json'
        if not contract_file.is_file():
            pending.extend(plan.loc[plan.partition.eq(i),'inat_taxon_id'].astype(int).tolist())
            continue
        contract=_read_json(contract_file)
        if (contract.get('source_species_frame_sha256')!=frame_sha or
                contract.get('shard')!=i or contract.get('nshards')!=shards or
                contract.get('expected_source_species')!=expected or
                contract.get('n_source_species_in_shard')!=int(plan.partition.eq(i).sum()) or
                contract.get('photo_cap','INVALID') is not None or
                contract.get('outcome_blind') is not True):
            raise ValueError(f'Source shard {i} violates fresh frozen census')
        kind=contract.get('source_kind','inat_api_observations_cursor')
        if kind not in ('inat_api_observations_cursor','inaturalist_open_dataset_monthly_snapshot'):
            raise ValueError(f'Unknown photo opportunity provider for shard {i}')
        if source_kind is None:
            source_kind=kind
        elif kind!=source_kind:
            raise ValueError('Mixed API and monthly Open Data source contracts are not a single frozen census')
        if kind=='inaturalist_open_dataset_monthly_snapshot':
            hashes=contract.get('input_metadata_sha256')
            if (not isinstance(hashes,dict) or set(hashes)!={'taxa','observations','observers','photos'} or
                    any(not isinstance(v,str) or len(v)!=64 or
                        any(c not in '0123456789abcdef' for c in v) for v in hashes.values())):
                raise ValueError('Snapshot taxon, observation and photo byte hashes required')
            if contract.get('created_d2_field_semantics')!='immutable_source_snapshot_date_not_API_created_d2_query':
                raise ValueError('Monthly bulk source date is not API created_d2: missing distinction')
            if source_metadata_hashes is None:
                source_metadata_hashes=hashes
            elif hashes!=source_metadata_hashes:
                raise ValueError('Inconsistent source snapshot input hashes between shards')
            if expected == 42111:
                claim={key:contract.get(key) for key in
                       ('source_archive_key','source_archive_sha256','source_archive_pin_sha256')}
                if (claim['source_archive_key']!=
                        f"metadata/inaturalist-open-data-{now_compact_date(contract.get('created_d2'))}.tar.gz"
                        or any(not isinstance(claim[k],str) or len(claim[k])!=64 or
                               any(c not in '0123456789abcdef' for c in claim[k])
                               for k in ('source_archive_sha256','source_archive_pin_sha256'))):
                    raise ValueError('Full FCP monthly source lacks dated archive SHA pin')
                if source_archive_claim is None:
                    source_archive_claim=claim
                elif claim!=source_archive_claim:
                    raise ValueError('Different monthly source archive pins across shards')
        now=str(contract.get('created_d2'))
        if date is not None and now!=date:
            raise ValueError('Source metadata date changed between shards')
        date=now
        expected_ids=set(plan.loc[plan.partition.eq(i),'inat_taxon_id'].astype(int))
        observed_ids={int(p.name) for p in (folder/'species').glob('[0-9]*') if p.is_dir()} if (folder/'species').exists() else set()
        if observed_ids-expected_ids:
            raise ValueError(f'Shard {i} contains source taxon IDs not in its partition')
        for tax,name in plan.loc[plan.partition.eq(i),['inat_taxon_id','species']].itertuples(index=False,name=None):
            d=folder/'species'/str(int(tax))
            check=d/'checkpoint.json'
            if not check.exists():
                pending.append(int(tax));continue
            record=_read_json(check)
            if (record.get('taxon_id')!=int(tax) or record.get('species')!=str(name) or
                record.get('created_d2')!=date or
                record.get('photo_cap','INVALID') is not None or
                record.get('colour_pixels_opened') is not False):
                raise ValueError(f'Wrong metadata identity/outcome firewall for {tax}')
            if record.get('status')!='complete':
                pending.append(int(tax));continue
            manifest=d/'photo_manifest.ndjson'
            if not manifest.is_file():
                raise ValueError(f'Completed taxon {tax} has no photo manifest')
            if kind=='inaturalist_open_dataset_monthly_snapshot':
                if (record.get('source_kind')!=kind or
                        record.get('snapshot_source_sha256')!=source_metadata_hashes or
                        record.get('metadata_sha256')!=sha256_file(manifest)):
                    raise ValueError(f'Completed bulk snapshot manifest identity/hash mismatch: {tax}')
            dirs.append(d)
    audit={'status':'SOURCE_CENSUS_COMPLETE' if not pending else 'INCOMPLETE_METADATA_NOT_AUTHORIZED',
           'source_species':len(frame),'completed_species':len(dirs),
           'incomplete_species':len(pending),'sample_missing_taxa':pending[:30],
           'source_frame_sha256':frame_sha,'frozen_created_d2':date,
           'source_kind':source_kind,'source_input_sha256':source_metadata_hashes,
           'source_archive_claim':source_archive_claim,
           'source_metadata_only':True,'biological_images_opened':False,
           'historical_measured_photos_used':False,'no_photo_number_cap':True}
    return dirs,audit


def freeze_fresh(frame_path:Path, metadata_root:Path, output:Path, *, shards:int=128, expected:int=42111,
                 workers:int=128)->dict:
    dirs,audit=required_source_dirs(frame_path,metadata_root,shards=shards,expected=expected)
    if audit['status']!='SOURCE_CENSUS_COMPLETE':
        # Fail BEFORE producing worker packets or biological classification inputs.
        raise ValueError(f'Cannot open any pixels before all source metadata is complete: {audit}')
    if output.exists() and any(output.iterdir()):
        raise ValueError('Fresh study output directory must be empty; no accidental overwrite')
    builder=_builder()
    result=builder.freeze_inputs(dirs,output,partitions=workers,already_measured=None)
    if (result['n_original_photo_ids_preloaded']!=0 or
            result['n_new_metadata_photo_ids_already_in_original']!=0 or
            result['source_species_manifests']!=expected or
            result['photo_cap'] is not None):
        raise ValueError('Fresh study unexpectedly reused historical photo measurements')
    protected=['firewall_receipt.json','metadata_join_key.csv','acquisition_key.csv']
    for part in range(workers):
        protected += [f'worker_partition_{part:03}.csv',f'acquisition_partition_{part:03}.csv']
    frozen_files={p:sha256_file(output/p) for p in protected}
    source={**audit,'study_id':STUDY,'status':'FRESH_BLIND_METADATA_FROZEN',
            'frozen_file_sha256':frozen_files,
            'worker_packets':workers,'all_available_licensed_photo_ids':result['eligible_total_photo_ids'],
            'new_measurements_completed':0,'historical_measured_CSV_read':False,
            'frozen_roi_method_commit':SOURCE_METHOD_COMMIT,
            'frozen_roi_measurement_blob':SOURCE_MEASUREMENT_BLOB,
            'reused_old_colour_labels':False,
            'measurement_status':'NOT_STARTED_NO_IMAGE_PIXELS_OPENED',
            'blind_firewall_sha256':sha256_file(output/'firewall_receipt.json')}
    (output/'fresh_study_freeze.json').write_text(json.dumps(source,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    return source


def validate_method_attestation(manifest_path:Path)->dict:
    doc=_read_json(manifest_path)
    required={'source_commit':SOURCE_METHOD_COMMIT,
              'measurement_git_blob':SOURCE_MEASUREMENT_BLOB,
              'detector_sha256':SOURCE_DETECTOR_SHA256,
              'efficient_sam_encoder_sha256':SOURCE_ENCODER_SHA256,
              'efficient_sam_decoder_sha256':SOURCE_DECODER_SHA256,
              'location_blind':True,'no_old_measurement_csv_inputs':True,
              'persisted_image_pixels':False}
    for k,val in required.items():
        if doc.get(k)!=val:
            raise ValueError(f'Non-identical or non-blind FCP measurement runtime: {k}')
    return doc


def reassemble_fresh(firewall:Path, sealed_dir:Path, attestation:Path, outdir:Path)->dict:
    """One result per newly remeasured source photo; never backfill old colours."""
    freeze=_read_json(firewall/'fresh_study_freeze.json')
    prior=_read_json(firewall/'firewall_receipt.json')
    if freeze['status']!='FRESH_BLIND_METADATA_FROZEN' or freeze['study_id']!=STUDY:
        raise ValueError('Only fully frozen fresh study inputs may be assembled')
    if prior.get('n_original_photo_ids_preloaded',-1)!=0 or prior.get('n_new_metadata_photo_ids_already_in_original',-1)!=0:
        raise ValueError('Historic photo measurements appear in supposedly fresh study')
    if freeze['blind_firewall_sha256']!=sha256_file(firewall/'firewall_receipt.json'):
        raise ValueError('Source photo firewall changed after freeze')
    for relative,digest in freeze['frozen_file_sha256'].items():
        if sha256_file(firewall/relative)!=digest:
            raise ValueError(f'Fresh photo source/worker bytes drifted after freeze: {relative}')
    att=validate_method_attestation(attestation)
    worker_parts=sorted(firewall.glob('worker_partition_*.csv'))
    if len(worker_parts)!=int(freeze['worker_packets']):
        raise ValueError('Missing blind worker partition')
    terminal_parts=sorted(sealed_dir.glob('terminal_*.csv'))
    if not terminal_parts:
        raise ValueError('No new ROI-v4 terminal measurements provided')
    outdir.mkdir(parents=True,exist_ok=True)
    db_path=outdir/'fresh_assembly.sqlite'
    if db_path.exists() or (outdir/'newly_remeasured_photos.csv').exists():
        raise ValueError('Cannot overwrite fresh reassembly/duplicate this run')
    db=sqlite3.connect(db_path)
    complete=False
    try:
        db.execute('CREATE TABLE exp (measurement_id TEXT PRIMARY KEY, photo_id TEXT UNIQUE, taxon INTEGER, species TEXT, observation TEXT, observer TEXT, license TEXT, observed_on TEXT, latitude TEXT, longitude TEXT, accuracy TEXT, worker_seen INTEGER DEFAULT 0, attribution TEXT, worker_part INTEGER, worker_index INTEGER, worker_filename TEXT)')
        db.execute('CREATE TABLE obs (measurement_id TEXT PRIMARY KEY, json TEXT NOT NULL)')
        with (firewall/'metadata_join_key.csv').open('r',encoding='utf-8',newline='') as f:
            reader=csv.DictReader(f)
            if not set(ANNOTATION)|{'measurement_id'} <= set(reader.fieldnames or []):
                raise ValueError('Fresh metadata join key incomplete')
            for v in reader:
                mid=v['measurement_id']
                db.execute('INSERT INTO exp(measurement_id,photo_id,taxon,species,observation,observer,license,observed_on,latitude,longitude,accuracy,attribution) VALUES (?,?,?,?,?,?,?,?,?,?,?,?)',
                        (mid,v['photo_id'],int(v['inat_taxon_id']),v['species'],v['observation_id'],v['observer_id'],v['photo_license'],v['observed_on'],v['latitude'],v['longitude'],v.get('positional_accuracy_m',''),v['attribution']))
        expected_count=db.execute('SELECT COUNT(*) FROM exp').fetchone()[0]
        if expected_count!=int(freeze['all_available_licensed_photo_ids']):
            raise ValueError('All-photo source denominator does not match metadata join identity')
        worker_count_by_part={}
        for part,path in enumerate(worker_parts):
            n_part=0
            with path.open('r',encoding='utf-8',newline='') as f:
                reader=csv.DictReader(f)
                if tuple(reader.fieldnames or ())!=WORKER:
                    raise ValueError('Worker sees non-blind or changed metadata columns')
                for worker_index,row in enumerate(reader):
                    cur=db.execute('UPDATE exp SET worker_seen=worker_seen+1,worker_part=?,worker_index=?,worker_filename=? WHERE measurement_id=? AND license=?',
                                   (part,worker_index,row['image_filename'],row['measurement_id'],row['photo_license']))
                    if cur.rowcount!=1:
                        raise ValueError('Worker entry unknown or licence inconsistent with source')
                    n_part+=1
            worker_count_by_part[part]=n_part
        if db.execute('SELECT COUNT(*) FROM exp WHERE worker_seen != 1').fetchone()[0]!=0:
            raise ValueError('Worker packet coverage duplicates/missing')
        expected_terminal={'measurement_id','morph','measurement_status'}
        new_cols={f'flower_fraction_{c}' for c in BIO}
        unit_indices=set()
        expected_unit_counts={part:0 for part in worker_count_by_part}
        units_by_part={part:[] for part in worker_count_by_part}
        for path in terminal_parts:
            receipt_path=path.with_suffix('.json')
            if not receipt_path.is_file():
                raise ValueError('Fresh terminal batch lacks durable method/row receipt')
            receipt=_read_json(receipt_path)
            if (receipt.get('status')!='COMPLETE_WORK_UNIT' or
                    receipt.get('terminal_csv_sha256')!=sha256_file(path) or
                    receipt.get('image_bytes_persisted') is not False or
                    receipt.get('reused_historical_measurement_csv') is not False or
                    receipt.get('full_study_complete') is not False or
                    receipt.get('method')!=att):
                raise ValueError('Fresh terminal work unit lacks frozen method or photo checksum')
            part,batch=receipt.get('partition'),receipt.get('batch')
            if (type(part) is not int or part not in worker_count_by_part or
                    type(batch) is not int or batch<0 or (part,batch) in unit_indices):
                raise ValueError('Fresh work unit partition/batch duplicated or invalid')
            unit_indices.add((part,batch))
            batch_size=receipt.get('compute_batch_size_not_photo_cap')
            claimed_batches=receipt.get('number_batches_for_this_partition')
            if type(batch_size) is not int or batch_size<1 or type(claimed_batches) is not int:
                raise ValueError('Invalid work-unit batching contract')
            from math import ceil
            if claimed_batches!=ceil(worker_count_by_part[part]/batch_size):
                raise ValueError('Work-unit batch count inconsistent with photo census')
            low=batch*batch_size
            high=min((batch+1)*batch_size,worker_count_by_part[part])
            if high<=low:
                raise ValueError('Incomplete fresh ROI-v4 work unit: unexpected empty batch')
            expected_source_rows=list(db.execute(
                'SELECT measurement_id,worker_filename,license FROM exp WHERE worker_part=? AND worker_index>=? AND worker_index<? ORDER BY worker_index',
                (part,low,high)))
            if len(expected_source_rows)!=high-low:
                raise ValueError('Frozen source worker sequence did not match expected batch')
            batch_hash=hashlib.sha256(
                pd.DataFrame(expected_source_rows,columns=list(WORKER)).to_csv(index=False).encode()).hexdigest()
            if receipt.get('selected_source_worker_sha256')!=batch_hash:
                raise ValueError('Terminal batch worker subset SHA256 does not match frozen source order')
            if receipt.get('n_source_photo_ids')!=len(expected_source_rows):
                raise ValueError('Incomplete fresh ROI-v4 terminal batch relative to frozen source')
            expected_unit_photo_ids={row[0] for row in expected_source_rows}
            units_by_part[part].append((batch,batch_size))
            seen_in_file=0
            with path.open('r',encoding='utf-8',newline='') as f:
                reader=csv.DictReader(f)
                if expected_terminal-set(reader.fieldnames or ()) or new_cols-set(reader.fieldnames or ()):
                    raise ValueError('Missing original ROI-v4 terminal morphology / nine palette outputs')
                if {'species','observer_id','photo_id','latitude','longitude','photo_url_large'} & set(reader.fieldnames or ()):
                    raise ValueError('ROI-v4 terminal results leaked identity or geography')
                for row in reader:
                    if row['measurement_id'] not in expected_unit_photo_ids:
                        raise ValueError('Terminal photo measurement does not belong to its frozen source batch')
                    status=row['measurement_status']
                    is_class=status=='classified_four_state_morph'
                    if is_class!=(row['morph'] in MORPHS):
                        raise ValueError('Photo classification status inconsistent with frozen colour states')
                    if not is_class and row['morph']!='mixed_uncertain':
                        raise ValueError('Unclassifiable image must remain terminal mixed_uncertain')
                    if is_class:
                        vec=np.array([float(row[f'flower_fraction_{c}']) for c in BIO],float)
                        if not np.isfinite(vec).all() or (vec<-1e-10).any() or abs(float(vec.sum())-1)>1e-7:
                            raise ValueError('Classifiable photo has invalid original nine-colour palette')
                    db.execute('INSERT INTO obs(measurement_id,json) VALUES (?,?)',
                               (row['measurement_id'],json.dumps(row,ensure_ascii=False)))
                    seen_in_file+=1
            if seen_in_file!=receipt.get('n_terminal_photo_ids') or seen_in_file!=receipt.get('n_source_photo_ids'):
                raise ValueError('Incomplete fresh ROI-v4 terminal batch relative to source manifest')
            expected_unit_counts[part]+=seen_in_file
        for part,n in worker_count_by_part.items():
            sizes={size for _,size in units_by_part[part]}
            from math import ceil
            target_batch_count=(ceil(n/next(iter(sizes))) if len(sizes)==1 else 0)
            if (expected_unit_counts[part]!=n or
                    len(sizes)!=(1 if n else 0) or
                    set(batch for batch,_ in units_by_part[part])!=set(range(target_batch_count))):
                raise ValueError('Incomplete fresh ROI-v4 partitions/batches')
        n_terminal=db.execute('SELECT COUNT(*) FROM obs').fetchone()[0]
        if n_terminal!=expected_count:
            raise ValueError(f'Incomplete fresh ROI-v4 measurement, terminal={n_terminal} expected={expected_count}')
        if db.execute('SELECT COUNT(*) FROM obs LEFT JOIN exp USING(measurement_id) WHERE exp.measurement_id IS NULL').fetchone()[0]:
            raise ValueError('Measurement not present in frozen photo opportunities')
        dest=outdir/'newly_remeasured_photos.csv'
        count=0;classified=0
        with dest.open('w',encoding='utf-8',newline='') as f:
            header=['inat_taxon_id','species','photo_id','observation_id','observer_id','observed_on','latitude','longitude','positional_accuracy_m','photo_license','attribution','morph','measurement_status','global_classifiable',*[f'flower_fraction_{c}' for c in BIO]]
            writer=csv.DictWriter(f,fieldnames=header);writer.writeheader()
            for tax,sp,pid,obs,ob,day,lat,lon,accuracy,lic,attribution,raw in db.execute('SELECT exp.taxon,exp.species,exp.photo_id,exp.observation,exp.observer,exp.observed_on,exp.latitude,exp.longitude,exp.accuracy,exp.license,exp.attribution,obs.json FROM exp JOIN obs USING(measurement_id) ORDER BY exp.taxon,CAST(exp.photo_id AS INTEGER)'):
                measured=json.loads(raw)
                usable=measured['measurement_status']=='classified_four_state_morph'
                writer.writerow({'inat_taxon_id':tax,'species':sp,'photo_id':pid,'observation_id':obs,'observer_id':ob,
                                 'observed_on':day,'latitude':lat,'longitude':lon,'positional_accuracy_m':accuracy,'photo_license':lic,'attribution':attribution,
                                 'morph':measured['morph'],'measurement_status':measured['measurement_status'],
                                 'global_classifiable':str(usable).lower(),
                                 **{f'flower_fraction_{c}':measured[f'flower_fraction_{c}'] for c in BIO}})
                count+=1;classified+=int(usable)
        output={'schema':'fcp_fresh_complete_remeasurement_v1','status':'COMPLETE_FRESH_MEASURED_PHOTO_CENSUS',
                'study_id':STUDY,'n_frozen_licensed_photo_ids':expected_count,
                'n_terminal_measurements':count,'n_classifiable':classified,
                'n_unclassifiable':count-classified,'n_source_taxa':freeze['source_species'],
                'method_attestation_sha256':sha256_file(attestation),
                'source_firewall_sha256':freeze['blind_firewall_sha256'],
                'new_photo_table_sha256':sha256_file(dest),
                'historical_measurements_reused':False,'same_photo_IDs_as_earlier_may_have_been_remeasured':True,
                'not_independent_source_species_confirmatory_replication':True,
                'frozen_H2_null_computed':False,'environment_caused_FCP':False}
        (outdir/'fresh_complete_receipt.json').write_text(json.dumps(output,indent=2)+'\n',encoding='utf-8')
        complete=True
        return output
    finally:
        db.close()
        # Only a COMPLETE receipt makes the reassembly deliverable. An incomplete
        # partial CSV is deleted on error; source terminal partitions are untouched.
        if not complete:
            for p in (outdir/'newly_remeasured_photos.csv',outdir/'fresh_complete_receipt.json'):
                if p.exists():p.unlink()
        if db_path.exists():db_path.unlink()


def main(argv:list[str]|None=None)->int:
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest='command',required=True)
    a=sub.add_parser('freeze',help='Freeze only complete fresh source-photo census')
    a.add_argument('--source-species',required=True,type=Path)
    a.add_argument('--metadata-root',required=True,type=Path)
    a.add_argument('--outdir',required=True,type=Path)
    a.add_argument('--shards',type=int,default=128)
    a.add_argument('--worker-partitions',type=int,default=128)
    a.add_argument('--expected-species',type=int,default=42111)
    b=sub.add_parser('reassemble',help='No reused photo rows: all original ROI-v4 terminal outcomes required')
    b.add_argument('--firewall',required=True,type=Path)
    b.add_argument('--sealed-results',required=True,type=Path)
    b.add_argument('--method-attestation',required=True,type=Path)
    b.add_argument('--outdir',required=True,type=Path)
    opt=p.parse_args(argv)
    if opt.command=='freeze':
        result=freeze_fresh(opt.source_species,opt.metadata_root,opt.outdir,
                            shards=opt.shards,workers=opt.worker_partitions,expected=opt.expected_species)
    else:
        result=reassemble_fresh(opt.firewall,opt.sealed_results,opt.method_attestation,opt.outdir)
    print(json.dumps({k:v for k,v in result.items() if k in ('status','source_species','all_available_licensed_photo_ids','n_terminal_measurements','n_classifiable')},ensure_ascii=False))
    return 0


if __name__=='__main__':
    raise SystemExit(main())
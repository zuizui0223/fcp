#!/usr/bin/env python3
"""ONE-OFF real-iNaturalist fresh ROI-v4 runtime pilot (not a 42,111-taxon study).

Every selected photograph is a frozen source ID from a prior provider metadata
artifact; selection is species-equal and hash based, never outcome based. Photos
are downloaded with a polite delay, measured by the EXACT original ROI-v4 model
in a three-column blind worker process, and discarded after the subprocess.
Outputs are explicitly quarantined from all fresh H2/W, null, soil and manuscript
inference. An acquisition failure remains a terminal missing observation.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd
from PIL import Image

from fcp_pipeline.fcp_allphotos import sha256_file
from fcp_pipeline.fcp_fresh_study import (
    SOURCE_METHOD_COMMIT, SOURCE_MEASUREMENT_BLOB, BIO,
    SOURCE_ENCODER_SHA256, SOURCE_DECODER_SHA256, SOURCE_DETECTOR_SHA256,
)
from scripts.acquisition.run_fcp_fresh_image_partition import verify_frozen_method

SOURCE_FRAME_SHA = '5d871bd1f1190c68fd513bb527e6c93de1606dc35ff8cfff5f6187850dc382cc'
SOURCE_ARTIFACT = 11928883989
SOURCE_RUN = 38639092021
SOURCE_CREATED_DATE = '2026-10-14'
ORIGINAL_RIBES_MANIFEST_SHA = '78c7e533c19a9099dd1f2de96dc7864ee65253abfe17dbfa2c536252a1665403'
SALT = 'FCP_SEPARATE_REAL_PIXEL_TECHNICAL_PILOT_20261015_V1'
WORKER_FIELDS = ('measurement_id', 'image_filename', 'photo_license')
IMAGE_URL = re.compile(r'^/photos/([0-9]+)/large\.(jpg|jpeg|png)
MAX_BYTES = 12 * 1024 * 1024
MAX_PER_SPECIES = 24


def selected_metadata(metadata_root: Path, n_each: int) -> tuple[list[dict], dict]:
    if not 1 <= n_each <= MAX_PER_SPECIES:
        raise ValueError('Real-photo tech pilot per-species bound must be 1..24')
    contract = json.loads((metadata_root / 'source_contract.json').read_text(encoding='utf-8'))
    if (contract.get('source_species_frame_sha256') != SOURCE_FRAME_SHA or
            contract.get('created_d2') != SOURCE_CREATED_DATE or
            contract.get('photo_cap', 'INVALID') is not None or
            contract.get('outcome_blind') is not True or
            contract.get('expected_source_species') != 42111):
        raise ValueError('Pilot source must be the verified 42,111 opportunity-frame dated metadata, never historic labelled colours')
    if len(list(metadata_root.glob('species/[0-9]*/checkpoint.json'))) != 9:
        raise ValueError('Expected original nine attempted species; do not accept a changed provider artifact')
    photos=[]
    archives={}
    taxa=[]
    global_ids=set()
    for ck in sorted(metadata_root.glob('species/[0-9]*/checkpoint.json'), key=lambda x:int(x.parent.name)):
        status=json.loads(ck.read_text(encoding='utf-8'))
        tid=int(ck.parent.name)
        if (status.get('taxon_id') != tid or status.get('created_d2') != SOURCE_CREATED_DATE or
                status.get('colour_pixels_opened') is not False or status.get('photo_cap','INVALID') is not None):
            raise ValueError('Unfrozen or already-measured metadata source taxon')
        if status.get('status') != 'complete':
            continue
        manifest=ck.parent/'photo_manifest.ndjson'
        digest=sha256_file(manifest)
        if tid==47129 and digest!=ORIGINAL_RIBES_MANIFEST_SHA:
            raise ValueError('Previously independently verified Ribes source manifest byte identity changed')
        if status.get('metadata_sha256') not in (None, digest):
            raise ValueError('Completed source checkpoint disagrees with its photo manifest')
        species=[]
        seen={}
        conflicting_photo_ids=set()
        with manifest.open('r',encoding='utf-8') as f:
            for line in f:
                obj=json.loads(line)
                pid=str(obj['photo_id'])
                if (int(obj['inat_taxon_id'])!=tid or obj['species']!=status['species'] or
                        not pid.isdigit() or not obj.get('photo_license')):
                    raise ValueError('Original dated photo metadata provenance / license / creator changed')
                if pid in seen:
                    if seen[pid]!=obj:
                        # The same original image appeared in multiple observations.
                        # It is unsuitable for photo-site attribution. Quarantine
                        # its photo ID entirely; never count repeated observations
                        # as independent flowers. This exclusion uses metadata only.
                        conflicting_photo_ids.add(pid)
                    continue
                seen[pid]=obj
                species.append(obj)
        species=[v for v in species if str(v['photo_id']) not in conflicting_photo_ids]
        if not species:
            continue
        # Deterministic, species-equal, blind to image colour and geography.
        species.sort(key=lambda x:hashlib.sha256(f'{SALT}:{tid}:{x["photo_id"]}'.encode()).digest())
        chosen=species[:n_each]
        for obj in chosen:
            pid=str(obj['photo_id'])
            if pid in global_ids:
                raise ValueError('Provider photo ID appears in two purported species')
            global_ids.add(pid)
        photos.extend(chosen)
        archives[str(tid)]={'species':status['species'],'full_manifest_sha256':digest,
                            'complete_source_unique_unambiguous_photo_ids':len(species),
                            'duplicate_conflicting_original_photo_ids_quarantined':len(conflicting_photo_ids),
                            'pilot_selected':len(chosen)}
        taxa.append(tid)
    if len(taxa)!=4 or 47129 not in taxa:
        raise ValueError(f'Prior artifact audit showed exactly four completed taxa, got {taxa}')
    selected_ids=sorted((int(v['inat_taxon_id']),str(v['photo_id'])) for v in photos)
    selected_hash=hashlib.sha256(json.dumps(selected_ids,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    return photos,{'source_artifact_id':SOURCE_ARTIFACT,'source_run_id':SOURCE_RUN,
                   'source_snapshot_date':SOURCE_CREATED_DATE,'initial_frame_sha256':SOURCE_FRAME_SHA,
                   'completed_taxa_available':taxa,'species':archives,'selected_photo_id_pairs_sha256':selected_hash,
                   'selection_outcome_blind':True,'pilot_only_max_photos_per_species':n_each}


def original_photo_bytes(url:str,pid:str)->bytes:
    uri=urlparse(url)
    match=IMAGE_URL.fullmatch(uri.path)
    if (uri.scheme!='https' or uri.netloc!='inaturalist-open-data.s3.amazonaws.com' or
            uri.params or uri.query or uri.fragment or match is None or match.group(1)!=pid):
        raise ValueError('Real photo URL is not the frozen original iNaturalist Open Data S3 object')
    request=Request(url,headers={'User-Agent':'fcp-real-fresh-roi-technical-pilot/1.0 (github.com/zuizui0223/fcp)',
                                 'Accept':'image/*'})
    with urlopen(request,timeout=45) as f:
        raw=f.read(MAX_BYTES+1)
    if not raw or len(raw)>MAX_BYTES:
        raise ValueError('Provider image empty or beyond bounded temporary image memory')
    # Full decode check, not only filename or HTTP status.
    from io import BytesIO
    with Image.open(BytesIO(raw)) as image:
        image.verify()
    return raw


def run_pilot(meta:Path,outdir:Path,source:Path,detector:Path,sam:Path,
              *,n_each:int=24,delay_seconds:float=1.1)->dict:
    if delay_seconds<1.0:
        raise ValueError('At most one original photo request per second')
    if outdir.exists() and any(outdir.iterdir()):
        raise ValueError('One-shot real-photo output must be empty; never silently overwrite pilot')
    chosen, provenance=selected_metadata(meta,n_each)
    method=verify_frozen_method(source,detector,sam)
    outdir.mkdir(parents=True,exist_ok=True)
    failures={}
    hashes={}
    actual_worker=[]
    with tempfile.TemporaryDirectory(prefix='fcp_original_real_photo_pilot_') as scratch:
        temp=Path(scratch)
        images=temp/'images';images.mkdir()
        for i,row in enumerate(chosen):
            if i:time.sleep(delay_seconds)
            mid='FCPR-PILOT-'+hashlib.sha256(f'{SALT}:{row["photo_id"]}'.encode()).hexdigest().upper()[:24]
            name=mid+'.jpg'
            try:
                raw=original_photo_bytes(str(row['photo_url_large']),str(row['photo_id']))
                (images/name).write_bytes(raw)
                hashes[mid]=hashlib.sha256(raw).hexdigest()
                actual_worker.append({'measurement_id':mid,'image_filename':name,'photo_license':row['photo_license']})
            except Exception as exc:
                failures[mid]=f'{type(exc).__name__}:{str(exc)[:190]}'
            row['measurement_id']=mid
            print(f'original_photo_source_progress={i+1}/{len(chosen)} downloaded={len(actual_worker)} failures={len(failures)}',flush=True)
        raw_output=temp/'original_unmodified_roi_v4.csv'
        if actual_worker:
            worker=temp/'blind_worker.csv'
            with worker.open('w',encoding='utf-8',newline='') as f:
                w=csv.DictWriter(f,fieldnames=WORKER_FIELDS);w.writeheader();w.writerows(actual_worker)
            subprocess.run([sys.executable,str(source),'--worker-manifest',str(worker),
                            '--images-dir',str(images),'--detector-weight',str(detector),
                            '--efficient-sam-dir',str(sam),'--output-csv',str(raw_output)],
                           check=True,timeout=2400)
            measured=pd.read_csv(raw_output,dtype={'measurement_id':str}).fillna('')
            if (len(measured)!=len(actual_worker) or measured.measurement_id.duplicated().any() or
                set(measured.measurement_id)!={x['measurement_id'] for x in actual_worker}):
                raise ValueError('Unmodified original ROI-v4 returned incorrect terminal source photo identities')
        else:
            measured=pd.DataFrame(columns=['measurement_id','morph','measurement_status',*[f'flower_fraction_{x}' for x in BIO]])
        merged={r['measurement_id']:r for r in measured.to_dict(orient='records')}
        final=[]
        for source_row in chosen:
            mid=source_row['measurement_id']
            result=merged.get(mid)
            if result is None:
                if mid not in failures:
                    raise ValueError('Source photo omitted by original ROI-v4 evaluator')
                result={'measurement_id':mid,'morph':'mixed_uncertain',
                        'measurement_status':'image_acquisition_failed',
                        'failure_reasons':failures[mid],
                        **{f'flower_fraction_{x}':0.0 for x in BIO}}
            good=result['measurement_status']=='classified_four_state_morph'
            if good!=(result['morph'] in ('white','yellow_orange','red_pink','blue_purple')):
                raise ValueError('Frozen original ROI-v4 morphological terminal status contradiction')
            final.append({**{k:source_row.get(k) for k in (
                 'inat_taxon_id','species','observation_id','photo_id','observer_id',
                 'observed_on','latitude','longitude','positional_accuracy_m',
                 'photo_license','attribution')},
                 'source_url':source_row['photo_url_large'],
                 'technical_pilot_not_full_FCP_study':True,
                 **result,'original_image_bytes_sha256':hashes.get(mid,''),
                 'global_classifiable':bool(good)})
    data=pd.DataFrame(final)
    if len(data)!=len(chosen) or data.photo_id.astype(str).duplicated().any():
        raise ValueError('Pilot terminal photo ledger incomplete or duplicated')
    photo_csv=outdir/'real_original_roi_v4_new_photo_pilot.csv'
    data.to_csv(photo_csv,index=False)
    classified=data.measurement_status.eq('classified_four_state_morph')
    receipt={
        'schema':'fcp_real_original_roi_v4_quarantined_technical_pilot_v1',
        'status':'COMPLETE_REAL_ORIGINAL_PHOTO_TECHNICAL_PILOT_NOT_42111_CENSUS',
        'run_utc':datetime.now(timezone.utc).isoformat(),
        'source':provenance,
        'n_completed_source_species_sampled':len(provenance['completed_taxa_available']),
        'n_real_source_original_photo_ids':len(chosen),
        'n_real_original_image_bytes_downloaded':len(actual_worker),
        'n_real_original_images_processed_by_frozen_roi_v4':len(measured),
        'n_original_provider_download_failures':len(failures),
        'n_real_photo_colour_classified':int(classified.sum()),
        'n_roi_unclassifiable':int((~classified & ~data.measurement_status.eq('image_acquisition_failed')).sum()),
        'per_species_terminal_counts':{str(t):{
           'n_selected_photo_ids':int(len(group)),
           'n_downloaded':int(group.original_image_bytes_sha256.ne('').sum()),
           'n_classifiable':int(group.global_classifiable.sum()),
           'n_image_acquisition_failed':int(group.measurement_status.eq('image_acquisition_failed').sum()),
        } for t,group in data.groupby('inat_taxon_id')},
        'status_counts':{str(k):int(v) for k,v in data.measurement_status.value_counts().items()},
        'raw_source_URLs_only_not_old_photo_colours':True,
        'original_roi_v4_method_hashes':method,
        'real_fresh_photo_terminal_csv_sha256':sha256_file(photo_csv),
        'new_photographs_scored_for_biological_results':False,
        'full_source_census_complete':False,
        'eligible_for_main_fcp_W_null_or_climate_inference':False,
        'n_original_historical_colour_rows_loaded':0,
        'source_image_pixels_persisted':False,
        'selected_photo_cap_applies_only_to_quarantined_technical_pilot':True,
    }
    (outdir/'real_original_roi_v4_new_photo_pilot_receipt.json').write_text(
        json.dumps(receipt,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print('FCP_REAL_PHOTO_PILOT_RESULT '+json.dumps({
        'status':receipt['status'],'n_original_photo_ids':len(chosen),
        'n_real_images_downloaded':len(actual_worker),
        'n_real_images_roi_measured':len(measured),
        'n_classifiable':int(classified.sum()),
        'status_counts':receipt['status_counts'],
        'per_species':receipt['per_species_terminal_counts'],
        'terminal_csv_sha256':receipt['real_fresh_photo_terminal_csv_sha256']},sort_keys=True),flush=True)
    return receipt


def main(argv=None)->int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-metadata-root',type=Path,required=True)
    parser.add_argument('--output-root',type=Path,required=True)
    parser.add_argument('--original-measure-script',type=Path,required=True)
    parser.add_argument('--original-detector',type=Path,required=True)
    parser.add_argument('--original-efficient-sam',type=Path,required=True)
    parser.add_argument('--per-complete-species',type=int,default=24)
    parser.add_argument('--delay-seconds',type=float,default=1.1)
    a=parser.parse_args(argv)
    run_pilot(a.source_metadata_root,a.output_root,a.original_measure_script,
              a.original_detector,a.original_efficient_sam,
              n_each=a.per_complete_species,delay_seconds=a.delay_seconds)
    return 0

if __name__=='__main__':
    raise SystemExit(main())
, flags=re.IGNORECASE)
MAX_BYTES = 12 * 1024 * 1024
MAX_PER_SPECIES = 24


def selected_metadata(metadata_root: Path, n_each: int) -> tuple[list[dict], dict]:
    if not 1 <= n_each <= MAX_PER_SPECIES:
        raise ValueError('Real-photo tech pilot per-species bound must be 1..24')
    contract = json.loads((metadata_root / 'source_contract.json').read_text(encoding='utf-8'))
    if (contract.get('source_species_frame_sha256') != SOURCE_FRAME_SHA or
            contract.get('created_d2') != SOURCE_CREATED_DATE or
            contract.get('photo_cap', 'INVALID') is not None or
            contract.get('outcome_blind') is not True or
            contract.get('expected_source_species') != 42111):
        raise ValueError('Pilot source must be the verified 42,111 opportunity-frame dated metadata, never historic labelled colours')
    if len(list(metadata_root.glob('species/[0-9]*/checkpoint.json'))) != 9:
        raise ValueError('Expected original nine attempted species; do not accept a changed provider artifact')
    photos=[]
    archives={}
    taxa=[]
    global_ids=set()
    for ck in sorted(metadata_root.glob('species/[0-9]*/checkpoint.json'), key=lambda x:int(x.parent.name)):
        status=json.loads(ck.read_text(encoding='utf-8'))
        tid=int(ck.parent.name)
        if (status.get('taxon_id') != tid or status.get('created_d2') != SOURCE_CREATED_DATE or
                status.get('colour_pixels_opened') is not False or status.get('photo_cap','INVALID') is not None):
            raise ValueError('Unfrozen or already-measured metadata source taxon')
        if status.get('status') != 'complete':
            continue
        manifest=ck.parent/'photo_manifest.ndjson'
        digest=sha256_file(manifest)
        if tid==47129 and digest!=ORIGINAL_RIBES_MANIFEST_SHA:
            raise ValueError('Previously independently verified Ribes source manifest byte identity changed')
        if status.get('metadata_sha256') not in (None, digest):
            raise ValueError('Completed source checkpoint disagrees with its photo manifest')
        species=[]
        seen={}
        conflicting_photo_ids=set()
        with manifest.open('r',encoding='utf-8') as f:
            for line in f:
                obj=json.loads(line)
                pid=str(obj['photo_id'])
                if (int(obj['inat_taxon_id'])!=tid or obj['species']!=status['species'] or
                        not pid.isdigit() or not obj.get('photo_license')):
                    raise ValueError('Original dated photo metadata provenance / license / creator changed')
                if pid in seen:
                    if seen[pid]!=obj:
                        # The same original image appeared in multiple observations.
                        # It is unsuitable for photo-site attribution. Quarantine
                        # its photo ID entirely; never count repeated observations
                        # as independent flowers. This exclusion uses metadata only.
                        conflicting_photo_ids.add(pid)
                    continue
                seen[pid]=obj
                species.append(obj)
        species=[v for v in species if str(v['photo_id']) not in conflicting_photo_ids]
        if not species:
            continue
        # Deterministic, species-equal, blind to image colour and geography.
        species.sort(key=lambda x:hashlib.sha256(f'{SALT}:{tid}:{x["photo_id"]}'.encode()).digest())
        chosen=species[:n_each]
        for obj in chosen:
            pid=str(obj['photo_id'])
            if pid in global_ids:
                raise ValueError('Provider photo ID appears in two purported species')
            global_ids.add(pid)
        photos.extend(chosen)
        archives[str(tid)]={'species':status['species'],'full_manifest_sha256':digest,
                            'complete_source_unique_unambiguous_photo_ids':len(species),
                            'duplicate_conflicting_original_photo_ids_quarantined':len(conflicting_photo_ids),
                            'pilot_selected':len(chosen)}
        taxa.append(tid)
    if len(taxa)!=4 or 47129 not in taxa:
        raise ValueError(f'Prior artifact audit showed exactly four completed taxa, got {taxa}')
    selected_ids=sorted((int(v['inat_taxon_id']),str(v['photo_id'])) for v in photos)
    selected_hash=hashlib.sha256(json.dumps(selected_ids,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    return photos,{'source_artifact_id':SOURCE_ARTIFACT,'source_run_id':SOURCE_RUN,
                   'source_snapshot_date':SOURCE_CREATED_DATE,'initial_frame_sha256':SOURCE_FRAME_SHA,
                   'completed_taxa_available':taxa,'species':archives,'selected_photo_id_pairs_sha256':selected_hash,
                   'selection_outcome_blind':True,'pilot_only_max_photos_per_species':n_each}


def original_photo_bytes(url:str,pid:str)->bytes:
    uri=urlparse(url)
    match=IMAGE_URL.fullmatch(uri.path)
    if (uri.scheme!='https' or uri.netloc!='inaturalist-open-data.s3.amazonaws.com' or
            uri.params or uri.query or uri.fragment or match is None or match.group(1)!=pid):
        raise ValueError('Real photo URL is not the frozen original iNaturalist Open Data S3 object')
    request=Request(url,headers={'User-Agent':'fcp-real-fresh-roi-technical-pilot/1.0 (github.com/zuizui0223/fcp)',
                                 'Accept':'image/*'})
    with urlopen(request,timeout=45) as f:
        raw=f.read(MAX_BYTES+1)
    if not raw or len(raw)>MAX_BYTES:
        raise ValueError('Provider image empty or beyond bounded temporary image memory')
    # Full decode check, not only filename or HTTP status.
    from io import BytesIO
    with Image.open(BytesIO(raw)) as image:
        image.verify()
    return raw


def run_pilot(meta:Path,outdir:Path,source:Path,detector:Path,sam:Path,
              *,n_each:int=24,delay_seconds:float=1.1)->dict:
    if delay_seconds<1.0:
        raise ValueError('At most one original photo request per second')
    if outdir.exists() and any(outdir.iterdir()):
        raise ValueError('One-shot real-photo output must be empty; never silently overwrite pilot')
    chosen, provenance=selected_metadata(meta,n_each)
    method=verify_frozen_method(source,detector,sam)
    outdir.mkdir(parents=True,exist_ok=True)
    failures={}
    hashes={}
    actual_worker=[]
    with tempfile.TemporaryDirectory(prefix='fcp_original_real_photo_pilot_') as scratch:
        temp=Path(scratch)
        images=temp/'images';images.mkdir()
        for i,row in enumerate(chosen):
            if i:time.sleep(delay_seconds)
            mid='FCPR-PILOT-'+hashlib.sha256(f'{SALT}:{row["photo_id"]}'.encode()).hexdigest().upper()[:24]
            name=mid+'.jpg'
            try:
                raw=original_photo_bytes(str(row['photo_url_large']),str(row['photo_id']))
                (images/name).write_bytes(raw)
                hashes[mid]=hashlib.sha256(raw).hexdigest()
                actual_worker.append({'measurement_id':mid,'image_filename':name,'photo_license':row['photo_license']})
            except Exception as exc:
                failures[mid]=f'{type(exc).__name__}:{str(exc)[:190]}'
            row['measurement_id']=mid
            print(f'original_photo_source_progress={i+1}/{len(chosen)} downloaded={len(actual_worker)} failures={len(failures)}',flush=True)
        raw_output=temp/'original_unmodified_roi_v4.csv'
        if actual_worker:
            worker=temp/'blind_worker.csv'
            with worker.open('w',encoding='utf-8',newline='') as f:
                w=csv.DictWriter(f,fieldnames=WORKER_FIELDS);w.writeheader();w.writerows(actual_worker)
            subprocess.run([sys.executable,str(source),'--worker-manifest',str(worker),
                            '--images-dir',str(images),'--detector-weight',str(detector),
                            '--efficient-sam-dir',str(sam),'--output-csv',str(raw_output)],
                           check=True,timeout=2400)
            measured=pd.read_csv(raw_output,dtype={'measurement_id':str}).fillna('')
            if (len(measured)!=len(actual_worker) or measured.measurement_id.duplicated().any() or
                set(measured.measurement_id)!={x['measurement_id'] for x in actual_worker}):
                raise ValueError('Unmodified original ROI-v4 returned incorrect terminal source photo identities')
        else:
            measured=pd.DataFrame(columns=['measurement_id','morph','measurement_status',*[f'flower_fraction_{x}' for x in BIO]])
        merged={r['measurement_id']:r for r in measured.to_dict(orient='records')}
        final=[]
        for source_row in chosen:
            mid=source_row['measurement_id']
            result=merged.get(mid)
            if result is None:
                if mid not in failures:
                    raise ValueError('Source photo omitted by original ROI-v4 evaluator')
                result={'measurement_id':mid,'morph':'mixed_uncertain',
                        'measurement_status':'image_acquisition_failed',
                        'failure_reasons':failures[mid],
                        **{f'flower_fraction_{x}':0.0 for x in BIO}}
            good=result['measurement_status']=='classified_four_state_morph'
            if good!=(result['morph'] in ('white','yellow_orange','red_pink','blue_purple')):
                raise ValueError('Frozen original ROI-v4 morphological terminal status contradiction')
            final.append({**{k:source_row.get(k) for k in (
                 'inat_taxon_id','species','observation_id','photo_id','observer_id',
                 'observed_on','latitude','longitude','positional_accuracy_m',
                 'photo_license','attribution')},
                 'source_url':source_row['photo_url_large'],
                 'technical_pilot_not_full_FCP_study':True,
                 **result,'original_image_bytes_sha256':hashes.get(mid,''),
                 'global_classifiable':bool(good)})
    data=pd.DataFrame(final)
    if len(data)!=len(chosen) or data.photo_id.astype(str).duplicated().any():
        raise ValueError('Pilot terminal photo ledger incomplete or duplicated')
    photo_csv=outdir/'real_original_roi_v4_new_photo_pilot.csv'
    data.to_csv(photo_csv,index=False)
    classified=data.measurement_status.eq('classified_four_state_morph')
    receipt={
        'schema':'fcp_real_original_roi_v4_quarantined_technical_pilot_v1',
        'status':'COMPLETE_REAL_ORIGINAL_PHOTO_TECHNICAL_PILOT_NOT_42111_CENSUS',
        'run_utc':datetime.now(timezone.utc).isoformat(),
        'source':provenance,
        'n_completed_source_species_sampled':len(provenance['completed_taxa_available']),
        'n_real_source_original_photo_ids':len(chosen),
        'n_real_original_image_bytes_downloaded':len(actual_worker),
        'n_real_original_images_processed_by_frozen_roi_v4':len(measured),
        'n_original_provider_download_failures':len(failures),
        'n_real_photo_colour_classified':int(classified.sum()),
        'n_roi_unclassifiable':int((~classified & ~data.measurement_status.eq('image_acquisition_failed')).sum()),
        'per_species_terminal_counts':{str(t):{
           'n_selected_photo_ids':int(len(group)),
           'n_downloaded':int(group.original_image_bytes_sha256.ne('').sum()),
           'n_classifiable':int(group.global_classifiable.sum()),
           'n_image_acquisition_failed':int(group.measurement_status.eq('image_acquisition_failed').sum()),
        } for t,group in data.groupby('inat_taxon_id')},
        'status_counts':{str(k):int(v) for k,v in data.measurement_status.value_counts().items()},
        'raw_source_URLs_only_not_old_photo_colours':True,
        'original_roi_v4_method_hashes':method,
        'real_fresh_photo_terminal_csv_sha256':sha256_file(photo_csv),
        'new_photographs_scored_for_biological_results':False,
        'full_source_census_complete':False,
        'eligible_for_main_fcp_W_null_or_climate_inference':False,
        'n_original_historical_colour_rows_loaded':0,
        'source_image_pixels_persisted':False,
        'selected_photo_cap_applies_only_to_quarantined_technical_pilot':True,
    }
    (outdir/'real_original_roi_v4_new_photo_pilot_receipt.json').write_text(
        json.dumps(receipt,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print('FCP_REAL_PHOTO_PILOT_RESULT '+json.dumps({
        'status':receipt['status'],'n_original_photo_ids':len(chosen),
        'n_real_images_downloaded':len(actual_worker),
        'n_real_images_roi_measured':len(measured),
        'n_classifiable':int(classified.sum()),
        'status_counts':receipt['status_counts'],
        'per_species':receipt['per_species_terminal_counts'],
        'terminal_csv_sha256':receipt['real_fresh_photo_terminal_csv_sha256']},sort_keys=True),flush=True)
    return receipt


def main(argv=None)->int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-metadata-root',type=Path,required=True)
    parser.add_argument('--output-root',type=Path,required=True)
    parser.add_argument('--original-measure-script',type=Path,required=True)
    parser.add_argument('--original-detector',type=Path,required=True)
    parser.add_argument('--original-efficient-sam',type=Path,required=True)
    parser.add_argument('--per-complete-species',type=int,default=24)
    parser.add_argument('--delay-seconds',type=float,default=1.1)
    a=parser.parse_args(argv)
    run_pilot(a.source_metadata_root,a.output_root,a.original_measure_script,
              a.original_detector,a.original_efficient_sam,
              n_each=a.per_complete_species,delay_seconds=a.delay_seconds)
    return 0

if __name__=='__main__':
    raise SystemExit(main())

#!/usr/bin/env python3
"""Fresh-only uncapped image processing, one checkpointed work-unit at a time.

Inputs are NEW photo URLs, never historical measurement rows. All measurement
runs call the original frozen ROI-v4+EfficientSAM+fixed-palette evaluator.
`batch_size` limits temporary memory, NOT species inclusion or final photos.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import re
import time
import subprocess
import sys
import tempfile
from io import BytesIO
from pathlib import Path
from typing import Callable
from urllib.request import Request, urlopen
from urllib.parse import urlsplit
from urllib.error import HTTPError, URLError

import numpy as np
import pandas as pd
from PIL import Image

from fcp_pipeline.fcp_fresh_study import (
    SOURCE_METHOD_COMMIT,SOURCE_MEASUREMENT_BLOB,SOURCE_DETECTOR_SHA256,
    SOURCE_ENCODER_SHA256,SOURCE_DECODER_SHA256,BIO,MORPHS,WORKER,
    validate_method_attestation,
)
from fcp_pipeline.fcp_allphotos import sha256_file

ACQ=('measurement_id','image_filename','photo_url_large','photo_license')
MAX_IMAGE_BYTES=12*1024*1024
DEFAULT_PROVIDER_PAUSE_SECONDS=1.1
_ORIGINAL_OBJECT=re.compile(r'^/photos/[0-9]+/large\.(?:jpg|jpeg|png)$',re.I)


def _git_blob(path:Path)->str:
    raw=path.read_bytes()
    return hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\x00'+raw).hexdigest()


def image_loader(url:str)->bytes:
    """Download original source image with bounded memory and transient retry."""
    u=urlsplit(url)
    if (u.scheme!='https' or u.netloc!='inaturalist-open-data.s3.amazonaws.com'
            or u.query or u.fragment or not _ORIGINAL_OBJECT.fullmatch(u.path)):
        raise ValueError('Photo URL is not a frozen original provider S3 object')
    req=Request(url,headers={'User-Agent':'fcp-fresh-photo-measurement/1.0 (github.com/zuizui0223/fcp)',
                            'Accept':'image/*'})
    for attempt in range(3):
        try:
            with urlopen(req,timeout=45) as response:
                buf=response.read(MAX_IMAGE_BYTES+1)
            if not buf or len(buf)>MAX_IMAGE_BYTES:
                raise ValueError('Original photo empty or exceeds 12MiB bounded acquisition')
            with Image.open(BytesIO(buf)) as img:
                img.verify()
            return buf
        except HTTPError as exc:
            if exc.code not in (408,425,429,500,502,503,504) or attempt==2:
                raise
        except (URLError,TimeoutError,ConnectionError):
            if attempt==2:
                raise
        time.sleep(2**attempt)
    raise RuntimeError('Original photo retries exhausted')


def verify_frozen_method(source:Path, detector:Path, sam:Path)->dict:
    if _git_blob(source)!=SOURCE_MEASUREMENT_BLOB:
        raise ValueError('Original frozen ROI-v4 measurement source changed')
    for k,path,target in (
        ('detector_sha256',detector,SOURCE_DETECTOR_SHA256),
        ('efficient_sam_encoder_sha256',sam/'efficient_sam_vitt_encoder.onnx',SOURCE_ENCODER_SHA256),
        ('efficient_sam_decoder_sha256',sam/'efficient_sam_vitt_decoder.onnx',SOURCE_DECODER_SHA256)):
        if not path.is_file() or sha256_file(path)!=target:
            raise ValueError(f'Frozen original runtime byte hash mismatch: {k}')
    return {'source_commit':SOURCE_METHOD_COMMIT,
            'measurement_git_blob':SOURCE_MEASUREMENT_BLOB,
            'detector_sha256':SOURCE_DETECTOR_SHA256,
            'efficient_sam_encoder_sha256':SOURCE_ENCODER_SHA256,
            'efficient_sam_decoder_sha256':SOURCE_DECODER_SHA256,
            'location_blind':True,'no_old_measurement_csv_inputs':True,
            'persisted_image_pixels':False}


def preflight(firewall:Path,part:int,batch:int,batch_size:int):
    if not (part>=0 and batch>=0 and batch_size>0):
        raise ValueError('Positive bounded compute unit required')
    frozen=json.loads((firewall/'fresh_study_freeze.json').read_text())
    if frozen['status']!='FRESH_BLIND_METADATA_FROZEN' or part>=frozen['worker_packets']:
        raise ValueError('Fresh source photo census is not authorized or partition wrong')
    if frozen['blind_firewall_sha256']!=sha256_file(firewall/'firewall_receipt.json'):
        raise ValueError('Fresh source firewall SHA changed')
    for name in (f'worker_partition_{part:03}.csv',f'acquisition_partition_{part:03}.csv'):
        if sha256_file(firewall/name)!=frozen['frozen_file_sha256'][name]:
            raise ValueError(f'Photo manifest changed after frozen source checkpoint: {name}')
    # The worker CSV has EXACTLY three blind fields, no coordinates, identity,
    # source URL, or environmental outcomes. URL lookup happens only here.
    workers=pd.read_csv(firewall/f'worker_partition_{part:03}.csv',dtype=str).fillna('')
    key=pd.read_csv(firewall/f'acquisition_partition_{part:03}.csv',dtype=str).fillna('')
    if tuple(workers.columns)!=WORKER or tuple(key.columns)!=ACQ:
        raise ValueError('Worker is non-blind or source acquisition key invalid')
    if workers.measurement_id.duplicated().any() or key.measurement_id.duplicated().any():
        raise ValueError('Worker/acquisition photo IDs must be unique')
    if set(workers.measurement_id)!=set(key.measurement_id):
        raise ValueError('Missing source URLs or source photo identity')
    if len(workers)==0:
        raise ValueError('This partition is empty; no image evaluation required')
    first=batch*batch_size
    if first>=len(workers):
        raise ValueError('Requested batch past last source photograph')
    sub=workers.iloc[first:first+batch_size].copy()
    acq=sub.merge(key,on=['measurement_id','image_filename','photo_license'],how='left',validate='one_to_one')
    if acq.photo_url_large.isna().any():
        raise ValueError('Unmatched provider image URL')
    return sub,acq,math.ceil(len(workers)/batch_size)


def validate_terminal(rows:pd.DataFrame, expected:list[str])->None:
    needed={'measurement_id','morph','measurement_status',*[f'flower_fraction_{x}' for x in BIO]}
    if needed-set(rows):
        raise ValueError('Original measurement schema missing frozen palette')
    if len(rows)!=len(expected) or rows.measurement_id.duplicated().any() or set(rows.measurement_id)!=set(expected):
        raise ValueError('Terminal results must include every fresh source image exactly once')
    if {'species','photo_id','observer_id','latitude','longitude','photo_url_large'} & set(rows):
        raise ValueError('Blind measurement process was exposed to source identity')
    for _,row in rows.iterrows():
        classified=row.measurement_status=='classified_four_state_morph'
        if classified!=(row.morph in MORPHS):
            raise ValueError('Original 4-state classification gate failed')
        if not classified and row.morph!='mixed_uncertain':
            raise ValueError('Classifiable failure cannot be replaced with new colour')
        if classified:
            z=pd.to_numeric(row[[f'flower_fraction_{x}' for x in BIO]],errors='raise').to_numpy(float)
            if not np.isfinite(z).all() or (z < -1e-12).any() or abs(z.sum()-1)>1e-7:
                raise ValueError('Frozen nine-colour composition invalid')


def evaluate_with_frozen_source(worker:Path,images:Path,source:Path,detector:Path,sam:Path,out:Path):
    subprocess.run([sys.executable,str(source),'--worker-manifest',str(worker),
                    '--images-dir',str(images),'--detector-weight',str(detector),
                    '--efficient-sam-dir',str(sam),'--output-csv',str(out)],check=True)
    return pd.read_csv(out,dtype={'measurement_id':str})


def run_batch(firewall:Path, part:int,batch:int,batch_size:int, outdir:Path,
              *, source:Path,detector:Path,sam:Path,
              loader:Callable[[str],bytes]|None=None,
              evaluator:Callable|None=None,
              _test_only_skip_method_check:bool=False)->dict:
    worker,acq,number_batches=preflight(firewall,part,batch,batch_size)
    if _test_only_skip_method_check:
        if evaluator is None:raise ValueError('Frozen source checks may only be bypassed for synthetic evaluator tests')
        fingerprint={'synthetic_test_only':True,'no_old_measurement_csv_inputs':True}
    else:
        fingerprint=verify_frozen_method(source,detector,sam)
    outdir.mkdir(parents=True,exist_ok=True)
    if not _test_only_skip_method_check:
        seal=outdir/'frozen_method_attestation.json'
        if seal.is_file():
            if json.loads(seal.read_text(encoding='utf-8'))!=fingerprint:
                raise ValueError('Original frozen method drifted after first image batch')
        else:
            seal.write_text(json.dumps(fingerprint,sort_keys=True,indent=2)+'\n',encoding='utf-8')
    stem=f'terminal_p{part:03}_b{batch:07}'
    result_file=outdir/(stem+'.csv')
    receipt_file=outdir/(stem+'.json')
    if result_file.exists() or receipt_file.exists():
        raise ValueError('Terminal work unit already exists; no silent overwrite/reclassification')
    failures=[]
    with tempfile.TemporaryDirectory(prefix='fcp_new_roi_v4_') as tmp:
        work=Path(tmp)
        images=work/'images';images.mkdir()
        okay=[]
        for index,r in enumerate(acq.itertuples(index=False)):
            if loader is None and index:
                time.sleep(DEFAULT_PROVIDER_PAUSE_SECONDS)
            try:
                payload=(loader or image_loader)(str(r.photo_url_large))
                if not payload:raise ValueError('empty image')
                (images/str(r.image_filename)).write_bytes(payload)
                okay.append({'measurement_id':r.measurement_id,'image_filename':r.image_filename,
                             'photo_license':r.photo_license})
            except Exception as exc:
                failures.append({'measurement_id':r.measurement_id,'morph':'mixed_uncertain',
                                 'measurement_status':'image_acquisition_failed',
                                 'failure_reasons':f'{type(exc).__name__}:{str(exc)[:160]}',
                                 **{f'flower_fraction_{c}':0.0 for c in BIO}})
        measured=pd.DataFrame(columns=['measurement_id','morph','measurement_status',*[f'flower_fraction_{c}' for c in BIO]])
        if okay:
            blind=pd.DataFrame(okay,columns=list(WORKER))
            target=work/'blind_worker.csv';blind.to_csv(target,index=False)
            measured=(evaluator or evaluate_with_frozen_source)(target,images,source,detector,sam,work/'blind_results.csv')
        results=pd.concat([measured,pd.DataFrame(failures)],ignore_index=True,sort=False)
        validate_terminal(results,worker.measurement_id.astype(str).tolist())
        results=results.sort_values('measurement_id',kind='stable')
        results.to_csv(result_file,index=False)
        record={'schema':'fcp_fresh_blind_terminal_partition_v1','status':'COMPLETE_WORK_UNIT',
                'partition':part,'batch':batch,'compute_batch_size_not_photo_cap':batch_size,
                'number_batches_for_this_partition':number_batches,
                'n_source_photo_ids':len(worker),'n_terminal_photo_ids':len(results),
                'n_classifiable':int(results.measurement_status.eq('classified_four_state_morph').sum()),
                'n_acquisition_failed':len(failures),
                'original_provider_rate_limit_seconds_between_photos':(
                    DEFAULT_PROVIDER_PAUSE_SECONDS if loader is None else None),
                'bounded_provider_original_image_bytes':MAX_IMAGE_BYTES,
                'photo_failures_preserved_not_replaced':True,
                'selected_source_worker_sha256':hashlib.sha256(worker.to_csv(index=False).encode()).hexdigest(),
                'terminal_csv_sha256':sha256_file(result_file),
                'method':fingerprint,
                'image_bytes_persisted':False,'reused_historical_measurement_csv':False,
                'full_study_complete':False}
        receipt_file.write_text(json.dumps(record,indent=2)+'\n')
        return record


def main(argv:list[str]|None=None)->int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--firewall',type=Path,required=True)
    p.add_argument('--partition',type=int,required=True)
    p.add_argument('--batch',type=int,required=True)
    p.add_argument('--batch-size',type=int,default=64)
    p.add_argument('--original-measure-script',type=Path,required=True)
    p.add_argument('--original-detector',type=Path,required=True)
    p.add_argument('--original-efficient-sam',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(argv)
    result=run_batch(a.firewall,a.partition,a.batch,a.batch_size,a.output,
                     source=a.original_measure_script,detector=a.original_detector,
                     sam=a.original_efficient_sam)
    print(json.dumps({k:result[k] for k in ('status','partition','batch','n_source_photo_ids','n_classifiable','n_acquisition_failed')}))
    return 0


if __name__=='__main__':
    raise SystemExit(main())
#!/usr/bin/env python3
"""Extract and verify a prehashed *dated* iNaturalist Open Data metadata bundle.

NO old photo colour CSVs or raw image pixels are read. A full fixed-date tar.gz
must be downloaded to an adequately provisioned disk before this script runs.
The source tar sha256 is a mandatory, externally recorded acquisition digest.
"""
from __future__ import annotations
import argparse
import gzip
import hashlib
import json
import re
import shutil
import tarfile
from datetime import date
from pathlib import Path
from typing import BinaryIO

TABLES = ('taxa', 'observations', 'observers', 'photos')
SOURCE_KEY = re.compile(r'^metadata/inaturalist-open-data-(\d{8})\.tar\.gz$')
CHUNK = 1 << 20


def sha256_file(path: Path) -> str:
    h=hashlib.sha256()
    with path.open('rb') as f:
        for part in iter(lambda: f.read(CHUNK), b''):
            h.update(part)
    return h.hexdigest()


def unpack_dated_metadata(archive:Path, outdir:Path, *, s3_key:str,
                          expected_archive_sha256:str) -> dict:
    if not SOURCE_KEY.fullmatch(s3_key):
        raise ValueError('Only a dated iNaturalist metadata archive is admissible (never latest)')
    day=SOURCE_KEY.fullmatch(s3_key).group(1)
    snapshot_date=date.fromisoformat(f'{day[:4]}-{day[4:6]}-{day[6:]}').isoformat()
    if not (len(expected_archive_sha256)==64 and
            all(x in '0123456789abcdef' for x in expected_archive_sha256)):
        raise ValueError('Pre-acquisition archive SHA256 must be explicit')
    actual=sha256_file(archive)
    if actual!=expected_archive_sha256:
        raise ValueError('Dated source archive SHA256 mismatch')
    if outdir.exists() and any(outdir.iterdir()):
        raise ValueError('Snapshot extraction output is not empty; cannot mix source archives')
    outdir.mkdir(parents=True,exist_ok=True)
    source_members={}
    result={}
    # Streaming gzip/tar means member bytes are never silently loaded in memory.
    with tarfile.open(archive,'r|gz') as reader:
        for member in reader:
            if member.isdir():
                continue
            if not member.isfile():
                raise ValueError('Source tar contains a non-regular file')
            name=Path(member.name).name
            stem=None
            raw_is_gz=False
            for table in TABLES:
                if name==f'{table}.csv.gz':
                    stem=table;raw_is_gz=True;break
                if name==f'{table}.csv':
                    stem=table;break
            if stem is None:
                # Other provider tables (projects, observation projects) are not FCP inputs.
                continue
            if stem in result:
                raise ValueError(f'Original snapshot contains duplicated metadata table: {stem}')
            # The real dated provider bundle has a top-level directory named
            # inaturalist-open-data-YYYYMMDD (confirmed from S3 range probe).
            allowed_prefix=f'inaturalist-open-data-{day}/'
            if (('/' in member.name and not member.name.startswith(('./','metadata/',allowed_prefix)))
                    or '..' in Path(member.name).parts):
                raise ValueError('Unexpected provider member path; snapshot layout not recognized')
            source=reader.extractfile(member)
            if source is None:
                raise ValueError('Cannot read source table from dated tar')
            target=outdir/f'{stem}.csv.gz'
            with target.open('wb') as destination:
                if raw_is_gz:
                    shutil.copyfileobj(source,destination,CHUNK)
                else:
                    # tar.gz bundles may contain UNCOMPRESSED CSV tables.
                    # Recompress deterministically to avoid huge multi-100GB plain CSV output.
                    with gzip.GzipFile(fileobj=destination,mode='wb',filename='',mtime=0) as zipped:
                        shutil.copyfileobj(source,zipped,CHUNK)
            with gzip.open(target,'rb') as check:
                header=check.readline(10000)
            if b'\t' not in header or not header.endswith((b'\n',b'\r\n')):
                raise ValueError(f'{stem}: extracted table is not official tab-separated metadata')
            result[stem]={'path':str(target.resolve()),'sha256':sha256_file(target),
                          'bytes':target.stat().st_size,'member_name':member.name,
                          'original_member_gz':raw_is_gz,
                          'header_fields':header.decode('utf-8').rstrip('\r\n').split('\t')}
    missing=set(TABLES)-set(result)
    if missing:
        raise ValueError(f'Incomplete original snapshot bundle missing tables: {sorted(missing)}')
    receipt={'schema':'fcp_dated_inaturalist_metadata_bundle_pin_v1',
             'status':'FOUR_ORIGINAL_METADATA_TABLES_EXTRACTED_AND_PINNED',
             'source_kind':'inaturalist_open_dataset_monthly_snapshot',
             's3_object_key':s3_key,'snapshot_date':snapshot_date,
             'full_bundle_sha256':actual,'full_bundle_bytes':archive.stat().st_size,
             'metadata_files_sha256':{k:result[k]['sha256'] for k in TABLES},
             'tables':result,'source_images_opened':False,
             'historical_photo_colours_read':False,
             'biological_photo_measurement_authorized':False,
             'hashes_recorded_before_any_FCP_outcome':True}
    (outdir/'fcp_source_archive_pin.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n')
    return receipt


def main(argv=None)->int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--archive',type=Path,required=True)
    p.add_argument('--s3-key',required=True)
    p.add_argument('--archive-sha256',required=True)
    p.add_argument('--outdir',type=Path,required=True)
    a=p.parse_args(argv)
    receipt=unpack_dated_metadata(a.archive,a.outdir,s3_key=a.s3_key,
                                  expected_archive_sha256=a.archive_sha256)
    print(json.dumps({'status':receipt['status'],'archive_sha256':receipt['full_bundle_sha256'],
                      'snapshot_date':receipt['snapshot_date'],
                      'metadata_files_sha256':receipt['metadata_files_sha256']},sort_keys=True))
    return 0

if __name__=='__main__':
    raise SystemExit(main())

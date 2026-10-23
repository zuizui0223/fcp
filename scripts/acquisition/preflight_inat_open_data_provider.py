#!/usr/bin/env python3
"""Small, real S3 metadata format probe. NOT a complete FCP census.

Only HEAD metadata and the first 128KiB of each mutable 'latest' gzip object.
Never opens images or creates a biological measurement dataset.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
from urllib.request import Request, urlopen
from urllib.parse import urlencode
from urllib.error import HTTPError
import tarfile
from xml.etree import ElementTree
import zlib

BUCKET = 'https://inaturalist-open-data.s3.amazonaws.com/'
EXPECTED = {
    'taxa': {'taxon_id', 'ancestry', 'name', 'rank'},
    'observations': {'observation_uuid', 'taxon_id', 'observer_id', 'latitude', 'longitude',
                     'positional_accuracy', 'quality_grade', 'observed_on'},
    'photos': {'photo_id', 'observation_uuid', 'observer_id', 'extension', 'license'},
}
MAX_PROBE_BYTES = 131072


def probe_object(name: str) -> dict:
    url = f'{BUCKET}{name}.csv.gz'
    with urlopen(Request(url, method='HEAD'), timeout=25) as response:
        if response.status != 200:
            raise RuntimeError(f'{name}: S3 HEAD returned {response.status}')
        headers = response.headers
        size = int(headers.get('Content-Length', '-1'))
        if size <= MAX_PROBE_BYTES:
            raise RuntimeError(f'{name}: unlikely real public Open Data object size: {size}')
        rec = {'object_key': f'{name}.csv.gz', 's3_head_content_length': size,
               's3_last_modified': headers.get('Last-Modified'),
               's3_etag': headers.get('ETag'),
               's3_version_id': headers.get('x-amz-version-id')}
    request = Request(url, headers={'Range': f'bytes=0-{MAX_PROBE_BYTES-1}'})
    with urlopen(request, timeout=25) as response:
        if response.status != 206 or not response.headers.get('Content-Range', '').startswith('bytes 0-'):
            raise RuntimeError(f'{name}: S3 ignored small range limit; refusing large download')
        sample = response.read(MAX_PROBE_BYTES + 1)
        if len(sample) > MAX_PROBE_BYTES:
            raise RuntimeError(f'{name}: sample exceeded cap')
    try:
        preview = zlib.decompressobj(16 + zlib.MAX_WBITS).decompress(sample)
        header = preview.split(b'\n', 1)[0].decode('utf-8').rstrip('\r')
    except (zlib.error, UnicodeDecodeError) as exc:
        raise RuntimeError(f'{name}: truncated gzip did not expose original metadata header') from exc
    fields = [s.strip() for s in header.split('\t')]
    if not EXPECTED[name].issubset(fields):
        raise ValueError(f'{name}: real provider schema mismatch: absent {sorted(EXPECTED[name] - set(fields))}; fields={fields}')
    return {**rec, 'header_fields': fields, 'required_fields_present': True,
            'gzip_sample_only_bytes': len(sample)}


def list_dated_archives() -> dict:
    # Listing is OPTIONAL. Some public S3 buckets allow individual downloads but forbid ListBucket.
    url = BUCKET + '?' + urlencode({'list-type': '2',
                                   'prefix': 'metadata/inaturalist-open-data-',
                                   'max-keys': '1000'})
    try:
        with urlopen(Request(url), timeout=25) as response:
            data = response.read(256_001)
            if len(data) > 256000 or response.status != 200:
                return {'status': 'list_unavailable_or_too_large'}
        tree = ElementTree.fromstring(data)
        items = []
        for node in tree.iter():
            if node.tag.endswith('Contents'):
                values = {child.tag.split('}')[-1]: child.text for child in node}
                key = values.get('Key', '')
                if key.startswith('metadata/inaturalist-open-data-') and key.endswith('.tar.gz'):
                    items.append({'key': key, 'size': int(values['Size'])})
        dated = [x for x in items if 'latest' not in x['key']]
        return {'status': 'listed', 'total_candidates': len(items),
                'dated_snapshots': sorted(dated, key=lambda q: q['key'])[-20:]}
    except (OSError, HTTPError, ValueError, ElementTree.ParseError) as exc:
        return {'status': 'list_unavailable', 'reason': type(exc).__name__}


def probe_dated_bundle_prefix() -> dict:
    """Read <=256KiB of a dated archive to identify its real tar member schema."""
    key='metadata/inaturalist-open-data-20260927.tar.gz'
    url=BUCKET+key
    limit=262144
    with urlopen(Request(url,method='HEAD'),timeout=25) as response:
        if response.status!=200:
            raise RuntimeError('S3 dated archive HEAD failed')
        source_size=int(response.headers.get('Content-Length','-1'))
        etag=response.headers.get('ETag')
    with urlopen(Request(url,headers={'Range':f'bytes=0-{limit-1}'}),timeout=25) as response:
        if response.status!=206:
            raise RuntimeError('Dated metadata archive did not honor HTTP range')
        content=response.read(limit+1)
        if len(content)>limit:
            raise RuntimeError('Dated archive prefix exceeded exact safety budget')
    # The first gzip record yields a plain POSIX tar header; this is NOT a
    # checksum of the full dated snapshot or a complete source download.
    prefix=zlib.decompressobj(16+zlib.MAX_WBITS).decompress(content)
    if len(prefix)<512:
        raise RuntimeError('Tar prefix too short for first metadata member')
    member=tarfile.TarInfo.frombuf(prefix[:512],encoding='utf-8',errors='surrogateescape')
    return {'s3_key':key,'source_archive_bytes':source_size,
            'etag_not_a_sha256':etag,'first_member_name':member.name,
            'first_member_size':member.size,'first_member_type':member.type.decode('latin1'),
            'partial_range_bytes':len(content),'full_archive_sha256_verified':False}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args(argv)
    receipt = {'schema': 'fcp_inat_real_provider_schema_preflight_v1',
               'checked_at_utc': datetime.now(timezone.utc).isoformat(),
               'status': 'PROVIDER_SCHEMA_ONLY_NOT_A_FROZEN_SNAPSHOT',
               'source': 'public iNaturalist Open Data mutable CSV root',
               'measured_photo_rows': 0, 'used_prior_colour_labels': False,
               'fully_downloaded_any_source': False, 'snapshot_sha256_pinned': False,
               'objects': {name: probe_object(name) for name in EXPECTED},
               'dated_archive_availability': list_dated_archives(),
               'dated_bundle_prefix': probe_dated_bundle_prefix()}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(receipt, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(receipt, ensure_ascii=False, sort_keys=True))
    return 0

if __name__ == '__main__':
    raise SystemExit(main())

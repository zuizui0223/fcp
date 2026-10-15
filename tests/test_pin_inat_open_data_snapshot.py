"""Tests only synthetic miniature tar sources: never real biological image files."""
from __future__ import annotations
import gzip
import hashlib
import io
import json
import tarfile
from pathlib import Path

import pytest
from scripts.acquisition.pin_inat_open_data_snapshot import unpack_dated_metadata,sha256_file


def archive(tmp:Path,*,as_gz:bool=True,missing:str|None=None)->Path:
    target=tmp/'provider_20260927.tar.gz'
    fields={'taxa':'taxon_id\tancestry\trank\tname\n1\t\tspecies\tA\n',
            'observations':'observation_uuid\ttaxon_id\n1\t1\n',
            'observers':'observer_id\tlogin\tname\n1\ta\tAlice\n',
            'photos':'photo_id\tobservation_uuid\n12\t1\n'}
    with tarfile.open(target,'w:gz') as writer:
        for name,text in fields.items():
            if name==missing:continue
            payload=gzip.compress(text.encode(),mtime=0) if as_gz else text.encode()
            path=f'inaturalist-open-data-20260927/{name}.csv'+('.gz' if as_gz else '')
            info=tarfile.TarInfo(path);info.size=len(payload)
            writer.addfile(info,io.BytesIO(payload))
    return target


@pytest.mark.parametrize('original_members_gzip',[True,False])
def test_all_four_tables_from_dated_original_archive(tmp_path,original_members_gzip):
    p=archive(tmp_path,as_gz=original_members_gzip)
    out=tmp_path/'out'
    r=unpack_dated_metadata(p,out,s3_key='metadata/inaturalist-open-data-20260927.tar.gz',
                            expected_archive_sha256=sha256_file(p))
    assert r['snapshot_date']=='2026-09-27'
    assert set(r['metadata_files_sha256'])=={'taxa','observations','observers','photos'}
    assert all((out/f'{x}.csv.gz').exists() for x in r['metadata_files_sha256'])
    assert all(x['original_member_gz']==original_members_gzip for x in r['tables'].values())
    assert r['historical_photo_colours_read'] is False
    assert json.loads((out/'fcp_source_archive_pin.json').read_text())==r
    with pytest.raises(ValueError,match='not empty'):
        unpack_dated_metadata(p,out,s3_key='metadata/inaturalist-open-data-20260927.tar.gz',
                              expected_archive_sha256=sha256_file(p))


def test_reject_wrong_sha_and_mutable_latest_and_missing_metadata(tmp_path):
    p=archive(tmp_path)
    with pytest.raises(ValueError,match='SHA256 mismatch'):
        unpack_dated_metadata(p,tmp_path/'wrong',s3_key='metadata/inaturalist-open-data-20260927.tar.gz',
                              expected_archive_sha256='0'*64)
    with pytest.raises(ValueError,match='dated'):
        unpack_dated_metadata(p,tmp_path/'latest',s3_key='metadata/inaturalist-open-data-latest.tar.gz',
                              expected_archive_sha256=sha256_file(p))
    q=archive(tmp_path,missing='observers')
    with pytest.raises(ValueError,match='missing tables'):
        unpack_dated_metadata(q,tmp_path/'incomplete',s3_key='metadata/inaturalist-open-data-20260927.tar.gz',
                              expected_archive_sha256=sha256_file(q))

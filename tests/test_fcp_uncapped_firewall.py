"""Synthetic tests of blind acquisition-data separation and complete manifest gates."""
import csv
import importlib.util
import json
from pathlib import Path
import pytest

SCRIPT=Path(__file__).resolve().parents[1]/'scripts/acquisition/build_fcp_uncapped_blind_batches.py'
spec=importlib.util.spec_from_file_location('allphotos_firewall',SCRIPT)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def make_dir(root,count):
    d=root/'species';d.mkdir()
    (d/'checkpoint.json').write_text(json.dumps(dict(status='complete',taxon_id=19,
            photo_cap=None,colour_pixels_opened=False)))
    with (d/'photo_manifest.ndjson').open('w') as f:
        for i in range(count):
            f.write(json.dumps(dict(inat_taxon_id=19,species='S',observation_id=i//3+1,
                        photo_id=i+1,photo_url_large=f'https://example.com/{i}/large.jpg',
                        photo_license='cc-by',observer_id=i//12+1,latitude=35.0,
                        longitude=137.0,observed_on='2024-07-01'))+'\n')
    return d


def test_every_photo_included_and_workers_blind(tmp_path):
    d=make_dir(tmp_path,117)
    res=m.freeze_inputs([d],tmp_path/'out',partitions=7)
    assert res['eligible_total_photo_ids']==117
    allworkers=[]
    for path in (tmp_path/'out').glob('worker_partition_*.csv'):
        with path.open() as f:
            a=list(csv.DictReader(f))
            assert len(a)==0 or set(a[0])==set(m.WORKER)
            allworkers.extend(a)
    assert len(allworkers)==117
    assert len({x['measurement_id'] for x in allworkers})==117
    assert all('https://' not in json.dumps(x) for x in allworkers)
    assert 'photo_url_large' in (tmp_path/'out'/'acquisition_key.csv').read_text().splitlines()[0]
    assert 'species' in (tmp_path/'out'/'metadata_join_key.csv').read_text().splitlines()[0]


def test_incomplete_manifest_cannot_freeze(tmp_path):
    d=make_dir(tmp_path,3)
    src=d/'checkpoint.json';x=json.loads(src.read_text());x['status']='in_progress';src.write_text(json.dumps(x))
    with pytest.raises(ValueError,match='Partial/nonblind'):
        m.freeze_inputs([d],tmp_path/'out')


def test_idempotent_checkpoint_duplicate_and_conflicting_duplicates(tmp_path):
    d=make_dir(tmp_path,3)
    manifest=d/'photo_manifest.ndjson'
    line=manifest.read_text().splitlines()[0]
    manifest.write_text(manifest.read_text()+line+'\n')
    result=m.freeze_inputs([d],tmp_path/'out')
    assert result['eligible_total_photo_ids']==3
    assert result['duplicate_retry_rows_reconciled']==1
    entry=json.loads(line)
    entry['longitude']=130.0
    manifest.write_text(manifest.read_text()+json.dumps(entry)+'\n')
    with pytest.raises(ValueError,match='conflicting source metadata'):
        m.freeze_inputs([d],tmp_path/'conflict')


def test_previously_measured_image_excluded_and_source_checked(tmp_path, monkeypatch):
    import hashlib
    d=make_dir(tmp_path,4)
    source=tmp_path/'previous.csv'
    source.write_text('photo_id,inat_taxon_id\n1,19\n2,19\n')
    checksum=hashlib.sha256(source.read_bytes()).hexdigest()
    monkeypatch.setitem(m.LEGACY_PHOTO_SHA256,'discovery',checksum)
    receipt=m.freeze_inputs([d],tmp_path/'out',already_measured={'discovery':source})
    assert receipt['eligible_total_photo_ids']==2
    assert receipt['n_original_photo_ids_preloaded']==2
    assert receipt['n_new_metadata_photo_ids_already_in_original']==2
    with pytest.raises(ValueError,match='SHA256 mismatch'):
        source.write_text(source.read_text()+'3,19\n')
        m.freeze_inputs([d],tmp_path/'bad',already_measured={'discovery':source})

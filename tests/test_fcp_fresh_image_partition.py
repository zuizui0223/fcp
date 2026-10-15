"""Synthetic byte-transfer tests; no iNaturalist requests or real flower results."""
from __future__ import annotations

import csv
from io import BytesIO
from pathlib import Path

import pandas as pd
import pytest
from PIL import Image

from fcp_pipeline.fcp_fresh_study import BIO,freeze_fresh
from test_fcp_fresh_study import setup_source

# Production code is in a script to keep source URL access outside flower workers.
import importlib.util
SCRIPT=Path(__file__).resolve().parents[1]/'scripts/acquisition/run_fcp_fresh_image_partition.py'
spec=importlib.util.spec_from_file_location('fcp_fresh_image_partition',SCRIPT)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def fake_png()->bytes:
    a=BytesIO();Image.new('RGB',(40,40),(255,255,255)).save(a,format='PNG');return a.getvalue()


def fake_measure(blind,images,source,detector,sam,out):
    rows=list(csv.DictReader(blind.open()))
    assert rows and set(rows[0])=={'measurement_id','image_filename','photo_license'}
    assert all((images/r['image_filename']).exists() for r in rows)
    return pd.DataFrame([{'measurement_id':r['measurement_id'],'morph':'white',
                          'measurement_status':'classified_four_state_morph',
                          **{f'flower_fraction_{c}':int(c=='white') for c in BIO}}
                         for r in rows])


def test_batched_fresh_measurement_uses_every_photo_and_cannot_claim_complete(tmp_path):
    frame,src=setup_source(tmp_path)
    firewall=tmp_path/'frozen';freeze_fresh(frame,src,firewall,shards=2,workers=2,expected=4)
    key=pd.read_csv(firewall/'acquisition_partition_000.csv',dtype=str)
    assert len(key)>1
    results=tmp_path/'sealed'
    for batch in range(len(key)):
        receipt=m.run_batch(firewall,0,batch,1,results,source=tmp_path/'test.py',
                            detector=tmp_path/'model.pt',sam=tmp_path/'sam',
                            loader=lambda _:fake_png(),evaluator=fake_measure,
                            _test_only_skip_method_check=True)
        assert receipt['n_source_photo_ids']==1
        assert receipt['full_study_complete'] is False
        assert receipt['reused_historical_measurement_csv'] is False
        assert receipt['number_batches_for_this_partition']==len(key)
        assert receipt['method']['synthetic_test_only'] is True
    assert len(list(results.glob('terminal_*.csv')))==len(key)
    with pytest.raises(ValueError,match='already exists'):
        m.run_batch(firewall,0,0,1,results,source=tmp_path/'x',detector=tmp_path/'y',sam=tmp_path/'z',
                    loader=lambda _:fake_png(),evaluator=fake_measure,_test_only_skip_method_check=True)
    with pytest.raises(ValueError,match='past last source photograph'):
        m.preflight(firewall,0,len(key),1)


def test_failed_download_is_terminal_and_unselected_photo_stays_future_work(tmp_path):
    frame,src=setup_source(tmp_path)
    firewall=tmp_path/'firewall';freeze_fresh(frame,src,firewall,shards=2,workers=1,expected=4)
    calls=[]
    def loader(url):
        calls.append(url)
        if len(calls)==2:
            raise OSError('synthetic 429')
        return fake_png()
    receipt=m.run_batch(firewall,0,0,3,tmp_path/'out',source=tmp_path/'x',detector=tmp_path/'y',sam=tmp_path/'z',
                        loader=loader,evaluator=fake_measure,_test_only_skip_method_check=True)
    assert receipt['n_source_photo_ids']==3
    assert receipt['n_terminal_photo_ids']==3
    assert receipt['n_acquisition_failed']==1
    assert receipt['n_classifiable']==2
    result=pd.read_csv(next((tmp_path/'out').glob('terminal_*.csv')))
    assert result.loc[result.morph.eq('mixed_uncertain'),'measurement_status'].eq('image_acquisition_failed').all()
    assert receipt['full_study_complete'] is False
    assert receipt['image_bytes_persisted'] is False


def test_production_method_refuses_unauthorized_weights_before_images(tmp_path):
    frame,src=setup_source(tmp_path)
    firewall=tmp_path/'freeze';freeze_fresh(frame,src,firewall,shards=2,workers=1,expected=4)
    with pytest.raises(ValueError,match='Original frozen ROI-v4 measurement source changed'):
        source=tmp_path/'changed.py';source.write_text('print("untrusted")')
        m.run_batch(firewall,0,0,1,tmp_path/'out',source=source,detector=tmp_path/'detector.pt',sam=tmp_path/'sam',
                    loader=lambda _:fake_png(),evaluator=fake_measure)
    assert not list((tmp_path/'out').glob('terminal_*.csv'))


def test_synthetic_batch_cannot_pass_biological_assembly(tmp_path):
    from fcp_pipeline.fcp_fresh_study import reassemble_fresh
    from test_fcp_fresh_study import attestation
    frame,src=setup_source(tmp_path)
    firewall=tmp_path/'freeze';freeze_fresh(frame,src,firewall,shards=2,workers=1,expected=4)
    m.run_batch(firewall,0,0,100,tmp_path/'fake',source=tmp_path/'x',detector=tmp_path/'y',sam=tmp_path/'z',
                loader=lambda _:fake_png(),evaluator=fake_measure,_test_only_skip_method_check=True)
    with pytest.raises(ValueError,match='frozen method'):
        reassemble_fresh(firewall,tmp_path/'fake',attestation(tmp_path),tmp_path/'reject')

def test_real_s3_loader_preserves_uppercase_extension_and_bans_foreign_urls(monkeypatch):
    import urllib.error
    payload=fake_png()
    calls=[]
    class FakeResponse:
        def __enter__(self):return self
        def __exit__(self,*args):return False
        def read(self,n):
            assert n==m.MAX_IMAGE_BYTES+1
            return payload
    def open_local_only(req,timeout):
        calls.append(req.full_url)
        assert timeout==45
        return FakeResponse()
    monkeypatch.setattr(m,'urlopen',open_local_only)
    allowed='https://inaturalist-open-data.s3.amazonaws.com/photos/1234/large.JPG'
    assert m.image_loader(allowed)==payload
    assert calls==[allowed]
    for url in ('http://inaturalist-open-data.s3.amazonaws.com/photos/1234/large.jpg',
                'https://other.example.org/photos/1234/large.jpg',
                'https://inaturalist-open-data.s3.amazonaws.com/photos/1234/large.jpg?resize=1',
                'https://inaturalist-open-data.s3.amazonaws.com/private/1234.jpg'):
        with pytest.raises(ValueError,match='frozen original provider'):
            m.image_loader(url)
    assert calls==[allowed]


def test_real_s3_loader_retries_transient_http_but_never_replaces_source_photo(monkeypatch):
    from urllib.error import HTTPError
    calls=[]
    sleeps=[]
    payload=fake_png()
    class FakeResponse:
        def __enter__(self):return self
        def __exit__(self,*args):return False
        def read(self,n):return payload
    def opener(req,timeout):
        calls.append(req.full_url)
        if len(calls)<3:
            raise HTTPError(req.full_url,429,'rate limited',None,None)
        return FakeResponse()
    monkeypatch.setattr(m,'urlopen',opener)
    monkeypatch.setattr(m.time,'sleep',sleeps.append)
    url='https://inaturalist-open-data.s3.amazonaws.com/photos/5555/large.jpeg'
    assert m.image_loader(url)==payload
    assert calls==[url]*3 and sleeps==[1,2]
    calls.clear();sleeps.clear()
    def missing(req,timeout):
        calls.append(req.full_url)
        raise HTTPError(req.full_url,404,'original image missing',None,None)
    monkeypatch.setattr(m,'urlopen',missing)
    with pytest.raises(HTTPError):
        m.image_loader(url)
    assert calls==[url] and not sleeps


def test_production_image_request_spacing_does_not_change_photo_coverage(tmp_path,monkeypatch):
    frame,src=setup_source(tmp_path)
    firewall=tmp_path/'freeze'
    freeze_fresh(frame,src,firewall,shards=2,workers=2,expected=4)
    delays=[]
    monkeypatch.setattr(m,'image_loader',lambda url:fake_png())
    monkeypatch.setattr(m.time,'sleep',delays.append)
    r=m.run_batch(firewall,0,0,3,tmp_path/'sealed',source=tmp_path/'x',
                  detector=tmp_path/'y',sam=tmp_path/'z',evaluator=fake_measure,
                  _test_only_skip_method_check=True)
    assert r['n_terminal_photo_ids']==3
    assert delays==[m.DEFAULT_PROVIDER_PAUSE_SECONDS]*2
    assert r['original_provider_rate_limit_seconds_between_photos']==m.DEFAULT_PROVIDER_PAUSE_SECONDS
    assert r['photo_failures_preserved_not_replaced'] is True

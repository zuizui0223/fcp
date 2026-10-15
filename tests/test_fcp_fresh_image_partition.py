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
            raise m.PermanentOriginalPhotoError('synthetic frozen-photo 404')
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
    monkeypatch.setattr(m,'_LAST_PROVIDER_REQUEST_MONOTONIC',None)
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
    monkeypatch.setattr(m,'_LAST_PROVIDER_REQUEST_MONOTONIC',None)
    clock=[0.0]
    monkeypatch.setattr(m.time,'monotonic',lambda:clock[0])
    def fake_sleep(seconds):
        sleeps.append(seconds)
        clock[0]+=seconds
    monkeypatch.setattr(m.time,'sleep',fake_sleep)
    url='https://inaturalist-open-data.s3.amazonaws.com/photos/5555/large.jpeg'
    assert m.image_loader(url)==payload
    assert calls==[url]*3
    assert sleeps[0]==1
    assert sleeps[1]==pytest.approx(m.DEFAULT_PROVIDER_PAUSE_SECONDS-1)
    assert sleeps[2]==2
    calls.clear();sleeps.clear()
    def missing(req,timeout):
        calls.append(req.full_url)
        raise HTTPError(req.full_url,404,'original image missing',None,None)
    monkeypatch.setattr(m,'urlopen',missing)
    with pytest.raises(HTTPError):
        m.image_loader(url)
    assert calls==[url] and sleeps==pytest.approx([m.DEFAULT_PROVIDER_PAUSE_SECONDS])


def test_real_provider_spacing_across_one_photo_batches(tmp_path,monkeypatch):
    frame,src=setup_source(tmp_path)
    firewall=tmp_path/'freeze'
    freeze_fresh(frame,src,firewall,shards=2,workers=2,expected=4)
    clock=[0.0]
    attempt_times=[]
    delays=[]
    monkeypatch.delenv('FCP_PROVIDER_GLOBAL_RATE_FILE',raising=False)
    monkeypatch.setattr(m,'_LAST_PROVIDER_REQUEST_MONOTONIC',None)
    monkeypatch.setattr(m.time,'monotonic',lambda:clock[0])
    def delay(seconds):
        assert seconds>=0
        delays.append(seconds)
        clock[0]+=seconds
    monkeypatch.setattr(m.time,'sleep',delay)
    payload=fake_png()
    class FakeResponse:
        def __enter__(self):return self
        def __exit__(self,*args):return False
        def read(self,n):return payload
    def opener(req,timeout):
        attempt_times.append(clock[0])
        return FakeResponse()
    monkeypatch.setattr(m,'urlopen',opener)
    # Synthetic source fixtures use non-provider dummy URLs. Route only this
    # test's two images through the REAL S3 loader with a valid fake S3 key.
    real_s3_loader=m.image_loader
    monkeypatch.setattr(m,'image_loader',lambda _url:real_s3_loader(
        'https://inaturalist-open-data.s3.amazonaws.com/photos/17/large.JPG'))
    out=tmp_path/'sealed'
    first=m.run_batch(firewall,0,0,1,out,source=tmp_path/'x',
                  detector=tmp_path/'y',sam=tmp_path/'z',evaluator=fake_measure,
                  _test_only_skip_method_check=True)
    second=m.run_batch(firewall,0,1,1,out,source=tmp_path/'x',
                  detector=tmp_path/'y',sam=tmp_path/'z',evaluator=fake_measure,
                  _test_only_skip_method_check=True)
    assert first['n_terminal_photo_ids']==second['n_terminal_photo_ids']==1
    assert attempt_times==[0.0,m.DEFAULT_PROVIDER_PAUSE_SECONDS]
    assert delays==[m.DEFAULT_PROVIDER_PAUSE_SECONDS]
    assert first['photo_failures_preserved_not_replaced'] is True
    assert first['original_provider_rate_limit_seconds_between_photos']==m.DEFAULT_PROVIDER_PAUSE_SECONDS


def test_global_file_enforces_spacing_even_after_process_state_reset(tmp_path,monkeypatch):
    path=tmp_path/'shared-original-photo-requests.txt'
    clock=[100.0]
    calls=[]
    monkeypatch.setenv('FCP_PROVIDER_GLOBAL_RATE_FILE',str(path))
    monkeypatch.setattr(m.time,'time',lambda:clock[0])
    monkeypatch.setattr(m.time,'sleep',lambda d:clock.__setitem__(0,clock[0]+d))
    monkeypatch.setattr(m,'_LAST_PROVIDER_REQUEST_MONOTONIC',None)
    class FakeResponse:
        def __enter__(self):return self
        def __exit__(self,*args):return False
        def read(self,n):return fake_png()
    def opener(req,timeout):
        calls.append(clock[0])
        return FakeResponse()
    monkeypatch.setattr(m,'urlopen',opener)
    url='https://inaturalist-open-data.s3.amazonaws.com/photos/17/large.jpg'
    m.image_loader(url)
    # Simulates starting another worker process with no in-memory request history.
    monkeypatch.setattr(m,'_LAST_PROVIDER_REQUEST_MONOTONIC',None)
    m.image_loader(url)
    assert calls==[100.0,101.1]
    assert float(path.read_text())==101.1


def test_broken_global_rate_state_is_fatal_not_ordinary_image_failure(tmp_path,monkeypatch):
    monkeypatch.setenv('FCP_PROVIDER_GLOBAL_RATE_FILE','relative-not-trusted')
    url='https://inaturalist-open-data.s3.amazonaws.com/photos/17/large.jpg'
    with pytest.raises(m.ProviderRateGateError,match='absolute'):
        m.image_loader(url)


def test_persistent_rate_limit_aborts_original_photo_work_unit_without_sealing(tmp_path,monkeypatch):
    """A 429 exhaustion is provider throttling, never missing flower phenotype."""
    from urllib.error import HTTPError
    frame,src=setup_source(tmp_path)
    firewall=tmp_path/'frozen'
    freeze_fresh(frame,src,firewall,shards=2,workers=2,expected=4)
    calls=[]
    def denied(req,timeout):
        calls.append(req.full_url)
        raise HTTPError(req.full_url,429,'provider busy',None,None)
    monkeypatch.setattr(m,'urlopen',denied)
    monkeypatch.setattr(m,'_pace_provider_request',lambda:None)
    monkeypatch.setattr(m.time,'sleep',lambda _:None)
    frozen_url='https://inaturalist-open-data.s3.amazonaws.com/photos/17/large.JPG'
    with pytest.raises(m.ProviderTransientAccessError,match='HTTP 429'):
        m.image_loader(frozen_url)
    assert calls==[frozen_url]*3
    calls.clear()
    target=tmp_path/'unsealed'
    with pytest.raises(m.ProviderTransientAccessError,match='HTTP 429'):
        m.run_batch(firewall,0,0,2,target,
                    source=tmp_path/'mock_source.py',detector=tmp_path/'mock_detector.pt',
                    sam=tmp_path/'mock_sam',loader=lambda _:m.image_loader(frozen_url),
                    evaluator=fake_measure,_test_only_skip_method_check=True)
    assert not list(target.glob('terminal_*.csv'))
    assert not list(target.glob('terminal_*.json'))


def test_provider_forbidden_and_network_timeouts_fail_closed(monkeypatch):
    from urllib.error import HTTPError,URLError
    url='https://inaturalist-open-data.s3.amazonaws.com/photos/18/large.jpg'
    monkeypatch.setattr(m,'_pace_provider_request',lambda:None)
    monkeypatch.setattr(m.time,'sleep',lambda _:None)
    calls=[]
    def forbidden(req,timeout):
        calls.append(req.full_url)
        raise HTTPError(req.full_url,403,'access denied',None,None)
    monkeypatch.setattr(m,'urlopen',forbidden)
    with pytest.raises(m.ProviderTransientAccessError,match='HTTP 403'):
        m.image_loader(url)
    assert len(calls)==1
    def unavailable(req,timeout):
        calls.append(req.full_url)
        raise URLError('temporary network outage')
    monkeypatch.setattr(m,'urlopen',unavailable)
    with pytest.raises(m.ProviderTransientAccessError,match='connection unavailable'):
        m.image_loader(url)
    assert len(calls)==4  # One 403 plus three bounded network retries


def test_permanently_missing_original_can_still_be_terminal_without_substitution(tmp_path):
    """404 is source-object missingness, not a global S3 throttle."""
    from urllib.error import HTTPError
    frame,src=setup_source(tmp_path)
    firewall=tmp_path/'frozen'
    freeze_fresh(frame,src,firewall,shards=2,workers=2,expected=4)
    def missing(_url):
        raise m.PermanentOriginalPhotoError('frozen original provider photo HTTP 404')
    target=tmp_path/'sealed'
    result=m.run_batch(firewall,0,0,1,target,source=tmp_path/'source.py',
                       detector=tmp_path/'detector.pt',sam=tmp_path/'sam',
                       loader=missing,evaluator=fake_measure,
                       _test_only_skip_method_check=True)
    assert result['n_terminal_photo_ids']==1
    assert result['n_acquisition_failed']==1
    assert result['n_classifiable']==0
    rows=pd.read_csv(next(target.glob('terminal_*.csv')))
    assert rows.measurement_status.tolist()==['image_acquisition_failed']
    assert rows.morph.tolist()==['mixed_uncertain']



def test_unexpected_local_image_write_error_cannot_become_a_sealed_photo_failure(tmp_path, monkeypatch):
    frame,src=setup_source(tmp_path)
    firewall=tmp_path/'freeze'
    freeze_fresh(frame,src,firewall,shards=2,workers=2,expected=4)
    orig_write_bytes=Path.write_bytes
    def fail_only_original_image(self, payload):
        if self.name.startswith('FCPR-') and self.suffix=='.jpg':
            raise OSError('simulated filesystem full while writing frozen photo')
        return orig_write_bytes(self,payload)
    monkeypatch.setattr(Path,'write_bytes',fail_only_original_image)
    output=tmp_path/'output'
    with pytest.raises(OSError,match='filesystem full'):
        m.run_batch(firewall,0,0,1,output,source=tmp_path/'m.py',
                    detector=tmp_path/'d.pt',sam=tmp_path/'sam',
                    loader=lambda _:fake_png(),evaluator=fake_measure,
                    _test_only_skip_method_check=True)
    assert not list(output.glob('terminal_*.csv'))
    assert not list(output.glob('terminal_*.json'))


def test_real_original_s3_404_yields_explicit_photo_terminal_but_418_aborts(monkeypatch):
    from urllib.error import HTTPError
    monkeypatch.setattr(m,'_pace_provider_request',lambda:None)
    url='https://inaturalist-open-data.s3.amazonaws.com/photos/15/large.JPG'
    def not_found(req,timeout):
        raise HTTPError(req.full_url,404,'not found',None,None)
    monkeypatch.setattr(m,'urlopen',not_found)
    with pytest.raises(m.PermanentOriginalPhotoError,match='HTTP 404'):
        m.image_loader(url)
    def unknown(req,timeout):
        raise HTTPError(req.full_url,418,'unexpected response',None,None)
    monkeypatch.setattr(m,'urlopen',unknown)
    with pytest.raises(m.ProviderTransientAccessError,match='Unexpected provider HTTP 418'):
        m.image_loader(url)


def test_undecodable_frozen_photo_is_explicit_not_a_local_disk_error(monkeypatch):
    monkeypatch.setattr(m,'_pace_provider_request',lambda:None)
    class FakeResponse:
        def __enter__(self):return self
        def __exit__(self,*args):return False
        def read(self,n):return b'not-a-jpeg'
    monkeypatch.setattr(m,'urlopen',lambda req,timeout:FakeResponse())
    url='https://inaturalist-open-data.s3.amazonaws.com/photos/15/large.jpg'
    with pytest.raises(m.PermanentOriginalPhotoError,match='undecodable original photo'):
        m.image_loader(url)


def test_invalid_source_photo_url_is_not_a_per_photo_missing_terminal(tmp_path):
    frame,src=setup_source(tmp_path)
    firewall=tmp_path/'freeze'
    freeze_fresh(frame,src,firewall,shards=2,workers=2,expected=4)
    # This test source deliberately has example.org rather than the
    # only authorized provider. No attempt to download any URL may occur.
    output=tmp_path/'bad_source'
    with pytest.raises(ValueError,match='not a frozen original provider'):
        m.run_batch(firewall,0,0,1,output,
                    source=tmp_path/'m.py',detector=tmp_path/'d.pt',
                    sam=tmp_path/'sam',loader=m.image_loader,
                    evaluator=fake_measure,_test_only_skip_method_check=True)
    assert not list(output.glob('terminal_*.csv'))

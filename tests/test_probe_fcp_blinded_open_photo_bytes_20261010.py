"""Never conflate retrieved pixels with original source equivalence or flower traits."""
import hashlib
import importlib.util
from io import BytesIO
from pathlib import Path
import pytest
from PIL import Image

PATH=Path(__file__).resolve().parents[1]/"scripts/analysis/probe_fcp_blinded_open_photo_bytes_20261010.py"
SPEC=importlib.util.spec_from_file_location("pilot",PATH)
m=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(m)

def jpeg():
    import random
    i=Image.frombytes("RGB",(96,96),random.Random(42).randbytes(96*96*3));b=BytesIO()
    i.save(b,"JPEG",quality=94)
    v=b.getvalue()
    assert len(v)>1024
    return v

class FakeResponse:
    def __init__(self,url,body):self.url=url;self.body=body
    def __enter__(self):return self
    def __exit__(self,*args):return False
    def geturl(self):return self.url
    def read(self,n):return self.body[:n]

def test_candidate_urls_are_constrained_to_public_open_data_and_same_photo_id():
    vals=m.candidate_urls("230101")
    assert len(vals)==6
    assert all(v.startswith(m.HOST+"/photos/230101/") for v in vals)
    with pytest.raises(ValueError):m.candidate_urls("https://example.org/something")

def test_byte_identical_image_is_verified_without_saving_image(tmp_path):
    data=jpeg();h=hashlib.sha256(data).hexdigest()
    fn=lambda req,timeout:FakeResponse(req.full_url,data)
    got,errors=m.acquire_one("230101",h,fn,sleep=lambda _:None)
    assert got["status"]=="byte_match_verified"
    assert got["new_bytes_sha256"]==h
    assert not errors
    assert not list(tmp_path.glob("*.jpg"))

def test_mismatch_is_not_renamed_as_historical_photo():
    data=jpeg();fn=lambda req,timeout:FakeResponse(req.full_url,data)
    got,_=m.acquire_one("230101","0"*64,fn,sleep=lambda _:None)
    assert got["status"]=="photo_id_retrieved_bytes_differ"
    assert got["new_bytes_sha256"]!= "0"*64

def test_redirect_to_nonmatching_host_fails_closed():
    data=jpeg()
    fn=lambda req,timeout:FakeResponse("https://host.invalid/photos/230101/large.jpeg",data)
    got,_=m.acquire_one("230101",hashlib.sha256(data).hexdigest(),fn,sleep=lambda _:None)
    assert got["status"]=="unavailable_on_open_data_host"

def test_nonimage_data_is_not_accepted():
    data=b"X"*5000;fn=lambda req,timeout:FakeResponse(req.full_url,data)
    got,_=m.acquire_one("230101",hashlib.sha256(data).hexdigest(),fn,sleep=lambda _:None)
    assert got["status"]=="unavailable_on_open_data_host"

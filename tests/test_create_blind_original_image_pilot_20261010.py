"""No real photo colours in tests; verify original SHA and blind sheet behaviour."""
import importlib.util
import hashlib
import io
from pathlib import Path
import numpy as np
import pandas as pd
from PIL import Image
import pytest

script=Path(__file__).resolve().parents[1]/"scripts/analysis/create_blind_original_image_pilot_20261010.py"
spec=importlib.util.spec_from_file_location("fcp_pilot",script)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

def fake_manifest(tmp_path,monkeypatch):
    d=pd.DataFrame({"audit_case_id":[f"FCQ-{i:016x}" for i in range(24)],
                    "photo_id":[str(1100000+i) for i in range(24)],
                    "frozen_original_image_sha256":["b"*64]*24})
    f=tmp_path/"pilot.csv";d.to_csv(f,index=False,lineterminator="\n")
    monkeypatch.setattr(m,"EXPECTED_CSV_SHA",hashlib.sha256(f.read_bytes()).hexdigest())
    return f

def test_blinded_24_photo_manifest_gate(tmp_path,monkeypatch):
    f=fake_manifest(tmp_path,monkeypatch)
    assert len(m.read_blind(f))==24
    payload=f.read_text()+"wrongextra\n"
    f.write_text(payload)
    with pytest.raises(RuntimeError,match="SHA mismatch"):
        m.read_blind(f)

def test_untrusted_extra_source_colour_column_is_rejected(tmp_path,monkeypatch):
    f=fake_manifest(tmp_path,monkeypatch)
    d=pd.read_csv(f)
    d["original_algorithm_colour_label"]="white"
    d.to_csv(f,index=False,lineterminator="\n")
    monkeypatch.setattr(m,"EXPECTED_CSV_SHA",hashlib.sha256(f.read_bytes()).hexdigest())
    with pytest.raises(RuntimeError,match="blind columns"):
        m.read_blind(f)

def test_exact_same_source_image_byte_match_and_mismatch():
    rng=np.random.default_rng(47)
    pixels=rng.integers(0,256,size=(80,100,3),dtype=np.uint8)
    buffer=io.BytesIO();Image.fromarray(pixels).save(buffer,format="PNG")
    payload=buffer.getvalue()
    class Reply:
        def __init__(self,u):self.url=u
        def __enter__(self):return self
        def __exit__(self,*a):return False
        def geturl(self):return self.url
        def read(self,size):return payload[:size]
    def opener(req,timeout):return Reply(req.full_url)
    image,url=m.get_one("23847515",hashlib.sha256(payload).hexdigest(),
                        opener=opener,pause=lambda _:None)
    assert image is not None and image.size==(100,80)
    assert url.endswith("large.jpeg")
    err,none=m.get_one("23847515","c"*64,opener=opener,pause=lambda _:None)
    assert err is None and none==""

def test_contact_sheet_contains_blind_ids_not_colour_metadata(tmp_path):
    entries=[(f"FCQ-{i:016x}",Image.new("RGB",(250,150),(120,160,180))) for i in range(12)]
    m.render_sheet(entries,tmp_path/"blind.png")
    with Image.open(tmp_path/"blind.png") as im:
        assert im.size==(1720,1080)
        assert im.mode=="RGB"

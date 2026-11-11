"""405-case SHA/source blinding tests; no internet during synthetic tests."""
import hashlib
import importlib.util
from pathlib import Path
import pandas as pd
import pytest

PATH=Path(__file__).resolve().parents[1]/"scripts/analysis/probe_fcp_original405_blinded_image_bytes_20261010.py"
SPEC=importlib.util.spec_from_file_location("full405",PATH)
m=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(m)

def make_files(tmp_path,*,leak=False):
    b=pd.DataFrame({"audit_case_id":["FCQ-a","FCQ-b"],"photo_id":["11","22"]})
    if leak:b["original_algorithm_colour_label"]=["white","red_pink"]
    k=pd.DataFrame({"audit_case_id":["FCQ-a","FCQ-b"],"photo_id":["11","22"],
                    "image_sha256":["a"*64,"b"*64],"original_algorithm_colour_label":["white","red_pink"]})
    bp=tmp_path/"blind.csv";kp=tmp_path/"key.csv"
    b.to_csv(bp,index=False);k.to_csv(kp,index=False)
    return bp,kp

def patch_expected(monkeypatch,b,k):
    monkeypatch.setattr(m,"N_CASES",2)
    monkeypatch.setattr(m,"EXPECTED_BLIND_SHA",hashlib.sha256(b.read_bytes()).hexdigest())
    monkeypatch.setattr(m,"EXPECTED_KEY_SHA",hashlib.sha256(k.read_bytes()).hexdigest())

def test_sealed_key_not_used_as_image_colour_labels(tmp_path,monkeypatch):
    b,k=make_files(tmp_path);patch_expected(monkeypatch,b,k)
    x=m.manifest(b,k)
    assert set(x.columns)=={"audit_case_id","photo_id","image_sha256"}
    assert len(x)==2

def test_fail_closed_on_reviewer_colour_leak(tmp_path,monkeypatch):
    b,k=make_files(tmp_path,leak=True);patch_expected(monkeypatch,b,k)
    with pytest.raises(RuntimeError,match="Blinding"):m.manifest(b,k)

def test_sha_drift_fails_closed(tmp_path,monkeypatch):
    b,k=make_files(tmp_path);patch_expected(monkeypatch,b,k)
    b.write_text(b.read_text().replace("FCQ-a","FCQ-modified"))
    with pytest.raises(RuntimeError,match="SHA"):m.manifest(b,k)

def test_unrecovered_count_not_converted_to_colour_error(tmp_path,monkeypatch):
    b,k=make_files(tmp_path);patch_expected(monkeypatch,b,k)
    def fake(pid,h):
        if pid=="11":return {"status":"byte_match_verified"},[]
        return {"status":"unavailable_on_open_data_host"},[]
    z=m.run(b,k,tmp_path/"out",acquirer=fake,sleep=lambda _:None)
    assert z["n_exact_historical_image_bytes"]==1
    assert z["n_unavailable_on_open_data_host"]==1
    assert z["label_reannotation_completed"] is False

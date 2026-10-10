"""Second-stage original-target-taxon disclosure test, without source colour leakage."""
from pathlib import Path
from io import BytesIO
from zipfile import ZipFile,ZIP_DEFLATED
import hashlib
import importlib.util
import pandas as pd
import pytest

f=Path(__file__).resolve().parents[1]/"scripts/analysis/build_fcp_stage2_taxon_only_review_20261010.py"
s=importlib.util.spec_from_file_location("taxon_gate",f)
m=importlib.util.module_from_spec(s);s.loader.exec_module(m)


def make_source(tmp_path,monkeypatch):
    n=405
    k=pd.DataFrame({
        "audit_case_id":[f"FCQ-{i:016x}" for i in range(n)],
        "photo_id":[str(500000+i) for i in range(n)],
        "species":[f"Species_{i%169}" for i in range(n)],
        "original_algorithm_colour_label":["white"]*n,
        "kind":["high_leverage"]*n,
        "cohort":["third"]*n,
        "latitude":["21.01"]*n,
        "longitude":["132.01"]*n,
    })
    b=k[["audit_case_id","photo_id"]]
    key=k.to_csv(index=False,lineterminator="\n").encode()
    blind=b.to_csv(index=False,lineterminator="\n").encode()
    monkeypatch.setattr(m,"KEY_SHA",hashlib.sha256(key).hexdigest())
    monkeypatch.setattr(m,"BLIND_SHA",hashlib.sha256(blind).hexdigest())
    archive=tmp_path/"source.zip"
    with ZipFile(archive,"w",compression=ZIP_DEFLATED) as z:
        z.writestr("photo_review/UNBLINDING_KEY_do_not_show_reviewers.csv",key)
        z.writestr("photo_review/BLINDED_reannotation_photo_queue.csv",blind)
    return archive


def test_all_405_target_taxa_recovered_with_colour_and_source_group_hidden(tmp_path,monkeypatch):
    source=make_source(tmp_path,monkeypatch)
    out=tmp_path/"STAGE2_OPEN_AFTER_FIRSTPASS.csv"
    r=m.build(source,out)
    assert r["n_cases"]==405
    assert r["n_distinct_target_taxa"]==169
    assert r["source_colour_and_design_groups_sealed"] is True
    assert r["human_reviews_claimed_by_this_file"]==0
    d=pd.read_csv(out,dtype=str)
    assert list(d.columns)==["audit_case_id","photo_id","focal_taxon_name_no_prior_colour"]
    assert not {"original_algorithm_colour_label","kind","cohort","latitude"}.intersection(d.columns)


def test_changed_source_key_sha_rejected(tmp_path,monkeypatch):
    source=make_source(tmp_path,monkeypatch)
    monkeypatch.setattr(m,"KEY_SHA","b"*64)
    with pytest.raises(RuntimeError,match="SHA seal mismatch"):
        m.build(source,tmp_path/"must_not_exist.csv")
    assert not (tmp_path/"must_not_exist.csv").exists()

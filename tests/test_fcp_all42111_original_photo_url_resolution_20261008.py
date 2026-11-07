"""Source batch observations/URLs from all FCP 42111 historical photo IDs, not pixels."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd
import pytest

PATH=Path(__file__).resolve().parents[1]/"scripts/acquisition/resolve_fcp_all42111_existing_photo_metadata_20261008.py"
sp=importlib.util.spec_from_file_location("resolve42111",PATH)
m=importlib.util.module_from_spec(sp)
sp.loader.exec_module(m)


def source():
    return pd.DataFrame([
        {"inat_taxon_id":1,"species":"Iris albus","observation_id":101,"photo_id":1001,"after_observer_cap":0},
        {"inat_taxon_id":2,"species":"Viola amara","observation_id":202,"photo_id":2002,"after_observer_cap":30},
        {"inat_taxon_id":3,"species":"Salvia bravo","observation_id":303,"photo_id":3003,"after_observer_cap":100}
    ])


def response(photo3=True,licensed3=True):
    return {"results":[
      {"id":101,"taxon":{"id":1},"photos":[{"id":1001,"url":"https://example.org/p/1001/square.jpg","license_code":"cc-by"}]},
      {"id":202,"taxon":{"id":2},"photos":[{"id":999,"url":"https://example.org/p/999/medium.jpg","license_code":"cc-by"}]},
      {"id":303,"taxon":{"id":30},"photos":[{"id":3003,"url":"https://example.org/p/3003/small.png","license_code":"cc-by" if licensed3 else "all-rights-reserved"}] if photo3 else []},
    ]}


def test_selected_exact_photo_not_any_other_image_and_retains_capacity_zero():
    out=m.decode_batch(source(),response())
    assert len(out)==3
    a,b,c=out
    assert a["photo_url_status"]=="VALID_SOURCE_PHOTO_URL_AND_LICENSE"
    assert a["photo_url_large"]=="https://example.org/p/1001/large.jpg"
    assert a["historical_source_photo_capacity"]==0
    assert a["source_taxon_identity_match"] is True
    assert b["photo_url_status"]=="ORIGINAL_PHOTO_ID_MISSING"
    assert b["photo_url_large"]==""
    assert c["photo_url_status"]=="TAXON_IDENTITY_MISMATCH"
    assert c["source_taxon_identity_match"] is False
    assert c["photo_url_large"]==""


def test_relicensing_blocks_old_image_and_does_not_relabel_colour():
    out=m.decode_batch(source(),response(licensed3=False))
    assert out[2]["photo_url_status"]=="TAXON_IDENTITY_MISMATCH"
    assert out[2]["photo_url_large"]==""
    assert "flower_colour" not in out[2]
    assert "morph" not in out[2]


def test_missing_observation_becomes_missing_state_not_monochromatic():
    r=response()
    r["results"]=r["results"][1:]
    out=m.decode_batch(source(),r)
    assert out[0]["photo_url_status"]=="OBSERVATION_ID_NOT_RETURNED"
    assert out[0]["photo_url_large"]==""


def test_duplicate_observation_id_cannot_be_accepted():
    r=response()
    r["results"].append(r["results"][0])
    with pytest.raises(ValueError,match="Duplicate observation"):
        m.decode_batch(source(),r)


def test_batch_identity_and_api_rate_policy_are_fixed(monkeypatch):
    monkeypatch.setattr(m,"EXPECTED_SPECIES",3)
    monkeypatch.setattr(m,"IDS_PER_REQUEST",2)
    calls=[]
    def fake(ids):
        calls.append(ids)
        valid=[r for r in response()["results"] if int(r["id"]) in ids]
        return {"results":valid}
    rows,summary=m.resolve_all(source(),fake,pause=lambda _:None)
    assert calls==[[101,202],[303]]
    assert len(rows)==3
    assert summary["n_observation_batches"]==2
    assert summary["n_url_and_license_verified"]==1
    assert summary["n_url_and_license_verified_original_capacity_zero"]==1
    assert summary["n_name_id_mismatches_in_source_taxonomy"]==1
    assert rows.inat_taxon_id.nunique()==3


def test_failed_api_batch_preserves_all_original_source_taxa(monkeypatch):
    monkeypatch.setattr(m,"EXPECTED_SPECIES",3)
    monkeypatch.setattr(m,"IDS_PER_REQUEST",2)
    def bad(ids):
        if ids[0]==101:raise ConnectionError("synthetic 429")
        return {"results":[x for x in response()["results"] if int(x["id"]) in ids]}
    rows,summary=m.resolve_all(source(),bad,pause=lambda _:None)
    assert rows.inat_taxon_id.nunique()==3
    assert summary["n_failed_batches"]==1
    assert summary["species_by_current_original_photo_URL_status"]["SOURCE_API_QUERY_FAILURE"]==2
    assert summary["n_url_and_license_verified"]==0


def test_NO_API_interaction_in_photo_identity_validation():
    assert m.EXPECTED_SPECIES==42111
    assert m.IDS_PER_REQUEST==100
    assert m.REQUEST_START_MIN_INTERVAL_SECONDS>=1
    assert m.RETRIES<=1
    assert "urlopen" not in m.original_photo_manifest.__code__.co_names


def test_path_endpoint_over_30_id_recovery_uses_search_query(monkeypatch):
    """Never submit 100 comma IDs to the broken /observations/{id} route."""
    from urllib.parse import parse_qs, urlparse
    from contextlib import contextmanager

    seen = []
    class Response:
        def read(self):
            return b'{"results":[]}'
    @contextmanager
    def fake_urlopen(request,timeout):
        seen.append(request.full_url)
        yield Response()

    monkeypatch.setattr(m, "urlopen", fake_urlopen)
    j=m.fetch_api_batch(list(range(100000,100100)),rate_sleep=lambda _:None)
    assert j == {"results":[]}
    assert len(seen)==1
    parsed=urlparse(seen[0])
    assert parsed.path == "/v1/observations"
    qs=parse_qs(parsed.query)
    assert qs["per_page"] == ["200"]
    assert len(qs["id"][0].split(","))==100
    assert "/observations/100000" not in seen[0]


def test_taxon_mismatch_has_no_verifiable_photo_url():
    g=source()
    out=m.decode_batch(g,response())
    c=out[2]
    assert c["photo_url_status"]=="TAXON_IDENTITY_MISMATCH"
    assert c["photo_url_large"] == ""

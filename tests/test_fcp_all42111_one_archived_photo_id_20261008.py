"""Entire 42111 global species photo identity sample is metadata only."""
import importlib.util
from pathlib import Path
import pandas as pd
import pytest

SRC=Path(__file__).resolve().parents[1]/"scripts/analysis/freeze_fcp_all42111_one_archived_photo_id_20261008.py"
sp=importlib.util.spec_from_file_location("archived_image_identity",SRC)
m=importlib.util.module_from_spec(sp)
sp.loader.exec_module(m)


def frames():
    species=pd.DataFrame({
        "inat_taxon_id":[1,2,3],
        "species":["Iris albus","Viola amara","Salvia bravo"],
        "after_observer_cap":[0,5,100]
    })
    pics=pd.DataFrame({
        "inat_taxon_id":[1,2,2,3,3],
        "observation_id":[101,101,202,303,304],
        "photo_id":[1101,1101,2202,3303,4404],
        "cell_id":[0,0,1,12,14],
        "species":["Iris albus","Viola amara","Viola amara",
                   "Salvia bravo","Salvia bravo"],
        "source_scan":["v1","v1","v2","v2","v1"]
    })
    return species,pics


def test_full_species_photo_identity_even_if_capacity_scan_zero(monkeypatch):
    monkeypatch.setattr(m,"N_WORLD_SPECIES",3)
    s,p=frames()
    final,report=m.freeze_photo_identifiers(s,p)
    assert len(final)==3
    assert final.inat_taxon_id.nunique()==3
    assert final.photo_id.nunique()==3
    assert final.observation_id.nunique()==3
    assert report["n_still_no_eligible_photo_by_later_capacity_scan"]==1
    assert report["n_original_archived_photo_ids_even_when_later_capacity_zero"]==1
    assert not final.white_or_pigment_state_measured.any()
    assert not final.original_photo_URL_verified.any()
    assert int(final.loc[final.inat_taxon_id==2,"photo_id"].iloc[0])==2202


def test_duplicated_archive_rows_are_deduplicated(monkeypatch):
    monkeypatch.setattr(m,"N_WORLD_SPECIES",3)
    s,p=frames()
    p=pd.concat([p,p.iloc[[0,1]]],ignore_index=True)
    x,report=m.freeze_photo_identifiers(s,p)
    assert len(x)==3
    assert report["n_distinct_archived_photo_IDs"]==3


def test_missing_one_species_photo_causes_hard_stop(monkeypatch):
    monkeypatch.setattr(m,"N_WORLD_SPECIES",3)
    s,p=frames()
    p=p.loc[p.inat_taxon_id!=3]
    with pytest.raises(RuntimeError,match="not every 42111|n_missing"):
        m.freeze_photo_identifiers(s,p)


def test_genuinely_impossible_global_unique_photo_id_is_rejected(monkeypatch):
    monkeypatch.setattr(m,"N_WORLD_SPECIES",3)
    s,p=frames()
    p=p.loc[p.photo_id!=2202]
    with pytest.raises(RuntimeError,match="Cannot make all historical"):
        m.freeze_photo_identifiers(s,p)


def test_invalid_geographic_cell_photo_is_rejected(monkeypatch):
    monkeypatch.setattr(m,"N_WORLD_SPECIES",3)
    s,p=frames()
    p.loc[0,"cell_id"]=162
    with pytest.raises(RuntimeError,match="fixed global equal-area"):
        m.freeze_photo_identifiers(s,p)


def test_image_files_not_touched_by_photo_id_selection():
    assert "requests" not in SRC.read_text()
    assert "photo_url_large" not in SRC.read_text()
    assert len(m.SALT)>25

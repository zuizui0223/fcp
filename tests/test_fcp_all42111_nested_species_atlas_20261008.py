"""Protect all 42,111 taxa and distinguish colour-photo opportunity from FCP."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

SOURCE=Path(__file__).resolve().parents[1]/"scripts/analysis/build_fcp_all42111_nested_species_atlas_20261008.py"
sp=importlib.util.spec_from_file_location("all42111_fcp",SOURCE)
m=importlib.util.module_from_spec(sp)
sp.loader.exec_module(m)


def source_frames():
    # Exact historic capacity frequencies, with all taxa, NOT just U100.
    levels=[(0,789),(1,7378),(2,9332),(5,6311),(10,5316),
            (20,2655),(30,1577),(40,1121),(50,862),(60,1212),
            (80,828),(100,4730)]
    depths=np.concatenate([np.full(n,value,dtype=int) for value,n in levels])
    assert len(depths)==42111
    ids=np.arange(1,42112,dtype=int)
    sp=pd.DataFrame({
        "inat_taxon_id":ids,"species":[f"Genus{i} example" for i in ids]
    })
    cap=pd.DataFrame({"inat_taxon_id":ids,
                      "after_observer_cap":depths,
                      "maximum_span_km":np.zeros(len(ids)),
                      "request_error":[""]*len(ids)})
    links1=pd.DataFrame({"inat_taxon_id":ids,"cell_id":ids%162})
    links2=pd.DataFrame({"inat_taxon_id":ids,"cell_id":ids%162})
    return sp,cap,links1,links2


def save_sources(tmp_path,frames):
    paths=[tmp_path/f"source_{i}.csv.gz" for i in range(4)]
    for p,data in zip(paths,frames):
        data.to_csv(p,index=False,compression="gzip")
    return paths


def test_exact_original_all42111_availability_thresholds():
    frames=source_frames()
    a=frames[1]
    assert len(a)==42111
    for min_n,total in m.CAPACITY_EXPECTED.items():
        assert int(a.after_observer_cap.ge(min_n).sum())==total
    assert int(a.after_observer_cap.eq(0).sum())==789
    assert m.availability_band(0)=="0_no_metadata_eligible_photos"
    assert m.availability_band(1)=="1_photo"
    assert m.availability_band(100)=="100_plus_photos"


def test_all_42111_retained_even_without_any_photo(tmp_path,monkeypatch):
    inputs=save_sources(tmp_path,source_frames())
    monkeypatch.setattr(m,"checked",lambda path,fingerprint:None)
    x,info=m.load_universe(*inputs)
    assert len(x)==42111
    assert x.inat_taxon_id.nunique()==42111
    assert int(x.observed_no_eligible_photo_metadata.sum())==789
    assert int(x.potential_40plus_ITV_measurement.sum())==8753
    assert int(x.potential_100plus_high_depth_geography.sum())==4730
    assert not x.genetic_FCP_polymorphism_validated.any()
    assert not x.white_achromatic_state_validated.any()
    assert x.atlas_label_status.eq("NOT_MEASURED_FROM_THIS_METADATA").all()
    assert info["photo_depth_threshold_counts"]["20"]==12985
    assert info["n_world_equal_area_cells"]==162
    assert info["species_geographic_cell_counts"]["n_2plus_cells"]==0
    for quota in m.PHOTO_QUOTAS:
        bound=int(x[f"max_photo_rows_at_cap_{quota}"].sum())
        assert info["photographic_preprocessing_plan_by_max_per_species"][str(quota)]["sum_of_possible_nested_photo_rows"]==bound
        assert bound>=0
        assert int(x[f"max_photo_rows_at_cap_{quota}"].max())<=quota
    assert x.loc[x.after_observer_cap==0,"source_only_descriptive_photo_possible"].eq(False).all()


def test_taxon_and_cell_index_uniqueness_fails_closed(tmp_path,monkeypatch):
    sp,cap,a,b=source_frames()
    sp.loc[0,"inat_taxon_id"]=sp.loc[1,"inat_taxon_id"]
    monkeypatch.setattr(m,"checked",lambda path,fingerprint:None)
    with pytest.raises(RuntimeError,match="Repeated or missing original taxon"):
        m.load_universe(*save_sources(tmp_path,(sp,cap,a,b)))


def test_empty_cell_geography_is_not_accepted(tmp_path,monkeypatch):
    sp,cap,a,b=source_frames()
    a=a.loc[a.inat_taxon_id!=1]
    b=b.loc[b.inat_taxon_id!=1]
    monkeypatch.setattr(m,"checked",lambda path,fingerprint:None)
    with pytest.raises(RuntimeError,match="no historically recorded geographic cell"):
        m.load_universe(*save_sources(tmp_path,(sp,cap,a,b)))


def test_wrong_historical_source_bytes_are_rejected_before_data_usage(tmp_path):
    p=tmp_path/"wrong.csv"
    p.write_text("inat_taxon_id,species\n1,Wrong example\n")
    with pytest.raises(RuntimeError,match="SHA256 mismatch"):
        m.checked(p,m.SPECIES_SHA)


def test_one_photo_never_misreported_as_polymorphism(tmp_path,monkeypatch):
    monkeypatch.setattr(m,"checked",lambda path,fingerprint:None)
    x,_=m.load_universe(*save_sources(tmp_path,source_frames()))
    only=x.loc[x.after_observer_cap==1]
    assert len(only)==7378
    assert only.potential_visible_colour_variation_screen.eq(False).all()
    assert only.potential_40plus_ITV_measurement.eq(False).all()
    assert only.source_only_descriptive_photo_possible.eq(True).all()

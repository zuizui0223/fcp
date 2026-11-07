"""Protect whole-4,730 study's new 2,730 species metadata acquisition."""
import importlib.util
from pathlib import Path
import sys

import pandas as pd
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts/acquisition"))
a=importlib.util.spec_from_file_location("all4730_acquire",ROOT/"scripts/acquisition/fcp_all4730_new2730_metadata_shard_20261008.py")
m=importlib.util.module_from_spec(a)
a.loader.exec_module(m)
s=importlib.util.spec_from_file_location("all4730_seal",ROOT/"scripts/analysis/seal_fcp_all4730_new2730_metadata_20261008.py")
n=importlib.util.module_from_spec(s)
s.loader.exec_module(n)


def test_frozen_2730_photo_taxa_partition_exactly_once():
    d=ROOT/"results/fcp_next2000_prospective_preopening_20261008"
    pool=m.read_allocation(d/"prospective_selected_2000_species.tsv",
                           d/"prospective_reserve_730_species.tsv")
    assert len(pool)==2730
    assert pool.inat_taxon_id.nunique()==2730
    shards=[m.shard_frame(pool,j) for j in range(20)]
    ids=pd.concat(shards).prospective_rank.astype(int)
    assert ids.nunique()==2730
    assert sorted(ids)==list(range(1,2731))
    assert all(len(z) in (136,137) for z in shards)
    assert pool.groupby("allocation_group").size().to_dict()=={
        "new_primary":2000,"new_replication":730}


def test_shard_budget_cannot_be_optimized_after_access():
    pool=pd.DataFrame({"prospective_rank":[1,2,3]})
    with pytest.raises(ValueError,match="Shard identity"):
        m.shard_frame(pool,0,12)
    with pytest.raises(ValueError,match="Shard identity"):
        m.shard_frame(pool,30,20)


def test_metadata_query_keeps_no_flower_colour_columns():
    from fcp_pipeline.random_photo_h9_pool import h9_query_for_species
    params=h9_query_for_species(12345,per_page=m.PER_PAGE,
         maximum_positional_accuracy_m=m.MAX_ACCURACY_M,
         allowed_photo_licenses=m.ALLOWED_LICENSES)
    assert params["per_page"]==200
    assert params["page"]==1
    assert params["order_by"]=="random"
    assert params["quality_grade"]=="research"
    assert params["term_id"]==12 and params["term_value_id"]==13
    assert "morph" not in params and "pigment" not in params
    assert "latitude" not in params and "longitude" not in params


def mock_audit():
    records=[]
    for i in range(1,2731):
        records.append({"inat_taxon_id":i,"prospective_rank":i,
                        "allocation_group":"new_primary" if i<=2000 else "new_replication",
                        "retained":0,"request_error":""})
    return pd.DataFrame(records)


def test_complete_attempt_denominator_with_zero_photos_is_not_a_biological_null():
    a=mock_audit()
    f=pd.DataFrame()
    r=[{"request_errors":0}]*20
    x=n.verify_unique_opportunities(a,f,r)
    assert x["n_species_attempted"]==2730
    assert x["n_primary_full100"]==0 and x["n_replication_full100"]==0
    assert x["n_metadata_photo_records"]==0
    assert not x["met_capacity_and_transport_gate"]


def test_duplicate_or_missing_taxon_causes_hard_failure():
    a=mock_audit()
    a.loc[0,"inat_taxon_id"]=a.loc[1,"inat_taxon_id"]
    with pytest.raises(RuntimeError,match="exactly once"):
        n.verify_unique_opportunities(a,pd.DataFrame(),[{"request_errors":0}]*20)


def test_old_photo_exclusion_inputs_cannot_be_omitted():
    assert len(m.OLD_EXCLUSIONS)==5
    assert m.PREVIOUS_EXCLUSION_UNION==228461
    for value in m.OLD_EXCLUSIONS.values():
        assert len(value)==64
    with pytest.raises(RuntimeError,match="unavailable"):
        m.exclusions(Path("/tmp/fcp_test_dir_that_does_not_exist"))


def test_all_successful_photo_metadata_gate_constrains_primary_and_replication():
    a=mock_audit()
    a.loc[a.prospective_rank<=1100,"retained"]=100
    a.loc[a.prospective_rank>2000,"retained"]=100
    x=n.verify_unique_opportunities(a,pd.DataFrame(),[{"request_errors":0}]*20)
    assert x["n_primary_full100"]==1100
    assert x["n_replication_full100"]==730
    assert x["met_capacity_and_transport_gate"] is True

"""Small synthetic *identity* tests. No florals or ecological claims."""
import pandas as pd
import pytest
from scripts.analysis import fcp_u100_4730_five_cohort_source_id_audit_20261017 as m


def make_df(taxon:int,first_photo:int,first_obs:int):
    return pd.DataFrame({
        "photo_id":list(range(first_photo,first_photo+100)),
        "observation_id":list(range(first_obs,first_obs+100)),
        "inat_taxon_id":[taxon]*100,
        "species":[f"Synthetic_{taxon}"]*100,
    })


def five_synthetic(monkeypatch):
    data={}
    for i,name in enumerate(m.COHORTS):
        df=make_df(i+1,i*1000+10,i*1000+100010)
        if name=="new_U100":
            # One photo and one observation reused in discovery.
            df.loc[0,"photo_id"]=10
            df.loc[0,"observation_id"]=100010
        data[name]=df
    audit={}
    for key,df in data.items():
        audit[key]=m.checked_source(df,key,n=100,taxa=1)
    monkeypatch.setattr(m,"FULL_U100",7)
    monkeypatch.setattr(m,"AUTHORIZED",5)
    monkeypatch.setattr(m,"UNMET",2)
    monkeypatch.setattr(m,"FULL_PHOTOS",500)
    monkeypatch.setattr(m,"FINAL_EFFECTIVE_TAXA",4)
    monkeypatch.setattr(m,"FINAL_EFFECTIVE_PHOTOS",400)
    monkeypatch.setattr(m,"EXCLUDED_NEW_TAXA",frozenset({5}))
    monkeypatch.setattr(m,"PINNED_OLD_DISCOVERY_PHOTO_COLLISIONS",frozenset({10}))
    monkeypatch.setattr(m,"PINNED_OLD_DISCOVERY_OBSERVATION_COLLISIONS",frozenset({100010}))
    monkeypatch.setattr(m,"ACTUAL",{k:{"n":100,"taxa":1} for k in m.COHORTS})
    full=pd.DataFrame({
        "inat_taxon_id":list(range(1,8)),
        "species":[f"Synthetic_{i}" for i in range(1,8)],
        "final_status":["original_149900_photo_measured"]*3+
                       ["wave4_49800_photo_id_authorized",
                        "new_metadata_100_of_100_photo_IDs"]+
                       ["new_H9_metadata_under_100_or_request_error"]*2,
    })
    return audit,full


def test_real_five_source_collision_rule_keeps_four_complete_100photo_species(monkeypatch):
    audit,full=five_synthetic(monkeypatch)
    r=m.inspect_all(audit,full)
    assert r["raw_selected_original_source_photo_ID_occurrences"]==500
    assert r["globally_distinct_original_source_photo_IDs"]==499
    assert r["technically_excluded_entire_new_100photo_taxa"]==[5]
    assert r["original_source_100photo_taxa_retain_for_analysis"]==4
    assert r["original_source_100photo_rows_retain_for_analysis"]==400
    assert r["not_photo_colour_or_ecological_statistical_inference"] is True


def test_reject_new_unexpected_reused_photo_even_when_same_total_taxa(monkeypatch):
    audit,full=five_synthetic(monkeypatch)
    audit["new_U100"]["photographs"].add(11)
    audit["new_U100"]["photo_to_taxon"][11]=5
    with pytest.raises(RuntimeError,match="differs from frozen evidence"):
        m.inspect_all(audit,full)


def test_reject_unexpected_original_validation_photo_reuse(monkeypatch):
    audit,full=five_synthetic(monkeypatch)
    audit["validation"]["photographs"].add(10)
    with pytest.raises(RuntimeError,match="earlier frozen"):
        m.inspect_all(audit,full)


def test_reject_missing_original_100photo_u100_taxonomy_status(monkeypatch):
    audit,full=five_synthetic(monkeypatch)
    full.loc[6,"final_status"]="new_metadata_100_of_100_photo_IDs"
    with pytest.raises(RuntimeError,match="status|ledger"):
        m.inspect_all(audit,full)


def test_reject_duplicate_source_photo_id_within_original_cohort():
    df=make_df(1,100,10000)
    df.loc[0,"photo_id"]=df.loc[1,"photo_id"]
    with pytest.raises(RuntimeError,match="row count/ID"):
        m.checked_source(df,"discovery",n=100,taxa=1)


def test_reject_under_100_original_species_and_unexpected_observation():
    df=make_df(1,100,10000)
    with pytest.raises(RuntimeError,match="100 photos"):
        m.checked_source(df.iloc[:99],"discovery",n=99,taxa=1)
    df.loc[0,"observation_id"]=df.loc[1,"observation_id"]
    with pytest.raises(RuntimeError,match="duplicated original iNaturalist observation"):
        m.checked_source(df,"discovery",n=100,taxa=1)


def test_reject_reused_photo_with_unmatched_old_new_taxon_id(monkeypatch):
    audit,full=five_synthetic(monkeypatch)
    audit["new_U100"]["photo_to_taxon"][10]=900
    with pytest.raises(RuntimeError,match="not correctly identified"):
        m.inspect_all(audit,full)


def test_source_photo_taxon_must_match_original_u100_binomial(monkeypatch):
    audit,full=five_synthetic(monkeypatch)
    full.loc[4,"species"]="different_biological_name"
    with pytest.raises(RuntimeError,match="taxon names differ"):
        m.inspect_all(audit,full)

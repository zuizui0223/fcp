"""Synthetic 100-photo source-ID assertions; no actual images or p-values."""
import pandas as pd
import pytest
from scripts.analysis import fcp_u100_4730_five_cohort_source_id_audit_20261017 as m


def test_input_cohort_exact_100_unique_photo_ids_and_taxa():
    df=pd.DataFrame({
        "photo_id":list(range(1,101))+list(range(101,201)),
        "observation_id":list(range(1001,1101))+list(range(1101,1201)),
        "inat_taxon_id":[11]*100+[22]*100,
        "species":["Aster alpha"]*100+["Aster beta"]*100,
    })
    valid=m.checked_source(df,"synthetic",n=200,taxa=2)
    assert valid["photo_rows"]==200
    assert valid["taxa"]==2
    assert len(valid["photographs"])==200
    assert len(valid["observations"])==200
    df.loc[0,"photo_id"]=df.loc[1,"photo_id"]
    with pytest.raises(RuntimeError,match="photographic row count/ID identity drift"):
        m.checked_source(df,"synthetic",n=200,taxa=2)


def fixture(monkeypatch):
    groups=m.COHORTS
    monkeypatch.setattr(m,"FULL_U100",11)
    monkeypatch.setattr(m,"AUTHORIZED",10)
    monkeypatch.setattr(m,"UNMET",1)
    monkeypatch.setattr(m,"FULL_PHOTOS",1000)
    monkeypatch.setattr(m,"ACTUAL",{g:{"n":200,"taxa":2} for g in groups})
    seen={}
    manifest=[]
    for i,g in enumerate(groups):
        taxa=[i*2+1,i*2+2]
        ids=list(range(200*i+1,200*i+201))
        table=pd.DataFrame({
            "photo_id":ids,
            "observation_id":[j+10000 for j in ids],
            "inat_taxon_id":[taxa[0]]*100+[taxa[1]]*100,
            "species":[f"Testus species{i*2+1}"]*100+[f"Testus species{i*2+2}"]*100,
        })
        seen[g]=m.checked_source(table,g,n=200,taxa=2)
        states=(["original_149900_photo_measured"] if i<3
                else ["wave4_49800_photo_id_authorized"] if i==3
                else ["new_metadata_100_of_100_photo_IDs"])
        manifest.extend({"inat_taxon_id":tid,"species":f"Testus species{tid}",
                         "final_status":states[0]} for tid in taxa)
    manifest.append({"inat_taxon_id":11,"species":"Testus species11",
                     "final_status":"new_H9_metadata_under_100_or_request_error"})
    return seen,pd.DataFrame(manifest)


def test_five_cohort_real_original_denominators_reconcile_by_unique_id(monkeypatch):
    d,full=fixture(monkeypatch)
    audit=m.inspect_all(d,full)
    assert audit["actual_original_source_photo_IDs"]==1000
    assert audit["all_five_cohort_original_complete100_photo_taxa"]==10
    assert audit["original_eligible_taxa_below_historical_fixed100"]==1
    assert audit["cross_five_cohort_real_photo_ID_overlaps"]==0
    assert audit["not_biological_inference"] is True


def test_cross_cohort_photo_reuse_cannot_pass(monkeypatch):
    d,full=fixture(monkeypatch)
    d["wave4"]["photographs"].remove(601)
    d["wave4"]["photographs"].add(1)
    with pytest.raises(RuntimeError,match="reused across"):
        m.inspect_all(d,full)


def test_cross_cohort_species_or_taxa_reuse_cannot_pass(monkeypatch):
    d,full=fixture(monkeypatch)
    d["new_U100"]["taxon_ids"].remove(9)
    d["new_U100"]["taxon_ids"].add(1)
    with pytest.raises(RuntimeError,match="reused across"):
        m.inspect_all(d,full)


def test_cross_cohort_observation_reuse_cannot_pass(monkeypatch):
    d,full=fixture(monkeypatch)
    d["new_U100"]["observations"].remove(10801)
    d["new_U100"]["observations"].add(10001)
    with pytest.raises(RuntimeError,match="observation reused"):
        m.inspect_all(d,full)


def test_original_u100_source_name_change_cannot_pass(monkeypatch):
    d,full=fixture(monkeypatch)
    full.loc[full.inat_taxon_id.eq(5),"species"]="Aster invented"
    with pytest.raises(RuntimeError,match="taxon-name relation"):
        m.inspect_all(d,full)


def test_original_u100_under100_taxon_must_not_be_imputed(monkeypatch):
    d,full=fixture(monkeypatch)
    full.loc[full.inat_taxon_id.eq(10),"final_status"]="new_H9_metadata_under_100_or_request_error"
    full.loc[full.inat_taxon_id.eq(11),"final_status"]="new_metadata_100_of_100_photo_IDs"
    with pytest.raises(RuntimeError,match="original less-than-100 photo taxon imputed"):
        m.inspect_all(d,full)

"""Fail-closed 60 Validation gap-two species metadata tests: zero live API calls."""
from __future__ import annotations
from types import SimpleNamespace
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"/"acquisition"))
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"/"analysis"))
from audit_fcp_validation_gap2_live_metadata_pilot_20261009 import (
    anchor_context_gap2,check_queue,decide_species
)

def requirement(year,required,found,error=False,trunc=False):
    return {"target_year":year,"required":required,"found":found,
            "api_error":error,"pagination_incomplete":trunc}

def test_metadata_query_envelope_matches_aggregation_schema():
    row=requirement(2021,1,1)
    assert "target_year" in row and "year" not in row
    assert decide_species([row,requirement(2023,1,1)])=="POSSIBLE_COMPLETE_METADATA_ONLY"

def test_both_years_need_two_distinct_observer_candidates():
    assert decide_species([requirement(2021,1,1),requirement(2023,1,1)])=="POSSIBLE_COMPLETE_METADATA_ONLY"
    assert decide_species([requirement(2021,1,1),requirement(2023,1,0)])=="NO_COMPLETE_METADATA_CANDIDATE_FIRST_PAGES"

def test_two_observers_needed_in_one_year():
    assert decide_species([requirement(2023,2,2)])=="POSSIBLE_COMPLETE_METADATA_ONLY"
    assert decide_species([requirement(2023,2,1)])=="NO_COMPLETE_METADATA_CANDIDATE_FIRST_PAGES"

def test_errors_and_paginated_unknown_do_not_turn_into_negative_or_positive():
    assert decide_species([requirement(2021,1,0,error=True),requirement(2023,1,1)])=="API_ERROR_UNRESOLVED"
    assert decide_species([requirement(2021,1,0,trunc=True),requirement(2023,1,1)])=="PAGINATION_INCOMPLETE_UNKNOWN"

def test_rejects_duplicate_year_or_missing_requirement():
    with pytest.raises(ValueError,match="Repeated"):
        decide_species([requirement(2021,1,1),requirement(2021,1,1)])
    with pytest.raises(ValueError,match="frozen two-slot"):
        decide_species([requirement(2021,1,1)])

def frame():
    records=[]
    for i in range(60):
        n=1 if i%2 else 2
        for year in ([2020,2022] if n==2 else [2022]):
            records.append({
                "cohort":"validation","priority_tier":2,
                "inat_taxon_id":f"{1000+i}","species":f"Taxon name{i}",
                "calendar_month":6,"target_year":year,
                "source_anchor_photo_id":f"photo{i}",
                "existing_distinct_observers_in_target_year":1 if n==2 else 0,
                "required_additional_distinct_observer_photos":1 if n==2 else 2,
                "live_metadata_status":"NOT_CHECKED"})
    return pd.DataFrame(records)

def test_queue_60_unique_taxa_120_observer_slots():
    x=check_queue(frame())
    assert x.inat_taxon_id.nunique()==60
    assert len(x)==90
    assert x.required_additional_distinct_observer_photos.sum()==120

def test_queue_rejects_status_leak_and_duplicate_year():
    x=frame()
    x.loc[x.index[0],"live_metadata_status"]="PHOTO_OPENED"
    with pytest.raises(ValueError,match="outcomes"):
        check_queue(x)
    y=frame()
    y.loc[y.index[0],"target_year"]=y.loc[y.index[1],"target_year"]
    with pytest.raises(ValueError,match="Repeat query"):
        check_queue(y)

def test_site_anchor_year_observer_guard():
    source=pd.DataFrame([
        {"inat_taxon_id":1100,"photo_id":"a","latitude":40.0,"longitude":-105.0,
         "month":6,"year":2020,"observer":"123"},
        {"inat_taxon_id":1100,"photo_id":"b","latitude":40.001,"longitude":-105.001,
         "month":6,"year":2022,"observer":"234"}
    ])
    row=SimpleNamespace(inat_taxon_id=1100,source_anchor_photo_id="a",
                        target_year=2020,calendar_month=6,
                        existing_distinct_observers_in_target_year=1,
                        required_additional_distinct_observer_photos=1)
    assert anchor_context_gap2(source,row)[2]=={"123"}
    row.required_additional_distinct_observer_photos=2
    with pytest.raises(RuntimeError,match="gap mismatch"):
        anchor_context_gap2(source,row)

"""Synthetic queue audit preserves source denominators, year/morph blindness and HOLD."""
import copy
import sys
from pathlib import Path
import pandas as pd
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"/"analysis"))
from build_fcp_replenishment_metadata_query_queue_20261008 import (
    build_queue, EXPECTED_FROZEN_GAP, REQUIRED_ORIGINAL
)


def source():
    rows=[]
    obs={}
    tid=0
    for c,n in REQUIRED_ORIGINAL.items():
        gaps=EXPECTED_FROZEN_GAP[c]
        obs[c]={"minimal_additional_observer_photo_slots_histogram":{**gaps,"3":0,"4":0}}
        values=([0]*gaps["0"]+[1]*gaps["1"]+[2]*gaps["2"]+
                [None]*(n-sum(gaps.values())))
        for gap in values:
            tid+=1
            rows.append({
                "cohort":c, "inat_taxon_id":str(tid), "species":f"Plant species{tid}",
                "two_year_month_already_observed": gap is not None,
                "minimum_additional_observer_photo_slots":gap,
                "candidate_years":"2020,2024" if gap is not None else None,
                "candidate_month":5 if gap is not None else None,
                "anchor_photo_id":f"photo{tid}" if gap is not None else None,
                "existing_observer_counts":({0:"2,2",1:"1,2",2:"1,1"}[gap]
                                              if gap is not None else None),
                "original_anchor_eligible":gap==0,
            })
    return pd.DataFrame(rows),{
        "schema":"fcp_same_siteyear_same_month_replenishment_gap_v1",
        "original_all_cohorts_gate_hold":True,"cohorts":obs
    }


def test_historical_global_query_queue_covers_all_gaps():
    d,r=source()
    result,queue=build_queue(d,r)
    assert len(d)==1109
    assert result["n_gap1_species"]==92
    assert result["n_gap2_species"]==171
    assert result["minimum_hypothetical_new_photo_slots_for_all_three_cohorts"]==50
    assert result["cohorts"]["validation"]["required_fraction_of_gap1_candidates_successful_if_gap2_unused"]==pytest.approx(20/22)
    assert len(queue)==92+2*171
    assert queue["required_additional_distinct_observer_photos"].sum()==92+2*171
    assert queue.live_metadata_status.eq("NOT_CHECKED").all()
    assert not {"morph","white","pigment","color","colour"} & set(queue)


def test_rejects_unexpected_photo_colour_input():
    d,r=source()
    d["morph"]="red_pink"
    with pytest.raises(ValueError,match="Photo-colour"):
        build_queue(d,r)


def test_rejects_missing_species():
    d,r=source()
    with pytest.raises(ValueError,match="Incomplete"):
        build_queue(d.iloc[1:],r)


def test_rejects_changed_observer_gap():
    d,r=source()
    k=d.index[d.minimum_additional_observer_photo_slots==1][0]
    d.loc[k,"existing_observer_counts"]="1,1"
    with pytest.raises(ValueError,match="disagree"):
        build_queue(d,r)


def test_rejects_qualified_photo_with_nonzero_gap():
    d,r=source()
    k=d.index[d.minimum_additional_observer_photo_slots==1][0]
    d.loc[k,"original_anchor_eligible"]=True
    with pytest.raises(ValueError,match="Nonzero"):
        build_queue(d,r)


def test_rejects_original_hold_change():
    d,r=source()
    r["original_all_cohorts_gate_hold"]=False
    with pytest.raises(ValueError,match="HOLD"):
        build_queue(d,r)

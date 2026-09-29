from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[1]
DESIGN=ROOT/"data"/"spatiotemporal_bridge_opportunity_design_v0_1.json"
RESULT=ROOT/"results"/"spatiotemporal_bridge_opportunity_v0_1"/"result.json"


def test_spatiotemporal_bridge_hold_is_frozen():
    d=json.loads(DESIGN.read_text())
    r=json.loads(RESULT.read_text())

    assert d["status"]=="FROZEN_BEFORE_SPATIOTEMPORAL_BRIDGE_OPPORTUNITY_AUDIT"
    assert d["qualification"]["minimum_unused_species_per_clade"]==5
    assert d["qualification"]["minimum_qualifying_clades"]==12
    assert d["qualification"]["minimum_total_selected_species"]==60
    assert d["qualification"]["maximum_selected_species_per_clade"]==10

    assert r["status"]=="HOLD_INSUFFICIENT_FRESH_CLADE_COVERAGE"
    assert r["outcome_blind"] is True
    assert r["forbidden_outcomes_opened"] is False
    assert r["bridge_outcome_computed"] is False

    assert r["frame"]["candidate_species"]==3230
    assert r["frame"]["prior_third_selected_species"]==500
    assert r["frame"]["unused_u100_species"]==2730
    assert r["frame"]["target_temporal_clades"]==28
    assert r["frame"]["matched_unused_species"]==127

    q=r["qualification"]
    assert q["minimum_unused_species_per_clade"]==5
    assert q["minimum_qualifying_clades"]==12
    assert q["observed_n_qualifying_clades"]==10
    assert q["gate_pass"] is False
    assert len(q["observed_qualifying_clades"])==10

    assert r["decision"]["existing_sparse_bridge_correlation"]=="DO_NOT_RUN"
    assert r["decision"]["new_bridge_selection"]=="DO_NOT_SELECT"
    assert r["decision"]["threshold_relaxation"]=="PROHIBITED"
    assert r["decision"]["family_level_or_synonym_broadening"]=="PROHIBITED_POST_HOC"

    assert r["science_changes"]["fcp_current_manuscript"] is False
    assert r["science_changes"]["chun_el_v0_3"] is False

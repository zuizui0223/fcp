from pathlib import Path
import csv
import json

ROOT=Path(__file__).resolve().parents[1]
DESIGN=ROOT/"data"/"general_spatiotemporal_turnover_design_v0_1.json"
MATRIX=ROOT/"data"/"general_spatiotemporal_turnover_eligibility_matrix_v0_1.csv"
AUDIT=ROOT/"data"/"general_spatiotemporal_turnover_source_audit_v0_1.json"
SEARCH=ROOT/"data"/"general_spatiotemporal_turnover_source_search_registry_v0_1.csv"


def test_general_turnover_design_is_frozen_and_non_flower_specific():
    d=json.loads(DESIGN.read_text())
    assert d["status"]=="FROZEN_BEFORE_GENERAL_SPATIOTEMPORAL_TURNOVER_SOURCE_AUDIT"
    assert set(d["response_vector"])=={
        "trait_time","trait_space","interaction_time","interaction_space"
    }
    assert set(d["predictors"])=={
        "context_heterogeneity","specialization","dispersal_mobility"
    }
    assert d["fixed_strata"]["interaction_type"]==[
        "mutualism","antagonism","mixed_or_other"
    ]
    assert d["programme_gate"]["trait_domain_min_independent_systems"]==12
    assert d["programme_gate"]["interaction_domain_min_independent_systems"]==12
    assert d["programme_gate"]["trait_domain_min_major_taxonomic_groups"]==3
    assert d["programme_gate"]["interaction_domain_min_interaction_types"]==2
    assert d["programme_gate"]["full_general_law_requires_both_domains"] is True
    assert d["relation_to_existing_programmes"]["active_FCP_manuscript_changed"] is False
    assert d["relation_to_existing_programmes"]["frozen_CHUN_EL_v0_3_changed"] is False


def test_interaction_space_requires_rewiring_decomposition():
    d=json.loads(DESIGN.read_text())
    r=d["response_vector"]["interaction_space"]
    assert r["required_decomposition"]==[
        "species_turnover_component","rewiring_among_cooccurring_species"
    ]
    assert r["primary_component"]=="rewiring_among_cooccurring_species"


def test_no_tier_a_source_is_promoted_yet():
    rows=list(csv.DictReader(MATRIX.open()))
    assert len(rows)>=7
    forbidden={"PASS","QUALIFIED","TIER_A","TIER_A_FULL_LINKED"}
    assert all(row["terminal_status"] not in forbidden for row in rows)
    audit=json.loads(AUDIT.read_text())
    assert audit["current_candidates"]["tier_A_promoted"]==0
    assert audit["current_candidates"]["decision"]=="NO_TIER_A_FULL_LINKED_SYSTEM_YET"
    assert audit["outcome_opening_authorized"] is False


def test_finite_source_search_is_exactly_six_classes():
    rows=list(csv.DictReader(SEARCH.open()))
    assert len(rows)==6
    assert [int(r["search_order"]) for r in rows]==[1,2,3,4,5,6]
    assert all(r["status"]=="PENDING_SOURCE_ONLY_AUDIT" for r in rows)
    assert all("stop class after" in r["stop_rule"] for r in rows)


def test_context_heterogeneity_is_only_shared_directional_accelerator():
    d=json.loads(DESIGN.read_text())
    p=d["predictors"]["context_heterogeneity"]["primary_prediction"]
    assert p=={
        "lambda_trait_time":"+",
        "lambda_trait_space":"+",
        "lambda_interaction_time":"+",
        "lambda_interaction_space":"+"
    }


def test_specialization_and_mobility_are_decouplers_not_rescues():
    d=json.loads(DESIGN.read_text())
    s=d["predictors"]["specialization"]["primary_prediction"]
    m=d["predictors"]["dispersal_mobility"]["primary_prediction"]
    assert s["lambda_interaction_time"]=="-"
    assert s["lambda_interaction_space_rewiring"]=="-"
    assert m["lambda_trait_space"]=="-"
    assert m["lambda_interaction_space"]=="-"
    assert m["lambda_trait_time"]=="no directional primary prediction"
    assert m["lambda_interaction_time"]=="no directional primary prediction"

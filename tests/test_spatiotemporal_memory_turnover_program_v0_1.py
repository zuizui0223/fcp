from pathlib import Path
import csv
import json

ROOT=Path(__file__).resolve().parents[1]
DESIGN=ROOT/"data"/"spatiotemporal_memory_turnover_design_v0_1.json"
MATRIX=ROOT/"data"/"spatiotemporal_memory_turnover_source_matrix_v0_1.csv"
PRED=ROOT/"data"/"spatiotemporal_memory_turnover_predictor_matrix_v0_2.csv"


def load():
    return json.loads(DESIGN.read_text())


def test_four_response_cells_share_one_effect_convention():
    d=load()
    assert d["status"]=="FROZEN_BEFORE_GENERALIZED_TRAIT_OR_INTERACTION_OUTCOME_OPENING"
    assert d["canonical_response"]["name"]=="null_centered_turnover_rho"
    assert set(d["response_matrix"])=={
        "trait_time","trait_space","interaction_time","interaction_space"
    }
    assert "faster turnover" in d["canonical_response"]["direction"]


def test_interaction_space_is_rewiring_not_species_turnover():
    d=load()
    x=d["interaction_decomposition"]
    assert x["required"] is True
    assert x["components"]==[
        "species_turnover",
        "rewiring_among_shared_species",
        "total_interaction_turnover",
    ]
    assert x["primary_for_cross_axis_coupling"]=="rewiring_among_shared_species"


def test_only_three_biological_predictor_axes_are_frozen():
    d=load()
    assert set(d["biological_predictors"])=={
        "context_heterogeneity",
        "ecological_specialization",
        "dispersal_mobility",
    }
    assert d["biological_predictors"]["context_heterogeneity"]["role"]=="primary common driver"
    assert d["biological_predictors"]["ecological_specialization"]["role"]=="primary decoupling candidate"
    assert d["biological_predictors"]["dispersal_mobility"]["role"]=="secondary decoupling candidate"


def test_time_space_pooling_requires_independent_replication():
    d=load()
    h=d["primary_hypotheses"]["H1_common_coupling"]
    assert h["strata"]==["traits","interactions"]
    assert ">=12 independent system units" in h["pooling_rule"]


def test_eligibility_and_informativeness_gates_are_frozen():
    d=load()
    t=d["eligibility_gates"]["trait_system"]
    i=d["eligibility_gates"]["interaction_system"]
    q=d["eligibility_gates"]["informativeness"]
    assert t["temporal_min_species"]==20
    assert t["spatial_min_species_with_repeated_measurements"]==8
    assert t["spatial_min_georeferenced_records_per_species"]==20
    assert t["spatial_min_unique_cells_per_species"]==5
    assert i["temporal_min_focal_taxa"]==20
    assert i["spatial_min_networks"]==5
    assert i["spatial_min_georeferenced_networks"]==5
    assert i["spatial_min_repeated_focal_taxa"]==10
    assert q["required"] is True
    assert q["no_posthoc_relaxation"] is True
    assert ">=80% directional recovery" in q["rule"]


def test_finite_source_family_is_locked_before_outcomes():
    d=load()
    assert d["source_order"]["trait_primary"]==["BIEN_4_2","AUSTRAITS_OPEN_RELEASE"]
    assert d["source_order"]["interaction_primary"]==["MANGAL_CURATED_NETWORKS"]
    assert d["source_order"]["interaction_external"]==["GLOBI_STABLE_ZENODO"]
    assert d["finite_family_rule"]["max_primary_trait_sources"]==2
    assert d["finite_family_rule"]["max_primary_interaction_sources"]==2
    assert d["finite_family_rule"]["source_switch_forbidden_after"]=="biological turnover outcome opened"


def test_source_matrix_keeps_holds_and_anchors_separate():
    rows={r["source_id"]:r for r in csv.DictReader(MATRIX.open())}
    assert rows["TRY_V7_CURRENT_DOWNLOAD_HOLD"]["selection_priority"]=="HOLD"
    assert rows["MANGAL_CURATED_NETWORKS"]["selection_priority"]=="P1"
    assert rows["GLOBI_STABLE_ZENODO"]["selection_priority"]=="P2"
    assert rows["CHUN_FCP_FLOWER_COLOUR"]["admission_role"]=="proof_of_concept_only"
    assert rows["IWE_INTERACTION_WINDOW"]["admission_role"]=="conceptual_interaction_anchor"


def test_current_papers_are_not_modified_by_programme_design():
    d=load()
    assert d["current_paper_science_changed"] is False
    assert d["motivating_anchors"]["use_as_generalization_evidence"] is False


def test_predictor_matrix_is_complete_and_preoutcome():
    rows=list(csv.DictReader(PRED.open()))
    assert len(rows)==12
    keys={(r["response_id"],r["predictor_id"]) for r in rows}
    assert keys=={(a,b) for a in ("TT","TS","IT","IS") for b in ("CH","ES","DM")}
    expected={
        ("TT","CH"):"positive",("TS","CH"):"positive",("IT","CH"):"positive",("IS","CH"):"positive",
        ("TT","ES"):"negative",("TS","ES"):"positive",("IT","ES"):"negative",("IS","ES"):"positive",
        ("TT","DM"):"no_directional_primary",("TS","DM"):"negative",
        ("IT","DM"):"no_directional_primary",("IS","DM"):"negative",
    }
    for r in rows:
        assert r["status"]=="PREOUTCOME_FROZEN"
        assert r["expected_direction"]==expected[(r["response_id"],r["predictor_id"])]

from pathlib import Path
import csv
import json

ROOT=Path(__file__).resolve().parents[1]
CONTRACT=ROOT/"data"/"spatiotemporal_predictor_independence_contract_v0_1.json"
MATRIX=ROOT/"data"/"spatiotemporal_predictor_source_matrix_v0_1.csv"
SOURCES=ROOT/"data"/"spatiotemporal_predictor_source_verification_v0_1.json"
TUNDRA=ROOT/"data"/"spatiotemporal_trait_geometry_calibration_tundra_v0_1.json"
NAMESPACE=ROOT/"results"/"globi_specialization_namespace_independence_v0_1"/"result.json"
CITATION=ROOT/"results"/"globi_citation_source_exclusion_v0_1"/"result.json"


def test_predictor_family_and_independence_contract():
    d=json.loads(CONTRACT.read_text())
    assert d["status"]=="FROZEN_BEFORE_GENERALIZED_TURNOVER_OUTCOME_OPENING"
    assert set(d["predictor_families"])=={
        "context_heterogeneity",
        "ecological_specialization",
        "dispersal_mobility",
    }
    assert d["predictor_families"]["context_heterogeneity"]["current_status"]=="ELIGIBLE_SOURCE_ROUTE"
    assert d["predictor_families"]["ecological_specialization"]["current_status"]=="ELIGIBLE_IN_PRINCIPLE_REQUIRES_GLOBI_SOURCE_EXCLUSION_AUDIT"
    assert d["predictor_families"]["dispersal_mobility"]["current_status"]=="HOLD_PRIMARY_PLANT_PANEL"
    assert d["predictor_missingness"]["no_imputation_from_response"] is True
    assert d["predictor_missingness"]["no_source_switch_after_response_opening"] is True
    assert d["current_science_changed"] is False


def test_predictor_source_matrix_is_complete():
    rows=list(csv.DictReader(MATRIX.open()))
    assert len(rows)==12
    keys={(r["response_id"],r["predictor_id"]) for r in rows}
    assert keys=={(a,b) for a in ("TT","TS","IT","IS") for b in ("CH","ES","DM")}

    expected={
        ("TT","CH"):"positive",("TS","CH"):"positive",
        ("IT","CH"):"positive",("IS","CH"):"positive",
        ("TT","ES"):"negative",("TS","ES"):"positive",
        ("IT","ES"):"negative",("IS","ES"):"positive",
        ("TT","DM"):"no_directional_primary",("TS","DM"):"negative",
        ("IT","DM"):"no_directional_primary",("IS","DM"):"negative",
    }
    for r in rows:
        assert r["expected_direction"]==expected[(r["response_id"],r["predictor_id"])]

    by={(r["response_id"],r["predictor_id"]):r for r in rows}
    assert by[("IT","ES")]["current_status"]=="HOLD_PENDING_SOURCE_EXCLUSION_AUDIT"
    assert by[("IS","ES")]["current_status"]=="HOLD_PENDING_SOURCE_EXCLUSION_AUDIT"
    assert by[("TT","DM")]["current_status"]=="HOLD_NO_PINNED_DIRECT_SOURCE"
    assert by[("TS","DM")]["current_status"]=="HOLD_NO_PINNED_DIRECT_SOURCE"


def test_source_metadata_and_calibration_boundaries():
    s=json.loads(SOURCES.read_text())
    assert s["status"]=="PREDICTOR_SOURCE_METADATA_VERIFIED_PREOUTCOME"
    assert s["sources"]["CHELSA_V2_1"]["version"]=="2.1"
    assert s["sources"]["GLOBI"]["stable_archive_doi"]=="10.5281/zenodo.3950589"
    assert s["sources"]["TRY_V7"]["current_status"]=="HOLD_DOWNLOAD_ACCESS_NOTICE"
    assert s["generalized_turnover_outcomes_opened"] is False

    t=json.loads(TUNDRA.read_text())
    assert t["status"]=="TRAIT_GEOMETRY_CALIBRATION_SOURCE_QUALIFIED"
    assert t["role"]=="CALIBRATION_ONLY_NOT_PRIMARY_GENERALIZATION_SOURCE"
    assert t["physical_support"]["rows"]==91970
    assert t["physical_support"]["traits_passing_paired_source_gate"]==14
    assert t["biological_turnover_outcomes_opened"] is False
    assert "not automatically" in t["predictor_boundary"]


def test_globi_source_independence_preflights_are_bounded():
    n=json.loads(NAMESPACE.read_text())
    assert n["status"]=="GLOBI_SPECIALIZATION_NAMESPACE_PREFLIGHT_COMPLETE"
    assert n["mangal_like_namespaces"]==["globalbioticinteractions/mangal"]
    assert n["interaction_rows_opened"] is False
    assert n["taxon_identities_opened"] is False
    assert n["specialization_values_computed"] is False

    x=json.loads(CITATION.read_text())
    assert x["status"]=="GLOBI_CITATION_EXCLUSION_MECHANISM_VERIFIED_FOR_CURRENT_MANGAL_CALIBRATION_DOIS"
    assert x["all_current_response_dois_found"] is True
    assert x["response_doi_match_counts"]["hadfield_2014"] >= 1
    assert x["response_doi_match_counts"]["havens_1992"] >= 1
    assert x["response_doi_match_counts"]["ricciardi_2010"] >= 1
    assert x["interaction_rows_opened"] is False
    assert x["taxon_identities_opened"] is False
    assert x["partner_entropy_computed"] is False

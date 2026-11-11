from pathlib import Path
import csv
import json

ROOT=Path(__file__).resolve().parents[1]
DESIGN=ROOT/"data"/"spatiotemporal_memory_turnover_design_v0_1.json"
MATRIX=ROOT/"data"/"spatiotemporal_memory_turnover_predictor_matrix_v0_2.csv"
GLOBI_DESIGN=ROOT/"data"/"globi_spatiotemporal_metadata_preflight_design_v0_1.json"
GLOBI_RESULT=ROOT/"results"/"globi_spatiotemporal_metadata_preflight_v0_1"/"result.json"


def test_predictor_matrix_is_complete_and_preoutcome():
    d=json.loads(DESIGN.read_text())
    rows=list(csv.DictReader(MATRIX.open()))
    assert len(rows)==12
    assert {(r["response_id"],r["predictor_id"]) for r in rows} == {
        (resp,pred)
        for resp in ("TT","TS","IT","IS")
        for pred in ("CH","ES","DM")
    }
    assert all(r["status"]=="PREOUTCOME_FROZEN" for r in rows)

    directions={(r["response_id"],r["predictor_id"]):r["expected_direction"] for r in rows}
    for resp in ("TT","TS","IT","IS"):
        assert directions[(resp,"CH")]=="positive"
    assert directions[("TT","ES")]=="negative"
    assert directions[("TS","ES")]=="positive"
    assert directions[("IT","ES")]=="negative"
    assert directions[("IS","ES")]=="positive"
    assert directions[("TS","DM")]=="negative"
    assert directions[("IS","DM")]=="negative"
    assert directions[("TT","DM")]=="no_directional_primary"
    assert directions[("IT","DM")]=="no_directional_primary"

    assert d["interaction_decomposition"]["primary_for_cross_axis_coupling"]=="rewiring_among_shared_species"
    assert d["current_paper_science_changed"] is False


def test_globi_preflight_is_schema_hold_not_biological_result():
    d=json.loads(GLOBI_DESIGN.read_text())
    r=json.loads(GLOBI_RESULT.read_text())
    assert d["status"]=="FROZEN_BEFORE_GLOBI_DATASET_INDEX_OPENING"
    assert d["source"]["interaction_rows_open_authorized"] is False
    assert r["status"]=="HOLD_GLOBI_DATASET_INDEX_LACKS_GEOMETRY_TIME_METADATA"
    assert r["dataset_index"]["rows"]==458
    assert r["dataset_index"]["header"]==["namespace"]
    assert r["dataset_index"]["geometry_time_gate_pass"] is False
    assert r["firewall"]["interaction_rows_opened"] is False
    assert r["firewall"]["taxon_identities_opened"] is False
    assert r["firewall"]["biological_outcome_computed"] is False
    assert r["decision"]["integrated_interaction_rows_may_be_opened_to_rescue_gate"] is False

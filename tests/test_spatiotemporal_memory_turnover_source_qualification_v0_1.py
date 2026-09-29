from pathlib import Path
import csv
import json

ROOT=Path(__file__).resolve().parents[1]
NODE=ROOT/"results"/"mangal_node_taxonomy_support_v0_1"/"result.json"
SEM=ROOT/"data"/"mangal_source_semantics_audit_v0_1.json"
MATRIX=ROOT/"data"/"spatiotemporal_memory_turnover_source_matrix_v0_1.csv"


def test_node_firewall_and_support_result():
    x=json.loads(NODE.read_text())
    assert x["status"]=="MANGAL_NODE_TAXONOMY_SUPPORT_AUDIT_COMPLETE"
    assert x["interaction_edges_opened"] is False
    assert x["trait_values_opened"] is False
    assert x["environment_values_opened"] is False
    assert x["necessary_support_passing_dataset_ids"]==[15,72,75]
    assert x["full_interaction_spatial_admission_authorized"] is False


def test_source_semantics_do_not_promote_fresh_mangal_system():
    x=json.loads(SEM.read_text())
    assert x["status"]=="MANGAL_SOURCE_SEMANTICS_AUDIT_COMPLETE_NO_FRESH_FULL_PAIR_PROMOTED"
    by={r["dataset_id"]:r for r in x["candidates"]}
    assert by[15]["terminal_status"]=="HOLD_REWIRING_STRUCTURALLY_UNIDENTIFIABLE"
    assert by[72]["terminal_status"]=="RETROSPECTIVE_CALIBRATION_CANDIDATE"
    assert by[75]["terminal_status"]=="HOLD_TEMPORAL_GUILD_SUPPORT_LT20"
    assert x["fresh_full_time_space_interaction_systems_promoted"]==0
    assert x["next_source_under_finite_family"]=="GLOBI_STABLE_ZENODO"
    assert x["full_interaction_outcome_opening_authorized"] is False


def test_source_matrix_moves_to_globi_without_rescue():
    rows={r["source_id"]:r for r in csv.DictReader(MATRIX.open())}
    assert rows["MANGAL_CURATED_NETWORKS"]["metadata_status"]=="SOURCE_FAMILY_EXHAUSTED_NO_FRESH_FULL_PAIR"
    assert rows["MANGAL_CURATED_NETWORKS"]["admission_role"]=="retrospective_calibration_only"
    assert rows["GLOBI_STABLE_ZENODO"]["selection_priority"]=="P2"

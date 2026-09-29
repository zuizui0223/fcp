#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
DESIGN=ROOT/"data"/"spatiotemporal_memory_turnover_design_v0_1.json"
MATRIX=ROOT/"data"/"spatiotemporal_memory_turnover_predictor_matrix_v0_2.csv"
SOURCES=ROOT/"data"/"spatiotemporal_memory_turnover_source_matrix_v0_1.csv"
VERIFY=ROOT/"data"/"spatiotemporal_memory_turnover_source_verification_v0_1.json"
MANGAL_NODE=ROOT/"results"/"mangal_node_taxonomy_support_v0_1"/"result.json"
MANGAL_SEM=ROOT/"data"/"mangal_source_semantics_audit_v0_1.json"


def rows(path:Path):
    with path.open(newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))


def main()->int:
    d=json.loads(DESIGN.read_text())
    m=rows(MATRIX)
    s={r["source_id"]:r for r in rows(SOURCES)}
    v=json.loads(VERIFY.read_text())
    node=json.loads(MANGAL_NODE.read_text())
    sem=json.loads(MANGAL_SEM.read_text())

    assert d["status"]=="FROZEN_BEFORE_GENERALIZED_TRAIT_OR_INTERACTION_OUTCOME_OPENING"
    assert d["canonical_response"]["name"]=="null_centered_turnover_rho"
    assert set(d["response_matrix"])=={
        "trait_time","trait_space","interaction_time","interaction_space"
    }
    assert set(d["biological_predictors"])=={
        "context_heterogeneity","ecological_specialization","dispersal_mobility"
    }
    assert d["interaction_decomposition"]["required"] is True
    assert d["interaction_decomposition"]["primary_for_cross_axis_coupling"]=="rewiring_among_shared_species"
    assert d["motivating_anchors"]["use_as_generalization_evidence"] is False
    assert d["current_paper_science_changed"] is False

    assert len(m)==12
    keys={(r["response_id"],r["predictor_id"]) for r in m}
    assert keys=={(r,p) for r in ("TT","TS","IT","IS") for p in ("CH","ES","DM")}
    expect={
        ("TT","CH"):"positive",("TS","CH"):"positive",("IT","CH"):"positive",("IS","CH"):"positive",
        ("TT","ES"):"negative",("TS","ES"):"positive",("IT","ES"):"negative",("IS","ES"):"positive",
        ("TT","DM"):"no_directional_primary",("TS","DM"):"negative",
        ("IT","DM"):"no_directional_primary",("IS","DM"):"negative",
    }
    for r in m:
        assert r["expected_direction"]==expect[(r["response_id"],r["predictor_id"])]
        assert r["status"]=="PREOUTCOME_FROZEN"

    assert s["MANGAL_CURATED_NETWORKS"]["metadata_status"]=="SOURCE_FAMILY_EXHAUSTED_NO_FRESH_FULL_PAIR"
    assert s["MANGAL_CURATED_NETWORKS"]["admission_role"]=="retrospective_calibration_only"
    assert s["GLOBI_STABLE_ZENODO"]["selection_priority"]=="P2"
    assert s["TRY_V7_CURRENT_DOWNLOAD_HOLD"]["selection_priority"]=="HOLD"
    assert s["CHUN_FCP_FLOWER_COLOUR"]["admission_role"]=="proof_of_concept_only"

    assert v["programme_outcome_opening_authorized"] is False
    assert node["interaction_edges_opened"] is False
    assert node["full_interaction_spatial_admission_authorized"] is False
    assert sem["fresh_full_time_space_interaction_systems_promoted"]==0
    assert sem["full_interaction_outcome_opening_authorized"] is False
    assert sem["next_source_under_finite_family"]=="GLOBI_STABLE_ZENODO"

    print(json.dumps({
        "status":"SPATIOTEMPORAL_MEMORY_TURNOVER_CANONICAL_DESIGN_VALID",
        "responses":4,
        "predictors":3,
        "cells":12,
        "canonical_effect":"null_centered_turnover_rho",
        "interaction_rewiring_decomposition_required":True,
        "mangal_fresh_full_pair_promoted":0,
        "next_interaction_source":"GLOBI_STABLE_ZENODO",
        "generalized_biological_outcomes_opened":False
    },indent=2))
    return 0


if __name__=="__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
DESIGN=ROOT/"development"/"spatiotemporal_turnover"/"design_v0_1.json"
MATRIX=ROOT/"development"/"spatiotemporal_turnover"/"eligibility_matrix_v0_1.csv"
SOURCES=ROOT/"development"/"spatiotemporal_turnover"/"source_registry_v0_1.csv"
PRIOR=ROOT/"development"/"spatiotemporal_turnover"/"prior_art_registry_v0_1.csv"


def read_csv(path:Path):
    with path.open(newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))


def main()->int:
    d=json.loads(DESIGN.read_text())
    rows=read_csv(MATRIX)
    sources=read_csv(SOURCES)
    prior=read_csv(PRIOR)

    assert d["status"]=="FROZEN_GENERAL_SPATIOTEMPORAL_TURNOVER_DESIGN_PREOUTCOME"
    assert d["programme_role"]=="CROSS_PROGRAM_DEVELOPMENT_LANE_NOT_CURRENT_FCP_MANUSCRIPT"
    assert d["canonical_direction"].startswith("All four primary responses")
    assert set(x["id"] for x in d["responses"].values())=={"TT","TS","IT","IS"}
    assert set(x["id"] for x in d["predictors"].values())=={"AB","BS","MO"}
    assert all(x["flower_colour_required"] is False for x in d["responses"].values())
    assert d["responses"]["interaction_space"]["primary_component"].startswith("beta_os")
    assert d["responses"]["interaction_space"]["required_diagnostic"].startswith("beta_st")
    assert d["governance"]["current_FCP_manuscript_changed"] is False
    assert d["governance"]["frozen_CHUN_EL_v0_3_changed"] is False
    assert d["governance"]["biological_outcomes_opened_by_this_design"] is False

    assert len(rows)==12
    keys={(r["response_id"],r["predictor_id"]) for r in rows}
    assert len(keys)==12
    assert keys=={(r,p) for r in ("TT","TS","IT","IS") for p in ("AB","BS","MO")}

    for pred_key,pred in d["predictors"].items():
        pid=pred["id"]
        expected=pred["primary_prediction"]
        for rid, direction in expected.items():
            row=next(x for x in rows if x["response_id"]==rid and x["predictor_id"]==pid)
            assert row["expected_direction"]==direction

    src={x["source_id"]:x for x in sources}
    assert src["TRY_V7"]["current_status"]=="HOLD_PROVIDER_DOWNLOAD_2026_09_23"
    assert src["GLOBI_STABLE"]["outcome_use_policy"]=="SOURCE_AUDIT_REQUIRED_BEFORE_OUTCOME"
    assert src["FCP_DISTTRAIT"]["outcome_use_policy"]=="DO_NOT_COUNT_AS_GENERAL_REPLICATION"
    assert src["CHUN_TEMPORAL"]["outcome_use_policy"]=="DO_NOT_COUNT_AS_GENERAL_REPLICATION"

    ids={x["prior_id"] for x in prior}
    assert {"POISOT_2012","FRUND_2021","ROHR_BASCOMPTE_2014","WARD_2026","MCLEOD_2020"}.issubset(ids)

    print(json.dumps({
        "status":"SPATIOTEMPORAL_TURNOVER_DESIGN_V0_1_VALID",
        "responses":4,
        "predictors":3,
        "cells":12,
        "flower_colour_required":False,
        "biological_outcomes_opened":False,
        "interaction_rewiring_decomposition_required":True
    },indent=2))
    return 0


if __name__=="__main__":
    raise SystemExit(main())

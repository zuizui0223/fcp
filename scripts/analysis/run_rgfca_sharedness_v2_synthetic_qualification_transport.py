#!/usr/bin/env python3
"""Transport-geometry compatibility wrapper for frozen RGFCA-v2 synthetic runner.

Scientific mapping is unchanged. This wrapper only teaches the runner the canonical
status/key names emitted by the successful transport-only 150x300 geometry recovery.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

import run_rgfca_sharedness_v2_synthetic_qualification as base

ROOT = Path(__file__).resolve().parents[2]
ACT = ROOT / "docs/supporting/rgfca_sharedness_v2_geometry_activation_amendment_v1.json"


def load_transport_inputs():
    contract=json.loads(base.CONTRACT.read_text()); mapping=json.loads(base.MAPPING.read_text())
    gate=json.loads(base.GATE_AMEND.read_text()); missing=json.loads(base.MISSING_RULE.read_text())
    result=json.loads(base.GEOMETRY_RESULT.read_text()); act=json.loads(ACT.read_text())
    if contract["status"]!="frozen_before_high_depth_capacity_pilot_final_result_and_before_any_v2_synthetic_outcome": raise RuntimeError("parent contract drift")
    if mapping["status"]!="frozen_after_complete_150x300_geometry_before_any_v2_synthetic_score": raise RuntimeError("mapping drift")
    if gate["status"]!="frozen_before_any_v2_synthetic_world_score" or missing["status"]!="frozen_before_any_v2_synthetic_world_score" or act["status"]!="frozen_before_any_v2_synthetic_world_score": raise RuntimeError("pre-outcome amendments missing")
    req=act["required"]
    if result.get("protocol")!=act["required_result_protocol"] or result.get("status")!=act["required_result_status"]: raise RuntimeError("transport geometry status drift")
    for k,v in req.items():
        if result.get(k)!=v: raise RuntimeError(f"transport geometry gate drift: {k}")
    frame=pd.read_csv(base.GEOMETRY)
    required={"species_order","sample_role","photo_order","latitude","longitude"}
    if not required.issubset(frame.columns): raise RuntimeError(f"geometry missing columns {sorted(required-set(frame.columns))}")
    if len(frame)!=45000 or frame["species_order"].nunique()!=150: raise RuntimeError("geometry row/species census drift")
    frame=frame.sort_values(["species_order","photo_order"],kind="mergesort").reset_index(drop=True)
    counts=frame.groupby("species_order",sort=True).size()
    if len(counts)!=150 or not counts.eq(300).all(): raise RuntimeError("not exactly 300 rows/species")
    roles=frame.groupby("species_order",sort=True)["sample_role"].agg(lambda x: tuple(pd.unique(x.astype(str))))
    if any(len(x)!=1 for x in roles): raise RuntimeError("species role not unique")
    train=np.array([int(i) for i,x in roles.items() if x[0]=="training"],dtype=int)
    test=np.array([int(i) for i,x in roles.items() if x[0]=="evaluation"],dtype=int)
    if len(train)!=75 or len(test)!=75 or set(train)&set(test): raise RuntimeError("75/75 split drift")
    if sorted(np.r_[train,test].tolist())!=list(range(150)): raise RuntimeError("species_order incomplete")
    lat=pd.to_numeric(frame["latitude"],errors="raise").to_numpy(float).reshape(150,300)
    lon=pd.to_numeric(frame["longitude"],errors="raise").to_numpy(float).reshape(150,300)
    xyz=base.latlon_xyz(lat.ravel(),lon.ravel()).reshape(150,300,3)
    if not np.isfinite(xyz).all() or not np.allclose(np.linalg.norm(xyz,axis=2),1.0,atol=1e-12): raise RuntimeError("invalid geometry")
    return contract,mapping,frame,xyz,train,test

base.load_inputs=load_transport_inputs

if __name__=="__main__":
    p=argparse.ArgumentParser(); p.add_argument("--stage",choices=["calibration","evaluation"],required=True); p.add_argument("--arm-index",type=int,required=True); p.add_argument("--output",type=Path,required=True)
    a=p.parse_args(); base.run(a.stage,a.arm_index,a.output)

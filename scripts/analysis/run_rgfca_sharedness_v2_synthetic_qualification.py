#!/usr/bin/env python3
"""Run one frozen RGFCA-v2 synthetic sharedness arm.

No observed flower/background colour, image pixel, climate value, 6-species result,
or 34-species result is read. The only empirical input is the frozen 150 x 300
metadata geometry.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import logsumexp

from benchmark_rgfca_sharedness_v2_exact_kernel import exact_axis_scores_vectorized

ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "docs/supporting/rgfca_sharedness_v2_synthetic_qualification_contract_v1.json"
MAPPING = ROOT / "docs/supporting/rgfca_sharedness_v2_synthetic_technical_mapping_v1.json"
GATE_AMEND = ROOT / "docs/supporting/rgfca_sharedness_v2_synthetic_gate_amendment_v1.json"
MISSING_RULE = ROOT / "docs/supporting/rgfca_sharedness_v2_synthetic_missing_axis_rule_v1.json"
GEOMETRY_RESULT = ROOT / "docs/supporting/rgfca_sharedness_v2_geometry_result_v1.json"
GEOMETRY = ROOT / "data/frozen/rgfca_sharedness_v2_geometry_photos_v1.csv"
FRACTIONS = np.array([0.0, 0.1, 0.25, 0.5, 0.75, 1.0], dtype=float)
RETENTIONS = (0.4, 0.6, 0.8)
REPS = 250


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def seed_for(root: int, *parts: object) -> int:
    raw = "|".join(map(str, (root,) + parts)).encode("utf-8")
    return int.from_bytes(hashlib.sha256(raw).digest()[:8], "little")


def unit_vectors(rng: np.random.Generator, shape: tuple[int, ...]) -> np.ndarray:
    x = rng.normal(size=shape)
    n = np.linalg.norm(x, axis=-1, keepdims=True)
    if np.any(n <= 1e-12):
        raise RuntimeError("unit-vector draw failure")
    return x / n


def latlon_xyz(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    latr = np.deg2rad(lat)
    lonr = np.deg2rad(lon)
    c = np.cos(latr)
    return np.c_[c * np.cos(lonr), c * np.sin(lonr), np.sin(latr)]


def load_inputs() -> tuple[dict, dict, pd.DataFrame, np.ndarray, np.ndarray, np.ndarray]:
    contract = json.loads(CONTRACT.read_text())
    mapping = json.loads(MAPPING.read_text())
    gate = json.loads(GATE_AMEND.read_text())
    missing = json.loads(MISSING_RULE.read_text())
    result = json.loads(GEOMETRY_RESULT.read_text())
    if contract["status"] != "frozen_before_high_depth_capacity_pilot_final_result_and_before_any_v2_synthetic_outcome":
        raise RuntimeError("parent contract drift")
    if mapping["status"] != "frozen_after_complete_150x300_geometry_before_any_v2_synthetic_score":
        raise RuntimeError("technical mapping drift")
    if gate["status"] != "frozen_before_any_v2_synthetic_world_score" or missing["status"] != "frozen_before_any_v2_synthetic_world_score":
        raise RuntimeError("pre-outcome amendments missing")
    if result.get("status") != "complete_frozen_metadata_geometry_150x300" or result.get("complete_fixed_frame") is not True:
        raise RuntimeError("150x300 geometry not complete")
    if [int(result[k]) for k in ("species", "training_species", "evaluation_species", "target_photos_per_species", "retained_photos")] != [150, 75, 75, 300, 45000]:
        raise RuntimeError("geometry census drift")
    frame = pd.read_csv(GEOMETRY)
    required = {"species_order", "sample_role", "photo_order", "latitude", "longitude"}
    if not required.issubset(frame.columns):
        raise RuntimeError(f"geometry missing columns {sorted(required-set(frame.columns))}")
    if len(frame) != 45000 or frame["species_order"].nunique() != 150:
        raise RuntimeError("geometry row/species census drift")
    frame = frame.sort_values(["species_order", "photo_order"], kind="mergesort").reset_index(drop=True)
    counts = frame.groupby("species_order", sort=True).size()
    if len(counts) != 150 or not counts.eq(300).all():
        raise RuntimeError("not exactly 300 rows/species")
    roles = frame.groupby("species_order", sort=True)["sample_role"].agg(lambda x: tuple(pd.unique(x.astype(str))))
    if any(len(x) != 1 for x in roles):
        raise RuntimeError("species role not unique")
    train_ids = np.array([int(i) for i, x in roles.items() if x[0] == "training"], dtype=int)
    test_ids = np.array([int(i) for i, x in roles.items() if x[0] == "evaluation"], dtype=int)
    if len(train_ids) != 75 or len(test_ids) != 75 or set(train_ids) & set(test_ids):
        raise RuntimeError("75/75 species split drift")
    if sorted(np.r_[train_ids, test_ids].tolist()) != list(range(150)):
        raise RuntimeError("species_order must be complete 0..149")
    lat = pd.to_numeric(frame["latitude"], errors="raise").to_numpy(float).reshape(150, 300)
    lon = pd.to_numeric(frame["longitude"], errors="raise").to_numpy(float).reshape(150, 300)
    xyz = latlon_xyz(lat.ravel(), lon.ravel()).reshape(150, 300, 3)
    if not np.isfinite(xyz).all() or not np.allclose(np.linalg.norm(xyz, axis=2), 1.0, atol=1e-12):
        raise RuntimeError("invalid geometry coordinates")
    return contract, mapping, frame, xyz, train_ids, test_ids


def nuisance_arms(mapping: dict) -> list[dict]:
    out = []
    for s in mapping["nuisance_scenarios"]:
        for r in RETENTIONS:
            row = dict(s)
            row.update(kind="nuisance", retention=float(r))
            out.append(row)
    return out


def positive_scenarios(mapping: dict) -> list[dict]:
    grid = mapping["positive_scenarios"]["cartesian_primary_grid"]
    out = []
    for a in grid["amplitude"]:
        for f in grid["shared_fraction"]:
            for sd in grid["species_threshold_sd"]:
                out.append({"id":"partially_shared_boundaries", "kind":"positive", "amplitude":float(a), "shared_fraction":float(f), "threshold_sd":float(sd)})
    for x in mapping["positive_scenarios"]["diagnostic_extra"]:
        out.append({"id":"partially_shared_boundaries", "kind":"positive", "amplitude":float(x["amplitude"]), "shared_fraction":float(x["shared_fraction"]), "threshold_sd":float(x["species_threshold_sd"])})
    if len(out) != 19:
        raise RuntimeError("positive scenario census drift")
    return out


def evaluation_arms(mapping: dict) -> list[dict]:
    out = nuisance_arms(mapping)
    for s in positive_scenarios(mapping):
        for r in RETENTIONS:
            row = dict(s)
            row["retention"] = float(r)
            out.append(row)
    if len(out) != 78:
        raise RuntimeError("evaluation arm census drift")
    return out


def arm_id(a: dict) -> str:
    keys = [k for k in sorted(a) if k != "kind"]
    return "|".join(f"{k}={a[k]}" for k in keys)


def clustered_axes(rng: np.random.Generator, n: int, concentration: float) -> np.ndarray:
    center = unit_vectors(rng, (1, 3))[0]
    raw = concentration * center[None, :] + rng.normal(size=(n, 3))
    return raw / np.linalg.norm(raw, axis=1)[:, None]


def retention_indices(root: int, stage: str, rep: int, retention: float, xyz: np.ndarray, coordinate_dependent: bool, slope: float = 0.0) -> np.ndarray:
    k = {0.4:120, 0.6:180, 0.8:240}[float(retention)]
    rng = np.random.default_rng(seed_for(root, "retention", stage, retention, rep))
    if not coordinate_dependent:
        priority = rng.random((150, 300))
    else:
        arng = np.random.default_rng(seed_for(root, "retention-coordinate", stage, rep))
        axis = unit_vectors(arng, (1, 3))[0]
        proj = np.einsum("snd,d->sn", xyz, axis)
        u = np.clip(rng.random((150,300)), 1e-12, 1-1e-12)
        gumbel = -np.log(-np.log(u))
        priority = -(float(slope) * proj + gumbel)  # lower is selected; equivalent to top-k utility
    order = np.argsort(priority, axis=1, kind="stable")[:, :k]
    return order.astype(np.int32)


def synthetic_world(root: int, stage: str, arm: dict, rep: int, xyz: np.ndarray) -> tuple[np.ndarray,np.ndarray,np.ndarray]:
    rng = np.random.default_rng(seed_for(root, "world", stage, arm_id(arm), rep))
    nsp, nph = 150, 300
    sid = str(arm["id"])
    normals = unit_vectors(rng, (nsp, 3))
    threshold_sd = float(arm.get("threshold_sd", 0.0))
    thresholds = rng.normal(0.0, threshold_sd, nsp) if threshold_sd > 0 else np.zeros(nsp)
    amp = float(arm.get("amplitude", arm.get("flower_amplitude", 0.0)))
    if sid == "directionally_clustered_nonshared":
        normals = clustered_axes(rng, nsp, float(arm["axis_concentration"]))
    if arm.get("kind") == "positive":
        common = unit_vectors(rng, (1,3))[0]
        chosen = rng.permutation(nsp)[:round(nsp * float(arm["shared_fraction"]))]
        normals[chosen] = common
    response_dir = unit_vectors(rng, (nsp, 3))
    fcent = rng.normal(size=(nsp,1,3))
    bcent = rng.normal(size=(nsp,1,3))
    matched = rng.normal(scale=0.5, size=(nsp,nph,3))
    fnoise = rng.normal(size=(nsp,nph,3))
    bnoise = rng.normal(size=(nsp,nph,3))
    if sid == "no_structure" or sid in {"matched_flower_background_geographic_gradient", "background_only_geographic_gradient"}:
        biological = np.zeros((nsp,nph,1))
    else:
        proj = np.einsum("snd,sd->sn", xyz, normals)
        biological = np.tanh((proj - thresholds[:,None]) / 0.15)[...,None]
    flower = fcent + amp * biological * response_dir[:,None,:] + fnoise + matched
    background = bcent + bnoise + matched
    if sid in {"matched_flower_background_geographic_gradient", "background_only_geographic_gradient"}:
        gaxis = unit_vectors(rng, (1,3))[0]
        gdir = unit_vectors(rng, (1,3))[0]
        grad = np.tanh(np.einsum("snd,d->sn", xyz, gaxis)/0.15)[...,None] * gdir[None,None,:]
        if sid == "matched_flower_background_geographic_gradient":
            scale = float(arm["matched_gradient_amplitude"])
            flower = flower + scale * grad
            background = background + scale * grad
        else:
            background = background + float(arm["background_gradient_amplitude"]) * grad
    coord = sid == "coordinate_dependent_retention_independent_boundary"
    idx = retention_indices(root, stage, rep, float(arm["retention"]), xyz, coord, float(arm.get("retention_logit_slope",0.0)))
    return flower, background, idx


def axis_log_evidence(xyz: np.ndarray, flower: np.ndarray, background: np.ndarray, idx: np.ndarray) -> np.ndarray:
    out = np.zeros((150,96), dtype=float)
    for s in range(150):
        take = idx[s]
        raw = exact_axis_scores_vectorized(xyz[s,take], flower[s,take], background[s,take])
        with np.errstate(invalid="ignore"):
            axis = np.nanmean(raw, axis=1)
        finite = np.isfinite(axis)
        if int(finite.sum()) < 2:
            continue
        centered = axis[finite] - axis[finite].mean()
        rms = float(np.sqrt(np.mean(centered**2)))
        if not np.isfinite(rms) or rms <= 1e-12:
            continue
        z = centered / rms
        lognorm = float(logsumexp(z) - np.log(len(z)))
        out[s,finite] = z - lognorm
        # non-finite axes stay neutral log LR = 0 by frozen rule.
    if not np.isfinite(out).all():
        raise RuntimeError("non-finite axis evidence")
    return out


def predictive_score(lr: np.ndarray, train_ids: np.ndarray, test_ids: np.ndarray) -> tuple[float,float,float,float]:
    tr = lr[train_ids]
    te = lr[test_ids]
    lf = np.full(len(FRACTIONS), -np.inf); np.log(FRACTIONS, out=lf, where=FRACTIONS>0)
    lnf = np.full(len(FRACTIONS), -np.inf); np.log(1-FRACTIONS, out=lnf, where=FRACTIONS<1)
    mix = np.logaddexp(lnf[:,None,None], lf[:,None,None] + tr[None,:,:])
    lp = mix.sum(axis=1)  # fraction x axis
    posterior = np.exp(lp - logsumexp(lp))
    weighted_axis = np.einsum("fk,f->k", posterior, FRACTIONS)
    independent = float(np.einsum("fk,f->", posterior, 1-FRACTIONS))
    ratio = independent + np.einsum("k,tk->t", weighted_axis, np.exp(te))
    score = float(np.mean(np.log(np.maximum(ratio, np.finfo(float).tiny))))
    mean_f = float(np.einsum("fk,f->", posterior, FRACTIONS))
    axis_mass = posterior.sum(axis=0)
    entropy = float(-np.sum(axis_mass*np.log(np.maximum(axis_mass,np.finfo(float).tiny)))/np.log(96))
    return score, mean_f, entropy, float(axis_mass.max())


def run(stage: str, arm_index: int, output: Path) -> None:
    _, mapping, _, xyz, train_ids, test_ids = load_inputs()
    root = int(mapping["execution"]["deterministic_seed_root"])
    arms = nuisance_arms(mapping) if stage == "calibration" else evaluation_arms(mapping)
    if not 0 <= arm_index < len(arms):
        raise ValueError("arm index out of range")
    arm = arms[arm_index]
    rows=[]
    t0 = __import__("time").perf_counter()
    for rep in range(REPS):
        flower, background, idx = synthetic_world(root, stage, arm, rep, xyz)
        lr = axis_log_evidence(xyz, flower, background, idx)
        score, mean_f, entropy, max_mass = predictive_score(lr, train_ids, test_ids)
        rows.append({"stage":stage,"arm_index":arm_index,"arm_id":arm_id(arm),"kind":arm["kind"],"scenario_id":arm["id"],"retention":arm["retention"],"replicate":rep,"world_score":score,"posterior_mean_sharing":mean_f,"axis_entropy":entropy,"max_axis_mass":max_mass,"amplitude":arm.get("amplitude",arm.get("flower_amplitude",np.nan)),"shared_fraction":arm.get("shared_fraction",0.0),"threshold_sd":arm.get("threshold_sd",0.0)})
        if (rep+1)%25==0:
            print(json.dumps({"stage":stage,"arm":arm_index,"completed":rep+1}), flush=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(output,index=False,lineterminator="\n")
    manifest={"protocol":mapping["protocol"],"status":"complete_v2_synthetic_arm","stage":stage,"arm_index":arm_index,"arm_id":arm_id(arm),"rows":len(rows),"seconds":float(__import__("time").perf_counter()-t0),"observed_flower_colour_opened":False,"observed_background_colour_opened":False,"image_pixels_opened":False,"lineage":{"contract_sha256":sha256_file(CONTRACT),"mapping_sha256":sha256_file(MAPPING),"gate_amendment_sha256":sha256_file(GATE_AMEND),"missing_axis_rule_sha256":sha256_file(MISSING_RULE),"geometry_result_sha256":sha256_file(GEOMETRY_RESULT),"geometry_csv_sha256":sha256_file(GEOMETRY),"csv_sha256":sha256_file(output)}}
    output.with_suffix(".json").write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")
    print(json.dumps(manifest,indent=2,sort_keys=True))


def main() -> int:
    p=argparse.ArgumentParser(); p.add_argument("--stage",choices=["calibration","evaluation"],required=True); p.add_argument("--arm-index",type=int,required=True); p.add_argument("--output",type=Path,required=True)
    a=p.parse_args(); run(a.stage,a.arm_index,a.output); return 0

if __name__=="__main__": raise SystemExit(main())

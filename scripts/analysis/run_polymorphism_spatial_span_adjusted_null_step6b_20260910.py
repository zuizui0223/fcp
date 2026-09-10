#!/usr/bin/env python3
"""Propagate the span-adjusted D-spatial statistic through frozen within-species spatial null arrays."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "polymorphism_spatial_span_adjusted_null_step6b_20260910"
OUT.mkdir(parents=True, exist_ok=True)

DISC_ATTR = ROOT / "results" / "polymorphism_species_attributes_step4_association_20260910" / "species_analysis_table.csv"
DISC_STEP6 = ROOT / "results" / "polymorphism_spatial_span_adjusted_step6_20260910" / "result.json"
RES_SPATIAL = ROOT / "results" / "polymorphism_spatial_reserve_step5b_20260909" / "reserve_species_membership_and_observed_metrics.csv"
RES_RAW = ROOT / "data" / "derived" / "rgfca_reserve_replication_measured_photos_v1.csv"

EXPECTED_DISC = 369
EXPECTED_RES = 363
EXPECTED_NULL = 999
EARTH_KM = 6371.0088
METRICS = ["primary", "observer_pair_exclusion", "calendar_quarter_stratification", "matched_background_differential"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--discovery-artifact-root", type=Path, required=True)
    p.add_argument("--reserve-artifact-root", type=Path, required=True)
    return p.parse_args()


def write_json(path: Path, obj: Any) -> None:
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def partial_spearman(d: np.ndarray, y: np.ndarray, span: np.ndarray, n: np.ndarray) -> float:
    vals = [np.asarray(v, float) for v in [d, y, span, n]]
    keep = np.logical_and.reduce([np.isfinite(v) for v in vals])
    d0, y0, s0, n0 = [v[keep] for v in vals]
    rd = stats.rankdata(d0).astype(float)
    ry = stats.rankdata(y0).astype(float)
    rs = stats.rankdata(s0).astype(float)
    rn = stats.rankdata(n0).astype(float)
    X = np.column_stack([np.ones(len(rd)), rs, rn])
    ed = rd - X @ np.linalg.lstsq(X, rd, rcond=None)[0]
    ey = ry - X @ np.linalg.lstsq(X, ry, rcond=None)[0]
    return float(np.corrcoef(ed, ey)[0, 1])


def max_great_circle_km(lat: np.ndarray, lon: np.ndarray) -> float:
    lat = np.radians(np.asarray(lat, float))
    lon = np.radians(np.asarray(lon, float))
    dlat = lat[:, None] - lat[None, :]
    dlon = lon[:, None] - lon[None, :]
    a = np.sin(dlat / 2.0) ** 2 + np.cos(lat[:, None]) * np.cos(lat[None, :]) * np.sin(dlon / 2.0) ** 2
    a = np.clip(a, 0.0, 1.0)
    return float(np.nanmax(2.0 * EARTH_KM * np.arcsin(np.sqrt(a))))


def reserve_spans(eligible: set[str]) -> pd.DataFrame:
    raw = pd.read_csv(RES_RAW, usecols=["species", "inat_taxon_id", "latitude", "longitude"])
    raw = raw[raw["species"].astype(str).isin(eligible)].copy()
    raw["latitude"] = pd.to_numeric(raw["latitude"], errors="coerce")
    raw["longitude"] = pd.to_numeric(raw["longitude"], errors="coerce")
    rows = []
    for (tid, sp), g in raw.groupby(["inat_taxon_id", "species"], sort=True):
        if len(g) != 100:
            raise RuntimeError(f"reserve raw geometry count drift for {sp}: {len(g)}")
        coords = g[["latitude", "longitude"]].to_numpy(float)
        if not np.isfinite(coords).all():
            raise RuntimeError(f"nonfinite reserve coordinates for {sp}")
        span = max_great_circle_km(coords[:, 0], coords[:, 1])
        rows.append({
            "inat_taxon_id": int(tid),
            "species": str(sp),
            "reserve_log1p_span": float(np.log1p(span)),
        })
    out = pd.DataFrame(rows)
    if len(out) != len(eligible):
        raise RuntimeError(f"reserve span count mismatch: {len(out)} vs {len(eligible)}")
    return out


def load_discovery_null(root: Path) -> tuple[pd.DataFrame, np.ndarray]:
    files = sorted(root.glob("global-rgfca-within-species-omnibus-shard-*/species_rho_permutations.csv"))
    if len(files) != 20:
        raise RuntimeError(f"expected 20 discovery shard files, found {len(files)}")
    long = pd.concat(
        [pd.read_csv(p, usecols=["species", "permutation_index", "rho"]) for p in files],
        ignore_index=True,
    )
    if long.duplicated(["species", "permutation_index"]).any():
        raise RuntimeError("duplicate discovery species/permutation rows")
    wide = long.pivot(index="species", columns="permutation_index", values="rho").sort_index()
    if len(wide) != EXPECTED_DISC or -1 not in wide.columns:
        raise RuntimeError(f"discovery wide shape/index mismatch: {wide.shape}")
    null_cols = list(range(EXPECTED_NULL))
    if not all(x in wide.columns for x in null_cols):
        raise RuntimeError("discovery null indices incomplete")
    observed = wide[-1].rename("spatial_observed_rho").reset_index()
    null = wide[null_cols].to_numpy(float)
    return observed, null


def load_reserve_null(root: Path) -> tuple[pd.DataFrame, np.ndarray]:
    shard_dirs = sorted([p for p in root.glob("reserve-inference-shard-*") if p.is_dir()])
    if len(shard_dirs) != 20:
        raise RuntimeError(f"expected 20 reserve shard dirs, found {len(shard_dirs)}")
    frames = []
    null_by_taxon: dict[int, np.ndarray] = {}
    for d in shard_dirs:
        receipt = json.loads((d / "receipt.json").read_text(encoding="utf-8"))
        if receipt.get("metrics") != METRICS or int(receipt.get("permutations", -1)) != EXPECTED_NULL:
            raise RuntimeError(f"reserve receipt drift: {d}")
        sp = pd.read_csv(d / "species.csv")
        z = np.load(d / "null.npz")
        tids = z["taxon_ids"].astype(int)
        arr = z["null"].astype(float)
        if arr.shape != (len(tids), len(METRICS), EXPECTED_NULL):
            raise RuntimeError(f"reserve null shape drift: {arr.shape}")
        frames.append(sp)
        for i, tid in enumerate(tids):
            if int(tid) in null_by_taxon:
                raise RuntimeError(f"duplicate reserve taxon {tid}")
            null_by_taxon[int(tid)] = arr[i]
    species = pd.concat(frames, ignore_index=True).sort_values("inat_taxon_id").reset_index(drop=True)
    if len(species) != EXPECTED_RES or species["inat_taxon_id"].duplicated().any():
        raise RuntimeError(f"reserve observed count drift: {len(species)}")
    tids = species["inat_taxon_id"].astype(int).to_numpy()
    if set(tids) != set(null_by_taxon):
        raise RuntimeError("reserve observed/null taxon identity mismatch")
    null = np.stack([null_by_taxon[int(t)] for t in tids], axis=0)
    return species, null


def null_test(d: np.ndarray, observed_y: np.ndarray, null_y: np.ndarray, span: np.ndarray, n: np.ndarray) -> dict[str, Any]:
    observed = partial_spearman(d, observed_y, span, n)
    vals = np.array([partial_spearman(d, null_y[:, k], span, n) for k in range(EXPECTED_NULL)], float)
    p_upper = float((1 + np.sum(vals >= observed - 1e-15)) / (EXPECTED_NULL + 1))
    return {
        "observed_partial_spearman": observed,
        "p_upper_geometry_preserving_spatial_null": p_upper,
        "null_mean": float(vals.mean()),
        "null_q025": float(np.quantile(vals, 0.025)),
        "null_q975": float(np.quantile(vals, 0.975)),
        "positive_supported_at_0_05": bool(observed > 0 and p_upper < 0.05),
        "null_values": vals,
    }


def main() -> None:
    args = parse_args()
    step6 = json.loads(DISC_STEP6.read_text(encoding="utf-8"))

    disc_attr = pd.read_csv(DISC_ATTR)
    disc_obs, disc_null = load_discovery_null(args.discovery_artifact_root)
    disc = disc_obs.merge(
        disc_attr[["species", "D", "n_classifiable", "log1p_span_primary"]],
        on="species",
        how="inner",
        validate="one_to_one",
    ).sort_values("species").reset_index(drop=True)
    if len(disc) != EXPECTED_DISC:
        raise RuntimeError("discovery join lost species")
    # load_discovery_null is species-sorted; preserve the same null order after the merge.
    disc_names = disc_obs.sort_values("species")["species"].astype(str).to_numpy()
    if not np.array_equal(disc["species"].astype(str).to_numpy(), disc_names):
        raise RuntimeError("discovery null order mismatch")

    Dd = disc["D"].to_numpy(float)
    nd = disc["n_classifiable"].to_numpy(float)
    sd = disc["log1p_span_primary"].to_numpy(float)
    yd = disc["spatial_observed_rho"].to_numpy(float)
    Dud = Dd * nd / (nd - 1.0)
    discovery = null_test(Dd, yd, disc_null, sd, nd)
    discovery_unb = null_test(Dud, yd, disc_null, sd, nd)

    if abs(discovery["observed_partial_spearman"] - float(step6["discovery"]["continuous_partial"]["partial_spearman"])) > 1e-12:
        raise RuntimeError("discovery Step6 partial anchor mismatch")

    res_metrics = pd.read_csv(RES_SPATIAL).sort_values("inat_taxon_id").reset_index(drop=True)
    rspan = reserve_spans(set(res_metrics["species"].astype(str)))
    res_species, res_null = load_reserve_null(args.reserve_artifact_root)
    reserve = res_species.merge(
        res_metrics[["inat_taxon_id", "species", "D", "n_classifiable", "rho_primary", "rho_matched_background_differential"]],
        on=["inat_taxon_id", "species"],
        how="inner",
        validate="one_to_one",
        suffixes=("_artifact", ""),
    ).merge(rspan, on=["inat_taxon_id", "species"], how="inner", validate="one_to_one")
    reserve = reserve.sort_values("inat_taxon_id").reset_index(drop=True)
    if len(reserve) != EXPECTED_RES:
        raise RuntimeError("reserve join lost species")
    if not np.array_equal(reserve["inat_taxon_id"].to_numpy(int), res_species["inat_taxon_id"].to_numpy(int)):
        raise RuntimeError("reserve null order mismatch")

    Dr = reserve["D"].to_numpy(float)
    nr = reserve["n_classifiable"].to_numpy(float)
    sr = reserve["reserve_log1p_span"].to_numpy(float)
    Dur = Dr * nr / (nr - 1.0)
    yr = reserve["rho_primary"].to_numpy(float)
    ybg = reserve["rho_matched_background_differential"].to_numpy(float)

    reserve_primary = null_test(Dr, yr, res_null[:, 0, :], sr, nr)
    reserve_primary_unb = null_test(Dur, yr, res_null[:, 0, :], sr, nr)
    reserve_bg = null_test(Dr, ybg, res_null[:, 3, :], sr, nr)
    reserve_bg_unb = null_test(Dur, ybg, res_null[:, 3, :], sr, nr)

    if abs(reserve_primary["observed_partial_spearman"] - float(step6["reserve_primary"]["continuous_partial"]["partial_spearman"])) > 1e-12:
        raise RuntimeError("reserve primary Step6 partial anchor mismatch")
    if abs(reserve_bg["observed_partial_spearman"] - float(step6["reserve_matched_background_differential"]["continuous_partial"]["partial_spearman"])) > 1e-12:
        raise RuntimeError("reserve background Step6 partial anchor mismatch")

    primary_survives_both = bool(discovery["positive_supported_at_0_05"] and reserve_primary["positive_supported_at_0_05"])
    bg_survives = bool(reserve_bg["positive_supported_at_0_05"])

    result = {
        "analysis": "polymorphism_spatial_span_adjusted_null_step6b",
        "date_jst": "2026-09-10",
        "protocol": "docs/POLYMORPHISM_SPATIAL_SPAN_ADJUSTED_NULL_STEP6B_PROTOCOL_20260910.md",
        "inferential_role": "post_step4_geometry_preserving_span_robustness",
        "upstream_discovery_run_id": 34088925008,
        "upstream_reserve_run_id": 34178957447,
        "spatial_null_permutations_reused": EXPECTED_NULL,
        "new_spatial_permutations_generated": False,
        "discovery_primary": {k: v for k, v in discovery.items() if k != "null_values"},
        "discovery_D_unbiased": {k: v for k, v in discovery_unb.items() if k != "null_values"},
        "reserve_primary": {k: v for k, v in reserve_primary.items() if k != "null_values"},
        "reserve_primary_D_unbiased": {k: v for k, v in reserve_primary_unb.items() if k != "null_values"},
        "reserve_matched_background_differential": {k: v for k, v in reserve_bg.items() if k != "null_values"},
        "reserve_matched_background_D_unbiased": {k: v for k, v in reserve_bg_unb.items() if k != "null_values"},
        "decision": {
            "adjusted_primary_gradient_supported_discovery_and_reserve": primary_survives_both,
            "adjusted_flower_specific_reserve_gradient_supported": bg_survives,
            "sampled_span_opportunity_sufficient_explanation": False if primary_survives_both else None,
            "claim_boundary": "If supported, sampled span and classifiable-photo count are insufficient to reproduce the D-spatial association under the original geometry-preserving within-species spatial null. This is not evidence of adaptation, causal range-size effects, population-genetic differentiation, or a shared boundary.",
        },
    }

    pd.DataFrame({
        "permutation_index": np.arange(EXPECTED_NULL),
        "discovery_partial_null": discovery["null_values"],
        "discovery_D_unbiased_partial_null": discovery_unb["null_values"],
    }).to_csv(OUT / "discovery_adjusted_partial_null.csv", index=False)
    pd.DataFrame({
        "permutation_index": np.arange(EXPECTED_NULL),
        "reserve_primary_partial_null": reserve_primary["null_values"],
        "reserve_primary_D_unbiased_partial_null": reserve_primary_unb["null_values"],
        "reserve_background_partial_null": reserve_bg["null_values"],
        "reserve_background_D_unbiased_partial_null": reserve_bg_unb["null_values"],
    }).to_csv(OUT / "reserve_adjusted_partial_null.csv", index=False)
    write_json(OUT / "result.json", result)

    lines = [
        "# Polymorphism spatial organization — Step 6b geometry-preserving adjusted-null result",
        "",
        f"**Adjusted primary gradient supported in discovery + reserve: `{primary_survives_both}`.**",
        f"**Adjusted reserve flower-minus-background gradient supported: `{bg_survives}`.**",
        "",
        "## D partial Spearman | sampled span + n_classifiable",
        "",
        f"- discovery: partial rho = **{discovery['observed_partial_spearman']:.6f}**, geometry-preserving p_upper = **{discovery['p_upper_geometry_preserving_spatial_null']:.6g}**",
        f"- reserve primary: partial rho = **{reserve_primary['observed_partial_spearman']:.6f}**, geometry-preserving p_upper = **{reserve_primary['p_upper_geometry_preserving_spatial_null']:.6g}**",
        f"- reserve matched-background differential: partial rho = **{reserve_bg['observed_partial_spearman']:.6f}**, geometry-preserving p_upper = **{reserve_bg['p_upper_geometry_preserving_spatial_null']:.6g}**",
        "",
        "## D_unbiased sensitivity",
        "",
        f"- discovery: partial rho = **{discovery_unb['observed_partial_spearman']:.6f}**, p_upper = **{discovery_unb['p_upper_geometry_preserving_spatial_null']:.6g}**",
        f"- reserve primary: partial rho = **{reserve_primary_unb['observed_partial_spearman']:.6f}**, p_upper = **{reserve_primary_unb['p_upper_geometry_preserving_spatial_null']:.6g}**",
        f"- reserve matched-background differential: partial rho = **{reserve_bg_unb['observed_partial_spearman']:.6f}**, p_upper = **{reserve_bg_unb['p_upper_geometry_preserving_spatial_null']:.6g}**",
        "",
        "## Interpretation boundary",
        "",
        "These P values reuse the original within-species spatial randomizations, so each species' coordinate geometry remains fixed. Support means sampled geographic opportunity plus classifiable-photo count is not sufficient to reproduce the observed polymorphism-spatial relationship. It does not establish adaptation, selection, causal range-size effects, population-genetic differentiation, or a shared global boundary.",
    ]
    (OUT / "RESULT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(json.dumps(result["decision"], indent=2))
    print((OUT / "RESULT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()

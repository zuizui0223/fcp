#!/usr/bin/env python3
"""Compute exact species-level four-state ambiguity bounds and stress-test FCP genus/spatial claims."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "polymorphism_ambiguity_bounds_step9_20260910"
OUT.mkdir(parents=True, exist_ok=True)

DISC_RAW = ROOT / "data" / "derived" / "global_monte_carlo_measured_photos_v1.csv"
RES_RAW = ROOT / "data" / "derived" / "rgfca_reserve_replication_measured_photos_v1.csv"
DISC_ATTR = ROOT / "results" / "polymorphism_species_attributes_step4_association_20260910" / "species_analysis_table.csv"
RES_ATTR = ROOT / "results" / "polymorphism_spatial_span_adjusted_step6_20260910" / "reserve_species_span_adjusted_input.csv"
STEP7_SHARED = ROOT / "results" / "polymorphism_genus_replication_step7_20260910" / "shared_repeated_genus_means.csv"

MORPHS = ["white", "yellow_orange", "red_pink", "blue_purple"]
N_PERM = 20_000
SPATIAL_NULL = 999
SEED = 20260921
EXPECTED_DISC = 369
EXPECTED_RES = 363
RES_METRICS = ["primary", "observer_pair_exclusion", "calendar_quarter_stratification", "matched_background_differential"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--discovery-artifact-root", type=Path, required=True)
    p.add_argument("--reserve-artifact-root", type=Path, required=True)
    return p.parse_args()


def write_json(path: Path, obj: Any) -> None:
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def d_from_counts(counts: np.ndarray) -> float:
    counts = np.asarray(counts, float)
    total = float(counts.sum())
    if total <= 0:
        return math.nan
    p = counts / total
    return float(1.0 - np.sum(p * p))


def completed_bounds(counts: np.ndarray, ambiguous: int) -> tuple[float, float, np.ndarray, np.ndarray]:
    c = np.asarray(counts, int).copy()
    A = int(ambiguous)
    if A < 0:
        raise RuntimeError("negative ambiguous count")
    if A == 0:
        d = d_from_counts(c)
        return d, d, c.copy(), c.copy()

    cmin = c.copy()
    cmin[int(np.argmax(cmin))] += A
    dmin = d_from_counts(cmin)

    cmax = c.copy()
    for _ in range(A):
        cmax[int(np.argmin(cmax))] += 1
    dmax = d_from_counts(cmax)
    if dmin > dmax + 1e-15:
        raise RuntimeError(f"invalid completion interval: {dmin} > {dmax}")
    return dmin, dmax, cmin, cmax


def build_endpoint_panel(path: Path, eligible: pd.DataFrame) -> pd.DataFrame:
    usecols = ["species", "inat_taxon_id", "morph", "measurement_status"]
    raw = pd.read_csv(path, usecols=usecols, low_memory=False)
    raw["species"] = raw["species"].astype(str)
    if not raw.groupby("species").size().eq(100).all():
        sizes = raw.groupby("species").size()
        raise RuntimeError(f"fixed 100-photo denominator failed in {path}: {sizes.min()}-{sizes.max()}")

    wanted = set(eligible["species"].astype(str))
    raw = raw[raw["species"].isin(wanted)].copy()
    rows: list[dict[str, Any]] = []
    for species, g in raw.groupby("species", sort=True):
        status = g["measurement_status"].astype(str)
        classifiable = status.eq("classified_four_state_morph")
        counts = np.array([int(np.sum(classifiable & g["morph"].astype(str).eq(m))) for m in MORPHS], int)
        n = int(counts.sum())
        A = int(status.eq("not_evaluable_ambiguous_palette_composition").sum())
        T = int(status.eq("not_evaluable_roi_or_flip_gate").sum())
        B = int(status.eq("not_evaluable_no_biological_palette_mass").sum())
        if n + A + T + B != 100:
            raise RuntimeError(f"status partition mismatch for {species}: {n}+{A}+{T}+{B}")
        d_obs = d_from_counts(counts)
        dmin, dmax, cmin, cmax = completed_bounds(counts, A)
        row: dict[str, Any] = {
            "species": species,
            "inat_taxon_id_endpoint": int(g["inat_taxon_id"].iloc[0]),
            "n_classifiable_reconstructed": n,
            "ambiguous_count": A,
            "technical_failure_count": T,
            "no_biological_palette_count": B,
            "ambiguous_palette_rate": A / 100.0,
            "technical_failure_rate": T / 100.0,
            "D_observed_reconstructed": d_obs,
            "D_min4": dmin,
            "D_max4": dmax,
            "D_interval_width": dmax - dmin,
        }
        for i, m in enumerate(MORPHS):
            row[f"count_{m}"] = int(counts[i])
            row[f"Dmin_count_{m}"] = int(cmin[i])
            row[f"Dmax_count_{m}"] = int(cmax[i])
        rows.append(row)
    panel = pd.DataFrame(rows)
    if len(panel) != len(eligible) or panel["species"].nunique() != len(eligible):
        raise RuntimeError(f"eligible endpoint panel mismatch: {len(panel)} != {len(eligible)}")
    return panel


def add_genus(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["genus"] = out["species"].astype(str).str.strip().str.split().str[0]
    return out


def repeated_group_pairs(labels: pd.Series) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict[str, Any]]:
    lab = labels.fillna("").astype(str).to_numpy()
    groups: dict[str, np.ndarray] = {}
    for name in sorted(set(lab)):
        if not name:
            continue
        idx = np.flatnonzero(lab == name)
        if len(idx) >= 2:
            groups[name] = idx
    species_idx = sorted({int(i) for idx in groups.values() for i in idx})
    remap = {old: new for new, old in enumerate(species_idx)}
    pi: list[int] = []
    pj: list[int] = []
    ww: list[float] = []
    ng = len(groups)
    for _name, idx in groups.items():
        npairs = len(idx) * (len(idx) - 1) // 2
        w = 1.0 / (ng * npairs)
        for a in range(len(idx) - 1):
            for b in range(a + 1, len(idx)):
                pi.append(remap[int(idx[a])])
                pj.append(remap[int(idx[b])])
                ww.append(w)
    meta = {
        "repeated_groups": int(ng),
        "species_in_repeated_groups": int(len(species_idx)),
        "gate_pass": bool(ng >= 10 and len(species_idx) >= 30),
    }
    return np.asarray(species_idx, int), np.asarray(pi, int), np.asarray(pj, int), np.asarray(ww, float), meta


def clustering_test(values: np.ndarray, labels: pd.Series, *, seed: int) -> dict[str, Any]:
    idx, pi, pj, weights, meta = repeated_group_pairs(labels)
    res = dict(meta)
    if not res["gate_pass"]:
        res.update({"tested": False, "clustering_gain": None, "p_lower": 1.0})
        return res
    z = np.asarray(values, float)[idx]
    if not np.isfinite(z).all() or not np.isclose(weights.sum(), 1.0):
        raise RuntimeError("invalid genus clustering input")
    observed = float(np.sum(np.abs(z[pi] - z[pj]) * weights))
    rng = np.random.default_rng(seed)
    null = np.empty(N_PERM, float)
    for start in range(0, N_PERM, 500):
        stop = min(start + 500, N_PERM)
        perms = np.stack([rng.permutation(len(z)) for _ in range(stop - start)])
        zz = z[perms]
        null[start:stop] = np.sum(np.abs(zz[:, pi] - zz[:, pj]) * weights[None, :], axis=1)
    nm = float(null.mean())
    gain = float(1.0 - observed / nm) if nm > 0 else math.nan
    p = float((1 + np.sum(null <= observed + 1e-15)) / (N_PERM + 1))
    res.update({
        "tested": True,
        "W_observed": observed,
        "null_mean_W": nm,
        "null_q025": float(np.quantile(null, 0.025)),
        "null_q975": float(np.quantile(null, 0.975)),
        "clustering_gain": gain,
        "p_lower": p,
        "supported_at_0_05": bool(gain > 0 and p < 0.05),
    })
    return res


def rank_residual(values: np.ndarray, controls: list[np.ndarray]) -> np.ndarray:
    y = stats.rankdata(np.asarray(values, float)).astype(float)
    cols = [np.ones(len(y), float)]
    for v in controls:
        v = np.asarray(v, float)
        if not np.isfinite(v).all():
            raise RuntimeError("nonfinite rank-control input")
        cols.append(stats.rankdata(v).astype(float))
    X = np.column_stack(cols)
    return y - X @ np.linalg.lstsq(X, y, rcond=None)[0]


def partial_rank(x: np.ndarray, y: np.ndarray, controls: list[np.ndarray]) -> float:
    ex = rank_residual(x, controls)
    ey = rank_residual(y, controls)
    return float(np.corrcoef(ex, ey)[0, 1])


def genus_endpoint_block(df: pd.DataFrame, *, span_col: str, seed_offset: int) -> dict[str, Any]:
    span = df[span_col].to_numpy(float)
    tech = df["technical_failure_rate"].to_numpy(float)
    out: dict[str, Any] = {}
    for j, col in enumerate(["D_min4", "D_max4"]):
        resid = rank_residual(df[col].to_numpy(float), [span, tech])
        out[col] = clustering_test(resid, df["genus"], seed=SEED + seed_offset + j)
    return out


def load_discovery_spatial_null(root: Path) -> tuple[pd.DataFrame, np.ndarray]:
    files = sorted(root.glob("global-rgfca-within-species-omnibus-shard-*/species_rho_permutations.csv"))
    if len(files) != 20:
        raise RuntimeError(f"expected 20 discovery spatial shards, found {len(files)}")
    long = pd.concat([pd.read_csv(p, usecols=["species", "permutation_index", "rho"]) for p in files], ignore_index=True)
    wide = long.pivot(index="species", columns="permutation_index", values="rho").sort_index()
    if len(wide) != EXPECTED_DISC or -1 not in wide.columns or not all(k in wide.columns for k in range(SPATIAL_NULL)):
        raise RuntimeError(f"discovery spatial null drift: {wide.shape}")
    return wide[-1].rename("spatial_observed_rho").reset_index(), wide[list(range(SPATIAL_NULL))].to_numpy(float)


def load_reserve_spatial_null(root: Path) -> tuple[pd.DataFrame, np.ndarray]:
    dirs = sorted([p for p in root.glob("reserve-inference-shard-*") if p.is_dir()])
    if len(dirs) != 20:
        raise RuntimeError(f"expected 20 reserve spatial shards, found {len(dirs)}")
    frames = []
    null_by_taxon: dict[int, np.ndarray] = {}
    for d in dirs:
        receipt = json.loads((d / "receipt.json").read_text(encoding="utf-8"))
        if receipt.get("metrics") != RES_METRICS or int(receipt.get("permutations", -1)) != SPATIAL_NULL:
            raise RuntimeError(f"reserve spatial receipt drift: {d}")
        sp = pd.read_csv(d / "species.csv")
        z = np.load(d / "null.npz")
        tids = z["taxon_ids"].astype(int)
        arr = z["null"].astype(float)
        if arr.shape != (len(tids), len(RES_METRICS), SPATIAL_NULL):
            raise RuntimeError(f"reserve null shape drift: {arr.shape}")
        frames.append(sp)
        for i, tid in enumerate(tids):
            null_by_taxon[int(tid)] = arr[i]
    species = pd.concat(frames, ignore_index=True).sort_values("inat_taxon_id").reset_index(drop=True)
    if len(species) != EXPECTED_RES or species["inat_taxon_id"].duplicated().any():
        raise RuntimeError("reserve spatial species drift")
    tids = species["inat_taxon_id"].astype(int).to_numpy()
    if set(tids) != set(null_by_taxon):
        raise RuntimeError("reserve spatial null identity mismatch")
    return species, np.stack([null_by_taxon[int(t)] for t in tids], axis=0)


def spatial_endpoint_test(d: np.ndarray, observed_y: np.ndarray, null_y: np.ndarray, span: np.ndarray, tech: np.ndarray) -> dict[str, Any]:
    observed = partial_rank(d, observed_y, [span, tech])
    vals = np.array([partial_rank(d, null_y[:, k], [span, tech]) for k in range(SPATIAL_NULL)], float)
    p = float((1 + np.sum(vals >= observed - 1e-15)) / (SPATIAL_NULL + 1))
    return {
        "observed_partial_spearman": observed,
        "p_upper_geometry_preserving_spatial_null": p,
        "null_mean": float(vals.mean()),
        "null_q025": float(np.quantile(vals, 0.025)),
        "null_q975": float(np.quantile(vals, 0.975)),
        "positive_supported_at_0_05": bool(observed > 0 and p < 0.05),
    }


def spatial_endpoint_block(df: pd.DataFrame, observed_y: np.ndarray, null_y: np.ndarray, *, span_col: str) -> dict[str, Any]:
    span = df[span_col].to_numpy(float)
    tech = df["technical_failure_rate"].to_numpy(float)
    return {
        col: spatial_endpoint_test(df[col].to_numpy(float), observed_y, null_y, span, tech)
        for col in ["D_min4", "D_max4"]
    }


def ambiguity_description(df: pd.DataFrame) -> dict[str, Any]:
    width = df["D_interval_width"].to_numpy(float)
    amb = df["ambiguous_palette_rate"].to_numpy(float)
    return {
        "n_species": int(len(df)),
        "fraction_A_zero": float(np.mean(df["ambiguous_count"].to_numpy(int) == 0)),
        "median_D_min4": float(df["D_min4"].median()),
        "median_D_observed": float(df["D"].median()),
        "median_D_max4": float(df["D_max4"].median()),
        "median_interval_width": float(np.median(width)),
        "q95_interval_width": float(np.quantile(width, 0.95)),
        "rho_observed_D_D_min4": float(stats.spearmanr(df["D"], df["D_min4"]).statistic),
        "rho_observed_D_D_max4": float(stats.spearmanr(df["D"], df["D_max4"]).statistic),
        "rho_interval_width_ambiguous_rate": float(stats.spearmanr(width, amb).statistic),
    }


def spearman_perm(x: np.ndarray, y: np.ndarray, *, seed: int) -> dict[str, Any]:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    rx = stats.rankdata(x).astype(float); rx -= rx.mean()
    ry = stats.rankdata(y).astype(float); ry -= ry.mean()
    denom = float(np.linalg.norm(rx) * np.linalg.norm(ry))
    observed = float(rx @ ry / denom)
    rng = np.random.default_rng(seed)
    null = np.empty(N_PERM, float)
    for i in range(N_PERM):
        null[i] = float(rx @ ry[rng.permutation(len(ry))] / denom)
    p = float((1 + np.sum(np.abs(null) >= abs(observed) - 1e-15)) / (N_PERM + 1))
    return {"n": int(len(x)), "rho": observed, "p_two_sided": p}


def shared_endpoint_concordance(discovery: pd.DataFrame, reserve: pd.DataFrame) -> dict[str, Any]:
    shared = pd.read_csv(STEP7_SHARED)
    if len(shared) != 23:
        raise RuntimeError(f"shared-genus anchor drift: {len(shared)}")
    genera = shared["genus"].astype(str).tolist()
    result: dict[str, Any] = {"shared_genera": 23}
    for j, col in enumerate(["D_min4", "D_max4"]):
        dx = []
        rx = []
        for genus in genera:
            dg = discovery[discovery["genus"].eq(genus)]
            rg = reserve[reserve["genus"].eq(genus)]
            if len(dg) < 2 or len(rg) < 2:
                raise RuntimeError(f"shared repeated genus lost: {genus}")
            dx.append(float(dg[col].mean()))
            rx.append(float(rg[col].mean()))
        result[col] = spearman_perm(np.asarray(dx), np.asarray(rx), seed=SEED + 800 + j)
    return result


def main() -> None:
    args = parse_args()
    disc_attr = add_genus(pd.read_csv(DISC_ATTR)).sort_values("species").reset_index(drop=True)
    res_attr = add_genus(pd.read_csv(RES_ATTR)).sort_values("species").reset_index(drop=True)
    if len(disc_attr) != EXPECTED_DISC or len(res_attr) != EXPECTED_RES:
        raise RuntimeError("eligible species count drift")

    dp = build_endpoint_panel(DISC_RAW, disc_attr)
    rp = build_endpoint_panel(RES_RAW, res_attr)
    discovery = disc_attr.merge(dp, on="species", validate="one_to_one").sort_values("species").reset_index(drop=True)
    reserve = res_attr.merge(rp, on="species", validate="one_to_one").sort_values("species").reset_index(drop=True)

    if np.max(np.abs(discovery["D"].to_numpy(float) - discovery["D_observed_reconstructed"].to_numpy(float))) > 1e-12:
        raise RuntimeError("discovery observed D reconstruction mismatch")
    if np.max(np.abs(reserve["D"].to_numpy(float) - reserve["D_observed_reconstructed"].to_numpy(float))) > 1e-12:
        raise RuntimeError("reserve observed D reconstruction mismatch")
    if not np.array_equal(discovery["n_classifiable"].astype(int), discovery["n_classifiable_reconstructed"].astype(int)):
        raise RuntimeError("discovery classifiable-count mismatch")
    if not np.array_equal(reserve["n_classifiable"].astype(int), reserve["n_classifiable_reconstructed"].astype(int)):
        raise RuntimeError("reserve classifiable-count mismatch")

    disc_desc = ambiguity_description(discovery)
    res_desc = ambiguity_description(reserve)
    disc_genus = genus_endpoint_block(discovery, span_col="log1p_span_primary", seed_offset=100)
    res_genus = genus_endpoint_block(reserve, span_col="reserve_log1p_span", seed_offset=200)

    disc_obs, disc_null = load_discovery_spatial_null(args.discovery_artifact_root)
    dspace = discovery.sort_values("species").reset_index(drop=True)
    disc_obs = disc_obs.sort_values("species").reset_index(drop=True)
    if not np.array_equal(dspace["species"].astype(str).to_numpy(), disc_obs["species"].astype(str).to_numpy()):
        raise RuntimeError("discovery spatial/endpoints order mismatch")
    disc_spatial = spatial_endpoint_block(
        dspace,
        disc_obs["spatial_observed_rho"].to_numpy(float),
        disc_null,
        span_col="log1p_span_primary",
    )

    res_obs, res_null = load_reserve_spatial_null(args.reserve_artifact_root)
    rspace = reserve.sort_values("inat_taxon_id_endpoint").reset_index(drop=True)
    if not np.array_equal(rspace["inat_taxon_id_endpoint"].astype(int).to_numpy(), res_obs["inat_taxon_id"].astype(int).to_numpy()):
        raise RuntimeError("reserve spatial/endpoints taxon order mismatch")
    res_primary = spatial_endpoint_block(
        rspace,
        res_obs["rho_primary"].to_numpy(float),
        res_null[:, 0, :],
        span_col="reserve_log1p_span",
    )
    res_background = spatial_endpoint_block(
        rspace,
        res_obs["rho_matched_background_differential"].to_numpy(float),
        res_null[:, 3, :],
        span_col="reserve_log1p_span",
    )

    concord = shared_endpoint_concordance(discovery, reserve)

    genus_endpoint_robust = bool(all(res_genus[c]["supported_at_0_05"] for c in ["D_min4", "D_max4"]))
    reserve_primary_endpoint_robust = bool(all(res_primary[c]["positive_supported_at_0_05"] for c in ["D_min4", "D_max4"]))
    reserve_bg_endpoint_robust = bool(all(res_background[c]["positive_supported_at_0_05"] for c in ["D_min4", "D_max4"]))

    result = {
        "analysis": "polymorphism_ambiguity_bounds_step9",
        "date_jst": "2026-09-10",
        "protocol": "docs/POLYMORPHISM_AMBIGUITY_BOUNDS_STEP9_PROTOCOL_20260910.md",
        "inferential_role": "uniform_endpoint_stress_test_under_four_state_ambiguity_completion",
        "new_image_acquisition": False,
        "new_spatial_permutations_generated": False,
        "spatial_null_permutations_reused": SPATIAL_NULL,
        "genus_permutations": N_PERM,
        "discovery_description": disc_desc,
        "reserve_description": res_desc,
        "discovery_genus": disc_genus,
        "reserve_genus": res_genus,
        "discovery_spatial": disc_spatial,
        "reserve_spatial_primary": res_primary,
        "reserve_spatial_matched_background_differential": res_background,
        "shared_23_genus_endpoint_concordance": concord,
        "decision": {
            "reserve_genus_uniform_endpoint_robust": genus_endpoint_robust,
            "reserve_primary_spatial_uniform_endpoint_robust": reserve_primary_endpoint_robust,
            "reserve_flower_minus_background_spatial_uniform_endpoint_robust": reserve_bg_endpoint_robust,
            "claim_boundary": "These are uniform D_min4/D_max4 completion stress tests. They do not exhaust arbitrary species-specific latent allocations and do not imply ambiguous rows are true four-state morphs.",
        },
    }

    discovery.to_csv(OUT / "discovery_ambiguity_endpoint_panel.csv", index=False)
    reserve.to_csv(OUT / "reserve_ambiguity_endpoint_panel.csv", index=False)
    write_json(OUT / "result.json", result)

    def gline(label: str, x: dict[str, Any]) -> str:
        return f"- {label}: gain **{x['clustering_gain']:.6f}**, p_lower **{x['p_lower']:.6g}**"
    def sline(label: str, x: dict[str, Any]) -> str:
        return f"- {label}: partial rho **{x['observed_partial_spearman']:.6f}**, p_upper **{x['p_upper_geometry_preserving_spatial_null']:.6g}**"

    lines = [
        "# Polymorphism ambiguity bounds — Step 9 result",
        "",
        f"**Reserve genus uniform-endpoint robust: `{genus_endpoint_robust}`.**",
        f"**Reserve primary spatial uniform-endpoint robust: `{reserve_primary_endpoint_robust}`.**",
        f"**Reserve flower-minus-background spatial uniform-endpoint robust: `{reserve_bg_endpoint_robust}`.**",
        "",
        "## Ambiguity interval width",
        "",
        f"- discovery median width **{disc_desc['median_interval_width']:.6f}**, 95th percentile **{disc_desc['q95_interval_width']:.6f}**, A=0 fraction **{disc_desc['fraction_A_zero']:.3f}**",
        f"- reserve median width **{res_desc['median_interval_width']:.6f}**, 95th percentile **{res_desc['q95_interval_width']:.6f}**, A=0 fraction **{res_desc['fraction_A_zero']:.3f}**",
        f"- discovery rho(observed D, D_min4) **{disc_desc['rho_observed_D_D_min4']:.6f}**, rho(observed D, D_max4) **{disc_desc['rho_observed_D_D_max4']:.6f}**",
        f"- reserve rho(observed D, D_min4) **{res_desc['rho_observed_D_D_min4']:.6f}**, rho(observed D, D_max4) **{res_desc['rho_observed_D_D_max4']:.6f}**",
        "",
        "## Genus clustering after span + technical-failure control",
        "",
        gline("discovery D_min4", disc_genus["D_min4"]),
        gline("discovery D_max4", disc_genus["D_max4"]),
        gline("reserve D_min4", res_genus["D_min4"]),
        gline("reserve D_max4", res_genus["D_max4"]),
        "",
        "## Geometry-preserving spatial endpoint tests",
        "",
        sline("discovery primary D_min4", disc_spatial["D_min4"]),
        sline("discovery primary D_max4", disc_spatial["D_max4"]),
        sline("reserve primary D_min4", res_primary["D_min4"]),
        sline("reserve primary D_max4", res_primary["D_max4"]),
        sline("reserve flower-minus-background D_min4", res_background["D_min4"]),
        sline("reserve flower-minus-background D_max4", res_background["D_max4"]),
        "",
        "## Shared 23 repeated genera",
        "",
        f"- D_min4 genus-mean concordance rho **{concord['D_min4']['rho']:.6f}**, p **{concord['D_min4']['p_two_sided']:.6g}**",
        f"- D_max4 genus-mean concordance rho **{concord['D_max4']['rho']:.6f}**, p **{concord['D_max4']['p_two_sided']:.6g}**",
        "",
        "## Interpretation boundary",
        "",
        "D_min4 and D_max4 are exact species-level diversity endpoints only under the assumption that every ambiguous-palette row belongs to one of the four existing morph classes. Downstream endpoint tests apply each endpoint rule uniformly across species; they do not exhaust arbitrary species-specific allocations or prove ambiguous rows are biological morphs.",
    ]
    (OUT / "RESULT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(result["decision"], indent=2))
    print((OUT / "RESULT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()

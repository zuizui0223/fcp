#!/usr/bin/env python3
"""Run the prospectively frozen polymorphism directionality Step 1 test.

The discrete polymorphism response is reconstructed in the 369-species discovery
cohort. M1/M2 are learned independently from the reserve cohort using the exact
zero-shift observer-equal recurrence reference procedure, then transferred to
flower-only discovery palette vectors without rotation.
"""
from __future__ import annotations

import importlib.util
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import linalg, stats

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "polymorphism_directionality_step1_20260909"
OUT.mkdir(parents=True, exist_ok=True)
DISCOVERY = ROOT / "data" / "derived" / "global_monte_carlo_measured_photos_v1.csv"
RESERVE = ROOT / "data" / "derived" / "rgfca_reserve_replication_measured_photos_v1.csv"
RECURRENCE = ROOT / "scripts" / "analysis" / "run_rgfca_global_signal_recurrence_20260909.py"
PALETTE = [
    "white", "yellow", "orange", "red", "pink", "magenta",
    "purple", "blue", "bronze", "green", "brown", "black",
]
MORPHS = ["white", "yellow_orange", "red_pink", "blue_purple"]
SEED = 20260909
N_LABEL_PERM = 20_000
N_BOOT = 5_000
N_ROT = 10_000


def load_recurrence_module():
    spec = importlib.util.spec_from_file_location("rgfca_recurrence_frozen", RECURRENCE)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load frozen recurrence module")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def bool_series(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False).astype(bool)
    return s.fillna("").astype(str).str.strip().str.lower().isin(["true", "1", "yes"])


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    return float(stats.spearmanr(a, b).statistic)


def learn_reference_axes() -> tuple[np.ndarray, dict]:
    mod = load_recurrence_module()
    df, sp, obs, values, x, y = mod.load_measurements(RESERVE)
    groups_for_shift, state_vectors, decompose = mod.make_engine(sp, obs, values, x, y)
    inv, eligible, xc, yc, state_sp = groups_for_shift(0.0, 0.0)
    state_ids, matrix = state_vectors(inv, eligible, bootstrap=False)
    n_species, n_states, eigenvalues, eigenvectors, explained, mad = decompose(
        state_ids, matrix, xc, yc, state_sp
    )
    axes = np.asarray(eigenvectors[:, :2], float)
    gram = axes.T @ axes
    if not np.allclose(gram, np.eye(2), atol=1e-10):
        raise RuntimeError(f"reference M1/M2 are not orthonormal: {gram}")
    return axes, {
        "reserve_rows_used_by_frozen_loader": int(len(df)),
        "reference_spatial_species": int(n_species),
        "reference_states": int(n_states),
        "M1_explained": float(explained[0]),
        "M2_explained": float(explained[1]),
        "M1_sum": float(axes[:, 0].sum()),
        "M2_sum": float(axes[:, 1].sum()),
    }


def load_discovery(axes: np.ndarray) -> tuple[pd.DataFrame, np.ndarray]:
    count_cols = [f"palette_count_{c}" for c in PALETTE]
    usecols = ["species", "morph", "global_classifiable"] + count_cols
    df = pd.read_csv(DISCOVERY, usecols=usecols)
    keep = bool_series(df["global_classifiable"]) & df["morph"].isin(MORPHS)
    df = df.loc[keep].copy()
    df["flower_total12"] = df[count_cols].sum(axis=1)
    if not df["flower_total12"].gt(0).all():
        raise RuntimeError("classifiable discovery rows contain non-positive flower total")
    flower = df[count_cols].to_numpy(float) / df["flower_total12"].to_numpy(float)[:, None]
    if not np.allclose(flower.sum(axis=1), 1.0, atol=1e-10):
        raise RuntimeError("flower proportions do not sum to one")
    df["M1_score"] = flower @ axes[:, 0]
    df["M2_score"] = flower @ axes[:, 1]
    return df, flower


def species_metrics(df: pd.DataFrame, flower: np.ndarray, axes: np.ndarray) -> pd.DataFrame:
    rows = []
    species = df["species"].astype(str).to_numpy()
    morph = df["morph"].astype(str).to_numpy()
    for sp in sorted(np.unique(species)):
        idx = np.flatnonzero(species == sp)
        if len(idx) < 40:
            continue
        counts = Counter(morph[idx])
        n = len(idx)
        p = np.array([counts.get(m, 0) / n for m in MORPHS], float)
        D = float(1.0 - np.sum(p * p))
        x = flower[idx]
        cov = np.cov(x, rowvar=False, ddof=1)
        total = float(np.trace(cov))
        if not np.isfinite(total) or total <= 0:
            continue
        v1 = float(axes[:, 0] @ cov @ axes[:, 0])
        v2 = float(axes[:, 1] @ cov @ axes[:, 1])
        score1 = df.iloc[idx]["M1_score"].to_numpy(float)
        score2 = df.iloc[idx]["M2_score"].to_numpy(float)
        rows.append({
            "species": sp,
            "n": n,
            "D": D,
            "second_fraction": float(np.sort(p)[-2]),
            "T": total,
            "R1": v1 / total,
            "R2": v2 / total,
            "R12": (v1 + v2) / total,
            "SD_M1": float(np.std(score1, ddof=1)),
            "SD_M2": float(np.std(score2, ddof=1)),
            "mean_M1": float(np.mean(score1)),
            "mean_M2": float(np.mean(score2)),
            **{f"p_{m}": float(p[i]) for i, m in enumerate(MORPHS)},
        })
    out = pd.DataFrame(rows).sort_values("species").reset_index(drop=True)
    if len(out) != 369:
        raise RuntimeError(f"fingerprint species mismatch: {len(out)} != 369")
    if abs(float(out["D"].max()) - 0.707645) > 1e-5:
        raise RuntimeError(f"D maximum fingerprint mismatch: {out['D'].max()}")
    s10 = float((out["second_fraction"] >= 0.10).mean())
    s20 = float((out["second_fraction"] >= 0.20).mean())
    if abs(s10 - 0.4661) > 5e-4 or abs(s20 - 0.2656) > 5e-4:
        raise RuntimeError(f"second-morph fingerprint mismatch: {s10}, {s20}")
    return out


def upper_label_permutation(d: np.ndarray, y: np.ndarray) -> tuple[float, np.ndarray]:
    rd = stats.rankdata(d).astype(float)
    ry = stats.rankdata(y).astype(float)
    rd -= rd.mean(); ry -= ry.mean()
    denom = float(np.linalg.norm(rd) * np.linalg.norm(ry))
    observed = float(rd @ ry / denom)
    rng = np.random.default_rng(SEED)
    null = np.empty(N_LABEL_PERM, float)
    for i in range(N_LABEL_PERM):
        null[i] = float(rd[rng.permutation(len(rd))] @ ry / denom)
    p = float((1 + np.sum(null >= observed)) / (N_LABEL_PERM + 1))
    return p, null


def bootstrap_rho(d: np.ndarray, y: np.ndarray) -> np.ndarray:
    rng = np.random.default_rng(SEED + 1)
    n = len(d)
    out = np.empty(N_BOOT, float)
    for i in range(N_BOOT):
        idx = rng.integers(0, n, size=n)
        out[i] = spearman(d[idx], y[idx])
    return out


def partial_spearman(d: np.ndarray, y: np.ndarray, total: np.ndarray) -> float:
    rd = stats.rankdata(d).astype(float)
    ry = stats.rankdata(y).astype(float)
    rt = stats.rankdata(np.log(total)).astype(float)
    X = np.column_stack([np.ones(len(rt)), rt])
    ed = rd - X @ np.linalg.lstsq(X, rd, rcond=None)[0]
    ey = ry - X @ np.linalg.lstsq(X, ry, rcond=None)[0]
    return float(np.corrcoef(ed, ey)[0, 1])


def random_rotation_null(d: np.ndarray, covs: np.ndarray, totals: np.ndarray) -> np.ndarray:
    # All centered compositions live in the zero-sum subspace.
    basis = linalg.null_space(np.ones((1, 12)))
    if basis.shape != (12, 11):
        raise RuntimeError(f"unexpected zero-sum basis shape {basis.shape}")
    rng = np.random.default_rng(SEED + 2)
    projections = np.empty((N_ROT, 12, 12), float)
    for b in range(N_ROT):
        q11, _ = np.linalg.qr(rng.normal(size=(11, 2)))
        a = basis @ q11[:, :2]
        projections[b] = a @ a.T

    rd = stats.rankdata(d).astype(float)
    rd -= rd.mean()
    norm_d = float(np.linalg.norm(rd))
    null = np.empty(N_ROT, float)
    chunk = 250
    for start in range(0, N_ROT, chunk):
        end = min(N_ROT, start + chunk)
        numer = np.einsum("bij,nij->bn", projections[start:end], covs, optimize=True)
        shares = numer / totals[None, :]
        ranks = stats.rankdata(shares, axis=1).astype(float)
        ranks -= ranks.mean(axis=1, keepdims=True)
        denom = norm_d * np.linalg.norm(ranks, axis=1)
        null[start:end] = (ranks @ rd) / denom
    return null


def main() -> None:
    axes, axis_meta = learn_reference_axes()
    discovery, flower = load_discovery(axes)
    sm = species_metrics(discovery, flower, axes)

    # Recover covariance matrices in the same species order for the rotation null.
    covs = []
    species_arr = discovery["species"].astype(str).to_numpy()
    for sp in sm["species"]:
        idx = np.flatnonzero(species_arr == sp)
        covs.append(np.cov(flower[idx], rowvar=False, ddof=1))
    covs = np.stack(covs)

    d = sm["D"].to_numpy(float)
    r12 = sm["R12"].to_numpy(float)
    totals = sm["T"].to_numpy(float)
    observed = spearman(d, r12)
    p_label, label_null = upper_label_permutation(d, r12)
    boots = bootstrap_rho(d, r12)
    partial = partial_spearman(d, r12, totals)
    rot_null = random_rotation_null(d, covs, totals)
    p_rotation = float((1 + np.sum(rot_null >= observed)) / (N_ROT + 1))

    diagnostics = {
        "rho_D_R1": spearman(d, sm["R1"].to_numpy(float)),
        "rho_D_R2": spearman(d, sm["R2"].to_numpy(float)),
        "rho_D_SD_M1": spearman(d, sm["SD_M1"].to_numpy(float)),
        "rho_D_SD_M2": spearman(d, sm["SD_M2"].to_numpy(float)),
        "rho_D_mean_M1": spearman(d, sm["mean_M1"].to_numpy(float)),
        "rho_D_mean_M2": spearman(d, sm["mean_M2"].to_numpy(float)),
        "rho_D_total_variance": spearman(d, totals),
    }
    primary_pass = bool(observed > 0 and p_label < 0.05)
    rotation_pass = bool(p_rotation < 0.05)
    survives = bool(primary_pass and rotation_pass)
    if survives:
        verdict = "CLAIM2_SURVIVES_DIRECTION_SPECIFIC"
    elif primary_pass:
        verdict = "CLAIM2_REJECT_OR_SUBSTANTIALLY_WEAKEN_GENERIC_ORIENTATION"
    else:
        verdict = "CLAIM2_REJECT_PRIMARY_NOT_SUPPORTED"

    result = {
        "analysis": "polymorphism_directionality_step1",
        "date_jst": "2026-09-09",
        "protocol": "docs/POLYMORPHISM_DIRECTIONALITY_STEP1_PROTOCOL_20260909.md",
        "new_image_acquisition": False,
        "discovery_species": int(len(sm)),
        "discovery_classifiable_rows": int(len(discovery)),
        "fingerprint": {
            "D_min": float(sm["D"].min()),
            "D_max": float(sm["D"].max()),
            "second_ge_0_10": float((sm["second_fraction"] >= 0.10).mean()),
            "second_ge_0_20": float((sm["second_fraction"] >= 0.20).mean()),
        },
        "axis_learning": axis_meta,
        "reference_loadings": {f"M{j+1}": {PALETTE[i]: float(axes[i, j]) for i in range(12)} for j in range(2)},
        "primary": {
            "rho_D_R12": observed,
            "label_permutation_upper_p": p_label,
            "bootstrap_95_percentile_ci": [float(np.quantile(boots, 0.025)), float(np.quantile(boots, 0.975))],
            "partial_spearman_controlling_log_T": partial,
        },
        "direction_specificity": {
            "random_2d_subspaces": N_ROT,
            "rotation_null_mean_rho": float(np.mean(rot_null)),
            "rotation_null_q025": float(np.quantile(rot_null, 0.025)),
            "rotation_null_q50": float(np.quantile(rot_null, 0.5)),
            "rotation_null_q975": float(np.quantile(rot_null, 0.975)),
            "rotation_null_max": float(np.max(rot_null)),
            "p_rotation_upper": p_rotation,
        },
        "diagnostics": diagnostics,
        "decision": {
            "primary_pass": primary_pass,
            "rotation_pass": rotation_pass,
            "claim2_survives": survives,
            "verdict": verdict,
        },
    }

    sm.to_csv(OUT / "species_metrics.csv", index=False)
    pd.DataFrame(axes, index=PALETTE, columns=["M1", "M2"]).to_csv(OUT / "reference_loadings.csv")
    pd.DataFrame({"rho": label_null}).to_csv(OUT / "label_permutation_null.csv", index=False)
    pd.DataFrame({"rho": boots}).to_csv(OUT / "species_bootstrap_rho.csv", index=False)
    pd.DataFrame({"rho": rot_null}).to_csv(OUT / "random_2d_rotation_null.csv", index=False)
    (OUT / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    md = [
        "# Polymorphism directionality Step 1 — result",
        "",
        f"**Frozen verdict: `{verdict}`.**",
        "",
        "## Primary total-variance-controlled test",
        "",
        f"- species: **{len(sm)}**",
        f"- rho(D, R12): **{observed:.6f}**",
        f"- upper-tail label-permutation p: **{p_label:.6g}**",
        f"- species-bootstrap 95% CI: **[{np.quantile(boots,0.025):.6f}, {np.quantile(boots,0.975):.6f}]**",
        f"- partial Spearman controlling rank(log T): **{partial:.6f}**",
        "",
        "## Direction specificity",
        "",
        f"- random 2D orientation null mean rho: **{np.mean(rot_null):.6f}**",
        f"- null 2.5–97.5%: **[{np.quantile(rot_null,0.025):.6f}, {np.quantile(rot_null,0.975):.6f}]**",
        f"- rotation-null p: **{p_rotation:.6g}**",
        "",
        "## Scale-sensitive diagnostics (not decision criteria)",
        "",
        f"- rho(D, SD(M1)): **{diagnostics['rho_D_SD_M1']:.6f}**",
        f"- rho(D, SD(M2)): **{diagnostics['rho_D_SD_M2']:.6f}**",
        f"- rho(D, mean M1): **{diagnostics['rho_D_mean_M1']:.6f}**",
        f"- rho(D, mean M2): **{diagnostics['rho_D_mean_M2']:.6f}**",
        f"- rho(D, total continuous variance T): **{diagnostics['rho_D_total_variance']:.6f}**",
        f"- rho(D, R1): **{diagnostics['rho_D_R1']:.6f}**",
        f"- rho(D, R2): **{diagnostics['rho_D_R2']:.6f}**",
        "",
        "## Frozen rule",
        "",
        f"- primary positive + p<0.05: **{primary_pass}**",
        f"- M1–M2 beats random 2D orientations at p<0.05: **{rotation_pass}**",
        f"- Claim 2 survives Step 1: **{survives}**",
        "",
        "No new image acquisition, axis rotation, threshold tuning or post-result species filtering was used.",
        "",
    ]
    (OUT / "RESULT.md").write_text("\n".join(md), encoding="utf-8")
    print("\n".join(md))


if __name__ == "__main__":
    main()

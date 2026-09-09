#!/usr/bin/env python3
"""Diagnose whether the frozen recurrent M3 association is flower or context driven."""
from __future__ import annotations

import importlib.util
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "polymorphism_m3_step3_20260909"
OUT.mkdir(parents=True, exist_ok=True)
DISCOVERY = ROOT / "data" / "derived" / "global_monte_carlo_measured_photos_v1.csv"
RESERVE = ROOT / "data" / "derived" / "rgfca_reserve_replication_measured_photos_v1.csv"
RECURRENCE = ROOT / "scripts" / "analysis" / "run_rgfca_global_signal_recurrence_20260909.py"
PALETTE = ["white", "yellow", "orange", "red", "pink", "magenta", "purple", "blue", "bronze", "green", "brown", "black"]
BIO = ["white", "yellow", "orange", "red", "pink", "magenta", "purple", "blue", "bronze"]
NUISANCE = ["green", "brown", "black"]
MORPHS = ["white", "yellow_orange", "red_pink", "blue_purple"]
N_PERM = 20_000
SEED = 20260909


def bool_series(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False).astype(bool)
    return s.fillna("").astype(str).str.strip().str.lower().isin(["true", "1", "yes"])


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    return float(stats.spearmanr(a, b).statistic)


def load_recurrence_module():
    spec = importlib.util.spec_from_file_location("rgfca_recurrence_frozen", RECURRENCE)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load frozen recurrence module")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def learn_m3() -> tuple[np.ndarray, dict]:
    mod = load_recurrence_module()
    df, sp, obs, values, x, y = mod.load_measurements(RESERVE)
    groups_for_shift, state_vectors, decompose = mod.make_engine(sp, obs, values, x, y)
    inv, eligible, xc, yc, state_sp = groups_for_shift(0.0, 0.0)
    state_ids, matrix = state_vectors(inv, eligible, bootstrap=False)
    n_species, n_states, _, eigenvectors, explained, _ = decompose(state_ids, matrix, xc, yc, state_sp)
    m3 = np.asarray(eigenvectors[:, 2], float)
    if not np.isclose(np.linalg.norm(m3), 1.0, atol=1e-10):
        raise RuntimeError("M3 is not unit length")
    return m3, {
        "reference_spatial_species": int(n_species),
        "reference_states": int(n_states),
        "M3_explained": float(explained[2]),
        "loadings": {PALETTE[i]: float(m3[i]) for i in range(len(PALETTE))},
    }


def discrete_D(df: pd.DataFrame, expected_n: int) -> pd.DataFrame:
    keep = bool_series(df["global_classifiable"]) & df["morph"].isin(MORPHS)
    rows = []
    for sp, g in df.loc[keep].groupby("species", sort=True):
        if len(g) < 40:
            continue
        cnt = Counter(g["morph"].astype(str))
        p = np.array([cnt.get(m, 0) / len(g) for m in MORPHS], float)
        rows.append({"species": str(sp), "D": float(1 - np.sum(p*p)), "n_classifiable": int(len(g))})
    out = pd.DataFrame(rows).sort_values("species").reset_index(drop=True)
    if len(out) != expected_n:
        raise RuntimeError(f"D frame mismatch {len(out)} != {expected_n}")
    return out


def rank_corr_perm(a: np.ndarray, b: np.ndarray, seed: int) -> tuple[float, float]:
    ra = stats.rankdata(a).astype(float); rb = stats.rankdata(b).astype(float)
    ra -= ra.mean(); rb -= rb.mean()
    denom = float(np.linalg.norm(ra) * np.linalg.norm(rb))
    obs = float(ra @ rb / denom)
    rng = np.random.default_rng(seed)
    ge = 0
    for _ in range(N_PERM):
        val = float(ra[rng.permutation(len(ra))] @ rb / denom)
        ge += int(val >= obs)
    return obs, float((1 + ge) / (N_PERM + 1))


def partial_rank_corr(a: np.ndarray, b: np.ndarray, control: np.ndarray) -> float:
    ra = stats.rankdata(a).astype(float)
    rb = stats.rankdata(b).astype(float)
    rc = stats.rankdata(control).astype(float)
    X = np.column_stack([np.ones(len(rc)), rc])
    ea = ra - X @ np.linalg.lstsq(X, ra, rcond=None)[0]
    eb = rb - X @ np.linalg.lstsq(X, rb, rcond=None)[0]
    return float(np.corrcoef(ea, eb)[0, 1])


def partial_perm_p(a: np.ndarray, b: np.ndarray, control: np.ndarray, seed: int) -> tuple[float, float]:
    ra = stats.rankdata(a).astype(float)
    rb = stats.rankdata(b).astype(float)
    rc = stats.rankdata(control).astype(float)
    X = np.column_stack([np.ones(len(rc)), rc])
    eb = rb - X @ np.linalg.lstsq(X, rb, rcond=None)[0]
    observed_ea = ra - X @ np.linalg.lstsq(X, ra, rcond=None)[0]
    observed = float(np.corrcoef(observed_ea, eb)[0, 1])
    rng = np.random.default_rng(seed)
    ge = 0
    for _ in range(N_PERM):
        rap = ra[rng.permutation(len(ra))]
        ea = rap - X @ np.linalg.lstsq(X, rap, rcond=None)[0]
        val = float(np.corrcoef(ea, eb)[0, 1])
        ge += int(val >= observed)
    return observed, float((1 + ge) / (N_PERM + 1))


def discovery_metrics(m3: np.ndarray) -> pd.DataFrame:
    count_cols = [f"palette_count_{c}" for c in PALETTE]
    use = ["species", "morph", "global_classifiable"] + count_cols
    df = pd.read_csv(DISCOVERY, usecols=use)
    dtab = discrete_D(df, 369)
    keep = bool_series(df["global_classifiable"]) & df["morph"].isin(MORPHS)
    work = df.loc[keep].copy()
    counts = work[count_cols].to_numpy(float)
    total = counts.sum(axis=1)
    if np.any(total <= 0):
        raise RuntimeError("classifiable discovery row has zero total palette mass")
    full = counts / total[:, None]
    bio_idx = [PALETTE.index(c) for c in BIO]
    nuisance_idx = [PALETTE.index(c) for c in NUISANCE]
    bio_total = counts[:, bio_idx].sum(axis=1)
    if np.any(bio_total <= 0):
        raise RuntimeError("classifiable discovery row has zero biological palette mass")
    bio12 = np.zeros_like(full)
    bio12[:, bio_idx] = counts[:, bio_idx] / bio_total[:, None]
    nuisance_fraction = counts[:, nuisance_idx].sum(axis=1) / total
    work["M3_full"] = full @ m3
    work["M3_bio9"] = bio12 @ m3
    work["nuisance_fraction"] = nuisance_fraction

    rows = []
    for sp, g in work.groupby("species", sort=True):
        rows.append({
            "species": str(sp),
            "mean_M3_full": float(g["M3_full"].mean()),
            "sd_M3_full": float(g["M3_full"].std(ddof=1)),
            "mean_M3_bio9": float(g["M3_bio9"].mean()),
            "mean_nuisance_fraction": float(g["nuisance_fraction"].mean()),
        })
    metrics = dtab.merge(pd.DataFrame(rows), on="species", how="inner", validate="one_to_one")
    if len(metrics) != 369:
        raise RuntimeError(f"discovery M3 metric frame mismatch: {len(metrics)}")
    return metrics


def reserve_metrics(m3: np.ndarray) -> pd.DataFrame:
    fc = [f"palette_count_{c}" for c in PALETTE]
    bc = [f"background_palette_count_{c}" for c in PALETTE]
    use = ["species", "morph", "global_classifiable"] + fc + bc
    df = pd.read_csv(RESERVE, usecols=use)
    dtab = discrete_D(df, 363)
    keep = bool_series(df["global_classifiable"]) & df["morph"].isin(MORPHS)
    work = df.loc[keep].copy()
    f = work[fc].to_numpy(float); b = work[bc].to_numpy(float)
    ft = f.sum(axis=1); bt = b.sum(axis=1)
    good = (ft > 0) & (bt > 0)
    work = work.loc[good].copy(); f = f[good]; b = b[good]; ft = ft[good]; bt = bt[good]
    F = f / ft[:, None]; B = b / bt[:, None]
    nidx = [PALETTE.index(c) for c in NUISANCE]
    work["F_M3"] = F @ m3
    work["B_M3"] = B @ m3
    work["Dsig_M3"] = (F - B) @ m3
    work["flower_nuisance"] = f[:, nidx].sum(axis=1) / ft
    rows = []
    for sp, g in work.groupby("species", sort=True):
        rows.append({
            "species": str(sp),
            "mean_F_M3": float(g["F_M3"].mean()),
            "mean_B_M3": float(g["B_M3"].mean()),
            "mean_Dsig_M3": float(g["Dsig_M3"].mean()),
            "mean_flower_nuisance": float(g["flower_nuisance"].mean()),
            "n_component_rows": int(len(g)),
        })
    metrics = dtab.merge(pd.DataFrame(rows), on="species", how="inner", validate="one_to_one")
    if len(metrics) != 363:
        raise RuntimeError(f"reserve M3 component frame mismatch: {len(metrics)}")
    return metrics


def main() -> None:
    m3, axis_meta = learn_m3()
    disc = discovery_metrics(m3)
    res = reserve_metrics(m3)

    d = disc["D"].to_numpy(float)
    rho_full, p_full = rank_corr_perm(d, disc["mean_M3_full"].to_numpy(float), SEED)
    rho_bio, p_bio = rank_corr_perm(d, disc["mean_M3_bio9"].to_numpy(float), SEED + 1)
    rho_partial, p_partial = partial_perm_p(d, disc["mean_M3_full"].to_numpy(float), disc["mean_nuisance_fraction"].to_numpy(float), SEED + 2)
    rho_sd = spearman(d, disc["sd_M3_full"].to_numpy(float))
    rho_nuis = spearman(d, disc["mean_nuisance_fraction"].to_numpy(float))

    rd = res["D"].to_numpy(float)
    rho_F = spearman(rd, res["mean_F_M3"].to_numpy(float))
    rho_B = spearman(rd, res["mean_B_M3"].to_numpy(float))
    rho_Dsig = spearman(rd, res["mean_Dsig_M3"].to_numpy(float))
    rho_Rnuis = spearman(rd, res["mean_flower_nuisance"].to_numpy(float))
    partial_FB = partial_rank_corr(rd, res["mean_F_M3"].to_numpy(float), res["mean_B_M3"].to_numpy(float))

    full_present = bool(rho_full > 0 and p_full < 0.05)
    bio_pass = bool(rho_bio > 0 and p_bio < 0.05)
    partial_pass = bool(rho_partial > 0 and p_partial < 0.05)
    if bio_pass and partial_pass and abs(rho_F) > abs(rho_B):
        verdict = "FLOWER_BIOLOGICAL_SUPPORTED"
    elif full_present and (not bio_pass or not partial_pass) and abs(rho_B) >= abs(rho_F):
        verdict = "BACKGROUND_CONTEXT_SUPPORTED"
    else:
        verdict = "MIXED_OR_UNRESOLVED"

    result = {
        "analysis": "polymorphism_m3_step3",
        "date_jst": "2026-09-09",
        "protocol": "docs/POLYMORPHISM_M3_STEP3_PROTOCOL_20260909.md",
        "axis": axis_meta,
        "discovery": {
            "n_species": int(len(disc)),
            "rho_D_mean_M3_full": rho_full,
            "p_mean_M3_full_upper": p_full,
            "rho_D_sd_M3_full": rho_sd,
            "rho_D_mean_M3_bio9": rho_bio,
            "p_mean_M3_bio9_upper": p_bio,
            "partial_rho_D_mean_M3_full_controlling_nuisance": rho_partial,
            "partial_p_upper": p_partial,
            "rho_D_mean_nuisance_fraction": rho_nuis,
        },
        "reserve_component_diagnostic": {
            "n_species": int(len(res)),
            "rho_D_mean_flower_M3": rho_F,
            "rho_D_mean_background_M3": rho_B,
            "rho_D_mean_flower_minus_background_M3": rho_Dsig,
            "rho_D_mean_flower_nuisance_fraction": rho_Rnuis,
            "partial_rho_D_flower_M3_controlling_background_M3": partial_FB,
        },
        "decision": {
            "full_discovery_association_present": full_present,
            "biological_only_pass": bio_pass,
            "nuisance_controlled_pass": partial_pass,
            "abs_background_ge_abs_flower_in_reserve": bool(abs(rho_B) >= abs(rho_F)),
            "verdict": verdict,
            "flowering_stage_claim_allowed": False,
        },
    }

    disc.to_csv(OUT / "discovery_species_m3_metrics.csv", index=False)
    res.to_csv(OUT / "reserve_species_m3_component_metrics.csv", index=False)
    pd.DataFrame({"palette": PALETTE, "M3": m3}).to_csv(OUT / "M3_loadings.csv", index=False)
    (OUT / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    lines = [
        "# Polymorphism M3 Step 3 — result",
        "",
        f"**Frozen diagnosis: `{verdict}`.**",
        "",
        "## Discovery (369 species; M3 learned independently in reserve)",
        "",
        f"- rho(D, mean M3 full-12): **{rho_full:.6f}**, permutation p **{p_full:.6g}**",
        f"- rho(D, SD M3 full-12): **{rho_sd:.6f}**",
        f"- rho(D, mean M3 after removing green/brown/black): **{rho_bio:.6f}**, p **{p_bio:.6g}**",
        f"- partial rho(D, mean M3 | flower nuisance fraction): **{rho_partial:.6f}**, p **{p_partial:.6g}**",
        f"- rho(D, flower nuisance fraction): **{rho_nuis:.6f}**",
        "",
        "## Reserve component diagnostic (363 species; diagnostic, not independent)",
        "",
        f"- rho(D, mean flower-only M3): **{rho_F:.6f}**",
        f"- rho(D, mean background-only M3): **{rho_B:.6f}**",
        f"- rho(D, mean flower-minus-background M3): **{rho_Dsig:.6f}**",
        f"- partial rho(D, flower M3 | background M3): **{partial_FB:.6f}**",
        f"- rho(D, flower nuisance fraction): **{rho_Rnuis:.6f}**",
        "",
        "## Interpretation boundary",
        "",
        f"- diagnosis: **{verdict}**",
        "- direct flowering-stage inference: **not allowed** (no frozen phenological-stage label)",
        "- Step 1 Claim 2 remains rejected regardless of this diagnosis.",
    ]
    (OUT / "RESULT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()

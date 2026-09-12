#!/usr/bin/env python3
"""Run the frozen H1 observer-disjoint four-state D reproducibility test."""
from __future__ import annotations

import hashlib
import importlib.util
import io
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "polymorphism_h1_observer_disjoint_d_20260913"
OUT.mkdir(parents=True, exist_ok=True)

PREFLIGHT_SCRIPT = ROOT / "scripts" / "analysis" / "run_polymorphism_h1_observer_split_preflight_20260913.py"
SOURCES = {
    "discovery": ROOT / "data" / "derived" / "global_monte_carlo_measured_photos_v1.csv",
    "reserve": ROOT / "data" / "derived" / "rgfca_reserve_replication_measured_photos_v1.csv",
}
MORPHS = ["white", "yellow_orange", "red_pink", "blue_purple"]
PRIMARY_MIN_HALF = 20
EXPECTED_SPECIES = {"discovery": 369, "reserve": 363}
EXPECTED_ASSIGNMENT_SHA256 = "f5e4333e596a69cdb032d840ec95eb02eadf1ec322658e2738e1b043fcde0cc1"
EXPECTED_PANEL_SHA256 = "bb30b6f73c6ecac06022c7ddd275cd01ea71d4bc40f54dd8889814819d1f92eb"
SEED = 20260913
N_BOOT = 5_000
N_PERM = 20_000


def load_preflight_module():
    spec = importlib.util.spec_from_file_location("h1_preflight_frozen", PREFLIGHT_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load frozen H1 preflight module")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def dataframe_csv_sha256(df: pd.DataFrame) -> str:
    buf = io.StringIO()
    df.to_csv(buf, index=False)
    return hashlib.sha256(buf.getvalue().encode("utf-8")).hexdigest()


def reconstruct_and_verify_preflight() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    mod = load_preflight_module()
    all_assign = []
    all_species = []
    cohorts: dict[str, Any] = {}
    for name, path in SOURCES.items():
        assign, panel, meta = mod.audit(name, path)
        all_assign.append(assign)
        all_species.append(panel)
        cohorts[name] = meta
    assignment = pd.concat(all_assign, ignore_index=True)
    panel = pd.concat(all_species, ignore_index=True)

    assign_sha = dataframe_csv_sha256(assignment)
    panel_sha = dataframe_csv_sha256(panel)
    if assign_sha != EXPECTED_ASSIGNMENT_SHA256:
        raise RuntimeError(f"preflight assignment hash mismatch: {assign_sha}")
    if panel_sha != EXPECTED_PANEL_SHA256:
        raise RuntimeError(f"preflight panel hash mismatch: {panel_sha}")

    # Frozen automatic rule must still select 20 before any morph column is opened.
    candidates = [20, 15, 10]
    selected = None
    for t in candidates:
        ok = all(int(cohorts[c]["eligible_by_min_half_threshold"][str(t)]) >= 100 for c in SOURCES)
        if ok:
            selected = t
            break
    if selected != PRIMARY_MIN_HALF:
        raise RuntimeError(f"frozen gate mismatch: selected {selected}, expected {PRIMARY_MIN_HALF}")
    for cohort, expected in EXPECTED_SPECIES.items():
        n = int(cohorts[cohort]["eligible_by_min_half_threshold"][str(PRIMARY_MIN_HALF)])
        if n != expected:
            raise RuntimeError(f"{cohort} primary gate fingerprint mismatch: {n} != {expected}")

    receipt = {
        "assignment_sha256": assign_sha,
        "panel_sha256": panel_sha,
        "selected_min_half": selected,
        "cohorts": cohorts,
        "morph_opened_during_preflight_verification": False,
        "D_computed_during_preflight_verification": False,
    }
    return assignment, panel, receipt


def bool_series(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False).astype(bool)
    return s.fillna("").astype(str).str.strip().str.lower().isin(["true", "1", "yes"])


def compute_half_metrics(cohort: str, path: Path, assignment: pd.DataFrame, panel: pd.DataFrame) -> pd.DataFrame:
    # Outcome opening begins here, after assignment hashes and the 20/half gate are verified.
    usecols = ["species", "observer_id", "global_classifiable", "morph"]
    df = pd.read_csv(path, usecols=usecols)
    df["species"] = df["species"].astype(str)
    df["observer_id"] = df["observer_id"].fillna("").astype(str)
    keep = bool_series(df["global_classifiable"]) & df["morph"].isin(MORPHS)
    df = df.loc[keep, ["species", "observer_id", "morph"]].copy()

    a = assignment.loc[assignment["cohort"] == cohort, ["species", "observer_id", "half"]].copy()
    a["species"] = a["species"].astype(str)
    a["observer_id"] = a["observer_id"].astype(str)
    if a.duplicated(["species", "observer_id"]).any():
        raise RuntimeError(f"{cohort}: duplicate assignment keys")
    df = df.merge(a, on=["species", "observer_id"], how="inner", validate="many_to_one")

    p = panel.loc[panel["cohort"] == cohort].copy()
    eligible = p.loc[p["n_min_half"] >= PRIMARY_MIN_HALF, "species"].astype(str)
    df = df[df["species"].isin(set(eligible))].copy()

    rows: list[dict[str, Any]] = []
    for species in sorted(eligible):
        row: dict[str, Any] = {"cohort": cohort, "species": species}
        for half in ["A", "B"]:
            x = df[(df["species"] == species) & (df["half"] == half)]
            n = int(len(x))
            if n < PRIMARY_MIN_HALF:
                raise RuntimeError(f"{cohort} {species} half {half}: n={n} below frozen gate")
            counts = x["morph"].value_counts()
            props = np.array([float(counts.get(m, 0)) / n for m in MORPHS], dtype=float)
            D = float(1.0 - np.sum(props * props))
            Dunb = float(D * n / (n - 1.0))
            row[f"n_{half}"] = n
            row[f"D_{half}"] = D
            row[f"D_unbiased_{half}"] = Dunb
            for i, morph in enumerate(MORPHS):
                row[f"p_{morph}_{half}"] = float(props[i])
        rows.append(row)
    out = pd.DataFrame(rows)
    if len(out) != EXPECTED_SPECIES[cohort]:
        raise RuntimeError(f"{cohort}: output species mismatch {len(out)} != {EXPECTED_SPECIES[cohort]}")
    return out


def spearman_fast(a: np.ndarray, b: np.ndarray) -> float:
    ra = stats.rankdata(a).astype(float)
    rb = stats.rankdata(b).astype(float)
    ra -= ra.mean()
    rb -= rb.mean()
    den = float(np.linalg.norm(ra) * np.linalg.norm(rb))
    if den == 0:
        return float("nan")
    return float(ra @ rb / den)


def lin_ccc(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    va = float(np.var(a, ddof=0))
    vb = float(np.var(b, ddof=0))
    cov = float(np.mean((a - a.mean()) * (b - b.mean())))
    den = va + vb + float((a.mean() - b.mean()) ** 2)
    return float(2.0 * cov / den) if den > 0 else float("nan")


def analyze_pair(a: np.ndarray, b: np.ndarray, seed: int) -> dict[str, Any]:
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if len(a) != len(b) or len(a) < 2 or not (np.isfinite(a).all() and np.isfinite(b).all()):
        raise RuntimeError("invalid H1 pair arrays")
    observed = spearman_fast(a, b)

    rng = np.random.default_rng(seed)
    boots = np.empty(N_BOOT, float)
    n = len(a)
    for i in range(N_BOOT):
        idx = rng.integers(0, n, size=n)
        boots[i] = spearman_fast(a[idx], b[idx])
    boots = boots[np.isfinite(boots)]
    if len(boots) < int(0.99 * N_BOOT):
        raise RuntimeError(f"too many nonfinite bootstrap correlations: {N_BOOT-len(boots)}")

    ra = stats.rankdata(a).astype(float); ra -= ra.mean()
    rb = stats.rankdata(b).astype(float); rb -= rb.mean()
    den = float(np.linalg.norm(ra) * np.linalg.norm(rb))
    null = np.empty(N_PERM, float)
    for i in range(N_PERM):
        null[i] = float(ra @ rb[rng.permutation(n)] / den)
    p = float((1 + np.sum(np.abs(null) >= abs(observed) - 1e-15)) / (N_PERM + 1))

    absdiff = np.abs(a - b)
    return {
        "n": int(n),
        "spearman_rho": float(observed),
        "bootstrap_replicates_finite": int(len(boots)),
        "bootstrap_95_percentile_ci": [
            float(np.quantile(boots, 0.025)),
            float(np.quantile(boots, 0.975)),
        ],
        "permutation_replicates": N_PERM,
        "permutation_two_sided_p": p,
        "lin_ccc": lin_ccc(a, b),
        "median_abs_difference": float(np.median(absdiff)),
        "q90_abs_difference": float(np.quantile(absdiff, 0.90)),
        "mean_A": float(np.mean(a)),
        "mean_B": float(np.mean(b)),
    }


def main() -> None:
    # Stage A: reconstruct and verify the frozen outcome-blind split before morph is opened.
    assignment, panel, preflight_receipt = reconstruct_and_verify_preflight()

    # Stage B: outcome opening.
    frames = []
    for cohort, path in SOURCES.items():
        frames.append(compute_half_metrics(cohort, path, assignment, panel))
    metrics = pd.concat(frames, ignore_index=True)

    analyses: dict[str, Any] = {}
    for j, cohort in enumerate(["discovery", "reserve"]):
        x = metrics[metrics["cohort"] == cohort].copy()
        analyses[cohort] = {
            "D": analyze_pair(x["D_A"].to_numpy(float), x["D_B"].to_numpy(float), SEED + 100 * j + 1),
            "D_unbiased": analyze_pair(
                x["D_unbiased_A"].to_numpy(float),
                x["D_unbiased_B"].to_numpy(float),
                SEED + 100 * j + 2,
            ),
        }

    reserve = analyses["reserve"]["D"]
    ci_low = float(reserve["bootstrap_95_percentile_ci"][0])
    pass_rho = bool(float(reserve["spearman_rho"]) >= 0.80)
    pass_ci = bool(ci_low > 0.70)
    pass_perm = bool(float(reserve["permutation_two_sided_p"]) < 0.001)
    supported = bool(pass_rho and pass_ci and pass_perm)
    verdict = "H1_OBSERVER_DISJOINT_D_REPRODUCIBLE" if supported else "H1_OBSERVER_DISJOINT_D_NOT_SUPPORTED"

    result = {
        "analysis": "polymorphism_h1_observer_disjoint_d",
        "date_jst": "2026-09-13",
        "protocol": "docs/POLYMORPHISM_H1_OBSERVER_DISJOINT_D_PROTOCOL_20260913.md",
        "preflight_freeze": "docs/POLYMORPHISM_H1_OBSERVER_SPLIT_PREFLIGHT_FREEZE_20260913.md",
        "new_image_acquisition": False,
        "historical_0_971_used_as_evidence": False,
        "primary_min_classifiable_per_half": PRIMARY_MIN_HALF,
        "four_states": MORPHS,
        "preflight_receipt": preflight_receipt,
        "analyses": analyses,
        "decision": {
            "primary_cohort": "reserve",
            "rho_threshold": 0.80,
            "bootstrap_lower_threshold_strictly_greater_than": 0.70,
            "permutation_p_threshold_strictly_less_than": 0.001,
            "pass_rho": pass_rho,
            "pass_bootstrap_lower": pass_ci,
            "pass_permutation": pass_perm,
            "supported": supported,
            "verdict": verdict,
        },
        "nonclaims": [
            "not an unbiased global prevalence estimate",
            "does not establish biological independence of photographs",
            "does not remove image-classification error",
            "does not identify ecological or evolutionary causes",
            "does not establish direction of colour evolution",
            "369/363 validation cohorts are not claimed representative of all 42111 species",
        ],
    }

    metrics.to_csv(OUT / "species_observer_disjoint_D.csv", index=False)
    (OUT / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    md = [
        "# H1 observer-disjoint D reproducibility — result",
        "",
        f"**Frozen verdict: `{verdict}`.**",
        "",
        f"Primary gate: >= {PRIMARY_MIN_HALF} classifiable rows in each observer-disjoint half.",
        "",
    ]
    for cohort in ["discovery", "reserve"]:
        r = analyses[cohort]["D"]
        md.extend([
            f"## {cohort}",
            "",
            f"- n species: **{r['n']}**",
            f"- Spearman rho(D_A, D_B): **{r['spearman_rho']:.6f}**",
            f"- bootstrap 95% CI: **[{r['bootstrap_95_percentile_ci'][0]:.6f}, {r['bootstrap_95_percentile_ci'][1]:.6f}]**",
            f"- permutation p: **{r['permutation_two_sided_p']:.6g}**",
            f"- Lin CCC: **{r['lin_ccc']:.6f}**",
            f"- median |D_A-D_B|: **{r['median_abs_difference']:.6f}**",
            f"- q90 |D_A-D_B|: **{r['q90_abs_difference']:.6f}**",
            "",
        ])
    (OUT / "RESULT.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

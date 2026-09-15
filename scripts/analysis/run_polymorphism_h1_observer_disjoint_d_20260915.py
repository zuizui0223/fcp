#!/usr/bin/env python3
"""Observer-disjoint reproducibility of the frozen four-state polymorphism score D.

Protocol is frozen in docs/POLYMORPHISM_H1_OBSERVER_DISJOINT_D_PROTOCOL_20260915.md.
Observer allocation never uses colour-state identity or D.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "polymorphism_h1_observer_disjoint_d_20260915"
OUT.mkdir(parents=True, exist_ok=True)

DISCOVERY = ROOT / "data" / "derived" / "global_monte_carlo_measured_photos_v1.csv"
RESERVE = ROOT / "data" / "derived" / "rgfca_reserve_replication_measured_photos_v1.csv"
PROTOCOL = "docs/POLYMORPHISM_H1_OBSERVER_DISJOINT_D_PROTOCOL_20260915.md"
MORPHS = ["white", "yellow_orange", "red_pink", "blue_purple"]
MIN_FULL = 40
MIN_HALF = 20
MIN_RESERVE_SPECIES = 100
RHO_GATE = 0.80
CCC_GATE = 0.80
N_BOOT = 10_000
N_PERM = 20_000
SEED = 20260915
SALT = "fcp-h1-observer-disjoint-d-v1"


def as_bool(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False).astype(bool)
    return s.fillna("").astype(str).str.strip().str.lower().isin(["true", "1", "yes"])


def stable_hash(cohort: str, species: str, observer: str) -> str:
    payload = f"{SALT}|{cohort}|{species}|{observer}".encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def split_observers(group: pd.DataFrame, cohort: str, species: str) -> dict[str, str]:
    """Assign whole observer groups using only group size + deterministic hash."""
    counts = group.groupby("observer_key", sort=False).size().rename("n").reset_index()
    counts["hash"] = counts["observer_key"].map(lambda x: stable_hash(cohort, species, x))
    counts = counts.sort_values(["n", "hash"], ascending=[False, True], kind="mergesort")
    total = {"A": 0, "B": 0}
    assignment: dict[str, str] = {}
    for row in counts.itertuples(index=False):
        if total["A"] < total["B"]:
            side = "A"
        elif total["B"] < total["A"]:
            side = "B"
        else:
            side = "A" if int(str(row.hash)[0], 16) % 2 == 0 else "B"
        assignment[str(row.observer_key)] = side
        total[side] += int(row.n)
    return assignment


def diversity(morphs: pd.Series) -> tuple[float, float]:
    counts = morphs.value_counts().reindex(MORPHS, fill_value=0).to_numpy(float)
    n = float(counts.sum())
    if n <= 1:
        return float("nan"), float("nan")
    p = counts / n
    d = float(1.0 - np.sum(p * p))
    return d, float(d * n / (n - 1.0))


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    rx = pd.Series(x).rank(method="average").to_numpy(float)
    ry = pd.Series(y).rank(method="average").to_numpy(float)
    if np.std(rx) == 0 or np.std(ry) == 0:
        return float("nan")
    return float(np.corrcoef(rx, ry)[0, 1])


def pearson(x: np.ndarray, y: np.ndarray) -> float:
    if np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def ccc(x: np.ndarray, y: np.ndarray) -> float:
    mx, my = float(np.mean(x)), float(np.mean(y))
    vx, vy = float(np.mean((x - mx) ** 2)), float(np.mean((y - my) ** 2))
    cov = float(np.mean((x - mx) * (y - my)))
    den = vx + vy + (mx - my) ** 2
    return float(2.0 * cov / den) if den > 0 else float("nan")


def bootstrap_pair(x: np.ndarray, y: np.ndarray, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    n = len(x)
    rho = np.empty(N_BOOT, float)
    concord = np.empty(N_BOOT, float)
    for i in range(N_BOOT):
        idx = rng.integers(0, n, size=n)
        rho[i] = spearman(x[idx], y[idx])
        concord[i] = ccc(x[idx], y[idx])
    rho = rho[np.isfinite(rho)]
    concord = concord[np.isfinite(concord)]
    return {
        "spearman_ci95": [float(np.quantile(rho, 0.025)), float(np.quantile(rho, 0.975))],
        "ccc_ci95": [float(np.quantile(concord, 0.025)), float(np.quantile(concord, 0.975))],
    }


def permutation_p(x: np.ndarray, y: np.ndarray, seed: int) -> float:
    observed = abs(spearman(x, y))
    rng = np.random.default_rng(seed)
    null = np.empty(N_PERM, float)
    for i in range(N_PERM):
        null[i] = abs(spearman(x, y[rng.permutation(len(y))]))
    return float((1 + np.sum(null >= observed - 1e-15)) / (N_PERM + 1))


def metric_summary(frame: pd.DataFrame, xcol: str, ycol: str, seed: int) -> dict:
    x = frame[xcol].to_numpy(float)
    y = frame[ycol].to_numpy(float)
    absdiff = np.abs(x - y)
    out = {
        "n": int(len(frame)),
        "spearman_rho": spearman(x, y),
        "spearman_permutation_p_two_sided": permutation_p(x, y, seed),
        "pearson_r": pearson(x, y),
        "lin_ccc": ccc(x, y),
        "mean_abs_difference": float(np.mean(absdiff)),
        "median_abs_difference": float(np.median(absdiff)),
    }
    out.update(bootstrap_pair(x, y, seed + 1000))
    return out


def load_and_split(path: Path, cohort: str, expected_full_n: int) -> tuple[pd.DataFrame, dict]:
    usecols = ["species", "morph", "global_classifiable", "observer_id"]
    df = pd.read_csv(path, usecols=usecols)
    classifiable = as_bool(df["global_classifiable"]) & df["morph"].isin(MORPHS) & df["species"].fillna("").astype(str).str.len().gt(0)
    admitted_all = df.loc[classifiable, usecols].copy()

    # Fingerprint the frozen full-depth cohort before observer missingness/splitting.
    full_counts = admitted_all.groupby("species").size()
    full_species = set(full_counts.index[full_counts >= MIN_FULL].astype(str))
    if len(full_species) != expected_full_n:
        raise RuntimeError(f"{cohort} full n>=40 fingerprint drift: {len(full_species)} != {expected_full_n}")

    admitted = admitted_all.loc[admitted_all["observer_id"].notna()].copy()
    admitted["species"] = admitted["species"].astype(str)
    admitted["observer_key"] = admitted["observer_id"].astype(str)

    rows: list[dict] = []
    reason_counts = {
        "full_n_ge_40": int(len(full_species)),
        "lt_2_observers": 0,
        "split_half_lt_20": 0,
        "evaluable": 0,
    }
    for species in sorted(full_species):
        g = admitted.loc[admitted["species"] == species].copy()
        if g["observer_key"].nunique() < 2:
            reason_counts["lt_2_observers"] += 1
            continue
        assignment = split_observers(g[["observer_key"]], cohort, species)
        g["split"] = g["observer_key"].map(assignment)
        ga = g.loc[g["split"] == "A"]
        gb = g.loc[g["split"] == "B"]
        if len(ga) < MIN_HALF or len(gb) < MIN_HALF:
            reason_counts["split_half_lt_20"] += 1
            continue
        da, dua = diversity(ga["morph"])
        db, dub = diversity(gb["morph"])
        rows.append({
            "cohort": cohort,
            "species": species,
            "n_A": int(len(ga)),
            "n_B": int(len(gb)),
            "n_observers_A": int(ga["observer_key"].nunique()),
            "n_observers_B": int(gb["observer_key"].nunique()),
            "D_A": da,
            "D_B": db,
            "D_unbiased_A": dua,
            "D_unbiased_B": dub,
        })
    out = pd.DataFrame(rows).sort_values("species").reset_index(drop=True)
    reason_counts["evaluable"] = int(len(out))
    reason_counts["fraction_of_full_depth_evaluable"] = float(len(out) / len(full_species))
    reason_counts["admitted_rows_all"] = int(len(admitted_all))
    reason_counts["admitted_rows_with_observer"] = int(len(admitted))
    reason_counts["rows_missing_observer"] = int(len(admitted_all) - len(admitted))
    return out, reason_counts


def main() -> None:
    disc, disc_gate = load_and_split(DISCOVERY, "discovery", 369)
    reserve, reserve_gate = load_and_split(RESERVE, "reserve", 363)
    if set(disc["species"]).intersection(set(reserve["species"])):
        raise RuntimeError("Discovery/reserve evaluable sets are not species-disjoint")

    disc_raw = metric_summary(disc, "D_A", "D_B", SEED + 1) if len(disc) else None
    res_raw = metric_summary(reserve, "D_A", "D_B", SEED + 2) if len(reserve) else None
    disc_unb = metric_summary(disc, "D_unbiased_A", "D_unbiased_B", SEED + 11) if len(disc) else None
    res_unb = metric_summary(reserve, "D_unbiased_A", "D_unbiased_B", SEED + 12) if len(reserve) else None

    if len(reserve) < MIN_RESERVE_SPECIES:
        verdict = "H1_OBSERVER_DISJOINT_NOT_EVALUABLE"
        supported = False
    else:
        supported = bool(res_raw["spearman_rho"] >= RHO_GATE and res_raw["lin_ccc"] >= CCC_GATE)
        verdict = (
            "H1_OBSERVER_DISJOINT_REPRODUCIBILITY_SUPPORTED"
            if supported
            else "H1_OBSERVER_DISJOINT_REPRODUCIBILITY_NOT_SUPPORTED"
        )

    combined = pd.concat([disc, reserve], ignore_index=True)
    combined.to_csv(OUT / "species_observer_disjoint_D.csv", index=False)

    result = {
        "analysis": "polymorphism_h1_observer_disjoint_d",
        "date_jst": "2026-09-15",
        "protocol": PROTOCOL,
        "new_image_acquisition": False,
        "historical_0_97_used_as_evidence": False,
        "split": {
            "salt": SALT,
            "observer_groups_kept_intact": True,
            "assignment_uses_morph_identity": False,
            "minimum_full_classifiable": MIN_FULL,
            "minimum_each_half": MIN_HALF,
        },
        "eligibility": {
            "discovery": disc_gate,
            "reserve": reserve_gate,
            "minimum_reserve_species": MIN_RESERVE_SPECIES,
        },
        "discovery_calibration": {
            "D": disc_raw,
            "D_unbiased": disc_unb,
        },
        "reserve_primary": {
            "D": res_raw,
            "D_unbiased": res_unb,
        },
        "decision_rule": {
            "reserve_spearman_min": RHO_GATE,
            "reserve_ccc_min": CCC_GATE,
            "reserve_n_min": MIN_RESERVE_SPECIES,
        },
        "decision": {
            "supported": supported,
            "verdict": verdict,
        },
        "hard_nonclaim": "Observer-disjoint reproducibility does not identify global prevalence, genetic morph discreteness, true population morph frequencies, or ignorable measurement missingness.",
    }
    (OUT / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    md = [
        "# H1 observer-disjoint D reproducibility — result",
        "",
        f"**Frozen-rule verdict: `{verdict}`.**",
        "",
        "## Reserve primary",
        "",
        f"- evaluable species: **{len(reserve)} / 363**",
    ]
    if res_raw is not None:
        md.extend([
            f"- Spearman rho(D_A, D_B): **{res_raw['spearman_rho']:.6f}**",
            f"- Lin CCC: **{res_raw['lin_ccc']:.6f}**",
            f"- 95% bootstrap rho interval: **[{res_raw['spearman_ci95'][0]:.6f}, {res_raw['spearman_ci95'][1]:.6f}]**",
            f"- 95% bootstrap CCC interval: **[{res_raw['ccc_ci95'][0]:.6f}, {res_raw['ccc_ci95'][1]:.6f}]**",
            f"- mean |D_A-D_B|: **{res_raw['mean_abs_difference']:.6f}**",
            f"- median |D_A-D_B|: **{res_raw['median_abs_difference']:.6f}**",
        ])
    md.extend([
        "",
        "## Discovery calibration",
        "",
        f"- evaluable species: **{len(disc)} / 369**",
    ])
    if disc_raw is not None:
        md.extend([
            f"- Spearman rho(D_A, D_B): **{disc_raw['spearman_rho']:.6f}**",
            f"- Lin CCC: **{disc_raw['lin_ccc']:.6f}**",
        ])
    md.extend([
        "",
        "Discovery cannot rescue a reserve failure. Species failing observer-split depth are not classified as biologically monomorphic.",
    ])
    (OUT / "README.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print(json.dumps(result["decision"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

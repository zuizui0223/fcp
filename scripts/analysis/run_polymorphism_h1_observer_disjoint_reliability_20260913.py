#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "polymorphism_h1_observer_disjoint_reliability_20260913"
DISC = ROOT / "data" / "derived" / "global_monte_carlo_measured_photos_v1.csv"
RES = ROOT / "data" / "derived" / "rgfca_reserve_replication_measured_photos_v1.csv"

MORPHS = ("white", "yellow_orange", "red_pink", "blue_purple")
MIN_FULL_CLASS = 40
MIN_HALF_PRIMARY = 20
MIN_HALF_SENS = 15
N_SPLITS = 200
BASE_SEED = 20260913
MIN_PAIRED_SPECIES = 100
MEDIAN_RHO_MIN = 2.0 / 3.0
Q05_RHO_MIN = 0.50


def as_bool(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False)
    return s.fillna("").astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def d_score(morphs: pd.Series) -> float:
    counts = morphs.value_counts().reindex(MORPHS, fill_value=0).to_numpy(dtype=float)
    n = counts.sum()
    if n <= 0:
        return float("nan")
    p = counts / n
    return float(1.0 - np.sum(p * p))


def full_d_frame(df: pd.DataFrame) -> pd.DataFrame:
    keep = as_bool(df["global_classifiable"]) & df["morph"].isin(MORPHS) & df["species"].fillna("").astype(str).str.len().gt(0)
    x = df.loc[keep, ["species", "morph"]].copy()
    tab = pd.crosstab(x["species"], pd.Categorical(x["morph"], categories=MORPHS))
    n = tab.sum(axis=1)
    tab = tab.loc[n >= MIN_FULL_CLASS]
    n = tab.sum(axis=1)
    p = tab.div(n, axis=0)
    d = 1.0 - (p * p).sum(axis=1)
    return pd.DataFrame({"species": tab.index.astype(str), "n_classifiable": n.astype(int).to_numpy(), "D": d.to_numpy(dtype=float)})


def stable_hash_int(*parts: object) -> int:
    text = " | ".join(str(x) for x in parts)
    return int(hashlib.sha256(text.encode("utf-8")).hexdigest(), 16)


def prepare_cohort(path: Path, cohort: str) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, pd.DataFrame]]:
    df = pd.read_csv(path, low_memory=False)
    required = {"species", "observer_id", "global_classifiable", "morph"}
    missing = sorted(required.difference(df.columns))
    if missing:
        raise RuntimeError(f"{cohort}: missing required columns {missing}")

    df = df.copy()
    df["species"] = df["species"].fillna("").astype(str).str.strip()
    df["observer_clean"] = df["observer_id"].fillna("").astype(str).str.strip()
    df["class4"] = as_bool(df["global_classifiable"]) & df["morph"].isin(MORPHS)
    df = df[df["species"].ne("")].copy()

    full = full_d_frame(df)

    known = df[df["observer_clean"].ne("")].copy()
    class_known = known[known["class4"]].copy()
    nclass = class_known.groupby("species").size().rename("n_classifiable_observer_known")
    nobs = known.groupby("species")["observer_clean"].nunique().rename("n_observers")
    nrows = known.groupby("species").size().rename("n_rows_observer_known")
    elig = pd.concat([nclass, nobs, nrows], axis=1).fillna(0).reset_index()
    for c in ["n_classifiable_observer_known", "n_observers", "n_rows_observer_known"]:
        elig[c] = elig[c].astype(int)
    elig["eligible"] = (elig["n_classifiable_observer_known"] >= MIN_FULL_CLASS) & (elig["n_observers"] >= 2)
    elig["cohort"] = cohort

    species_frames = {s: g.copy() for s, g in known.groupby("species", sort=False) if s in set(elig.loc[elig["eligible"], "species"])}
    return full, elig, species_frames


def observer_assignment(g: pd.DataFrame, species: str, seed: int) -> dict[str, int]:
    counts = g.groupby("observer_clean", sort=False).size().to_dict()
    items = []
    for obs, count in counts.items():
        h = stable_hash_int(seed, species, obs)
        items.append((str(obs), int(count), h))
    items.sort(key=lambda z: (-z[1], z[2]))

    totals = [0, 0]
    assignment: dict[str, int] = {}
    for obs, count, h in items:
        if totals[0] < totals[1]:
            side = 0
        elif totals[1] < totals[0]:
            side = 1
        else:
            side = h & 1
        assignment[obs] = side
        totals[side] += count

    if len(set(assignment.values())) != 2:
        raise RuntimeError(f"{species}: observer partition failed to create two sides")
    return assignment


def rank_average(x: np.ndarray) -> np.ndarray:
    return pd.Series(x).rank(method="average").to_numpy(dtype=float)


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 3:
        return float("nan")
    rx = rank_average(x)
    ry = rank_average(y)
    sx = float(np.std(rx))
    sy = float(np.std(ry))
    if sx == 0.0 or sy == 0.0:
        return float("nan")
    return float(np.corrcoef(rx, ry)[0, 1])


def ccc(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 2:
        return float("nan")
    mx, my = float(np.mean(x)), float(np.mean(y))
    vx, vy = float(np.var(x, ddof=1)), float(np.var(y, ddof=1))
    cov = float(np.cov(x, y, ddof=1)[0, 1])
    den = vx + vy + (mx - my) ** 2
    if den <= 0:
        return float("nan")
    return float(2.0 * cov / den)


def spearman_brown(rho: float) -> float:
    if not np.isfinite(rho) or rho <= -1.0:
        return float("nan")
    return float(2.0 * rho / (1.0 + rho))


def one_partition(species_frames: dict[str, pd.DataFrame], seed: int, min_half: int) -> dict[str, float]:
    xa: list[float] = []
    xb: list[float] = []
    for species, g in species_frames.items():
        assignment = observer_assignment(g, species, seed)
        side = g["observer_clean"].map(assignment).to_numpy(dtype=int)
        class4 = g["class4"].to_numpy(dtype=bool)
        morph = g["morph"].astype(str)

        mask_a = class4 & (side == 0)
        mask_b = class4 & (side == 1)
        na, nb = int(mask_a.sum()), int(mask_b.sum())
        if na < min_half or nb < min_half:
            continue
        da = d_score(morph.loc[g.index[mask_a]])
        db = d_score(morph.loc[g.index[mask_b]])
        if np.isfinite(da) and np.isfinite(db):
            xa.append(da)
            xb.append(db)

    x = np.asarray(xa, dtype=float)
    y = np.asarray(xb, dtype=float)
    rho = spearman(x, y)
    return {
        "paired_n": int(len(x)),
        "rho": rho,
        "spearman_brown": spearman_brown(rho),
        "ccc": ccc(x, y),
        "mae": float(np.mean(np.abs(x - y))) if len(x) else float("nan"),
        "bias_a_minus_b": float(np.mean(x - y)) if len(x) else float("nan"),
    }


def q(values: Iterable[float], p: float) -> float:
    x = np.asarray([v for v in values if np.isfinite(v)], dtype=float)
    return float(np.quantile(x, p)) if len(x) else float("nan")


def summarize(metrics: pd.DataFrame) -> dict:
    return {
        "partitions": int(len(metrics)),
        "partitions_with_defined_rho": int(metrics["rho"].notna().sum()),
        "paired_n_median": float(metrics["paired_n"].median()),
        "paired_n_min": int(metrics["paired_n"].min()),
        "paired_n_max": int(metrics["paired_n"].max()),
        "rho_median": float(metrics["rho"].median()),
        "rho_q05": q(metrics["rho"], 0.05),
        "rho_q95": q(metrics["rho"], 0.95),
        "spearman_brown_median": float(metrics["spearman_brown"].median()),
        "spearman_brown_q05": q(metrics["spearman_brown"], 0.05),
        "ccc_median": float(metrics["ccc"].median()),
        "ccc_q05": q(metrics["ccc"], 0.05),
        "ccc_q95": q(metrics["ccc"], 0.95),
        "mae_median": float(metrics["mae"].median()),
        "mae_q95": q(metrics["mae"], 0.95),
        "bias_median": float(metrics["bias_a_minus_b"].median()),
        "bias_abs_q95": q(metrics["bias_a_minus_b"].abs(), 0.95),
    }


def passes_primary(summary: dict) -> bool:
    return bool(
        summary["paired_n_median"] >= MIN_PAIRED_SPECIES
        and summary["rho_median"] >= MEDIAN_RHO_MIN
        and summary["rho_q05"] >= Q05_RHO_MIN
    )


def analyze_cohort(species_frames: dict[str, pd.DataFrame], cohort: str) -> tuple[pd.DataFrame, dict, dict]:
    rows = []
    for i in range(1, N_SPLITS + 1):
        seed = BASE_SEED + i
        for gate_name, min_half in [("primary20", MIN_HALF_PRIMARY), ("sensitivity15", MIN_HALF_SENS)]:
            m = one_partition(species_frames, seed, min_half)
            rows.append({"cohort": cohort, "split_index": i, "seed": seed, "gate": gate_name, "min_half_classifiable": min_half, **m})
    frame = pd.DataFrame(rows)
    primary = summarize(frame[frame["gate"] == "primary20"].copy())
    sensitivity = summarize(frame[frame["gate"] == "sensitivity15"].copy())
    return frame, primary, sensitivity


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    expected_sha = {
        "discovery": "ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4",
        "reserve": "0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6",
    }
    observed_sha = {"discovery": sha256_file(DISC), "reserve": sha256_file(RES)}
    if observed_sha != expected_sha:
        raise RuntimeError(f"Frozen measured-photo SHA256 drift: {observed_sha}")

    full_d, elig_d, frames_d = prepare_cohort(DISC, "discovery")
    full_r, elig_r, frames_r = prepare_cohort(RES, "reserve")
    if len(full_d) != 369:
        raise RuntimeError(f"Discovery full-D fingerprint drift: {len(full_d)}")
    if len(full_r) != 363:
        raise RuntimeError(f"Reserve full-D fingerprint drift: {len(full_r)}")
    if set(full_d["species"]).intersection(set(full_r["species"])):
        raise RuntimeError("Discovery and reserve cohorts are not species-disjoint")
    if abs(float(full_d["D"].max()) - 0.707645) > 1e-5:
        raise RuntimeError("Discovery D maximum fingerprint drift")

    metrics_d, primary_d, sens_d = analyze_cohort(frames_d, "discovery")
    metrics_r, primary_r, sens_r = analyze_cohort(frames_r, "reserve")
    metrics = pd.concat([metrics_d, metrics_r], ignore_index=True)

    reserve_pass = passes_primary(primary_r)
    discovery_pass = passes_primary(primary_d)
    if primary_r["paired_n_median"] < MIN_PAIRED_SPECIES:
        verdict = "H1_OBSERVER_DISJOINT_D_RELIABILITY_UNDERIDENTIFIED"
    elif reserve_pass:
        verdict = "H1_OBSERVER_DISJOINT_D_RELIABILITY_SUPPORTED"
    else:
        verdict = "H1_OBSERVER_DISJOINT_D_RELIABILITY_NOT_SUPPORTED"
    label = verdict + ("+DISCOVERY_CONSISTENT" if reserve_pass and discovery_pass else "")

    OUT.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(OUT / "partition_metrics.csv", index=False)
    pd.concat([elig_d, elig_r], ignore_index=True).to_csv(OUT / "species_observer_eligibility.csv", index=False)

    result = {
        "analysis": "polymorphism_h1_observer_disjoint_D_reliability",
        "date_jst": "2026-09-13",
        "protocol": "docs/POLYMORPHISM_H1_OBSERVER_DISJOINT_RELIABILITY_PROTOCOL_20260913.md",
        "role": "measurement_and_sampling_stability_not_biological_prevalence",
        "frozen_inputs_sha256": observed_sha,
        "four_states": list(MORPHS),
        "mixed_uncertain_promoted": False,
        "n_partitions": N_SPLITS,
        "partition_seed_base": BASE_SEED,
        "eligibility": {
            "full_D_fingerprint_discovery": int(len(full_d)),
            "full_D_fingerprint_reserve": int(len(full_r)),
            "observer_known_eligible_discovery": int(elig_d["eligible"].sum()),
            "observer_known_eligible_reserve": int(elig_r["eligible"].sum()),
            "minimum_full_classifiable_observer_known": MIN_FULL_CLASS,
            "minimum_distinct_observers": 2,
        },
        "primary_gate": {
            "minimum_classifiable_per_half": MIN_HALF_PRIMARY,
            "minimum_median_paired_species": MIN_PAIRED_SPECIES,
            "minimum_median_split_rho": MEDIAN_RHO_MIN,
            "minimum_q05_split_rho": Q05_RHO_MIN,
        },
        "discovery": {"primary20": primary_d, "sensitivity15": sens_d, "passes_primary_rule": discovery_pass},
        "reserve": {"primary20": primary_r, "sensitivity15": sens_r, "passes_primary_rule": reserve_pass},
        "decision": {
            "reserve_is_primary": True,
            "verdict": verdict,
            "label": label,
            "discovery_consistent": bool(discovery_pass),
            "sensitivity_can_rescue_primary": False,
        },
        "hard_nonclaims": [
            "not a global polymorphism prevalence estimate",
            "not proof of globally representative species sampling",
            "not validation of image-level state classification correctness",
            "not evidence against geographic morph-frequency structure",
            "not a causal ecological or phylogenetic result",
            "not an H2 geometry test",
        ],
    }
    (OUT / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    md = [
        "# H1 observer-disjoint D reliability",
        "",
        f"Verdict: **{label}**",
        "",
        "## Reserve primary",
        "",
        f"- observer-known eligible species before half-size gate: {result['eligibility']['observer_known_eligible_reserve']}",
        f"- median paired species per split: {primary_r['paired_n_median']:.0f}",
        f"- median split-half Spearman rho: {primary_r['rho_median']:.6f}",
        f"- 5th-95th percentile rho: {primary_r['rho_q05']:.6f} to {primary_r['rho_q95']:.6f}",
        f"- median Spearman-Brown projected full-estimate reliability: {primary_r['spearman_brown_median']:.6f}",
        f"- median CCC: {primary_r['ccc_median']:.6f}",
        f"- median absolute D difference: {primary_r['mae_median']:.6f}",
        "",
        "## Discovery diagnostic",
        "",
        f"- observer-known eligible species before half-size gate: {result['eligibility']['observer_known_eligible_discovery']}",
        f"- median paired species per split: {primary_d['paired_n_median']:.0f}",
        f"- median split-half Spearman rho: {primary_d['rho_median']:.6f}",
        f"- 5th-95th percentile rho: {primary_d['rho_q05']:.6f} to {primary_d['rho_q95']:.6f}",
        f"- median Spearman-Brown projected full-estimate reliability: {primary_d['spearman_brown_median']:.6f}",
        "",
        "The >=15-per-half sensitivity is diagnostic only and cannot rescue the primary gate.",
    ]
    (OUT / "RESULT.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print(json.dumps(result["decision"], indent=2, sort_keys=True))
    print(json.dumps({"discovery_primary20": primary_d, "reserve_primary20": primary_r}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

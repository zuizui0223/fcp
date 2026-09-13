#!/usr/bin/env python3
"""Canonical observer-disjoint reliability test for the current four-state D."""
from __future__ import annotations

import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "polymorphism_h1_observer_disjoint_D_reliability_20260913"
DISCOVERY = ROOT / "data" / "derived" / "global_monte_carlo_measured_photos_v1.csv"
RESERVE = ROOT / "data" / "derived" / "rgfca_reserve_replication_measured_photos_v1.csv"
MORPHS = ["white", "yellow_orange", "red_pink", "blue_purple"]
SALT = "h1-observer-disjoint-v1"
SEED = 20260913
N_BOOT = 10_000
FULL_MIN = 40
HALF_MIN = 20


def bool_series(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False).astype(bool)
    return s.fillna("").astype(str).str.strip().str.lower().isin(["true", "1", "yes"])


def valid_observer_series(s: pd.Series) -> pd.Series:
    z = s.fillna("").astype(str).str.strip()
    return z.ne("") & ~z.str.lower().isin(["nan", "none", "null", "na"])


def stable_tie(species: str, observer: str) -> str:
    return hashlib.sha256(f"{SALT}|{species}|{observer}".encode()).hexdigest()


def d_score(morphs) -> tuple[float, float, dict[str, int]]:
    vals = [str(x) for x in morphs]
    counts = Counter(vals)
    n = len(vals)
    if n == 0:
        return float("nan"), float("nan"), {m: 0 for m in MORPHS}
    p = np.array([counts.get(m, 0) / n for m in MORPHS], float)
    return float(1 - np.sum(p * p)), float(np.sort(p)[-2]), {m: int(counts.get(m, 0)) for m in MORPHS}


def assign_observers(species: str, df: pd.DataFrame):
    counts = df.groupby("observer_id", sort=False).size().to_dict()
    ordered = sorted(((str(o), int(n)) for o, n in counts.items()), key=lambda x: (-x[1], stable_tie(species, x[0])))
    assignment = {}
    na = nb = 0
    for obs, n in ordered:
        if na <= nb:
            assignment[obs] = "A"; na += n
        else:
            assignment[obs] = "B"; nb += n
    return assignment, na, nb


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 2:
        return float("nan")
    return float(stats.spearmanr(x, y).statistic)


def pearson(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 2 or np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def ccc(x: np.ndarray, y: np.ndarray) -> float:
    mx, my = float(np.mean(x)), float(np.mean(y))
    vx, vy = float(np.var(x)), float(np.var(y))
    cov = float(np.mean((x - mx) * (y - my)))
    den = vx + vy + (mx - my) ** 2
    return float(2 * cov / den) if den > 0 else float("nan")


def kappa_binary(a: np.ndarray, b: np.ndarray) -> dict:
    a = np.asarray(a, bool); b = np.asarray(b, bool)
    po = float(np.mean(a == b)); pa = float(np.mean(a)); pb = float(np.mean(b))
    pe = pa * pb + (1 - pa) * (1 - pb)
    return {"agreement": po, "kappa": float((po - pe) / (1 - pe)) if pe < 1 else float("nan"),
            "A_positive_fraction": pa, "B_positive_fraction": pb}


def bootstrap_spearman(x: np.ndarray, y: np.ndarray, seed: int):
    rng = np.random.default_rng(seed); n = len(x); vals = []
    for _ in range(N_BOOT):
        idx = rng.integers(0, n, n)
        r = spearman(x[idx], y[idx])
        if np.isfinite(r): vals.append(r)
    if not vals:
        return [float("nan"), float("nan")], 0
    a = np.asarray(vals)
    return [float(np.quantile(a, .025)), float(np.quantile(a, .975))], int(len(a))


def analyze(path: Path, cohort: str, expected_full: int, seed_offset: int):
    raw = pd.read_csv(path, usecols=["species", "observer_id", "morph", "global_classifiable"], low_memory=False)
    raw["species"] = raw["species"].fillna("").astype(str).str.strip()
    all_species = sorted(x for x in raw["species"].unique() if x)
    admitted = raw.loc[bool_series(raw["global_classifiable"]) & raw["morph"].isin(MORPHS)].copy()
    admitted["observer_valid"] = valid_observer_series(admitted["observer_id"])
    admitted.loc[admitted["observer_valid"], "observer_id"] = admitted.loc[admitted["observer_valid"], "observer_id"].astype(str).str.strip()

    rows = []; attrition = Counter(); full_depth_n = 0
    for sp in all_species:
        sdf = admitted.loc[admitted["species"] == sp].copy(); nfull = len(sdf)
        if nfull < FULL_MIN:
            attrition["full_n_classifiable_lt_40"] += 1; continue
        full_depth_n += 1
        Dfull, secondfull, fullcounts = d_score(sdf["morph"])
        valid = sdf.loc[sdf["observer_valid"]].copy()
        nmiss = nfull - len(valid); nobs = valid["observer_id"].nunique()
        if nobs < 2:
            attrition["fewer_than_2_valid_observers"] += 1; continue
        assignment, na0, nb0 = assign_observers(sp, valid)
        valid["split"] = valid["observer_id"].map(assignment)
        A = valid.loc[valid["split"] == "A"]; B = valid.loc[valid["split"] == "B"]
        na, nb = len(A), len(B)
        if (na, nb) != (na0, nb0): raise RuntimeError(f"split accounting mismatch: {sp}")
        if set(A["observer_id"]) & set(B["observer_id"]): raise RuntimeError(f"observer leakage: {sp}")
        if na < HALF_MIN or nb < HALF_MIN:
            attrition["observer_disjoint_half_lt_20"] += 1; continue
        DA, sA, cA = d_score(A["morph"]); DB, sB, cB = d_score(B["morph"])
        rows.append({"cohort": cohort, "species": sp, "n_full_classifiable": nfull, "n_missing_observer": nmiss,
                     "n_valid_observers": int(nobs), "n_A": na, "n_B": nb, "split_imbalance": abs(na-nb)/(na+nb),
                     "D_full": Dfull, "D_A": DA, "D_B": DB, "abs_D_difference": abs(DA-DB),
                     "signed_D_difference_A_minus_B": DA-DB, "second_fraction_full": secondfull,
                     "second_fraction_A": sA, "second_fraction_B": sB,
                     "H1_10_A": sA >= .10, "H1_10_B": sB >= .10, "H1_20_A": sA >= .20, "H1_20_B": sB >= .20,
                     **{f"full_{m}": fullcounts[m] for m in MORPHS}, **{f"A_{m}": cA[m] for m in MORPHS}, **{f"B_{m}": cB[m] for m in MORPHS}})
    if full_depth_n != expected_full: raise RuntimeError(f"{cohort} four-state fingerprint {full_depth_n} != {expected_full}")
    out = pd.DataFrame(rows).sort_values("species").reset_index(drop=True)
    if len(out) < 2: raise RuntimeError(f"{cohort}: reliability underidentified")
    x = out["D_A"].to_numpy(float); y = out["D_B"].to_numpy(float); ad = out["abs_D_difference"].to_numpy(float)
    rho = spearman(x, y); ci, nboot = bootstrap_spearman(x, y, SEED + seed_offset)
    slope, intercept = np.polyfit(x, y, 1) if np.std(x) > 0 else (float("nan"), float("nan"))
    metrics = {"cohort": cohort, "species_total_in_source": len(all_species), "full_n_classifiable_ge_40": full_depth_n,
               "reliability_eligible_n": len(out), "attrition": dict(sorted(attrition.items())),
               "spearman_D_A_D_B": rho, "spearman_bootstrap_95_percentile_ci": ci,
               "bootstrap_replicates_requested": N_BOOT, "bootstrap_replicates_finite": nboot,
               "pearson_D_A_D_B": pearson(x,y), "lin_ccc": ccc(x,y), "median_abs_D_difference": float(np.median(ad)),
               "q90_abs_D_difference": float(np.quantile(ad,.90)), "mean_signed_D_difference_A_minus_B": float(np.mean(x-y)),
               "ols_D_B_on_D_A": {"slope": float(slope), "intercept": float(intercept)},
               "H1_gate_stability_0_10": kappa_binary(out["H1_10_A"], out["H1_10_B"]),
               "H1_gate_stability_0_20": kappa_binary(out["H1_20_A"], out["H1_20_B"]),
               "discrepancy_opportunity_diagnostics": {"rho_absdiff_full_n_classifiable": spearman(ad,out["n_full_classifiable"].to_numpy(float)),
                   "rho_absdiff_n_valid_observers": spearman(ad,out["n_valid_observers"].to_numpy(float)),
                   "rho_absdiff_split_imbalance": spearman(ad,out["split_imbalance"].to_numpy(float))},
               "observer_leakage_count": 0}
    return out, metrics


def sha256_file(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1024*1024), b""): h.update(chunk)
    return h.hexdigest()


def main():
    expected={"discovery":"ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4",
              "reserve":"0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6"}
    observed={"discovery":sha256_file(DISCOVERY),"reserve":sha256_file(RESERVE)}
    if observed != expected: raise RuntimeError(f"source SHA drift: {observed}")
    dr, d = analyze(DISCOVERY,"discovery",369,0); rr, r = analyze(RESERVE,"reserve",363,1)
    lower=r["spearman_bootstrap_95_percentile_ci"][0]
    g1=bool(np.isfinite(r["spearman_D_A_D_B"]) and r["spearman_D_A_D_B"]>=.80)
    g2=bool(np.isfinite(lower) and lower>=.70); g3=bool(np.isfinite(r["lin_ccc"]) and r["lin_ccc"]>=.75)
    supported=g1 and g2 and g3
    verdict="H1_OBSERVER_DISJOINT_D_RELIABILITY_SUPPORTED" if supported else "H1_OBSERVER_DISJOINT_D_RELIABILITY_NOT_SUPPORTED"
    result={"analysis":"polymorphism_h1_observer_disjoint_D_reliability","date_jst":"2026-09-13",
            "protocol":"docs/POLYMORPHISM_H1_OBSERVER_DISJOINT_D_RELIABILITY_PROTOCOL_20260913.md",
            "historical_approx_0_971_used_as_evidence":False,"frozen_inputs_sha256":observed,"state_definition":MORPHS,
            "split":{"salt":SALT,"observer_disjoint":True,"assignment_uses_morph_identity":False,
                     "assignment_rule":"descending admitted-row count; tie SHA256(salt|species|observer); greedy to smaller half; A on tie",
                     "full_min":FULL_MIN,"half_min":HALF_MIN},"discovery":d,"reserve":r,
            "decision":{"decision_cohort":"reserve","rho_floor_0_80_pass":g1,"bootstrap_lower_floor_0_70_pass":g2,
                        "ccc_floor_0_75_pass":g3,"supported":supported,"verdict":verdict},
            "nonclaims":["not a global polymorphism prevalence estimate","not an untouched prospective biological replication",
                         "excluded species are not called monomorphic or unreliable","does not validate mixed_uncertain as a biological morph"]}
    OUT.mkdir(parents=True, exist_ok=True)
    pd.concat([dr,rr],ignore_index=True).to_csv(OUT/"species_split_metrics.csv",index=False)
    (OUT/"result.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__=="__main__": main()

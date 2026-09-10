#!/usr/bin/env python3
"""Decompose sampled-span and classification-yield adjustments for genus clustering."""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "polymorphism_genus_adjustment_decomposition_step7b_20260910"
OUT.mkdir(parents=True, exist_ok=True)

DISC = ROOT / "results" / "polymorphism_species_attributes_step4_association_20260910" / "species_analysis_table.csv"
RES = ROOT / "results" / "polymorphism_spatial_span_adjusted_step6_20260910" / "reserve_species_span_adjusted_input.csv"
STEP7 = ROOT / "results" / "polymorphism_genus_replication_step7_20260910" / "result.json"
STEP7_SHARED = ROOT / "results" / "polymorphism_genus_replication_step7_20260910" / "shared_repeated_genus_means.csv"
STEP2 = ROOT / "results" / "polymorphism_mixed_step2_20260909" / "result.json"

N_PERM = 20_000
SEED = 20260918
EXPECTED_DISCOVERY = 369
EXPECTED_RESERVE = 363


def write_json(path: Path, obj: Any) -> None:
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def add_genus(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    derived = out["species"].astype(str).str.strip().str.split().str[0]
    if "genus" in out.columns:
        stored = out["genus"].astype(str).str.strip()
        bad = stored.ne(derived)
        if bad.any():
            rows = out.loc[bad, ["species", "genus"]]
            raise RuntimeError(f"genus mismatch: {rows.head().to_dict('records')}")
    out["genus"] = derived
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
    if ng:
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
        "group_sizes": {k: int(len(v)) for k, v in groups.items()},
        "gate_pass": bool(ng >= 10 and len(species_idx) >= 30),
    }
    return np.asarray(species_idx, int), np.asarray(pi, int), np.asarray(pj, int), np.asarray(ww, float), meta


def cluster(values: np.ndarray, labels: pd.Series, *, seed: int) -> dict[str, Any]:
    idx, pi, pj, weights, meta = repeated_group_pairs(labels)
    res = dict(meta)
    if not res["gate_pass"]:
        res.update({"tested": False, "W_observed": None, "clustering_gain": None, "p_lower": 1.0})
        return res
    z = np.asarray(values, float)[idx]
    if not np.isfinite(z).all() or not np.isclose(weights.sum(), 1.0):
        raise RuntimeError("invalid genus clustering inputs")
    observed = float(np.sum(np.abs(z[pi] - z[pj]) * weights))
    rng = np.random.default_rng(seed)
    null = np.empty(N_PERM, float)
    for start in range(0, N_PERM, 500):
        stop = min(start + 500, N_PERM)
        perm = np.stack([rng.permutation(len(z)) for _ in range(stop - start)])
        zz = z[perm]
        null[start:stop] = np.sum(np.abs(zz[:, pi] - zz[:, pj]) * weights[None, :], axis=1)
    null_mean = float(null.mean())
    gain = float(1.0 - observed / null_mean) if null_mean > 0 else math.nan
    p = float((1 + np.sum(null <= observed + 1e-15)) / (N_PERM + 1))
    res.update({
        "tested": True,
        "W_observed": observed,
        "null_mean_W": null_mean,
        "null_q025": float(np.quantile(null, 0.025)),
        "null_q975": float(np.quantile(null, 0.975)),
        "clustering_gain": gain,
        "p_lower": p,
        "supported_at_0_05": bool(gain > 0 and p < 0.05),
    })
    return res


def residualize_rank(outcome: np.ndarray, controls: list[np.ndarray]) -> np.ndarray:
    y = stats.rankdata(np.asarray(outcome, float)).astype(float)
    cols = [np.ones(len(y), float)]
    for v in controls:
        v = np.asarray(v, float)
        if not np.isfinite(v).all():
            raise RuntimeError("nonfinite adjustment control")
        cols.append(stats.rankdata(v).astype(float))
    X = np.column_stack(cols)
    return y - X @ np.linalg.lstsq(X, y, rcond=None)[0]


def spearman_perm(x: np.ndarray, y: np.ndarray, *, seed: int) -> dict[str, Any]:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    rx = stats.rankdata(x).astype(float)
    ry = stats.rankdata(y).astype(float)
    rx -= rx.mean()
    ry -= ry.mean()
    denom = float(np.linalg.norm(rx) * np.linalg.norm(ry))
    observed = float(rx @ ry / denom)
    rng = np.random.default_rng(seed)
    null = np.empty(N_PERM, float)
    for i in range(N_PERM):
        null[i] = float(rx @ ry[rng.permutation(len(ry))] / denom)
    p = float((1 + np.sum(np.abs(null) >= abs(observed) - 1e-15)) / (N_PERM + 1))
    return {
        "n": int(len(x)),
        "rho": observed,
        "p_two_sided": p,
        "null_mean": float(null.mean()),
        "null_q025": float(np.quantile(null, 0.025)),
        "null_q975": float(np.quantile(null, 0.975)),
    }


def tranche_tests(df: pd.DataFrame, *, span_col: str, seed_offset: int) -> dict[str, Any]:
    D = df["D"].to_numpy(float)
    n = df["n_classifiable"].to_numpy(float)
    span = df[span_col].to_numpy(float)
    Dunb = D * n / (n - 1.0)
    genus = df["genus"]

    span_only_D = residualize_rank(D, [span])
    span_only_Dunb = residualize_rank(Dunb, [span])
    n_only_D = residualize_rank(D, [n])
    n_only_Dunb = residualize_rank(Dunb, [n])

    return {
        "raw_D_unbiased": cluster(Dunb, genus, seed=SEED + seed_offset + 1),
        "span_only_adjusted_D": cluster(span_only_D, genus, seed=SEED + seed_offset + 2),
        "span_only_adjusted_D_unbiased": cluster(span_only_Dunb, genus, seed=SEED + seed_offset + 3),
        "n_classifiable_only_adjusted_D": cluster(n_only_D, genus, seed=SEED + seed_offset + 4),
        "n_classifiable_only_adjusted_D_unbiased": cluster(n_only_Dunb, genus, seed=SEED + seed_offset + 5),
        "rho_D_n_classifiable": float(stats.spearmanr(D, n).statistic),
        "rho_D_unbiased_n_classifiable": float(stats.spearmanr(Dunb, n).statistic),
        "rho_D_span": float(stats.spearmanr(D, span).statistic),
        "rho_D_unbiased_span": float(stats.spearmanr(Dunb, span).statistic),
    }


def main() -> None:
    discovery = add_genus(pd.read_csv(DISC)).sort_values("species").reset_index(drop=True)
    reserve = add_genus(pd.read_csv(RES)).sort_values("species").reset_index(drop=True)
    step7 = json.loads(STEP7.read_text(encoding="utf-8"))
    step2 = json.loads(STEP2.read_text(encoding="utf-8"))
    shared = pd.read_csv(STEP7_SHARED)

    if len(discovery) != EXPECTED_DISCOVERY or discovery["species"].nunique() != EXPECTED_DISCOVERY:
        raise RuntimeError("discovery count drift")
    if len(reserve) != EXPECTED_RESERVE or reserve["species"].nunique() != EXPECTED_RESERVE:
        raise RuntimeError("reserve count drift")
    if int(step2["species_level"]["raw_rows_per_species_min"]) != 100 or int(step2["species_level"]["raw_rows_per_species_max"]) != 100:
        raise RuntimeError("discovery fixed 100-photo denominator not preserved")
    if not reserve["reserve_photo_rows"].eq(100).all():
        raise RuntimeError("reserve fixed 100-photo denominator not preserved")
    if len(shared) != 23:
        raise RuntimeError(f"Step7 shared-genus list drift: {len(shared)} != 23")

    disc_cov = repeated_group_pairs(discovery["genus"])[4]
    res_cov = repeated_group_pairs(reserve["genus"])[4]
    if disc_cov["repeated_groups"] != step7["discovery"]["raw_step4_seed_exact_reproduction"]["repeated_groups"]:
        raise RuntimeError("discovery repeated-genus coverage drift")
    if res_cov["repeated_groups"] != step7["reserve"]["raw"]["repeated_groups"]:
        raise RuntimeError("reserve repeated-genus coverage drift")

    discovery_tests = tranche_tests(discovery, span_col="log1p_span_primary", seed_offset=100)
    reserve_tests = tranche_tests(reserve, span_col="reserve_log1p_span", seed_offset=200)

    shared_names = shared["genus"].astype(str).tolist()
    d_unb = discovery.assign(D_unbiased=discovery["D"] * discovery["n_classifiable"] / (discovery["n_classifiable"] - 1.0))
    r_unb = reserve.assign(D_unbiased=reserve["D"] * reserve["n_classifiable"] / (reserve["n_classifiable"] - 1.0))
    rows = []
    for genus in shared_names:
        dg = d_unb[d_unb["genus"].eq(genus)]
        rg = r_unb[r_unb["genus"].eq(genus)]
        if len(dg) < 2 or len(rg) < 2:
            raise RuntimeError(f"shared genus lost repeated coverage: {genus}")
        rows.append({
            "genus": genus,
            "discovery_n": int(len(dg)),
            "reserve_n": int(len(rg)),
            "discovery_mean_D_unbiased": float(dg["D_unbiased"].mean()),
            "reserve_mean_D_unbiased": float(rg["D_unbiased"].mean()),
        })
    shared_unb = pd.DataFrame(rows)
    shared_unb_test = spearman_perm(
        shared_unb["discovery_mean_D_unbiased"].to_numpy(float),
        shared_unb["reserve_mean_D_unbiased"].to_numpy(float),
        seed=SEED + 300,
    )

    res_span_pass = bool(
        reserve_tests["span_only_adjusted_D"]["supported_at_0_05"]
        or reserve_tests["span_only_adjusted_D_unbiased"]["supported_at_0_05"]
    )
    res_n_pass = bool(
        reserve_tests["n_classifiable_only_adjusted_D"]["supported_at_0_05"]
        or reserve_tests["n_classifiable_only_adjusted_D_unbiased"]["supported_at_0_05"]
    )
    raw_unb_rep = bool(reserve_tests["raw_D_unbiased"]["supported_at_0_05"])

    if res_span_pass and not res_n_pass:
        classification = "SPAN_NOT_SUFFICIENT_BUT_SIGNAL_ENTANGLED_WITH_CLASSIFICATION_YIELD"
    elif res_span_pass and res_n_pass:
        classification = "GENUS_CLUSTERING_SURVIVES_SPAN_AND_N_SEPARATELY"
    elif not res_span_pass:
        classification = "SPAN_ADJUSTMENT_WEAKENS_GENUS_PROPENSITY_CLAIM"
    else:
        classification = "UNRESOLVED"

    result = {
        "analysis": "polymorphism_genus_adjustment_decomposition_step7b",
        "date_jst": "2026-09-10",
        "protocol": "docs/POLYMORPHISM_GENUS_ADJUSTMENT_DECOMPOSITION_STEP7B_PROTOCOL_20260910.md",
        "fixed_photo_denominator_verified": {"discovery": 100, "reserve": 100},
        "n_permutations_each_test": N_PERM,
        "discovery": discovery_tests,
        "reserve": reserve_tests,
        "step7_simultaneous_span_plus_n": {
            "discovery": step7["discovery"]["adjusted_span_and_n"],
            "reserve": step7["reserve"]["adjusted_span_and_n"],
        },
        "shared_23_genus_D_unbiased_concordance": shared_unb_test,
        "decision": {
            "reserve_raw_D_unbiased_genus_clustering_replicated": raw_unb_rep,
            "sampled_span_alone_insufficient_to_explain_reserve_genus_clustering": res_span_pass,
            "genus_clustering_survives_n_classifiable_only_adjustment_in_reserve": res_n_pass,
            "classification": classification,
            "claim_boundary": "Fixed 100-photo acquisition removes unequal raw effort, and D_unbiased addresses ordinary finite-sample bias, but non-random classification missingness remains unresolved; n_classifiable is a post-acquisition measurement-process variable and is not treated as a clean exogenous confounder.",
        },
    }

    shared_unb.to_csv(OUT / "shared_23_genus_D_unbiased_means.csv", index=False)
    write_json(OUT / "result.json", result)

    def line(label: str, x: dict[str, Any]) -> str:
        return f"- {label}: gain **{x['clustering_gain']:.6f}**, p_lower **{x['p_lower']:.6g}**"

    lines = [
        "# Polymorphism genus clustering — Step 7b adjustment decomposition result",
        "",
        f"**Decision: `{classification}`.**",
        f"**Reserve raw D_unbiased replication: `{raw_unb_rep}`.**",
        f"**Sampled span alone insufficient: `{res_span_pass}`.**",
        f"**Reserve genus clustering survives n_classifiable-only adjustment: `{res_n_pass}`.**",
        "",
        "## Discovery",
        "",
        line("raw D_unbiased", discovery_tests["raw_D_unbiased"]),
        line("span-only D", discovery_tests["span_only_adjusted_D"]),
        line("span-only D_unbiased", discovery_tests["span_only_adjusted_D_unbiased"]),
        line("n_classifiable-only D", discovery_tests["n_classifiable_only_adjusted_D"]),
        line("n_classifiable-only D_unbiased", discovery_tests["n_classifiable_only_adjusted_D_unbiased"]),
        "",
        "## Species-disjoint reserve",
        "",
        line("raw D_unbiased", reserve_tests["raw_D_unbiased"]),
        line("span-only D", reserve_tests["span_only_adjusted_D"]),
        line("span-only D_unbiased", reserve_tests["span_only_adjusted_D_unbiased"]),
        line("n_classifiable-only D", reserve_tests["n_classifiable_only_adjusted_D"]),
        line("n_classifiable-only D_unbiased", reserve_tests["n_classifiable_only_adjusted_D_unbiased"]),
        "",
        "## Fixed-denominator diagnostic",
        "",
        "- raw photo denominator per species: **100 in discovery and 100 in reserve**",
        f"- discovery rho(D, n_classifiable) = **{discovery_tests['rho_D_n_classifiable']:.6f}**",
        f"- reserve rho(D, n_classifiable) = **{reserve_tests['rho_D_n_classifiable']:.6f}**",
        f"- discovery rho(D, sampled span) = **{discovery_tests['rho_D_span']:.6f}**",
        f"- reserve rho(D, sampled span) = **{reserve_tests['rho_D_span']:.6f}**",
        "",
        "## Shared 23 repeated genera — D_unbiased",
        "",
        f"- Spearman rho = **{shared_unb_test['rho']:.6f}**, p_two_sided = **{shared_unb_test['p_two_sided']:.6g}**",
        "",
        "## Interpretation boundary",
        "",
        "A span-only surviving result rules out sampled geographic span as a sufficient explanation of the reserve genus pattern. Failure after n_classifiable adjustment is treated separately because n_classifiable is a post-acquisition classification-success variable under a fixed 100-photo denominator and is strongly entangled with mixed/uncertain colour measurement. Non-random classification missingness therefore remains a measurement limitation rather than a solved confound.",
    ]
    (OUT / "RESULT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(json.dumps(result["decision"], indent=2))
    print((OUT / "RESULT.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()

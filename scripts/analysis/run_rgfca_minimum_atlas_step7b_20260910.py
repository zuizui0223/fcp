#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DISC = ROOT / "data/derived/global_monte_carlo_measured_photos_v1.csv"
RES = ROOT / "data/derived/rgfca_reserve_replication_measured_photos_v1.csv"
OUT = ROOT / "results/rgfca_minimum_atlas_step7b_20260910"
OUT.mkdir(parents=True, exist_ok=True)
MORPHS = ["white", "yellow_orange", "red_pink", "blue_purple"]
MORPH_CODE = {m: i for i, m in enumerate(MORPHS)}
S_LEVELS = [25, 50, 100, 150, 250, 350]
N_LEVELS = [10, 20, 40, 60, 80, 100]
R = 200
SEED = 20260910
Q = np.arange(0.1, 1.0, 0.1)


def bseries(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False).astype(bool)
    return s.fillna("").astype(str).str.strip().str.lower().isin(["true", "1", "yes"])


def cell_id(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    la = np.clip(np.floor((lat + 90.0) / 10.0).astype(int), 0, 17)
    lo = np.clip(np.floor((lon + 180.0) / 10.0).astype(int), 0, 35)
    return la * 36 + lo


def load_pool(path: Path) -> dict[str, dict[str, np.ndarray]]:
    use = ["species", "morph", "global_classifiable", "latitude", "longitude"]
    df = pd.read_csv(path, usecols=use)
    df["_ok"] = bseries(df["global_classifiable"]) & df["morph"].isin(MORPHS) & df["latitude"].notna() & df["longitude"].notna()
    pools: dict[str, dict[str, np.ndarray]] = {}
    for sp, g in df.groupby("species", sort=True):
        if len(g) < 100:
            continue
        g = g.iloc[:100].copy()
        ok = g["_ok"].to_numpy(bool)
        mc = np.full(100, -1, dtype=np.int8)
        if ok.any():
            mc[ok] = np.array([MORPH_CODE[x] for x in g.loc[ok, "morph"].astype(str)], dtype=np.int8)
        cells = np.full(100, -1, dtype=np.int16)
        geo = g["latitude"].notna().to_numpy() & g["longitude"].notna().to_numpy()
        if geo.any():
            cells[geo] = cell_id(g.loc[geo, "latitude"].to_numpy(float), g.loc[geo, "longitude"].to_numpy(float)).astype(np.int16)
        pools[str(sp)] = {"morph": mc, "cell": cells}
    return pools


def d_value(morph: np.ndarray) -> float | None:
    x = morph[morph >= 0]
    if len(x) < 5:
        return None
    c = np.bincount(x, minlength=4).astype(float)
    p = c / c.sum()
    return float(1.0 - np.sum(p * p))


def second_fraction(morph: np.ndarray) -> float | None:
    x = morph[morph >= 0]
    if len(x) < 5:
        return None
    c = np.bincount(x, minlength=4).astype(float)
    p = c / c.sum()
    return float(np.sort(p)[-2])


def species_cell_props(morph: np.ndarray, cells: np.ndarray) -> dict[int, np.ndarray]:
    keep = (morph >= 0) & (cells >= 0)
    if not keep.any():
        return {}
    out: dict[int, np.ndarray] = {}
    for c in np.unique(cells[keep]):
        z = morph[keep & (cells == c)]
        cnt = np.bincount(z, minlength=4).astype(float)
        out[int(c)] = cnt / cnt.sum()
    return out


def reference(pools: dict[str, dict[str, np.ndarray]]) -> dict[str, Any]:
    cell_sum: dict[int, np.ndarray] = {}
    cell_n: dict[int, int] = {}
    ds: list[float] = []
    seconds: list[float] = []
    class_n: list[int] = []
    for sp in sorted(pools):
        x = pools[sp]["morph"]
        d = d_value(x)
        sec = second_fraction(x)
        if d is not None:
            ds.append(d); seconds.append(float(sec)); class_n.append(int((x >= 0).sum()))
        for c, p in species_cell_props(x, pools[sp]["cell"]).items():
            cell_sum[c] = cell_sum.get(c, np.zeros(4, float)) + p
            cell_n[c] = cell_n.get(c, 0) + 1
    cells = sorted(c for c, n in cell_n.items() if n >= 3)
    mat = np.stack([cell_sum[c] / cell_n[c] for c in cells])
    return {
        "cells": cells,
        "cell_index": {c: i for i, c in enumerate(cells)},
        "atlas": mat,
        "D": np.asarray(ds, float),
        "D_quantiles": np.quantile(np.asarray(ds, float), Q),
        "second_ge_0_10": float(np.mean(np.asarray(seconds) >= 0.10)),
        "n_species_D": int(len(ds)),
        "median_classifiable": float(np.median(class_n)),
    }


def metrics_from_acc(cell_sum: np.ndarray, cell_n: np.ndarray, ds: list[float], secs: list[float], ref: dict[str, Any], selected: int, photo_n: int, class_counts: list[int]) -> dict[str, Any]:
    present = cell_n > 0
    coverage = float(present.mean())
    if present.any():
        est = cell_sum[present] / cell_n[present, None]
        target = ref["atlas"][present]
        a = est.ravel(); b = target.ravel()
        denom = float(np.linalg.norm(a) * np.linalg.norm(b))
        cosine = float(a @ b / denom) if denom > 0 else np.nan
    else:
        cosine = np.nan
    if len(ds) >= 5:
        dq = np.quantile(np.asarray(ds, float), Q)
        d_mae = float(np.mean(np.abs(dq - ref["D_quantiles"])))
        mean_d = float(np.mean(ds)); med_d = float(np.median(ds))
        sec10 = float(np.mean(np.asarray(secs) >= 0.10))
    else:
        d_mae = np.nan; mean_d = np.nan; med_d = np.nan; sec10 = np.nan
    pass_one = bool(np.isfinite(cosine) and np.isfinite(d_mae) and coverage >= 0.80 and cosine >= 0.90 and d_mae <= 0.05)
    return {
        "species_count": selected,
        "raw_photos_per_species": photo_n,
        "nominal_budget": selected * photo_n,
        "atlas_cell_coverage": coverage,
        "atlas_composition_cosine": cosine,
        "D_quantile_MAE": d_mae,
        "mean_D": mean_d,
        "median_D": med_d,
        "second_ge_0_10_fraction": sec10,
        "median_classifiable_photos": float(np.median(class_counts)) if class_counts else np.nan,
        "fraction_species_ge5_classifiable": float(np.mean(np.asarray(class_counts) >= 5)) if class_counts else np.nan,
        "pass": pass_one,
    }


def run_tranche(pools: dict[str, dict[str, np.ndarray]], tranche: str, seed: int) -> tuple[pd.DataFrame, dict[str, Any]]:
    species = np.array(sorted(pools), dtype=object)
    if len(species) < 350:
        raise RuntimeError(f"{tranche}: only {len(species)} species with >=100 raw rows")
    ref = reference(pools)
    ref_cells = ref["cells"]; cindex = ref["cell_index"]
    rng = np.random.default_rng(seed)
    rows: list[dict[str, Any]] = []
    targets = set(S_LEVELS)
    for rep in range(R):
        sp_order = species[rng.permutation(len(species))]
        photo_orders = {sp: rng.permutation(100) for sp in sp_order[:max(S_LEVELS)]}
        for photo_n in N_LEVELS:
            sums = np.zeros((len(ref_cells), 4), float)
            nums = np.zeros(len(ref_cells), int)
            ds: list[float] = []; secs: list[float] = []; class_counts: list[int] = []
            for k, sp in enumerate(sp_order[:max(S_LEVELS)], start=1):
                idx = photo_orders[sp][:photo_n]
                morph = pools[sp]["morph"][idx]
                cells = pools[sp]["cell"][idx]
                class_counts.append(int((morph >= 0).sum()))
                d = d_value(morph); sec = second_fraction(morph)
                if d is not None:
                    ds.append(d); secs.append(float(sec))
                for c, p in species_cell_props(morph, cells).items():
                    ci = cindex.get(c)
                    if ci is not None:
                        sums[ci] += p; nums[ci] += 1
                if k in targets:
                    m = metrics_from_acc(sums, nums, ds, secs, ref, k, photo_n, class_counts)
                    m.update({"tranche": tranche, "realization": rep})
                    rows.append(m)
    rdf = pd.DataFrame(rows)
    summary_rows = []
    for (s, n), g in rdf.groupby(["species_count", "raw_photos_per_species"], sort=True):
        summary_rows.append({
            "tranche": tranche,
            "species_count": int(s),
            "raw_photos_per_species": int(n),
            "nominal_budget": int(s*n),
            "pass_rate": float(g["pass"].mean()),
            "median_coverage": float(g["atlas_cell_coverage"].median()),
            "q10_coverage": float(g["atlas_cell_coverage"].quantile(0.10)),
            "median_cosine": float(g["atlas_composition_cosine"].median()),
            "q10_cosine": float(g["atlas_composition_cosine"].quantile(0.10)),
            "median_D_quantile_MAE": float(g["D_quantile_MAE"].median()),
            "q90_D_quantile_MAE": float(g["D_quantile_MAE"].quantile(0.90)),
            "design_pass": bool(g["pass"].mean() >= 0.90),
        })
    sdf = pd.DataFrame(summary_rows)
    meta = {
        "species_universe": int(len(species)),
        "reference_cells_ge3_species": int(len(ref_cells)),
        "reference_D_species": int(ref["n_species_D"]),
        "reference_second_ge_0_10": float(ref["second_ge_0_10"]),
        "reference_median_classifiable": float(ref["median_classifiable"]),
        "passing_designs": sdf.loc[sdf["design_pass"], ["species_count", "raw_photos_per_species", "nominal_budget"]].to_dict("records"),
    }
    return rdf, {"summary": sdf, "meta": meta}


def main() -> None:
    dp = load_pool(DISC); rp = load_pool(RES)
    dr, d = run_tranche(dp, "discovery", SEED + 1)
    rr, r = run_tranche(rp, "reserve", SEED + 2)
    reps = pd.concat([dr, rr], ignore_index=True)
    sums = pd.concat([d["summary"], r["summary"]], ignore_index=True)
    reps.to_csv(OUT / "realization_metrics.csv.gz", index=False, compression="gzip")
    sums.to_csv(OUT / "design_summary.csv", index=False)
    common = sums.pivot_table(index=["species_count", "raw_photos_per_species", "nominal_budget"], columns="tranche", values="design_pass", aggfunc="first").reset_index()
    common["joint_pass"] = common.get("discovery", False).fillna(False).astype(bool) & common.get("reserve", False).fillna(False).astype(bool)
    passing = common.loc[common["joint_pass"]].sort_values(["nominal_budget", "species_count", "raw_photos_per_species"])
    if len(passing):
        z = passing.iloc[0]
        verdict = "COMMON_MINIMUM_FOUND"
        minimum = {"species_count": int(z.species_count), "raw_photos_per_species": int(z.raw_photos_per_species), "nominal_budget": int(z.nominal_budget)}
    else:
        verdict = "NO_COMMON_MINIMUM_WITHIN_GRID"
        minimum = None
    result = {
        "analysis": "rgfca_minimum_atlas_step7b",
        "protocol": "docs/RGFCA_MINIMUM_ATLAS_STEP7B_PROTOCOL_20260910.md",
        "realizations_per_design": R,
        "discovery": d["meta"],
        "reserve": r["meta"],
        "verdict": verdict,
        "minimum_common_design": minimum,
        "claim_boundary": "Minimum refers to recovery of the frozen species-equal 10-degree atlas plus D distribution under the predeclared thresholds; it is not a universal sampling law for all flower-colour questions.",
    }
    (OUT / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    md = ["# RGFCA Step 7B — minimum global atlas", "", f"**Verdict: `{verdict}`**", ""]
    if minimum:
        md += [f"- common minimum: **{minimum['species_count']} species x {minimum['raw_photos_per_species']} raw photos/species = {minimum['nominal_budget']} photos**", ""]
    for tranche, meta in [("discovery", d["meta"]), ("reserve", r["meta"])]:
        md += [f"## {tranche}", "", f"- species universe: **{meta['species_universe']}**", f"- reference atlas cells (>=3 species): **{meta['reference_cells_ge3_species']}**", f"- reference D species: **{meta['reference_D_species']}**", f"- passing designs: **{len(meta['passing_designs'])}**", ""]
    md += ["Primary pass requires >=90% of 200 realizations to achieve >=0.80 cell coverage, >=0.90 composition cosine, and <=0.05 D-quantile MAE in each tranche independently."]
    (OUT / "RESULT.md").write_text("\n".join(md) + "\n", encoding="utf-8")

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Audit FCP spatial holdout failures from saved species-level predictions.

Pure Python standard library; never refits, changes, or selects model coefficients.
Run from the repository root:
python scripts/analysis/diagnose_fcp_spatial_holdout_failure_20261019.py \
 --input results/fcp_environment_transport_20261019/spatial_species_metrics.csv \
 --outdir results/fcp_spatial_failure_diagnosis_20261019/python_replay
"""
import argparse
import csv
import hashlib
import json
import random
import statistics
from collections import defaultdict
from pathlib import Path

NUMERIC = (
    "buffer_km", "n_photos", "n_white", "logloss_gain", "brier_gain",
    "augmented_outside_train_feature_range_photos",
    "augmented_abs_train_z_gt_3_photos", "median_test_to_train_km"
)
FIELDS = ("inat_taxon_id", "species", "genus", "source_cohort", "comparison") + NUMERIC
GROUP_VARS = ("outside_frac", "extreme_frac", "n_photos", "balance")
SEED = 2026101920

def load(path):
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            missing = set(FIELDS) - set(r)
            if missing:
                raise ValueError(f"missing columns: {sorted(missing)}")
            for k in NUMERIC:
                r[k] = float(r[k])
            if r["n_photos"] <= 0 or not 0 <= r["n_white"] <= r["n_photos"]:
                raise ValueError("invalid class/photo denominator")
            r["outside_frac"] = r["augmented_outside_train_feature_range_photos"] / r["n_photos"]
            r["extreme_frac"] = r["augmented_abs_train_z_gt_3_photos"] / r["n_photos"]
            r["balance"] = 2 * min(r["n_white"], r["n_photos"]-r["n_white"])/r["n_photos"]
            rows.append(r)
    if len({(r["inat_taxon_id"], r["buffer_km"], r["comparison"]) for r in rows}) != len(rows):
        raise ValueError("duplicate species x buffer x comparison")
    return rows

def summary(rows):
    vals = sorted(r["logloss_gain"] for r in rows)
    if not vals:
        return {"n": 0}
    k10, k20 = int(len(vals)*.1), int(len(vals)*.2)
    return {
        "n": len(vals), "genera": len({r["genus"] for r in rows}),
        "n_photos": int(sum(r["n_photos"] for r in rows)),
        "mean_gain": statistics.fmean(vals),
        "median_gain": statistics.median(vals),
        "fraction_improved": sum(v > 0 for v in vals) / len(vals),
        "mean_excluding_worst_10pct": statistics.fmean(vals[k10:]),
        "mean_excluding_worst_20pct": statistics.fmean(vals[k20:]),
        "mean_brier_gain": statistics.fmean(r["brier_gain"] for r in rows)
    }

def percentile(values, proportion):
    z = sorted(values)
    return z[int((len(z)-1)*proportion)]

def quartile_group(rows, key):
    cuts = [percentile([r[key] for r in rows], p) for p in (.25, .5, .75)]
    groups = [[] for _ in range(4)]
    for r in rows:
        value = r[key]
        q = 0 if value <= cuts[0] else 1 if value <= cuts[1] else 2 if value <= cuts[2] else 3
        groups[q].append(r)
    return cuts, groups

def genus_bootstrap_q4_minus_q1(rows, qmap, repetitions, seed):
    by_genus = defaultdict(list)
    for r in rows:
        by_genus[r["genus"]].append(r)
    keys = sorted(by_genus)
    rng = random.Random(seed)
    draws = []
    for _ in range(repetitions):
        sums = [0.0, 0.0]
        n = [0, 0]
        for g in rng.choices(keys, k=len(keys)):
            for r in by_genus[g]:
                q = qmap[r["inat_taxon_id"]]
                if q in (0, 3):
                    idx = 0 if q == 0 else 1
                    sums[idx] += r["logloss_gain"]
                    n[idx] += 1
        if n[0] > 0 and n[1] > 0:
            draws.append(sums[1] / n[1] - sums[0] / n[0])
    return {
        "n_draws": len(draws),
        "difference_q4_minus_q1": None if not draws else statistics.fmean(draws),
        "genus_bootstrap_95pct": [percentile(draws, .025), percentile(draws, .975)] if draws else None,
        "method": "Fixed empirical quartile cuts; resample entire genus clusters. Exploratory, not cross-validation refitting."
    }

def write_csv(path, records):
    if not records:
        return
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--bootstrap", type=int, default=2000)
    a = p.parse_args()
    if a.bootstrap < 1:
        raise ValueError("bootstrap must be >= 1")
    all_rows = load(a.input)
    a.outdir.mkdir(parents=True, exist_ok=True)
    comparisons, quartiles, bootstrap = [], [], []
    grouped = defaultdict(list)
    for r in all_rows:
        grouped[(int(r["buffer_km"]), r["comparison"])].append(r)
    for (buffer_km, comparison), rows in sorted(grouped.items()):
        main_stats = {"buffer_km": buffer_km, "comparison": comparison, **summary(rows)}
        comparisons.append(main_stats)
        for var in GROUP_VARS:
            cuts, groups = quartile_group(rows, var)
            qmap = {r["inat_taxon_id"]: q for q, arr in enumerate(groups) for r in arr}
            for q, group in enumerate(groups):
                s = summary(group)
                quartiles.append({
                    "buffer_km": buffer_km, "comparison": comparison, "quartile_var": var,
                    "quartile": q + 1, "cut1": cuts[0], "cut2": cuts[1], "cut3": cuts[2],
                    **s
                })
            if var == "outside_frac" and all(groups[i] for i in (0,3)):
                bootstrap.append({
                    "buffer_km": buffer_km, "comparison": comparison,
                    **genus_bootstrap_q4_minus_q1(rows, qmap, a.bootstrap, SEED + buffer_km)
                })
    write_csv(a.outdir / "comparison_summary.csv", comparisons)
    write_csv(a.outdir / "quartile_summary.csv", quartiles)
    source_sha = hashlib.sha256(a.input.read_bytes()).hexdigest()
    (a.outdir / "receipt.json").write_text(json.dumps({
        "input_file": str(a.input), "input_sha256": source_sha,
        "input_rows": len(all_rows), "bootstrap_reps": a.bootstrap, "seed": SEED,
        "bootstrap_outside_q4_minus_q1": bootstrap,
        "classification": "posthoc internal diagnostic, no causal inference or model refit"
    }, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"rows": len(all_rows), "comparisons": comparisons, "bootstrap": bootstrap}, indent=2))

if __name__ == "__main__":
    main()

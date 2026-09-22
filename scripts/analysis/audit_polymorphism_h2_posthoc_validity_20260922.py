#!/usr/bin/env python3
"""Post-confirmatory validity diagnostics for the frozen third-cohort H2 result.

This script does not modify or rerun the frozen prospective decision. It:
1) decomposes the observed white-axis alignment relative to the frozen structured null;
2) stratifies species-level alignment by whether white is a primary/secondary coarse morph;
3) reruns a post hoc 299-world structured null starting from the coarse-gate species
   and reapplies the continuous-minor gate within every null world.

The measured CSV should be obtained from the immutable third-cohort artifact
(ID 10496492307, digest recorded in the audit result).
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np
import pandas as pd

BIO = ["white","yellow","orange","red","pink","magenta","purple","blue","bronze"]
FRACTIONS = [f"flower_fraction_{c}" for c in BIO]
MORPHS = ["white","yellow_orange","red_pink","blue_purple"]
EPS = 1e-12
MAX_KMEANS_ITER = 200


def _bool_series(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False).astype(bool)
    return s.fillna("").astype(str).str.strip().str.lower().isin(["true","1","yes"])


def _load_classifiable(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, usecols=["species","morph","global_classifiable"] + FRACTIONS)
    df["species"] = df["species"].astype(str)
    df["morph"] = df["morph"].fillna("").astype(str)
    x = df.loc[_bool_series(df["global_classifiable"]) & df["morph"].isin(MORPHS)].copy()
    p = x[FRACTIONS].to_numpy(float)
    mass = p.sum(axis=1)
    if np.any(~np.isfinite(p)) or np.any(mass <= 0):
        raise RuntimeError("invalid classifiable palette row")
    x.loc[:, FRACTIONS] = p / mass[:, None]
    return x


def _two_means(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    grand = x.mean(axis=0)
    i0 = int(np.argmax(np.sum((x-grand)**2, axis=1)))
    d0 = np.sum((x-x[i0])**2, axis=1)
    i1 = int(np.argmax(d0))
    if float(d0[i1]) <= EPS:
        labels = np.zeros(len(x), dtype=int)
        labels[len(x)//2:] = 1
    else:
        centers = np.vstack([x[i0],x[i1]])
        labels = np.full(len(x), -1, dtype=int)
        for _ in range(MAX_KMEANS_ITER):
            dist = np.sum((x[:,None,:]-centers[None,:,:])**2, axis=2)
            new = np.argmin(dist, axis=1).astype(int)
            if np.all(new == new[0]):
                only = int(new[0])
                new[int(np.argmax(dist[:,only]))] = 1-only
            if np.array_equal(new, labels):
                break
            labels = new
            for k in (0,1):
                centers[k] = x[labels == k].mean(axis=0)
        else:
            raise RuntimeError("two-means did not converge")
    return labels, np.bincount(labels, minlength=2)


def _coarse_second(g: pd.DataFrame) -> float:
    cnt = Counter(g["morph"].astype(str))
    order = sorted(MORPHS, key=lambda m: (-cnt.get(m,0), m))
    return float(cnt.get(order[1],0) / len(g))


def _prepare(work: pd.DataFrame, species: list[str]) -> dict:
    sub = work.loc[work["species"].isin(set(species))].reset_index(drop=True)
    sp = sub["species"].to_numpy(str)
    morph = sub["morph"].to_numpy(str)
    return {
        "p": sub[FRACTIONS].to_numpy(float),
        "species_idx": {s: np.flatnonzero(sp == s) for s in sorted(species)},
        "morph_idx": {m: np.flatnonzero(morph == m) for m in MORPHS},
    }


def _permute(prep: dict, rng: np.random.Generator) -> np.ndarray:
    src = prep["p"]
    out = src.copy()
    for m in MORPHS:
        idx = prep["morph_idx"][m]
        if len(idx) > 1:
            out[idx] = src[idx[rng.permutation(len(idx))]]
    return out


def _regated_units(p: np.ndarray, species_idx: dict[str,np.ndarray], threshold: float) -> np.ndarray:
    units = []
    for sp in sorted(species_idx):
        rows = p[species_idx[sp]]
        labels, counts = _two_means(np.sqrt(rows))
        if float(counts.min()/len(rows)) < threshold:
            continue
        if counts[0] > counts[1]:
            major, minor = 0, 1
        elif counts[1] > counts[0]:
            major, minor = 1, 0
        else:
            c0, c1 = rows[labels==0].mean(axis=0), rows[labels==1].mean(axis=0)
            major, minor = (0,1) if tuple(c0.tolist()) <= tuple(c1.tolist()) else (1,0)
        delta = rows[labels==minor].mean(axis=0) - rows[labels==major].mean(axis=0)
        norm = float(np.linalg.norm(delta))
        if np.isfinite(norm) and norm > EPS:
            units.append(delta/norm)
    return np.vstack(units)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--measured-csv", type=Path, required=True)
    ap.add_argument("--delta-csv", type=Path, required=True)
    ap.add_argument("--h2-result", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--null-reps", type=int, default=299)
    ap.add_argument("--seed", type=int, default=20260915)
    args = ap.parse_args()

    h2 = json.loads(args.h2_result.read_text())
    primary = h2["thresholds"]["primary_0_10"]
    observed = float(primary["observed_W"])
    q = np.array([1.0] + [-1.0/8.0]*8)
    q /= np.linalg.norm(q)

    d = pd.read_csv(args.delta_csv)
    delta = d[[f"delta_{c}" for c in BIO]].to_numpy(float)
    units = delta / np.linalg.norm(delta, axis=1)[:,None]
    wi = np.square(units @ q)
    white = d["H1_coarse_primary_morph"].eq("white") | d["H1_coarse_secondary_morph"].eq("white")

    work = _load_classifiable(args.measured_csv)
    coarse = [
        sp for sp,g in work.groupby("species", sort=True)
        if len(g) >= 40 and _coarse_second(g) >= 0.10
    ]
    prep = _prepare(work, coarse)
    rng = np.random.default_rng(args.seed)
    null, retained = [], []
    for _ in range(args.null_reps):
        u = _regated_units(_permute(prep,rng), prep["species_idx"], 0.10)
        retained.append(len(u))
        null.append(float(np.mean(np.square(u @ q))))
    null = np.asarray(null)
    retained = np.asarray(retained)

    out = {
        "role":"post_confirmatory_validity_diagnostic",
        "changes_frozen_h2_verdict":False,
        "frozen_h2":{
            "species":int(len(d)),
            "observed_W":observed,
            "isotropic_8d_expectation":0.125,
            "structured_null_median":float(primary["structured_null_summary"]["q50"]),
            "structured_null_upper_p":float(primary["structured_null_upper_p"]),
            "observed_to_null_median_ratio":float(primary["observed_to_null_median_ratio"]),
        },
        "white_involvement":{
            "white_in_primary_or_secondary_coarse_morph_n":int(white.sum()),
            "fraction":float(white.mean()),
            "mean_species_alignment_white_involved":float(wi[white].mean()),
            "mean_species_alignment_nonwhite_involved":float(wi[~white].mean()),
        },
        "gate_reapplied_null":{
            "starting_coarse_gate_species":int(len(coarse)),
            "null_replicates":int(args.null_reps),
            "seed":int(args.seed),
            "median_W":float(np.median(null)),
            "mean_W":float(np.mean(null)),
            "q025":float(np.quantile(null,.025)),
            "q975":float(np.quantile(null,.975)),
            "retained_species_median":float(np.median(retained)),
            "retained_species_min":int(retained.min()),
            "retained_species_max":int(retained.max()),
            "upper_p":float((1+np.sum(null>=observed))/(len(null)+1)),
        },
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir/"result_recomputed.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
    pd.DataFrame({"W":null,"retained_species":retained}).to_csv(args.output_dir/"gate_reapplied_null_299.csv",index=False)
    print(json.dumps(out,indent=2,sort_keys=True))


if __name__ == "__main__":
    main()

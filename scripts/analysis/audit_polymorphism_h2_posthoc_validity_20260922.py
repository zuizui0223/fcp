#!/usr/bin/env python3
"""Post-confirmatory validity audit for the third-cohort white-axis H2 result.

This script does not change the frozen prospective H2 verdict. It quantifies:
1) how much observed W is carried by species whose two most frequent coarse
   morphs include white;
2) a post hoc structured null that reapplies the continuous minor-cluster gate
   in every null world; and
3) numerical divergence between the frozen biological H2 implementation and
   the later generic disttrait implementation on the same real-data rows.

The third-cohort image pixels/background palette were not persisted, so direct
clip/exposure diagnostics cannot be reconstructed from the sealed artifact.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "packages" / "disttrait" / "src"))

from disttrait import one_vs_rest_contrast, structured_alignment_null, two_mode_axis  # noqa: E402

BIO = ["white", "yellow", "orange", "red", "pink", "magenta", "purple", "blue", "bronze"]
FRACTIONS = [f"flower_fraction_{c}" for c in BIO]
MORPHS = ["white", "yellow_orange", "red_pink", "blue_purple"]
MIN_CLASSIFIABLE = 40
THRESHOLD = 0.10
N_REAPPLIED_NULL = 299
SEED = 20260915
EPS = 1e-12


def _bool_series(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False).astype(bool)
    return s.fillna("").astype(str).str.strip().str.lower().isin(["true", "1", "yes"])


def load_classifiable(path: Path) -> pd.DataFrame:
    usecols = ["species", "morph", "global_classifiable"] + FRACTIONS
    df = pd.read_csv(path, usecols=usecols)
    keep = _bool_series(df["global_classifiable"]) & df["morph"].isin(MORPHS)
    work = df.loc[keep].copy()
    x = work[FRACTIONS].to_numpy(float)
    mass = x.sum(axis=1)
    if np.any(~np.isfinite(x)) or np.any(mass <= 0):
        raise RuntimeError("invalid classifiable palette rows")
    work.loc[:, FRACTIONS] = x / mass[:, None]
    work["species"] = work["species"].astype(str)
    return work


def coarse_gate(g: pd.DataFrame) -> tuple[float, str, str, dict[str, int]]:
    counts = Counter(g["morph"].astype(str))
    order = sorted(MORPHS, key=lambda m: (-counts.get(m, 0), m))
    primary, secondary = order[:2]
    return (
        counts.get(secondary, 0) / len(g),
        primary,
        secondary,
        {m: int(counts.get(m, 0)) for m in MORPHS},
    )


def deterministic_two_means(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray(x, dtype=float)
    grand = x.mean(axis=0)
    i0 = int(np.argmax(np.sum((x - grand) ** 2, axis=1)))
    d0 = np.sum((x - x[i0]) ** 2, axis=1)
    i1 = int(np.argmax(d0))
    if float(d0[i1]) <= EPS:
        labels = np.zeros(len(x), dtype=int)
        labels[len(x) // 2 :] = 1
    else:
        centers = np.vstack([x[i0], x[i1]])
        labels = np.full(len(x), -1, dtype=int)
        for _ in range(200):
            dist = np.sum((x[:, None, :] - centers[None, :, :]) ** 2, axis=2)
            new = np.argmin(dist, axis=1).astype(int)
            if np.all(new == new[0]):
                only = int(new[0])
                new[int(np.argmax(dist[:, only]))] = 1 - only
            if np.array_equal(new, labels):
                break
            labels = new
            for k in (0, 1):
                centers[k] = x[labels == k].mean(axis=0)
        else:
            raise RuntimeError("two-means did not converge")
    return labels, np.bincount(labels, minlength=2)


def frozen_unit(p: np.ndarray, *, reapply_minor_gate: bool) -> np.ndarray | None:
    labels, counts = deterministic_two_means(np.sqrt(p))
    if reapply_minor_gate and counts.min() / len(p) < THRESHOLD:
        return None
    if counts[0] > counts[1]:
        major, minor = 0, 1
    elif counts[1] > counts[0]:
        major, minor = 1, 0
    else:
        c0 = p[labels == 0].mean(axis=0)
        c1 = p[labels == 1].mean(axis=0)
        major, minor = (0, 1) if tuple(c0.tolist()) <= tuple(c1.tolist()) else (1, 0)
    delta = p[labels == minor].mean(axis=0) - p[labels == major].mean(axis=0)
    if not np.isclose(delta.sum(), 0.0, atol=1e-10):
        raise RuntimeError("Delta outside zero-sum subspace")
    norm = float(np.linalg.norm(delta))
    if not np.isfinite(norm) or norm <= EPS:
        return None
    return delta / norm


def plus_one_upper(observed: float, null: np.ndarray) -> float:
    return float((1 + np.sum(null >= observed)) / (len(null) + 1))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--measured-csv", type=Path, required=True)
    ap.add_argument("--frozen-delta-csv", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    work = load_classifiable(args.measured_csv)
    frozen = pd.read_csv(args.frozen_delta_csv)
    selected_species = frozen["species"].astype(str).tolist()
    q = one_vs_rest_contrast(9, focal_index=0)

    # Frozen observed species contributions and white-in-top-two decomposition.
    rows = []
    for _, row in frozen.iterrows():
        d = row[[f"delta_{c}" for c in BIO]].to_numpy(float)
        u = d / np.linalg.norm(d)
        rows.append(
            {
                "species": str(row["species"]),
                "w_i": float((u @ q) ** 2),
                "white_in_top_two": (
                    str(row["H1_coarse_primary_morph"]) == "white"
                    or str(row["H1_coarse_secondary_morph"]) == "white"
                ),
            }
        )
    contrib = pd.DataFrame(rows)
    observed = float(contrib["w_i"].mean())

    # Start at the observed coarse-state gate (185 species), then reapply the
    # continuous minor-cluster gate after each structured-null permutation.
    coarse_species = []
    for sp, g in work.groupby("species", sort=True):
        if len(g) < MIN_CLASSIFIABLE:
            continue
        second, *_ = coarse_gate(g)
        if second >= THRESHOLD:
            coarse_species.append(sp)

    sub = work.loc[work["species"].isin(coarse_species)].copy().reset_index(drop=True)
    compositions = sub[FRACTIONS].to_numpy(float)
    morphs = sub["morph"].astype(str).to_numpy()
    species = sub["species"].astype(str).to_numpy()
    species_indices = {sp: np.flatnonzero(species == sp) for sp in sorted(coarse_species)}
    morph_indices = {m: np.flatnonzero(morphs == m) for m in MORPHS}

    rng = np.random.default_rng(SEED)
    null_w = np.empty(N_REAPPLIED_NULL, dtype=float)
    null_n = np.empty(N_REAPPLIED_NULL, dtype=int)
    for b in range(N_REAPPLIED_NULL):
        perm = compositions.copy()
        for morph in MORPHS:
            idx = morph_indices[morph]
            if len(idx) > 1:
                perm[idx] = compositions[idx[rng.permutation(len(idx))]]
        units = []
        for sp in sorted(species_indices):
            u = frozen_unit(perm[species_indices[sp]], reapply_minor_gate=True)
            if u is not None:
                units.append(u)
        arr = np.vstack(units)
        null_w[b] = float(np.mean((arr @ q) ** 2))
        null_n[b] = len(arr)

    # disttrait real-data numerical audit.
    dt_axes = []
    direct_sensitive = []
    frozen_axis = {
        str(r["species"]): (
            r[[f"delta_{c}" for c in BIO]].to_numpy(float)
            / np.linalg.norm(r[[f"delta_{c}" for c in BIO]].to_numpy(float))
        )
        for _, r in frozen.iterrows()
    }
    for sp in selected_species:
        p = work.loc[work["species"] == sp, FRACTIONS].to_numpy(float)
        u = two_mode_axis(p).unit_axis
        dt_axes.append(u)
        cosine = float(abs(u @ frozen_axis[sp]))
        if cosine < 0.999999999:
            direct_sensitive.append({"species": sp, "abs_axis_cosine_to_frozen": cosine})
    dt_specieswise_w = float(np.mean((np.vstack(dt_axes) @ q) ** 2))

    selected = work.loc[work["species"].isin(selected_species)].copy().reset_index(drop=True)
    dt_structured = structured_alignment_null(
        selected[FRACTIONS].to_numpy(float),
        selected["species"].astype(str).to_numpy(),
        selected["morph"].astype(str).to_numpy(),
        q,
        n_permutations=1,
        seed=SEED,
        strata_order=MORPHS,
    )
    # Recover per-species axes along the structured-alignment observed path to
    # identify near-tie species responsible for numerical divergence.
    p = selected[FRACTIONS].to_numpy(float)
    p = p / p.sum(axis=1)[:, None]
    sp = selected["species"].astype(str).to_numpy()
    structured_sensitive = []
    for label in sorted(set(sp.tolist())):
        idx = np.flatnonzero(sp == label)
        u = two_mode_axis(p[idx]).unit_axis
        cosine = float(abs(u @ frozen_axis[label]))
        if cosine < 0.999999999:
            structured_sensitive.append({"species": label, "abs_axis_cosine_to_frozen": cosine})

    out = {
        "schema": "polymorphism_h2_posthoc_validity_diagnostic_v1",
        "date_jst": "2026-09-22",
        "role": "post_confirmatory_validity_and_robustness_diagnostic_only",
        "changes_frozen_h2_verdict": False,
        "white_in_top_two_coarse_morphs": {
            "species_with_white": int(contrib["white_in_top_two"].sum()),
            "species_total": int(len(contrib)),
            "fraction": float(contrib["white_in_top_two"].mean()),
            "mean_species_contribution_W_white_included": float(
                contrib.loc[contrib["white_in_top_two"], "w_i"].mean()
            ),
            "mean_species_contribution_W_white_not_included": float(
                contrib.loc[~contrib["white_in_top_two"], "w_i"].mean()
            ),
        },
        "gate_reapplied_structured_null": {
            "starting_species_after_coarse_gate": int(len(coarse_species)),
            "replicates": N_REAPPLIED_NULL,
            "seed": SEED,
            "null_mean": float(null_w.mean()),
            "null_median": float(np.median(null_w)),
            "null_q025": float(np.quantile(null_w, 0.025)),
            "null_q975": float(np.quantile(null_w, 0.975)),
            "retained_species_min": int(null_n.min()),
            "retained_species_median": float(np.median(null_n)),
            "retained_species_max": int(null_n.max()),
            "upper_p_plus_one": plus_one_upper(observed, null_w),
            "null_values_ge_observed": int(np.sum(null_w >= observed)),
        },
        "disttrait_numerical_audit": {
            "frozen_W": observed,
            "disttrait_specieswise_two_mode_W": dt_specieswise_w,
            "disttrait_structured_alignment_observed_W": float(dt_structured.observed),
            "specieswise_sensitive": direct_sensitive,
            "structured_path_sensitive": structured_sensitive,
            "interpretation": (
                "disttrait is a later generic implementation and is not a bitwise "
                "reproducer of the frozen biological H2 pipeline in near-tied "
                "two-means initializations"
            ),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(out, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

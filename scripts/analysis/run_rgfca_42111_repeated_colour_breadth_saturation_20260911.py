#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

BIOLOGICAL = ["white", "yellow_orange", "red_pink", "blue_purple"]
N_LEVELS = [100, 250, 500, 1000, 2000, 5000, 10000, 20000, 30000, 42111]
N_REP = 200
SEED0 = 20260911
L1_MAX = 0.05
MAX_PROP_DEV = 0.025
CLASSIFIABLE_DEV = 0.03
PASS_FRACTION = 0.95


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--measured", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args()


def metrics(q: pd.DataFrame, full_vec: np.ndarray, full_class_frac: float) -> dict[str, object]:
    classifiable = q["breadth_classifiable"].astype(bool).to_numpy()
    nclass = int(classifiable.sum())
    if nclass == 0:
        vec = np.zeros(4, dtype=float)
    else:
        morphs = q.loc[classifiable, "morph"].astype(str)
        vc = morphs.value_counts()
        vec = np.asarray([float(vc.get(m, 0) / nclass) for m in BIOLOGICAL], dtype=float)
    class_frac = float(nclass / len(q))
    l1 = float(np.abs(vec - full_vec).sum())
    max_dev = float(np.abs(vec - full_vec).max())
    denom = float(np.linalg.norm(vec) * np.linalg.norm(full_vec))
    cosine = float(np.dot(vec, full_vec) / denom) if denom > 0 else float("nan")
    cdev = abs(class_frac - full_class_frac)
    passed = bool(l1 <= L1_MAX and max_dev <= MAX_PROP_DEV and cdev <= CLASSIFIABLE_DEV)
    out: dict[str, object] = {
        "rows": int(len(q)),
        "classifiable": nclass,
        "classifiable_fraction": class_frac,
        "l1_to_full": l1,
        "cosine_to_full": cosine,
        "max_abs_colour_proportion_deviation": max_dev,
        "abs_classifiable_fraction_deviation": float(cdev),
        "simultaneous_pass": passed,
    }
    for m, v in zip(BIOLOGICAL, vec):
        out[f"p_{m}"] = float(v)
    return out


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    x = pd.read_csv(args.measured)
    if len(x) != 42111 or x["inat_taxon_id"].nunique() != 42111:
        raise RuntimeError("measured table is not exact 42,111 species")
    if "breadth_classifiable" not in x.columns:
        raise RuntimeError("measured table lacks frozen breadth_classifiable endpoint")
    x = x.sort_values("breadth_rank", kind="mergesort").reset_index(drop=True)
    if x["breadth_rank"].tolist() != list(range(1, 42112)):
        raise RuntimeError("breadth rank drift")

    full_class = x["breadth_classifiable"].astype(bool)
    nfull = int(full_class.sum())
    if nfull == 0:
        raise RuntimeError("complete experiment has zero classifiable anchors")
    vc = x.loc[full_class, "morph"].astype(str).value_counts()
    full_vec = np.asarray([float(vc.get(m, 0) / nfull) for m in BIOLOGICAL], dtype=float)
    full_class_frac = float(nfull / 42111)

    rows: list[dict[str, object]] = []
    nsp = len(x)
    for rep in range(N_REP):
        rng = np.random.default_rng(SEED0 + rep)
        perm = rng.permutation(nsp)
        for n in N_LEVELS:
            if n == 42111:
                q = x
            else:
                q = x.iloc[perm[:n]]
            r = metrics(q, full_vec, full_class_frac)
            r.update({"replicate": rep + 1, "seed": SEED0 + rep, "species_n": n})
            rows.append(r)
    raw = pd.DataFrame(rows)

    summaries: list[dict[str, object]] = []
    for n, g in raw.groupby("species_n", sort=True):
        pass_fraction = float(g["simultaneous_pass"].mean())
        summaries.append({
            "species_n": int(n),
            "replicates": int(len(g)),
            "simultaneous_pass_fraction": pass_fraction,
            "composition_saturated": bool(pass_fraction >= PASS_FRACTION),
            "l1_median": float(g["l1_to_full"].median()),
            "l1_q95": float(g["l1_to_full"].quantile(.95)),
            "max_abs_colour_dev_median": float(g["max_abs_colour_proportion_deviation"].median()),
            "max_abs_colour_dev_q95": float(g["max_abs_colour_proportion_deviation"].quantile(.95)),
            "classifiable_dev_median": float(g["abs_classifiable_fraction_deviation"].median()),
            "classifiable_dev_q95": float(g["abs_classifiable_fraction_deviation"].quantile(.95)),
            "cosine_median": float(g["cosine_to_full"].median()),
            "cosine_q05": float(g["cosine_to_full"].quantile(.05)),
        })
    summary = pd.DataFrame(summaries)
    eligible = summary.loc[(summary["species_n"] < 42111) & summary["composition_saturated"]]
    minimum = int(eligible["species_n"].min()) if len(eligible) else None
    verdict = "minimum_repeated_composition_dataset_identified" if minimum is not None else "no_subsampled_minimum_passed"

    raw.to_csv(args.output_dir / "repeated_colour_breadth_saturation_200x.csv.gz", index=False, compression="gzip", lineterminator="\n")
    summary.to_csv(args.output_dir / "repeated_colour_breadth_saturation_summary.csv", index=False, lineterminator="\n")
    result = {
        "analysis": "rgfca_42111_repeated_colour_breadth_saturation",
        "status": "complete_prespecified_repeated_saturation",
        "species_universe": 42111,
        "classifiable_full": nfull,
        "classifiable_fraction_full": full_class_frac,
        "full_four_state_composition": {m: float(v) for m, v in zip(BIOLOGICAL, full_vec)},
        "species_levels": N_LEVELS,
        "replicates_per_level": N_REP,
        "criterion": {
            "l1_max": L1_MAX,
            "maximum_absolute_colour_proportion_deviation": MAX_PROP_DEV,
            "absolute_classifiable_fraction_deviation": CLASSIFIABLE_DEV,
            "required_replicate_pass_fraction": PASS_FRACTION,
        },
        "verdict": verdict,
        "minimum_repeated_composition_species_n": minimum,
        "claim_boundary": "Minimum N applies only to the one-anchor species-equal observed-state composition, not within-species polymorphism or the geographic taxon-cell map."
    }
    (args.output_dir / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    min_text = f"{minimum:,}" if minimum is not None else "none below 42,111"
    (args.output_dir / "RESULT.md").write_text(
        "# RGFCA 42,111 — repeated colour-breadth saturation\n\n"
        f"- full classifiable anchors: **{nfull:,} / 42,111 ({full_class_frac:.3%})**\n"
        f"- repeated samples per N: **{N_REP}**\n"
        f"- minimum repeated composition dataset: **{min_text} species**\n"
        f"- verdict: **{verdict}**\n\n"
        "Pass requires >=95% of repeated samples to satisfy L1<=0.05, maximum single-colour deviation<=0.025, and classifiable-fraction deviation<=0.03 simultaneously.\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

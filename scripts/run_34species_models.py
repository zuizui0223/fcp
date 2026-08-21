#!/usr/bin/env python3
"""Run the frozen 34-species five-metric manuscript analysis."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from fcp_pipeline.constants import DEFAULT_PERMUTATIONS, DEFAULT_SEED, METRICS, MODEL_FORMULA
from fcp_pipeline.models import analyse_metrics
from fcp_pipeline.validation import validate_frozen_dataset, validate_model_results


def add_holm(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    out["wald_p_holm_five_metrics"] = multipletests(
        pd.to_numeric(out["wald_p_value_clustered"]), method="holm"
    )[1]
    out["permutation_p_holm_five_metrics"] = multipletests(
        pd.to_numeric(out["permutation_p_two_sided"]), method="holm"
    )[1]
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--permutations", type=int, default=DEFAULT_PERMUTATIONS)
    ap.add_argument("--seed", type=int, default=DEFAULT_SEED)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    data = validate_frozen_dataset(pd.read_csv(args.dataset))
    rng = np.random.default_rng(args.seed)

    rows, loo_rows = analyse_metrics(
        data, "baseline_unambiguous_34", METRICS, args.permutations, rng
    )
    results = add_holm(pd.DataFrame(rows))
    validate_model_results(results)
    loo = pd.DataFrame(loo_rows)

    results.to_csv(outdir / "environmental_niche_five_metric_models.csv", index=False)
    loo.to_csv(outdir / "environmental_niche_five_metric_leave_one_family_out.csv", index=False)
    data.sort_values("canonical_name").to_csv(
        outdir / "environmental_niche_five_metric_model_dataset.csv", index=False
    )

    manifest = {
        "status": "complete",
        "dataset_role": "frozen_34_species_baseline_only",
        "n_species": 34,
        "n_families": 25,
        "n_within": 20,
        "n_among": 14,
        "metrics": METRICS,
        "model_formula": MODEL_FORMULA,
        "covariance": "family-clustered sandwich",
        "permutations": args.permutations,
        "seed": args.seed,
        "multiplicity": "Holm correction across the five niche metrics within the frozen baseline",
        "results": results.to_dict("records"),
        "interpretation_guard": (
            "Five metrics are evaluated symmetrically. Occupied climatic breadth is not physiological "
            "tolerance; associations are exploratory and not causal."
        ),
    }
    (outdir / "environmental_niche_five_metric_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

BIOLOGICAL = ("blue_purple", "red_pink", "white", "yellow_orange")
N = 42111


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--measured", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    x = pd.read_csv(args.measured)
    if len(x) != N or x["inat_taxon_id"].nunique() != N:
        raise RuntimeError("measured table is not exact 42,111 species")
    classified = x["measurement_status"].astype(str).eq("classified_four_state_morph") & x["morph"].astype(str).isin(BIOLOGICAL)
    c = int(classified.sum())
    missing = N - c
    raw = {}
    bounds = {}
    rows = []
    for morph in BIOLOGICAL:
        k = int((classified & x["morph"].astype(str).eq(morph)).sum())
        raw[morph] = (k / c) if c else None
        lo = k / N
        hi = (k + missing) / N
        bounds[morph] = {"lower": lo, "upper": hi}
        rows.append({"morph": morph, "classified_count": k, "raw_classifiable_fraction": raw[morph], "full_denominator_lower": lo, "full_denominator_upper": hi})
    pd.DataFrame(rows).to_csv(args.output_dir / "colour_noassumption_bounds.csv", index=False, lineterminator="\n")
    result = {
        "analysis": "rgfca_42111_breadth_noassumption_bounds",
        "status": "complete_prespecified_noassumption_bounds",
        "species_denominator": N,
        "classifiable_species": c,
        "unclassified_or_missing_species": missing,
        "classifiable_fraction": c / N,
        "raw_classifiable_composition": raw,
        "full_denominator_noassumption_bounds": bounds,
        "bounds_assume_missing_at_random": False,
        "bounds_required_to_sum_to_one": False,
        "claim_boundary": "Bounds concern one frozen observed anchor per species; they do not identify modal species colour, polymorphism prevalence, D, or C*/S*."
    }
    (args.output_dir / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""External non-flower empirical transport using San Francisco street trees.

The analysis is deliberately bounded to a fixed public-domain subset stored in
packages/disttrait/fixtures/sf_street_trees_v1.csv. It tests whether the generic
continuous-trait spatial layer runs unchanged on a non-flower, non-FCP dataset.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from disttrait import (
    absolute_pairwise,
    continuous_spatial_permutation_null,
    distribution_spatial_association,
    great_circle_pairwise_km,
    species_equal_spatial_omnibus,
)


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_INPUT = ROOT / "packages" / "disttrait" / "fixtures" / "sf_street_trees_v1.csv"
DEFAULT_META = ROOT / "packages" / "disttrait" / "fixtures" / "sf_street_trees_v1_metadata.json"
N_PERMUTATIONS = 99
SEED = 20260919


def run(input_path: Path = DEFAULT_INPUT, metadata_path: Path = DEFAULT_META) -> dict:
    frame = pd.read_csv(input_path)
    required = {"tree_id", "taxon_label", "dbh", "latitude", "longitude"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"missing required columns: {missing}")
    if frame.isna().any().any():
        raise ValueError("fixture contains missing values")
    if (frame["dbh"] <= 0).any():
        raise ValueError("DBH must be positive")

    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if len(frame) != int(metadata["selection"]["total_rows"]):
        raise RuntimeError("fixture row-count drift")

    observed: list[float] = []
    nulls: list[np.ndarray] = []
    spreads: list[float] = []
    species_rows: list[dict[str, float | int | str]] = []

    for taxon, group in frame.groupby("taxon_label", sort=True):
        if len(group) != 80:
            raise RuntimeError(f"{taxon}: expected exactly 80 fixed observations")
        values = np.log(group["dbh"].to_numpy(float))
        rho, null = continuous_spatial_permutation_null(
            group["latitude"].to_numpy(float),
            group["longitude"].to_numpy(float),
            values,
            n_permutations=N_PERMUTATIONS,
            seed=SEED,
            key=str(taxon),
        )
        spread = float(np.std(values, ddof=1))
        observed.append(float(rho))
        nulls.append(null)
        spreads.append(spread)
        species_rows.append(
            {
                "taxon_label": str(taxon),
                "n": int(len(group)),
                "sd_log_dbh": spread,
                "spatial_rho": float(rho),
            }
        )

    observed_array = np.asarray(observed, dtype=float)
    null_array = np.vstack(nulls)
    conditioned = species_equal_spatial_omnibus(observed_array, null_array)
    spread_assoc = distribution_spatial_association(
        spreads,
        observed_array,
        spatial_null=null_array,
    )

    # Deliberately naive comparator: ignore taxon identity and correlate all
    # pairwise geographic distances with all pairwise absolute log-DBH differences.
    pooled_values = np.log(frame["dbh"].to_numpy(float))
    pooled_geo = great_circle_pairwise_km(
        frame["latitude"].to_numpy(float),
        frame["longitude"].to_numpy(float),
    )
    pooled_trait = absolute_pairwise(pooled_values)
    naive = spearmanr(pooled_geo, pooled_trait)

    return {
        "schema": "disttrait_sf_street_tree_empirical_transport_v1",
        "date_jst": "2026-09-19",
        "role": "external_nonflower_empirical_transport",
        "biological_claim": False,
        "source": metadata["source"],
        "selection": metadata["selection"],
        "trait": {
            "raw": "dbh",
            "transform": "natural log",
            "within_taxon_dissimilarity": "absolute difference in log(DBH)",
        },
        "analysis": {
            "taxa": int(len(species_rows)),
            "observations": int(len(frame)),
            "permutations_per_taxon": N_PERMUTATIONS,
            "species_conditioned_mean_rho": float(conditioned.observed),
            "species_conditioned_p_upper": float(conditioned.p_upper),
            "spread_spatial_rho": float(spread_assoc.observed),
            "spread_spatial_p_upper": float(spread_assoc.p_upper),
            "naive_pooled_rho": float(naive.statistic),
            "naive_pooled_p": float(naive.pvalue),
        },
        "taxa": species_rows,
        "claim_boundary": [
            "external empirical transport demonstration only",
            "taxon labels are inherited from the source inventory and are not taxonomically revalidated",
            "no ecological causal mechanism is inferred",
            "the naive pooled and species-conditioned statistics target different estimands",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_META)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    result = run(args.input, args.metadata)
    payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8")
    print(payload, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""External empirical transport on ShareTrait Gammarus metabolic-rate data."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from disttrait import continuous_spatial_permutation_null


ROOT = Path(__file__).resolve().parents[3]
FIXTURE = ROOT / "packages" / "disttrait" / "fixtures" / "sharetrait_gammarus_metabolic_v1.csv"
METADATA = ROOT / "packages" / "disttrait" / "fixtures" / "sharetrait_gammarus_metabolic_v1_metadata.json"
N_PERMUTATIONS = 999
SEED = 2026091909
KEY = "sharetrait|gammarus_insensibilis|metabolic_rate"


def _summary(x: np.ndarray) -> dict[str, float]:
    return {
        "mean": float(np.mean(x)),
        "sd": float(np.std(x, ddof=1)),
        "median": float(np.median(x)),
        "q025": float(np.quantile(x, 0.025)),
        "q975": float(np.quantile(x, 0.975)),
        "min": float(np.min(x)),
        "max": float(np.max(x)),
    }


def _upper_p(observed: float, null: np.ndarray) -> float:
    return float((1 + np.sum(null >= observed)) / (len(null) + 1))


def run() -> dict:
    frame = pd.read_csv(FIXTURE)
    metadata = json.loads(METADATA.read_text(encoding="utf-8"))

    required = {
        "species",
        "population_id",
        "site_id",
        "latitude",
        "longitude",
        "individual_id",
        "metabolic_rate",
        "metabolic_rate_unit",
    }
    missing = sorted(required - set(frame.columns))
    if missing:
        raise RuntimeError(f"fixture missing required columns: {missing}")
    if len(frame) != 375 or frame["individual_id"].nunique() != 375:
        raise RuntimeError("frozen 375-individual denominator drift")
    if frame["species"].nunique() != 1 or frame["species"].iloc[0] != "Gammarus_insensibilis":
        raise RuntimeError("species identity drift")
    if set(frame["metabolic_rate_unit"]) != {"Joule/day"}:
        raise RuntimeError("metabolic-rate unit drift")
    counts = frame.groupby("population_id").size().to_dict()
    if counts != {"TRAPOP076": 123, "TRAPOP077": 130, "TRAPOP078": 122}:
        raise RuntimeError(f"population denominator drift: {counts}")
    if frame["site_id"].nunique() != 3:
        raise RuntimeError("expected exactly three sampled sites")

    rate = frame["metabolic_rate"].to_numpy(float)
    if np.any(~np.isfinite(rate)) or np.any(rate <= 0):
        raise RuntimeError("metabolic-rate fixture contains nonpositive/nonfinite values")
    log_rate = np.log(rate)

    latitude = frame["latitude"].to_numpy(float)
    longitude = frame["longitude"].to_numpy(float)

    observed_log, null_log = continuous_spatial_permutation_null(
        latitude,
        longitude,
        log_rate,
        n_permutations=N_PERMUTATIONS,
        seed=SEED,
        key=KEY,
    )
    observed_raw, null_raw = continuous_spatial_permutation_null(
        latitude,
        longitude,
        rate,
        n_permutations=N_PERMUTATIONS,
        seed=SEED,
        key=KEY,
    )

    populations = []
    for pop, group in frame.groupby("population_id", sort=True):
        values = group["metabolic_rate"].to_numpy(float)
        populations.append(
            {
                "population_id": str(pop),
                "site_id": str(group["site_id"].iloc[0]),
                "location_label": str(group["location_label"].iloc[0]),
                "latitude": float(group["latitude"].iloc[0]),
                "longitude": float(group["longitude"].iloc[0]),
                "n": int(len(group)),
                "metabolic_rate": _summary(values),
                "log_metabolic_rate": _summary(np.log(values)),
            }
        )

    return {
        "schema": "disttrait_sharetrait_gammarus_empirical_v1",
        "date_jst": "2026-09-19",
        "role": "second_external_nonflower_empirical_transport",
        "source": metadata["source"],
        "fixture": {
            "path": str(FIXTURE.relative_to(ROOT)),
            "metadata": str(METADATA.relative_to(ROOT)),
            "rows": int(len(frame)),
            "individuals": int(frame["individual_id"].nunique()),
            "populations": int(frame["population_id"].nunique()),
            "sites": int(frame["site_id"].nunique()),
        },
        "trait": {
            "name": "standardized ShareTrait metabolic_rate value",
            "unit": "Joule/day",
            "primary_transform": "natural_log",
            "primary_pairwise_dissimilarity": "absolute_difference_in_log_metabolic_rate",
            "sensitivity_pairwise_dissimilarity": "absolute_difference_in_raw_metabolic_rate",
        },
        "primary_log_rate": {
            "observed_rho": float(observed_log),
            "p_upper": _upper_p(observed_log, null_log),
            "null_mean": float(np.mean(null_log)),
            "null_q025": float(np.quantile(null_log, 0.025)),
            "null_q975": float(np.quantile(null_log, 0.975)),
            "permutations": int(N_PERMUTATIONS),
        },
        "raw_rate_sensitivity": {
            "observed_rho": float(observed_raw),
            "p_upper": _upper_p(observed_raw, null_raw),
            "null_mean": float(np.mean(null_raw)),
            "null_q025": float(np.quantile(null_raw, 0.025)),
            "null_q975": float(np.quantile(null_raw, 0.975)),
            "permutations": int(N_PERMUTATIONS),
        },
        "population_descriptives": populations,
        "interpretation_boundary": {
            "supported_if_positive": "Metabolic-rate differences among the 375 observed individuals are geographically organized across the three sampled Gammarus populations under the frozen site coordinates.",
            "not_supported": [
                "a causal effect of latitude, temperature, climate, body mass or population history",
                "fine-scale within-site spatial organization because individuals inherit population-level site coordinates",
                "a cross-species comparative result because this fixture contains one species",
                "global transport validity for all continuous traits",
            ],
        },
    }


def main() -> int:
    print(json.dumps(run(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

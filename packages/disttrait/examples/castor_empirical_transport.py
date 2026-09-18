#!/usr/bin/env python3
"""External empirical transport of disttrait to continuous Ricinus traits.

The source dataset is not redistributed here. It is read from a commit-pinned
public GitHub URL maintained by the original study authors.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from disttrait import continuous_spatial_permutation_null


SOURCE_REPOSITORY = "The-Frederickson-Lab/castor_latitudinal_gradient"
SOURCE_COMMIT = "ff8cb26f83cd271c516e83510cd06c4b5563db22"
SOURCE_PATH = "morph_natint_20May24.csv"
SOURCE_GIT_BLOB = "c4dcca178d7c885585edbe0cdbeee89eef2a484a"
SOURCE_URL = (
    "https://raw.githubusercontent.com/"
    f"{SOURCE_REPOSITORY}/{SOURCE_COMMIT}/{SOURCE_PATH}"
)
N_PERMUTATIONS = 999
SEED = 20260919

TRAITS = {
    "leaf_area_proxy": "areaavg",
    "herbivory_fraction": "pherb",
    "leaf_blade_efn_density": "EFNLB",
}


def _one_trait(
    frame: pd.DataFrame,
    *,
    label: str,
    column: str,
) -> dict[str, float | int | str]:
    latitude = pd.to_numeric(frame["Latitude"], errors="coerce")
    longitude = pd.to_numeric(frame["Longitude"], errors="coerce")
    values = pd.to_numeric(frame[column], errors="coerce")
    keep = latitude.notna() & longitude.notna() & values.notna()
    work = pd.DataFrame(
        {
            "latitude": latitude[keep].astype(float),
            "longitude": longitude[keep].astype(float),
            "value": values[keep].astype(float),
        }
    )
    if len(work) < 20:
        raise RuntimeError(f"{label}: fewer than 20 complete georeferenced observations")

    observed, null = continuous_spatial_permutation_null(
        work["latitude"].to_numpy(float),
        work["longitude"].to_numpy(float),
        work["value"].to_numpy(float),
        n_permutations=N_PERMUTATIONS,
        seed=SEED,
        key=f"castor|{label}",
    )
    p_upper = float((1 + np.sum(null >= observed)) / (len(null) + 1))
    return {
        "source_column": column,
        "n_complete": int(len(work)),
        "observed_spatial_rho": float(observed),
        "p_upper": p_upper,
        "null_mean": float(np.mean(null)),
        "null_q025": float(np.quantile(null, 0.025)),
        "null_median": float(np.median(null)),
        "null_q975": float(np.quantile(null, 0.975)),
    }


def run_transport() -> dict:
    frame = pd.read_csv(SOURCE_URL)
    required = {"Plant ID", "Latitude", "Longitude", *TRAITS.values()}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise RuntimeError(f"source schema drift: missing columns {missing}")

    result = {
        "schema": "disttrait_castor_empirical_transport_v1",
        "source": {
            "repository": SOURCE_REPOSITORY,
            "commit": SOURCE_COMMIT,
            "path": SOURCE_PATH,
            "git_blob_sha": SOURCE_GIT_BLOB,
            "url": SOURCE_URL,
        },
        "dataset": {
            "taxon": "Ricinus communis",
            "rows": int(len(frame)),
            "plant_ids": int(frame["Plant ID"].astype(str).nunique()),
            "external_dataset_redistributed": False,
        },
        "analysis": {
            "trait_distance": "absolute pairwise difference",
            "spatial_statistic": "Spearman(pairwise great-circle distance, pairwise absolute trait difference)",
            "null": "complete-observation vertex permutation with geographic geometry fixed",
            "permutations": N_PERMUTATIONS,
            "seed": SEED,
        },
        "traits": {},
        "claim_boundary": [
            "external empirical transport demonstration, not a preregistered biological hypothesis test",
            "one plant species does not validate cross-species comparative performance",
            "source sampling design and observation biases are inherited from the external study",
            "spatial association does not identify an ecological mechanism",
        ],
    }
    for label, column in TRAITS.items():
        result["traits"][label] = _one_trait(frame, label=label, column=column)
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run_transport()
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

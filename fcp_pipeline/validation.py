"""Hard validation gates for the frozen 34-species manuscript pipeline."""
from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd

from .constants import EXPECTED_COUNTS, METRICS

FROZEN_DATASET_SHA256 = "bdc06dd671f41ce062ebf4ba687437909d9617b268657504c1c6c5e991d417ed"

REQUIRED_COLUMNS = {
    "canonical_name",
    "family",
    "spatial_scale",
    "classification_source",
    "n_climate_cells",
    *METRICS,
}


def validate_frozen_dataset(data: pd.DataFrame) -> pd.DataFrame:
    missing = sorted(REQUIRED_COLUMNS - set(data.columns))
    if missing:
        raise ValueError(f"Frozen dataset missing columns: {missing}")
    if data["canonical_name"].duplicated().any():
        dup = data.loc[data["canonical_name"].duplicated(), "canonical_name"].tolist()
        raise ValueError(f"Duplicate frozen species: {dup}")
    if len(data) != EXPECTED_COUNTS["species"]:
        raise ValueError(f"Expected 34 frozen species; found {len(data)}")
    if data["family"].nunique() != EXPECTED_COUNTS["families"]:
        raise ValueError(f"Expected 25 families; found {data['family'].nunique()}")
    within = int((data["spatial_scale"] == "within_population").sum())
    among = int((data["spatial_scale"] == "among_population").sum())
    if within != EXPECTED_COUNTS["within_population"] or among != EXPECTED_COUNTS["among_population"]:
        raise ValueError(f"Expected 20 within / 14 among; found {within} / {among}")
    if set(data["classification_source"].dropna()) != {"baseline_unambiguous"}:
        raise ValueError("Frozen paper dataset must contain baseline_unambiguous classifications only")
    if (pd.to_numeric(data["n_climate_cells"], errors="coerce") < 20).any():
        raise ValueError("Frozen paper dataset contains a species below the 20-cell threshold")
    return data.copy()


def validate_frozen_file(path: str | Path) -> pd.DataFrame:
    """Validate both exact bytes and biological/statistical invariants of the production freeze."""
    path = Path(path)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != FROZEN_DATASET_SHA256:
        raise ValueError(f"Frozen dataset SHA-256 changed: {digest} != {FROZEN_DATASET_SHA256}")
    data = validate_frozen_dataset(pd.read_csv(path))
    names = data["canonical_name"].astype(str).tolist()
    if names != sorted(names):
        raise ValueError("Production frozen dataset must remain in canonical_name ascending order")
    return data


def validate_model_results(results: pd.DataFrame) -> None:
    if len(results) != len(METRICS) or set(results["metric"]) != set(METRICS):
        raise ValueError("Exactly five symmetric niche metrics are required")
    if not (results["analysis_status"] == "complete").all():
        raise ValueError("One or more frozen five-metric models are incomplete")
    if not (pd.to_numeric(results["permutations_valid"]) == 9999).all():
        raise ValueError("Every frozen model must contain 9,999 valid permutations")

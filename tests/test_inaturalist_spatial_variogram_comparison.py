from __future__ import annotations

import csv
from pathlib import Path

import pytest

from scripts.data.compare_inaturalist_spatial_variograms import compare_variograms


FIELDS = [
    "canonical_name",
    "equal_pair_count_bin",
    "n_pairs",
    "distance_km_min",
    "distance_km_median",
    "distance_km_max",
    "standardized_colour_distance_median",
]


def write(path: Path, median: float) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerow(
            {
                "canonical_name": "Species one",
                "equal_pair_count_bin": 1,
                "n_pairs": 10,
                "distance_km_min": 1.0,
                "distance_km_median": 2.0,
                "distance_km_max": 3.0,
                "standardized_colour_distance_median": median,
            }
        )


def test_comparison_accepts_platform_roundoff(tmp_path: Path) -> None:
    expected = tmp_path / "expected.csv"
    observed = tmp_path / "observed.csv"
    write(expected, 1.1234567890123)
    write(observed, 1.1234567890124)
    compare_variograms(expected, observed)


def test_comparison_rejects_scientific_difference(tmp_path: Path) -> None:
    expected = tmp_path / "expected.csv"
    observed = tmp_path / "observed.csv"
    write(expected, 1.0)
    write(observed, 1.000001)
    with pytest.raises(RuntimeError, match="numerical mismatch"):
        compare_variograms(expected, observed)

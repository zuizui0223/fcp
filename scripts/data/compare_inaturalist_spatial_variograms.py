#!/usr/bin/env python3
"""Compare descriptive variograms across numerical platforms."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


IDENTITY_FIELDS = {
    "canonical_name",
    "equal_pair_count_bin",
    "n_pairs",
}
NUMERIC_FIELDS = {
    "distance_km_min",
    "distance_km_median",
    "distance_km_max",
    "standardized_colour_distance_median",
}
DEFAULT_ABSOLUTE_TOLERANCE = 1e-10


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), list(reader)


def compare_variograms(
    expected_path: Path,
    observed_path: Path,
    *,
    absolute_tolerance: float = DEFAULT_ABSOLUTE_TOLERANCE,
) -> None:
    expected_fields, expected = read_csv(expected_path)
    observed_fields, observed = read_csv(observed_path)
    required = IDENTITY_FIELDS | NUMERIC_FIELDS
    if expected_fields != observed_fields or set(expected_fields) != required:
        raise RuntimeError("variogram schema mismatch")
    if len(expected) != len(observed):
        raise RuntimeError("variogram row-count mismatch")
    for index, (first, second) in enumerate(zip(expected, observed), start=1):
        for field in IDENTITY_FIELDS:
            if first[field] != second[field]:
                raise RuntimeError(f"row {index} identity mismatch: {field}")
        for field in NUMERIC_FIELDS:
            if not math.isclose(
                float(first[field]),
                float(second[field]),
                rel_tol=0.0,
                abs_tol=absolute_tolerance,
            ):
                difference = abs(float(first[field]) - float(second[field]))
                raise RuntimeError(
                    f"row {index} numerical mismatch: {field}; absolute difference={difference}"
                )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("expected", type=Path)
    parser.add_argument("observed", type=Path)
    parser.add_argument(
        "--absolute-tolerance", type=float, default=DEFAULT_ABSOLUTE_TOLERANCE
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    compare_variograms(
        args.expected,
        args.observed,
        absolute_tolerance=args.absolute_tolerance,
    )
    print(
        f"variogram comparison passed with absolute_tolerance={args.absolute_tolerance:g}"
    )


if __name__ == "__main__":
    main()

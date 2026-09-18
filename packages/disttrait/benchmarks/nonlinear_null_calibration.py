#!/usr/bin/env python3
"""Targeted null-calibration audit for the nonlinear benchmark."""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from nonlinear_curvature_surface import (
    ALPHA,
    MISSING_FRACTIONS,
    REVERSAL_FRACTIONS,
    _cell_seed,
    run_world,
)


WORLDS_PER_CELL = 200


def run_calibration(*, worlds_per_cell: int = WORLDS_PER_CELL) -> dict:
    cells: list[dict[str, float | int]] = []
    raw_rows: list[dict[str, float]] = []

    for reversal_fraction in REVERSAL_FRACTIONS:
        for missing_fraction in MISSING_FRACTIONS:
            rows = [
                run_world(
                    seed=_cell_seed(0.0, reversal_fraction, missing_fraction, world),
                    effect_size=0.0,
                    reversal_fraction=reversal_fraction,
                    missing_fraction=missing_fraction,
                )
                for world in range(int(worlds_per_cell))
            ]
            frame = pd.DataFrame(rows)
            cells.append(
                {
                    "reversal_fraction": float(reversal_fraction),
                    "missing_fraction": float(missing_fraction),
                    "worlds": int(worlds_per_cell),
                    "equal_fpr": float(np.mean(frame["equal_matched_p"] < ALPHA)),
                    "linear_fpr": float(np.mean(frame["common_linear_p"] < ALPHA)),
                    "quadratic_fpr": float(np.mean(frame["common_quadratic_p"] < ALPHA)),
                    "meta_mean_fpr": float(np.mean(frame["meta_mean_p"] < ALPHA)),
                    "meta_heterogeneity_fpr": float(
                        np.mean(frame["meta_heterogeneity_p"] < ALPHA)
                    ),
                    "meta_omnibus_fpr": float(np.mean(frame["meta_omnibus_p"] < ALPHA)),
                }
            )
            for row in rows:
                raw_rows.append(
                    {
                        "reversal_fraction": float(reversal_fraction),
                        "missing_fraction": float(missing_fraction),
                        **row,
                    }
                )

    cells_frame = pd.DataFrame(cells)
    raw = pd.DataFrame(raw_rows)
    by_missingness = []
    for missing_fraction, group in raw.groupby("missing_fraction", sort=True):
        by_missingness.append(
            {
                "missing_fraction": float(missing_fraction),
                "worlds": int(len(group)),
                "equal_fpr": float(np.mean(group["equal_matched_p"] < ALPHA)),
                "linear_fpr": float(np.mean(group["common_linear_p"] < ALPHA)),
                "quadratic_fpr": float(np.mean(group["common_quadratic_p"] < ALPHA)),
                "meta_mean_fpr": float(np.mean(group["meta_mean_p"] < ALPHA)),
                "meta_heterogeneity_fpr": float(
                    np.mean(group["meta_heterogeneity_p"] < ALPHA)
                ),
                "meta_omnibus_fpr": float(np.mean(group["meta_omnibus_p"] < ALPHA)),
            }
        )

    return {
        "cells": cells,
        "by_missingness": by_missingness,
        "summary": {
            "worlds_per_cell": int(worlds_per_cell),
            "total_null_worlds": int(len(raw)),
            "max_equal_fpr": float(cells_frame["equal_fpr"].max()),
            "max_linear_fpr": float(cells_frame["linear_fpr"].max()),
            "max_quadratic_fpr": float(cells_frame["quadratic_fpr"].max()),
            "max_meta_mean_fpr": float(cells_frame["meta_mean_fpr"].max()),
            "max_meta_heterogeneity_fpr": float(
                cells_frame["meta_heterogeneity_fpr"].max()
            ),
            "max_meta_omnibus_fpr": float(cells_frame["meta_omnibus_fpr"].max()),
            "overall_equal_fpr": float(np.mean(raw["equal_matched_p"] < ALPHA)),
            "overall_meta_omnibus_fpr": float(
                np.mean(raw["meta_omnibus_p"] < ALPHA)
            ),
        },
    }


def main() -> int:
    print(json.dumps(run_calibration(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Summarize pre-colour RGFCA observation-bias control surfaces.

This script is descriptive and outcome blind. It combines the frozen target-group
sampling-availability surface with the frozen all-taxa platform-effort negative
control. It does not open candidate image pixels or flower-colour outcomes.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
TARGET = ROOT / "data/frozen/global_monte_carlo_sampling_availability_cells_v1.csv"
PLATFORM = ROOT / "data/frozen/global_monte_carlo_platform_effort_cells_v1.csv"
OUT = ROOT / "docs/supporting/rgfca_observation_bias_integrated_summary_v1.json"


def gini_nonnegative(values: np.ndarray) -> float:
    x = np.asarray(values, dtype=float)
    if np.any(x < 0):
        raise ValueError("Gini requires nonnegative values")
    if x.size == 0 or float(x.sum()) == 0.0:
        return float("nan")
    x = np.sort(x)
    n = x.size
    return float((2.0 * np.dot(np.arange(1, n + 1), x) / (n * x.sum())) - (n + 1) / n)


def qdict(values: pd.Series) -> dict[str, float]:
    v = pd.to_numeric(values, errors="coerce").dropna()
    qs = v.quantile([0.1, 0.25, 0.5, 0.75, 0.9])
    return {"q10": float(qs.loc[0.1]), "q25": float(qs.loc[0.25]), "median": float(qs.loc[0.5]), "q75": float(qs.loc[0.75]), "q90": float(qs.loc[0.9])}


def top_n_cells(frame: pd.DataFrame, value: str, n: int) -> list[int]:
    ranked = frame.sort_values([value, "cell_id"], ascending=[False, True], kind="mergesort")
    return ranked.head(n)["cell_id"].astype(int).tolist()


def main() -> int:
    target = pd.read_csv(TARGET)
    platform = pd.read_csv(PLATFORM)
    required_target = {"cell_id", "all_research_photo_records", "top_decile_all_photo_effort"}
    required_platform = {"cell_id", "all_taxa_research_photo_records", "target_group_all_research_photo_records"}
    if not required_target.issubset(target.columns):
        raise RuntimeError(f"target surface missing columns: {sorted(required_target - set(target.columns))}")
    if not required_platform.issubset(platform.columns):
        raise RuntimeError(f"platform surface missing columns: {sorted(required_platform - set(platform.columns))}")
    if len(target) != 162 or len(platform) != 162:
        raise RuntimeError("expected exactly 162 equal-area cells in both surfaces")
    if set(target.cell_id.astype(int)) != set(platform.cell_id.astype(int)):
        raise RuntimeError("cell-id sets differ between frozen surfaces")

    df = target[["cell_id", "all_research_photo_records", "top_decile_all_photo_effort"]].merge(
        platform[["cell_id", "all_taxa_research_photo_records", "target_group_all_research_photo_records"]],
        on="cell_id",
        how="inner",
        validate="one_to_one",
    )
    if not np.array_equal(
        df["all_research_photo_records"].astype(np.int64).to_numpy(),
        df["target_group_all_research_photo_records"].astype(np.int64).to_numpy(),
    ):
        raise RuntimeError("platform file target-group copy does not match frozen target-group surface")

    df["log_target"] = np.log1p(df["all_research_photo_records"].astype(float))
    df["log_platform"] = np.log1p(df["all_taxa_research_photo_records"].astype(float))
    df["target_fraction_of_platform"] = np.where(
        df["all_taxa_research_photo_records"] > 0,
        df["all_research_photo_records"] / df["all_taxa_research_photo_records"],
        np.nan,
    )

    X = np.column_stack([np.ones(len(df)), df["log_platform"].to_numpy(float)])
    y = df["log_target"].to_numpy(float)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    fitted = X @ beta
    residual = y - fitted
    ss_res = float(np.sum(residual**2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = float(1.0 - ss_res / ss_tot) if ss_tot > 0 else float("nan")
    df["log_target_given_platform_residual"] = residual

    frozen_target_top = sorted(df.loc[df["top_decile_all_photo_effort"].astype(bool), "cell_id"].astype(int).tolist())
    if len(frozen_target_top) != 17:
        raise RuntimeError(f"expected 17 frozen target-group high-effort cells, found {len(frozen_target_top)}")
    platform_top = sorted(top_n_cells(df, "all_taxa_research_photo_records", 17))
    target_set = set(frozen_target_top)
    platform_set = set(platform_top)
    intersection = sorted(target_set & platform_set)
    union = sorted(target_set | platform_set)

    target_total = int(df["all_research_photo_records"].sum())
    platform_total = int(df["all_taxa_research_photo_records"].sum())
    platform_top_mask = df["cell_id"].isin(platform_set)
    target_top_mask = df["cell_id"].isin(target_set)

    out = {
        "protocol": "rgfca-observation-bias-integrated-summary-v1",
        "status": "complete_precolour_integrated_observation_bias_summary",
        "candidate_image_pixels_opened": False,
        "flower_colour_used": False,
        "cells": int(len(df)),
        "global_counts": {
            "target_group_records": target_total,
            "all_taxa_platform_records": platform_total,
            "target_group_fraction_of_all_taxa_records": float(target_total / platform_total),
            "target_group_zero_cells": int((df["all_research_photo_records"] == 0).sum()),
            "platform_zero_cells": int((df["all_taxa_research_photo_records"] == 0).sum()),
        },
        "equal_area_concentration": {
            "target_group_gini_recomputed": gini_nonnegative(df["all_research_photo_records"].to_numpy(float)),
            "platform_gini_recomputed": gini_nonnegative(df["all_taxa_research_photo_records"].to_numpy(float)),
            "spearman_log_target_vs_log_platform": float(df[["log_target", "log_platform"]].corr(method="spearman").iloc[0, 1]),
            "pearson_log_target_vs_log_platform": float(df[["log_target", "log_platform"]].corr(method="pearson").iloc[0, 1]),
            "ols_log_target_on_log_platform": {
                "intercept": float(beta[0]),
                "slope": float(beta[1]),
                "r_squared_descriptive": r2,
                "note": "Descriptive equal-area-cell regression only; spatial autocorrelation is not used for inference here.",
            },
        },
        "top_decile_geometry": {
            "n_cells_each": 17,
            "frozen_target_group_top_cells": frozen_target_top,
            "platform_top_cells": platform_top,
            "intersection_cells": intersection,
            "intersection_count": len(intersection),
            "union_count": len(union),
            "jaccard": float(len(intersection) / len(union)),
            "target_records_share_inside_platform_top17": float(df.loc[platform_top_mask, "all_research_photo_records"].sum() / target_total),
            "platform_records_share_inside_frozen_target_top17": float(df.loc[target_top_mask, "all_taxa_research_photo_records"].sum() / platform_total),
        },
        "target_fraction_across_platform_nonzero_cells": qdict(df.loc[df["all_taxa_research_photo_records"] > 0, "target_fraction_of_platform"]),
        "platform_conditioned_target_residual": {
            "quantiles": qdict(df["log_target_given_platform_residual"]),
            "largest_positive_cell_ids": df.nlargest(10, "log_target_given_platform_residual")["cell_id"].astype(int).tolist(),
            "largest_negative_cell_ids": df.nsmallest(10, "log_target_given_platform_residual")["cell_id"].astype(int).tolist(),
            "interpretation": "Residuals are descriptive target-group over/under-representation relative to platform-wide record density, not corrected flower-colour weights.",
        },
        "interpretation": (
            "Platform-wide and target-group recording intensity are compared before any RGFCA colour outcome. "
            "High concordance indicates that generic platform activity is a major component of the target-group availability surface; "
            "remaining residual structure can also contain biological target-group availability and unmeasured observation processes."
        ),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

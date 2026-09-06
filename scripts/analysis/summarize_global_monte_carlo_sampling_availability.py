#!/usr/bin/env python3
"""Summarize the frozen RGFCA metadata-only sampling-availability surface."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
CELLS = ROOT / "data/frozen/global_monte_carlo_sampling_availability_cells_v1.csv"
MANIFEST = ROOT / "docs/supporting/global_monte_carlo_sampling_availability_manifest_v1.json"
OUT = ROOT / "docs/supporting/global_monte_carlo_sampling_availability_summary_v1.json"

LAYERS = [
    "all_research_photo_records",
    "license_eligible_photo_records",
    "flowering_annotated_eligible_records",
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def gini(values: np.ndarray) -> float:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x) & (x >= 0)]
    if len(x) == 0 or float(x.sum()) == 0.0:
        return float("nan")
    x = np.sort(x)
    n = len(x)
    return float((2.0 * np.sum((np.arange(1, n + 1)) * x) / (n * x.sum())) - (n + 1.0) / n)


def qdict(values: np.ndarray) -> dict[str, float | None]:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return {"q10": None, "q25": None, "median": None, "q75": None, "q90": None}
    qs = np.quantile(x, [0.10, 0.25, 0.50, 0.75, 0.90])
    return {k: float(v) for k, v in zip(["q10", "q25", "median", "q75", "q90"], qs)}


def spearman(x: pd.Series, y: pd.Series) -> float | None:
    pair = pd.concat([pd.to_numeric(x, errors="coerce"), pd.to_numeric(y, errors="coerce")], axis=1).dropna()
    if len(pair) < 3:
        return None
    r = pair.rank(method="average")
    value = float(r.iloc[:, 0].corr(r.iloc[:, 1]))
    return value if np.isfinite(value) else None


def main() -> None:
    if OUT.exists():
        raise RuntimeError("sampling-availability summary already exists; refusing overwrite")
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    if manifest.get("status") != "complete_metadata_only_sampling_availability_surface":
        raise RuntimeError("sampling-availability census is not complete")
    if manifest.get("candidate_image_pixels_opened") is not False or manifest.get("flower_colour_used") is not False:
        raise RuntimeError("sampling-availability lineage unexpectedly opened colour pixels")

    df = pd.read_csv(CELLS)
    if len(df) != 162 or df["cell_id"].nunique() != 162:
        raise RuntimeError("sampling-availability cell frame drift")
    if int(manifest.get("request_errors", -1)) != 0:
        raise RuntimeError("summary v1 expects complete error-free census")

    layer_summary: dict[str, object] = {}
    for layer in LAYERS:
        x = pd.to_numeric(df[layer], errors="raise").to_numpy(dtype=float)
        total = float(x.sum())
        top_mask = df["top_decile_all_photo_effort"].astype(bool).to_numpy()
        layer_summary[layer] = {
            "total_records": int(total),
            "zero_cells": int(np.count_nonzero(x == 0)),
            "nonzero_cells": int(np.count_nonzero(x > 0)),
            "maximum_cell_records": int(np.max(x)),
            "maximum_cell_share_of_global_records": (float(np.max(x) / total) if total > 0 else None),
            "gini_across_equal_area_cells": gini(x),
            "all_photo_effort_top_decile_cells_share_of_records": (float(x[top_mask].sum() / total) if total > 0 else None),
            "all_cells_quantiles": qdict(x),
            "nonzero_cell_quantiles": qdict(x[x > 0]),
        }

    licence_frac = pd.to_numeric(df["licence_eligible_fraction"], errors="coerce").to_numpy(dtype=float)
    flowering_frac = pd.to_numeric(df["flowering_annotation_fraction_given_licence"], errors="coerce").to_numpy(dtype=float)
    output = {
        "protocol": "rgfca-observation-bias-robustness-v1",
        "status": "complete_sampling_availability_bias_summary",
        "candidate_image_pixels_opened": False,
        "flower_colour_used": False,
        "cells": 162,
        "layer_summary": layer_summary,
        "selection_filter_summary": {
            "licence_eligible_fraction_across_defined_cells": qdict(licence_frac),
            "flowering_annotation_fraction_given_licence_across_defined_cells": qdict(flowering_frac),
            "spearman_log_all_vs_log_license": spearman(df["log1p_all_research_photo_records"], df["log1p_license_eligible_photo_records"]),
            "spearman_log_all_vs_log_flowering_annotated": spearman(df["log1p_all_research_photo_records"], df["log1p_flowering_annotated_eligible_records"]),
            "spearman_log_license_vs_log_flowering_annotated": spearman(df["log1p_license_eligible_photo_records"], df["log1p_flowering_annotated_eligible_records"]),
        },
        "robustness_consequence": {
            "primary_field_replaced_by_effort_correction": False,
            "top_decile_high_effort_cells_deleted_in_frozen_sensitivity": int(df["top_decile_all_photo_effort"].sum()),
            "availability_layers_used_as_pre_frozen_negative_controls": True,
            "same_species_conditioned_999_colour_null_required": True,
        },
        "interpretation": "The equal-area public-photo frame is expected to be strongly heterogeneous. These summaries quantify that heterogeneity before any RGFCA colour field and define robustness covariates; they are not estimates of true biological sampling probabilities.",
        "lineage": {
            "cell_surface_sha256": sha256_file(CELLS),
            "sampling_availability_manifest_sha256": sha256_file(MANIFEST),
        },
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()

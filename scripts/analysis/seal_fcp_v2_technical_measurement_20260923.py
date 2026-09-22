#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from fcp_pipeline.measurement_validity_v2 import seal_technical_table

EXPECTED_ROWS = 40_000
EXPECTED_PARTITIONS = 256
EXPECTED_HEAVY_ROWS = 8_000


def _q(values: pd.Series, q: float) -> float:
    x = pd.to_numeric(values, errors="coerce")
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return float("nan")
    return float(np.quantile(x.to_numpy(float), q))


def _bool_series(values: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(values):
        return values.fillna(False).astype(bool)
    return values.fillna("").astype(str).str.strip().str.casefold().isin(
        {"1", "true", "yes"}
    )


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--partition-results-dir", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args()

    technical_files = sorted(args.partition_results_dir.rglob("technical_b*_s*_p*.csv"))
    acquisition_files = sorted(args.partition_results_dir.rglob("acquisition_b*_s*_p*.csv"))
    if len(technical_files) != EXPECTED_PARTITIONS:
        raise RuntimeError(
            f"technical partition census incomplete: {len(technical_files)} != {EXPECTED_PARTITIONS}"
        )
    if len(acquisition_files) != EXPECTED_PARTITIONS:
        raise RuntimeError(
            f"acquisition partition census incomplete: {len(acquisition_files)} != {EXPECTED_PARTITIONS}"
        )

    acquisition = pd.concat(
        [pd.read_csv(path, dtype={"measurement_id": str}) for path in acquisition_files],
        ignore_index=True,
    )
    if len(acquisition) != EXPECTED_ROWS:
        raise RuntimeError(f"acquisition row census drift: {len(acquisition)}")
    if acquisition["measurement_id"].nunique() != EXPECTED_ROWS:
        raise RuntimeError("acquisition measurement IDs are not unique")

    technical_parts = []
    for path in technical_files:
        frame = pd.read_csv(path, dtype={"measurement_id": str})
        if len(frame):
            technical_parts.append(frame)
    technical = (
        pd.concat(technical_parts, ignore_index=True)
        if technical_parts
        else pd.DataFrame(columns=["measurement_id"])
    )
    if len(technical) and technical["measurement_id"].duplicated().any():
        raise RuntimeError("duplicate technical measurement IDs")
    if not set(technical.get("measurement_id", pd.Series(dtype=str))) <= set(
        acquisition["measurement_id"]
    ):
        raise RuntimeError("technical IDs outside acquisition census")

    merged = acquisition.merge(
        technical,
        on=["measurement_id", "source_image_sha256", "heavy_counterfactual"],
        how="left",
        validate="one_to_one",
        suffixes=("", "_technical"),
    )
    if len(merged) != EXPECTED_ROWS:
        raise RuntimeError("technical merge lost frozen rows")

    acquired = merged["acquisition_status"].fillna("").eq(
        "acquired_and_decode_verified"
    )
    available = merged["technical_status"].fillna("").eq(
        "technical_metrics_available"
    )
    merged["technical_metrics_available"] = available
    merged["heavy_counterfactual"] = _bool_series(merged["heavy_counterfactual"])

    if int(merged["heavy_counterfactual"].sum()) != EXPECTED_HEAVY_ROWS:
        raise RuntimeError("heavy-counterfactual denominator drift")

    flower_q95 = _q(merged.loc[available, "base_flower_near_clip_fraction"], 0.95)
    background_q95 = _q(
        merged.loc[available, "base_background_near_clip_fraction"], 0.95
    )
    separation_q05 = _q(
        merged.loc[available, "base_flower_background_lab_distance"], 0.05
    )
    if not all(np.isfinite(v) for v in (flower_q95, background_q95, separation_q05)):
        raise RuntimeError("technical stratum quantiles are not estimable")

    flower_threshold = float(max(0.01, flower_q95))
    background_threshold = float(max(0.01, background_q95))

    merged["stratum_high_flower_near_clip"] = (
        available
        & pd.to_numeric(
            merged["base_flower_near_clip_fraction"], errors="coerce"
        ).gt(flower_threshold)
    )
    merged["stratum_high_background_near_clip"] = (
        available
        & pd.to_numeric(
            merged["base_background_near_clip_fraction"], errors="coerce"
        ).gt(background_threshold)
    )
    merged["stratum_low_flower_background_lab_separation"] = (
        available
        & pd.to_numeric(
            merged["base_flower_background_lab_distance"], errors="coerce"
        ).lt(separation_q05)
    )

    flip_iou = pd.to_numeric(
        merged["base_horizontal_flip_mask_iou"], errors="coerce"
    )
    flip_de = pd.to_numeric(
        merged["base_horizontal_flip_colour_delta_e"], errors="coerce"
    )
    base_admitted = merged["base_roi_status"].fillna("").eq(
        "automated_colour_state_admitted"
    )
    merged["stratum_roi_unstable_fixed"] = available & (
        (~base_admitted) | flip_iou.lt(0.5) | flip_de.gt(5.0)
    )

    heavy_available = (
        available
        & merged["heavy_counterfactual"]
        & pd.to_numeric(merged["jitter_valid_variants"], errors="coerce").gt(0)
    )
    jitter_iou = pd.to_numeric(merged["jitter_mask_iou_min"], errors="coerce")
    jitter_de = pd.to_numeric(
        merged["jitter_colour_delta_e_max"], errors="coerce"
    )
    merged["stratum_heavy_prompt_unstable"] = heavy_available & (
        jitter_iou.lt(0.5) | jitter_de.gt(5.0)
    )

    # Strip acquisition bookkeeping that is not part of the sealed technical
    # estimand, but retain source SHA for the later biological-pass identity check.
    technical_drop = [
        c
        for c in (
            "acquisition_status",
            "image_bytes",
            "failure_reason",
            "failure_reason_technical",
        )
        if c in merged.columns
    ]
    sealed_table = merged.drop(columns=technical_drop).copy()

    seal = seal_technical_table(sealed_table)

    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    sealed_table.to_csv(
        out / "technical_table.csv.gz",
        index=False,
        compression={"method": "gzip", "compresslevel": 9, "mtime": 0},
        lineterminator="\n",
    )

    strata_cols = [
        "measurement_id",
        "technical_metrics_available",
        "heavy_counterfactual",
        "stratum_high_flower_near_clip",
        "stratum_high_background_near_clip",
        "stratum_low_flower_background_lab_separation",
        "stratum_roi_unstable_fixed",
        "stratum_heavy_prompt_unstable",
    ]
    sealed_table[strata_cols].to_csv(
        out / "technical_strata.csv",
        index=False,
        lineterminator="\n",
    )

    summary = {
        "schema": "fcp_v2_technical_seal_v1",
        "status": "RESPONSE_BLIND_TECHNICAL_SEAL_FROZEN",
        "rows": EXPECTED_ROWS,
        "partition_receipts": EXPECTED_PARTITIONS,
        "acquired_rows": int(acquired.sum()),
        "acquisition_failed_rows": int((~acquired).sum()),
        "technical_metrics_available_rows": int(available.sum()),
        "technical_metrics_unavailable_rows": int((~available).sum()),
        "heavy_counterfactual_rows": int(merged["heavy_counterfactual"].sum()),
        "technical_table_sha256": seal.table_sha256,
        "biological_outcomes_opened": False,
        "species_opened": False,
        "morph_opened": False,
        "q_white_or_W_opened": False,
        "technical_strata": {
            "high_flower_near_clip": {
                "threshold": flower_threshold,
                "q95": flower_q95,
                "rows": int(merged["stratum_high_flower_near_clip"].sum()),
            },
            "high_background_near_clip": {
                "threshold": background_threshold,
                "q95": background_q95,
                "rows": int(
                    merged["stratum_high_background_near_clip"].sum()
                ),
            },
            "low_flower_background_lab_separation": {
                "threshold_q05": separation_q05,
                "rows": int(
                    merged[
                        "stratum_low_flower_background_lab_separation"
                    ].sum()
                ),
            },
            "roi_unstable_fixed": {
                "flip_iou_min": 0.5,
                "flip_delta_e_max": 5.0,
                "rows": int(merged["stratum_roi_unstable_fixed"].sum()),
            },
            "heavy_prompt_unstable": {
                "jitter_iou_min": 0.5,
                "jitter_delta_e_max": 5.0,
                "evaluable_rows": int(heavy_available.sum()),
                "rows": int(merged["stratum_heavy_prompt_unstable"].sum()),
            },
        },
    }
    (out / "technical_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()

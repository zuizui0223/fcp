#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

EV_LEVELS = (-1.0, -0.5, 0.0, 0.5, 1.0)
HEAVY_EXTRA_CONDITIONS = 4 + 1 + 7  # full EV + neutral background + jitter
BASE_CONDITIONS = 1 + len(EV_LEVELS)  # baseline + fixed-mask EV


def _as_bool(value: object) -> bool:
    return str(value).strip().casefold() in {"1", "true", "yes"}


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--expected-worker-manifest", type=Path, required=True)
    p.add_argument("--acquisition-receipt", type=Path, required=True)
    p.add_argument("--biological-long", type=Path, required=True)
    p.add_argument("--output-terminal-receipt", type=Path, required=True)
    args = p.parse_args()

    expected = pd.read_csv(args.expected_worker_manifest, dtype=str).fillna("")
    receipt = pd.read_csv(args.acquisition_receipt, dtype=str).fillna("")
    biological = pd.read_csv(
        args.biological_long,
        dtype={"measurement_id": str},
    )

    required_expected = {
        "measurement_id",
        "heavy_counterfactual",
        "expected_source_sha256",
    }
    if not required_expected.issubset(expected.columns):
        raise RuntimeError("expected worker manifest schema drift")
    if expected["measurement_id"].nunique() != len(expected):
        raise RuntimeError("expected worker IDs are not unique")

    required_receipt = {
        "measurement_id",
        "biological_acquisition_status",
        "expected_source_sha256",
        "observed_source_sha256",
        "heavy_counterfactual",
        "failure_reason",
    }
    if not required_receipt.issubset(receipt.columns):
        raise RuntimeError("Pass-B acquisition receipt schema drift")
    if receipt["measurement_id"].nunique() != len(receipt):
        raise RuntimeError("Pass-B acquisition receipt IDs are not unique")
    if set(receipt["measurement_id"]) != set(expected["measurement_id"]):
        raise RuntimeError("Pass-B acquisition receipt coverage is incomplete")

    if len(biological):
        required_bio = {
            "measurement_id",
            "condition_id",
            "condition_family",
            "morph",
            "measurement_status",
            "counterfactual_status",
            "image_sha256",
        }
        if not required_bio.issubset(biological.columns):
            raise RuntimeError("Pass-B biological output schema drift")
        if biological.duplicated(["measurement_id", "condition_id"]).any():
            raise RuntimeError("duplicate Pass-B biological condition rows")
        if not set(biological["measurement_id"].astype(str)) <= set(
            expected["measurement_id"].astype(str)
        ):
            raise RuntimeError("biological output contains unexpected measurement IDs")

    joined = expected[
        ["measurement_id", "heavy_counterfactual", "expected_source_sha256"]
    ].merge(
        receipt[
            [
                "measurement_id",
                "biological_acquisition_status",
                "expected_source_sha256",
                "observed_source_sha256",
                "failure_reason",
            ]
        ],
        on=["measurement_id", "expected_source_sha256"],
        how="inner",
        validate="one_to_one",
    )
    if len(joined) != len(expected):
        raise RuntimeError("Pass-B terminal merge lost expected rows")

    bio_by_id = (
        biological.groupby("measurement_id", sort=False)
        if len(biological)
        else None
    )
    terminal_rows = []
    for row in joined.itertuples(index=False):
        measurement_id = str(row.measurement_id)
        heavy = _as_bool(row.heavy_counterfactual)
        source_status = str(row.biological_acquisition_status)
        observed_sha = str(row.observed_source_sha256)
        expected_sha = str(row.expected_source_sha256)

        if source_status == "acquired_sha_match":
            if observed_sha != expected_sha:
                raise RuntimeError(
                    "acquired_sha_match receipt contains unequal SHA256 values"
                )
            if bio_by_id is None or measurement_id not in bio_by_id.groups:
                raise RuntimeError(
                    "sha-matched Pass-B row lacks biological condition output"
                )
            group = bio_by_id.get_group(measurement_id)
            expected_conditions = (
                BASE_CONDITIONS
                + (HEAVY_EXTRA_CONDITIONS if heavy else 0)
            )
            if len(group) != expected_conditions:
                raise RuntimeError(
                    f"Pass-B condition coverage drift for {measurement_id}: "
                    f"{len(group)} != {expected_conditions}"
                )
            baseline = group.loc[group["condition_id"].eq("baseline")]
            if len(baseline) != 1:
                raise RuntimeError("Pass-B baseline condition coverage drift")
            baseline_row = baseline.iloc[0]
            terminal_status = "biological_measurement_complete"
            baseline_morph = str(baseline_row["morph"])
            baseline_measurement_status = str(
                baseline_row["measurement_status"]
            )
            condition_rows = int(len(group))
        elif source_status in {
            "image_acquisition_failed",
            "source_byte_drift",
        }:
            if bio_by_id is not None and measurement_id in bio_by_id.groups:
                raise RuntimeError(
                    "failed/drifted Pass-B source reached biological measurement"
                )
            terminal_status = source_status
            baseline_morph = "mixed_uncertain"
            baseline_measurement_status = source_status
            condition_rows = 0
        else:
            raise RuntimeError(
                f"unknown Pass-B source status: {source_status}"
            )

        terminal_rows.append(
            {
                "measurement_id": measurement_id,
                "pass_b_terminal_status": terminal_status,
                "heavy_counterfactual": heavy,
                "expected_source_sha256": expected_sha,
                "observed_source_sha256": observed_sha,
                "source_byte_match": bool(
                    source_status == "acquired_sha_match"
                ),
                "baseline_morph": baseline_morph,
                "baseline_measurement_status": baseline_measurement_status,
                "condition_rows": condition_rows,
                "failure_reason": str(row.failure_reason),
            }
        )

    terminal = pd.DataFrame(terminal_rows)
    if terminal["measurement_id"].nunique() != len(expected):
        raise RuntimeError("Pass-B terminal receipt ID coverage drift")

    args.output_terminal_receipt.parent.mkdir(
        parents=True,
        exist_ok=True,
    )
    terminal.to_csv(
        args.output_terminal_receipt,
        index=False,
        lineterminator="\n",
    )


if __name__ == "__main__":
    main()

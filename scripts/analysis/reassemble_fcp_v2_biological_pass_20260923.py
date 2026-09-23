#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

EXPECTED_ROWS = 40_000
EXPECTED_PARTITIONS = 256
EXPECTED_HEAVY = 8_000
NONHEAVY_CONDITIONS = 6
HEAVY_CONDITIONS = 18


def _bool_series(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False).astype(bool)
    return s.fillna("").astype(str).str.strip().str.casefold().isin(
        {"1", "true", "yes"}
    )


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--partition-results-dir", type=Path, required=True)
    p.add_argument("--sealed-join-key", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args()

    terminal_files = sorted(
        args.partition_results_dir.rglob("terminal_b*_s*_p*.csv")
    )
    biological_files = sorted(
        args.partition_results_dir.rglob("biological_b*_s*_p*.csv")
    )
    acquisition_files = sorted(
        args.partition_results_dir.rglob("acquisition_b*_s*_p*.csv")
    )
    for label, files in (
        ("terminal", terminal_files),
        ("biological", biological_files),
        ("acquisition", acquisition_files),
    ):
        if len(files) != EXPECTED_PARTITIONS:
            raise RuntimeError(
                f"{label} partition census incomplete: {len(files)} != {EXPECTED_PARTITIONS}"
            )

    terminals = pd.concat(
        [pd.read_csv(path, dtype={"measurement_id": str}) for path in terminal_files],
        ignore_index=True,
    )
    if len(terminals) != EXPECTED_ROWS:
        raise RuntimeError(
            f"Pass-B terminal row denominator drift: {len(terminals)}"
        )
    if terminals["measurement_id"].nunique() != EXPECTED_ROWS:
        raise RuntimeError("Pass-B terminal IDs are not unique")

    acquisitions = pd.concat(
        [pd.read_csv(path, dtype={"measurement_id": str}) for path in acquisition_files],
        ignore_index=True,
    )
    if len(acquisitions) != EXPECTED_ROWS:
        raise RuntimeError("Pass-B acquisition row denominator drift")
    if acquisitions["measurement_id"].nunique() != EXPECTED_ROWS:
        raise RuntimeError("Pass-B acquisition IDs are not unique")
    if set(acquisitions["measurement_id"]) != set(terminals["measurement_id"]):
        raise RuntimeError("Pass-B acquisition/terminal ID coverage differs")

    biological_parts = []
    for path in biological_files:
        frame = pd.read_csv(path, dtype={"measurement_id": str})
        if len(frame):
            biological_parts.append(frame)
    biological = (
        pd.concat(biological_parts, ignore_index=True)
        if biological_parts
        else pd.DataFrame()
    )
    if len(biological):
        if biological.duplicated(["measurement_id", "condition_id"]).any():
            raise RuntimeError("duplicate Pass-B biological condition rows")
        if not set(biological["measurement_id"]) <= set(terminals["measurement_id"]):
            raise RuntimeError("Pass-B biological IDs outside terminal denominator")

    join_key = pd.read_csv(
        args.sealed_join_key,
        dtype={"measurement_id": str},
    )
    if len(join_key) != EXPECTED_ROWS:
        raise RuntimeError("Pass-B sealed join-key denominator drift")
    if join_key["measurement_id"].nunique() != EXPECTED_ROWS:
        raise RuntimeError("Pass-B sealed join-key IDs are not unique")
    if set(join_key["measurement_id"]) != set(terminals["measurement_id"]):
        raise RuntimeError("Pass-B sealed join key differs from terminal IDs")
    if join_key["species"].nunique() != 400:
        raise RuntimeError("Pass-B sealed species denominator drift")

    terminal_status = terminals[
        "pass_b_terminal_status"
    ].astype(str).value_counts().to_dict()
    source_match = _bool_series(terminals["source_byte_match"])
    matched_ids = set(
        terminals.loc[source_match, "measurement_id"].astype(str)
    )
    if len(biological) and set(biological["measurement_id"]) != matched_ids:
        missing = len(matched_ids - set(biological["measurement_id"]))
        extra = len(set(biological["measurement_id"]) - matched_ids)
        raise RuntimeError(
            f"Pass-B biological source-match coverage drift: missing={missing} extra={extra}"
        )

    # Validate long-form condition coverage against the metadata-frozen heavy flag.
    heavy = _bool_series(terminals["heavy_counterfactual"])
    expected_condition_rows = (
        int((source_match & ~heavy).sum()) * NONHEAVY_CONDITIONS
        + int((source_match & heavy).sum()) * HEAVY_CONDITIONS
    )
    if len(biological) != expected_condition_rows:
        raise RuntimeError(
            f"Pass-B biological long denominator drift: "
            f"{len(biological)} != {expected_condition_rows}"
        )

    if len(biological):
        counts = biological.groupby("measurement_id", sort=False).size()
        expected_by_id = terminals.set_index("measurement_id")
        for mid, n in counts.items():
            row = expected_by_id.loc[mid]
            expected_n = (
                HEAVY_CONDITIONS
                if _bool_series(pd.Series([row["heavy_counterfactual"]])).iloc[0]
                else NONHEAVY_CONDITIONS
            )
            if int(n) != expected_n:
                raise RuntimeError(
                    f"Pass-B per-ID condition count drift for {mid}: {n} != {expected_n}"
                )
        baseline = biological.loc[
            biological["condition_id"].astype(str).eq("baseline")
        ]
        if baseline["measurement_id"].nunique() != int(source_match.sum()):
            raise RuntimeError("Pass-B baseline coverage differs from SHA-matched rows")

    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    terminals.to_csv(
        out / "terminal_receipt_40000.csv.gz",
        index=False,
        compression={"method": "gzip", "compresslevel": 9, "mtime": 0},
        lineterminator="\n",
    )
    acquisitions.to_csv(
        out / "acquisition_receipt_40000.csv.gz",
        index=False,
        compression={"method": "gzip", "compresslevel": 9, "mtime": 0},
        lineterminator="\n",
    )
    biological.to_csv(
        out / "biological_conditions_long.csv.gz",
        index=False,
        compression={"method": "gzip", "compresslevel": 9, "mtime": 0},
        lineterminator="\n",
    )

    summary = {
        "schema": "fcp_v2_biological_reassembly_v1",
        "status": "COMPLETE_PASS_B_TERMINAL_REASSEMBLY",
        "partition_receipts": EXPECTED_PARTITIONS,
        "terminal_rows": EXPECTED_ROWS,
        "species_sealed_until_complete_reassembly": 400,
        "heavy_counterfactual_rows_frozen": EXPECTED_HEAVY,
        "terminal_status_counts": {
            str(k): int(v) for k, v in terminal_status.items()
        },
        "source_byte_matched_rows": int(source_match.sum()),
        "source_byte_drift_rows": int(
            terminals["pass_b_terminal_status"].astype(str).eq(
                "source_byte_drift"
            ).sum()
        ),
        "image_acquisition_failed_rows": int(
            terminals["pass_b_terminal_status"].astype(str).eq(
                "image_acquisition_failed"
            ).sum()
        ),
        "biological_condition_rows": int(len(biological)),
        "baseline_rows": int(
            biological["condition_id"].astype(str).eq("baseline").sum()
            if len(biological)
            else 0
        ),
        "species_join_opened_after_complete_reassembly": True,
        "new_null_generated": False,
        "replacement_rows": 0,
    }
    (out / "reassembly_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

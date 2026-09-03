#!/usr/bin/env python3
"""Combine acquisition failures and blind ROI measurements into one terminal partition."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from fcp_pipeline.photo_first_measurement_execution import (
    BIOLOGICAL_MORPHS,
    validate_terminal_partition_results,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--selected-blind-manifest", type=Path, required=True)
    parser.add_argument("--acquisition-failures", type=Path, required=True)
    parser.add_argument("--measurement-results", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--semantic-shard", type=int, required=True)
    parser.add_argument("--compute-partition", type=int, required=True)
    return parser.parse_args()


def read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype={"measurement_id": str}).fillna("")


def main() -> int:
    args = parse_args()
    selected = read_csv(args.selected_blind_manifest)
    failures = read_csv(args.acquisition_failures)
    measured = read_csv(args.measurement_results)
    terminal = pd.concat([failures, measured], ignore_index=True, sort=False)
    expected = selected["measurement_id"].astype(str).tolist()
    validate_terminal_partition_results(terminal, expected)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    stem = f"partition_s{args.semantic_shard:02d}_p{args.compute_partition:02d}"
    result_path = args.output_dir / f"{stem}.csv"
    receipt_path = args.output_dir / f"{stem}.json"
    terminal = terminal.sort_values("measurement_id").reset_index(drop=True)
    terminal.to_csv(result_path, index=False, lineterminator="\n")
    receipt = {
        "status": "complete_random_photo_first_terminal_partition",
        "semantic_shard": int(args.semantic_shard),
        "compute_partition": int(args.compute_partition),
        "selected_rows": int(len(selected)),
        "terminal_rows": int(len(terminal)),
        "acquisition_failed_rows": int(
            terminal["measurement_status"].eq("image_acquisition_failed").sum()
        ),
        "classified_rows": int(terminal["morph"].isin(BIOLOGICAL_MORPHS).sum()),
        "mixed_uncertain_rows": int(terminal["morph"].eq("mixed_uncertain").sum()),
        "source_urls_present": False,
        "species_present": False,
        "coordinates_present": False,
    }
    receipt_path.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

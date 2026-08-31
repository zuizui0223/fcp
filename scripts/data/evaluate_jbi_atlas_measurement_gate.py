#!/usr/bin/env python3
"""Evaluate the frozen scale-out image-completeness gate without coordinates."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.atlas_measurement import (
    evaluate_scaleout_measurement_gate,
    validate_inference_contract,
    validate_measurement_result_rows,
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def parse_bool(value: str) -> bool:
    folded = value.strip().casefold()
    if folded == "true":
        return True
    if folded == "false":
        return False
    raise ValueError(f"background_features_available must be true/false, got {value!r}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--measurement-results", type=Path)
    inputs.add_argument("--measurement-results-dir", type=Path)
    parser.add_argument("--sealed-species-key", type=Path, required=True)
    parser.add_argument(
        "--inference-contract",
        type=Path,
        default=Path("docs/supporting/jbi_image_first_atlas_inference_contract_v3.json"),
    )
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    contract = json.loads(args.inference_contract.read_text(encoding="utf-8"))
    validate_inference_contract(contract)
    result_paths = (
        [args.measurement_results]
        if args.measurement_results is not None
        else sorted(args.measurement_results_dir.glob("measurement_shard_*.csv"))
    )
    if not result_paths:
        raise RuntimeError("no measurement result tables were found")
    results = [row for path in result_paths for row in read_csv(path)]
    for row in results:
        row["background_features_available"] = parse_bool(
            row["background_features_available"]
        )
    validate_measurement_result_rows(results)
    decision = evaluate_scaleout_measurement_gate(
        results,
        read_csv(args.sealed_species_key),
        contract["scaleout_measurement_gate"],
    )
    decision["protocol"] = contract["protocol"]
    decision["claim_ceiling"] = (
        "Measurement completeness only; coordinates remain unopened and no spatial, "
        "environmental or pollinator conclusion is allowed."
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(decision, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(decision, indent=2, sort_keys=True))
    if not decision["coordinate_join_permitted"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()

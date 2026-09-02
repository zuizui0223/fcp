#!/usr/bin/env python3
"""Evaluate terminal atlas v5 measurement completeness before opening coordinates.

A scientifically valid ``not_evaluable`` completeness outcome is a completed
analysis, not a technical execution failure.  This command therefore exits
zero after writing any valid gate decision.  Coordinate access remains governed
only by ``coordinate_join_permitted`` in the persisted decision JSON.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.atlas_measurement import (
    evaluate_scaleout_measurement_gate,
    validate_measurement_result_rows,
)
from fcp_pipeline.atlas_measurement_v5 import validate_measurement_execution_contract
from scripts.data.build_jbi_atlas_measurement_firewall_v5 import verify_repo_parent_blobs


DEFAULT_CONTRACT = ROOT / "docs/supporting/jbi_atlas_measurement_execution_contract_v5.json"
DEFAULT_INFERENCE = ROOT / "docs/supporting/jbi_image_first_atlas_inference_contract_v5.json"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    folded = str(value).strip().casefold()
    if folded == "true":
        return True
    if folded == "false":
        return False
    raise ValueError(f"expected boolean true/false, got {value!r}")


def load_complete_v5_measurement_bundle(
    directory: Path,
    *,
    contract: dict[str, Any],
    inference: dict[str, Any],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    paths = sorted(directory.glob("measurement_v5_shard_*.csv"))
    expected_shards = int(contract["technical_execution"]["measurement_shard_count"])
    if len(paths) != expected_shards:
        raise RuntimeError(
            f"v5 measurement bundle requires {expected_shards} shard tables, found {len(paths)}"
        )
    rows: list[dict[str, Any]] = []
    files: dict[str, str] = {}
    indexes: set[int] = set()
    model_ids: set[str] = set()
    weights: set[str] = set()
    roi_contracts: set[str] = set()
    for path in paths:
        try:
            index = int(path.stem.rsplit("_", 1)[1])
        except (IndexError, ValueError) as exc:
            raise RuntimeError(f"invalid v5 measurement shard filename: {path.name}") from exc
        manifest_path = path.with_suffix(".json")
        if not manifest_path.is_file():
            raise RuntimeError(f"missing v5 measurement shard manifest: {manifest_path.name}")
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        shard_rows = read_csv(path)
        if (
            manifest.get("status") != "complete_location_blind_roi_v4_measurement_v5_shard"
            or manifest.get("protocol") != contract["protocol"]
            or manifest.get("inference_version") != inference["version"]
            or manifest.get("superseded_v3_ordered_inference_used") is not False
            or manifest.get("shard_index") != index
            or manifest.get("shard_count") != expected_shards
            or manifest.get("coordinates_opened") is not False
            or manifest.get("taxon_names_opened") is not False
            or manifest.get("frozen_shard_denominator") != len(shard_rows)
            or manifest.get("terminal_records") != len(shard_rows)
            or manifest.get("result_sha256") != sha256(path)
        ):
            raise RuntimeError(f"v5 measurement shard evidence changed: {path.name}")
        if index in indexes:
            raise RuntimeError("v5 measurement shard index is duplicated")
        indexes.add(index)
        model_ids.add(str(manifest.get("model_id") or ""))
        weights.add(str(manifest.get("trained_weight_sha256") or ""))
        roi_contracts.add(str(manifest.get("roi_contract_sha256_lf_canonical_v1") or ""))
        for row in shard_rows:
            row["background_features_available"] = parse_bool(
                row["background_features_available"]
            )
        validate_measurement_result_rows(shard_rows)
        rows.extend(shard_rows)
        files[path.name] = sha256(path)
        files[manifest_path.name] = sha256(manifest_path)
    if indexes != set(range(expected_shards)):
        raise RuntimeError("v5 measurement shard set is incomplete")
    if any(
        len(values) != 1 or not next(iter(values))
        for values in (model_ids, weights, roi_contracts)
    ):
        raise RuntimeError("v5 measurement shards disagree on estimator identity")
    if len(rows) != 60000 or len({str(row["measurement_id"]) for row in rows}) != 60000:
        raise RuntimeError("v5 measurement bundle does not preserve the 60,000-record denominator")
    return rows, {
        "status": "pass_complete_location_blind_roi_v4_measurement_v5_bundle",
        "shard_count": expected_shards,
        "terminal_records": len(rows),
        "model_id": next(iter(model_ids)),
        "trained_weight_sha256": next(iter(weights)),
        "roi_contract_sha256_lf_canonical_v1": next(iter(roi_contracts)),
        "files": dict(sorted(files.items())),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--measurement-results-dir", type=Path, required=True)
    parser.add_argument("--sealed-species-key", type=Path, required=True)
    parser.add_argument("--measurement-contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--inference-v5", type=Path, default=DEFAULT_INFERENCE)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    contract = json.loads(args.measurement_contract.read_text(encoding="utf-8"))
    inference = json.loads(args.inference_v5.read_text(encoding="utf-8"))
    validate_measurement_execution_contract(contract, inference)
    verify_repo_parent_blobs(contract)
    results, bundle = load_complete_v5_measurement_bundle(
        args.measurement_results_dir,
        contract=contract,
        inference=inference,
    )
    species_key = read_csv(args.sealed_species_key)
    if len(species_key) != 60000:
        raise RuntimeError("sealed species key no longer has 60,000 rows")
    decision = evaluate_scaleout_measurement_gate(
        results,
        species_key,
        contract["location_blind_measurement"],
    )
    decision["protocol"] = contract["protocol"]
    decision["inference_version"] = inference["version"]
    decision["superseded_v3_ordered_inference_used"] = False
    decision["measurement_bundle"] = bundle
    decision["source_sha256"] = {
        "sealed_species_key": sha256(args.sealed_species_key),
        "measurement_contract": sha256(args.measurement_contract),
        "inference_v5": sha256(args.inference_v5),
    }
    decision["execution_status"] = "pass_measurement_gate_evaluation_completed"
    decision["scientific_outcome_is_process_failure"] = False
    decision["claim_ceiling"] = (
        "Measurement completeness classification only. A zero process exit means the "
        "frozen 60,000-record gate was evaluated successfully; it does not imply that "
        "coordinate_join_permitted is true and does not imply biological support."
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(decision, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(decision, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Evaluate the frozen scale-out image-completeness gate without coordinates."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.atlas_measurement import (
    EXECUTION_PROTOCOL,
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


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_complete_measurement_bundle(
    directory: Path,
) -> tuple[list[dict[str, str]], dict[str, object]]:
    """Load only a complete, self-hashed set of location-blind ROI-v4 shards."""

    paths = sorted(directory.glob("measurement_shard_*.csv"))
    if not paths:
        raise RuntimeError("no measurement result tables were found")
    rows: list[dict[str, str]] = []
    files: dict[str, str] = {}
    indexes: set[int] = set()
    shard_counts: set[int] = set()
    model_ids: set[str] = set()
    trained_weights: set[str] = set()
    roi_contracts: set[str] = set()
    for path in paths:
        try:
            index = int(path.stem.rsplit("_", 1)[1])
        except (IndexError, ValueError) as exc:
            raise RuntimeError(f"invalid measurement shard filename: {path.name}") from exc
        manifest_path = path.with_suffix(".json")
        if not manifest_path.is_file():
            raise RuntimeError(f"measurement shard manifest is missing: {manifest_path.name}")
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        shard_rows = read_csv(path)
        if (
            manifest.get("status")
            != "complete_location_blind_roi_v4_measurement_shard"
            or manifest.get("execution_protocol") != EXECUTION_PROTOCOL
            or manifest.get("shard_index") != index
            or manifest.get("coordinates_opened") is not False
            or manifest.get("taxon_names_opened") is not False
            or manifest.get("frozen_shard_denominator") != len(shard_rows)
            or manifest.get("terminal_records") != len(shard_rows)
            or manifest.get("result_sha256") != sha256(path)
        ):
            raise RuntimeError(f"measurement shard evidence changed: {path.name}")
        if index in indexes:
            raise RuntimeError("measurement shard index is duplicated")
        indexes.add(index)
        shard_counts.add(int(manifest.get("shard_count", 0)))
        model_ids.add(str(manifest.get("model_id") or ""))
        trained_weights.add(str(manifest.get("trained_weight_sha256") or ""))
        roi_contracts.add(str(manifest.get("roi_contract_sha256_lf_canonical_v1") or ""))
        files[path.name] = sha256(path)
        files[manifest_path.name] = sha256(manifest_path)
        rows.extend(shard_rows)
    if len(shard_counts) != 1:
        raise RuntimeError("measurement shards disagree on shard_count")
    shard_count = shard_counts.pop()
    if shard_count < 1 or indexes != set(range(shard_count)):
        raise RuntimeError("measurement shard set is incomplete")
    if any(len(values) != 1 or not next(iter(values)) for values in (
        model_ids,
        trained_weights,
        roi_contracts,
    )):
        raise RuntimeError("measurement shards disagree on estimator identity")
    return rows, {
        "status": "pass_complete_location_blind_roi_v4_measurement_bundle",
        "shard_count": shard_count,
        "model_id": next(iter(model_ids)),
        "trained_weight_sha256": next(iter(trained_weights)),
        "roi_contract_sha256_lf_canonical_v1": next(iter(roi_contracts)),
        "files": dict(sorted(files.items())),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--measurement-results-dir", type=Path, required=True)
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
    results, bundle = load_complete_measurement_bundle(args.measurement_results_dir)
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
    decision["measurement_bundle"] = bundle
    decision["source_sha256"] = {
        "sealed_species_key": sha256(args.sealed_species_key),
        "inference_contract": sha256(args.inference_contract),
    }
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

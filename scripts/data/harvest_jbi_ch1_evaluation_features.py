#!/usr/bin/env python3
"""Harvest and validate the frozen Chapter 1 evaluation feature shards.

This script is intentionally stricter than a simple concatenation.  It requires the
artifact rows to be exactly the predeclared 720 evaluation photograph IDs in the frozen
split, with 120 rows for each of six species and zero calibration overlap.  It never
selects features or observations from their spatial outcome.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable

import pandas as pd


EXPECTED_SPECIES = 6
EXPECTED_PER_SPECIES = 120
EXPECTED_TOTAL = EXPECTED_SPECIES * EXPECTED_PER_SPECIES
EXPECTED_SHARDS = EXPECTED_SPECIES * 6
EXPECTED_ROWS_PER_SHARD = 20

REQUIRED_TRUE = (
    "evaluation_row",
    "evaluation_feature_measurement",
)
REQUIRED_FALSE = (
    "calibration_only",
    "final_label",
)


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def canonical_row_hash(rows: Iterable[dict[str, Any]]) -> str:
    payload = "\n".join(
        json.dumps(row, sort_keys=True, ensure_ascii=False, separators=(",", ":"))
        for row in sorted(rows, key=lambda x: (str(x["species"]), str(x["photo_id"])))
    ).encode("utf-8")
    return sha256_bytes(payload)


def read_jsonl_files(root: Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    paths = sorted(root.rglob("*.jsonl"))
    if len(paths) != EXPECTED_SHARDS:
        raise ValueError(f"expected {EXPECTED_SHARDS} JSONL shards, found {len(paths)}")

    rows: list[dict[str, Any]] = []
    files: list[dict[str, Any]] = []
    for path in paths:
        raw = path.read_bytes()
        parsed = [
            json.loads(line)
            for line in raw.decode("utf-8-sig").splitlines()
            if line.strip()
        ]
        if len(parsed) != EXPECTED_ROWS_PER_SHARD:
            raise ValueError(
                f"{path}: expected {EXPECTED_ROWS_PER_SHARD} rows, found {len(parsed)}"
            )
        files.append(
            {
                "path": str(path.relative_to(root)),
                "rows": len(parsed),
                "sha256": sha256_bytes(raw),
            }
        )
        rows.extend(parsed)
    return rows, files


def validate_against_split(
    rows: list[dict[str, Any]],
    split: pd.DataFrame,
) -> dict[str, Any]:
    required_split = {"species", "photo_id", "split", "latitude", "longitude"}
    if not required_split.issubset(split.columns):
        raise ValueError(
            f"frozen split missing fields: {sorted(required_split - set(split.columns))}"
        )

    if len(rows) != EXPECTED_TOTAL:
        raise ValueError(f"expected {EXPECTED_TOTAL} artifact rows, found {len(rows)}")

    required_row = {"species", "photo_id", "blind_id", "feature_status"}
    for index, row in enumerate(rows):
        missing = required_row - set(row)
        if missing:
            raise ValueError(f"row {index} missing required fields: {sorted(missing)}")
        for field in REQUIRED_TRUE:
            if row.get(field) is not True:
                raise ValueError(f"row {index}: {field} must be true")
        for field in REQUIRED_FALSE:
            if row.get(field) is not False:
                raise ValueError(f"row {index}: {field} must be false")

    artifact_keys = [(str(r["species"]), str(r["photo_id"])) for r in rows]
    if len(set(artifact_keys)) != EXPECTED_TOTAL:
        duplicates = [k for k, n in Counter(artifact_keys).items() if n > 1]
        raise ValueError(f"evaluation artifacts contain duplicate species/photo IDs: {duplicates[:5]}")
    blind_ids = [str(r["blind_id"]) for r in rows]
    if len(set(blind_ids)) != EXPECTED_TOTAL:
        raise ValueError("evaluation artifacts contain duplicate blind IDs")

    split = split.copy()
    split["species"] = split["species"].astype(str)
    split["photo_id"] = split["photo_id"].astype(str)
    expected_eval = split.loc[split["split"].astype(str) == "evaluation"]
    calibration = split.loc[split["split"].astype(str) == "calibration"]
    if len(expected_eval) != EXPECTED_TOTAL:
        raise ValueError(f"frozen split has {len(expected_eval)} evaluation rows, expected 720")

    expected_keys = set(zip(expected_eval["species"], expected_eval["photo_id"], strict=True))
    observed_keys = set(artifact_keys)
    if observed_keys != expected_keys:
        missing = sorted(expected_keys - observed_keys)[:10]
        extra = sorted(observed_keys - expected_keys)[:10]
        raise ValueError(f"artifact IDs differ from frozen evaluation set; missing={missing}, extra={extra}")

    calibration_keys = set(zip(calibration["species"], calibration["photo_id"], strict=True))
    overlap = observed_keys & calibration_keys
    if overlap:
        raise ValueError(f"calibration ID leakage detected: {sorted(overlap)[:10]}")

    species_counts = Counter(str(r["species"]) for r in rows)
    if len(species_counts) != EXPECTED_SPECIES or set(species_counts.values()) != {EXPECTED_PER_SPECIES}:
        raise ValueError(f"evaluation rows must be 6 species x 120: {dict(species_counts)}")

    shard_pairs: Counter[tuple[str, int]] = Counter()
    for row in rows:
        if int(row.get("compute_shard_count", -1)) != 6:
            raise ValueError("every row must record compute_shard_count=6")
        shard_index = int(row.get("compute_shard_index", -1))
        if shard_index not in range(6):
            raise ValueError(f"invalid compute_shard_index={shard_index}")
        shard_pairs[(str(row["species"]), shard_index)] += 1
    if len(shard_pairs) != EXPECTED_SHARDS or set(shard_pairs.values()) != {EXPECTED_ROWS_PER_SHARD}:
        raise ValueError(f"expected 36 species/shard cells x 20 rows: {dict(shard_pairs)}")

    feature_status = Counter(str(r.get("feature_status", "missing")) for r in rows)
    methods = Counter(str(r.get("feature_method", r.get("localization_method", "missing"))) for r in rows)
    return {
        "total_rows": len(rows),
        "species_counts": dict(sorted(species_counts.items())),
        "shard_cells": len(shard_pairs),
        "rows_per_shard_values": sorted(set(shard_pairs.values())),
        "feature_status_counts": dict(sorted(feature_status.items())),
        "feature_method_counts": dict(sorted(methods.items())),
        "frozen_evaluation_id_match": True,
        "calibration_overlap_count": 0,
    }


def join_coordinates(rows: list[dict[str, Any]], split: pd.DataFrame) -> list[dict[str, Any]]:
    lookup = {
        (str(row.species), str(row.photo_id)): {
            "latitude": float(row.latitude),
            "longitude": float(row.longitude),
        }
        for row in split.itertuples(index=False)
        if str(row.split) == "evaluation"
    }
    joined: list[dict[str, Any]] = []
    for row in rows:
        key = (str(row["species"]), str(row["photo_id"]))
        if key not in lookup:
            raise ValueError(f"no frozen coordinate row for {key}")
        item = dict(row)
        # Coordinates come only from the frozen acquisition manifest, never from model output.
        item.update(lookup[key])
        joined.append(item)
    return sorted(joined, key=lambda x: (str(x["species"]), str(x["photo_id"])))


def flatten_for_csv(row: dict[str, Any]) -> dict[str, Any]:
    flat: dict[str, Any] = {}

    def visit(prefix: str, value: Any) -> None:
        if isinstance(value, dict):
            for key in sorted(value):
                visit(f"{prefix}.{key}" if prefix else str(key), value[key])
        elif isinstance(value, (list, tuple)):
            flat[prefix] = json.dumps(value, ensure_ascii=False, sort_keys=True)
        else:
            flat[prefix] = value

    visit("", row)
    return flat


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("artifact_root", type=Path)
    parser.add_argument(
        "--split",
        type=Path,
        default=Path("data/frozen/jbi_ch1_photo_split_v1.csv"),
    )
    parser.add_argument(
        "--output-jsonl",
        type=Path,
        default=Path("data/evaluation/jbi_ch1_evaluation_features_v1.jsonl"),
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=Path("data/evaluation/jbi_ch1_evaluation_features_v1.csv"),
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("docs/supporting/jbi_ch1_evaluation_feature_harvest_v1.json"),
    )
    args = parser.parse_args()

    rows, files = read_jsonl_files(args.artifact_root)
    split = pd.read_csv(args.split, dtype={"photo_id": str})
    validation = validate_against_split(rows, split)
    joined = join_coordinates(rows, split)

    args.output_jsonl.parent.mkdir(parents=True, exist_ok=True)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.parent.mkdir(parents=True, exist_ok=True)

    jsonl_text = "\n".join(
        json.dumps(row, sort_keys=True, ensure_ascii=False, separators=(",", ":"))
        for row in joined
    ) + "\n"
    args.output_jsonl.write_text(jsonl_text, encoding="utf-8")
    pd.DataFrame([flatten_for_csv(row) for row in joined]).to_csv(
        args.output_csv, index=False, lineterminator="\n"
    )

    manifest = {
        "protocol": "jbi-ch1-evaluation-feature-harvest-v1",
        "status": "pass",
        "source_workflow_run_id": 33281907575,
        "source_workflow_head_sha": "23260eb8a39c692c98cfc7717f10019978e9dd18",
        "split_path": str(args.split),
        "split_sha256": sha256_bytes(args.split.read_bytes()),
        "artifact_files": files,
        "artifact_file_count": len(files),
        "validation": validation,
        "canonical_row_sha256": canonical_row_hash(joined),
        "output_jsonl_sha256": sha256_bytes(args.output_jsonl.read_bytes()),
        "output_csv_sha256": sha256_bytes(args.output_csv.read_bytes()),
        "measurement_rules_changed_after_evaluation_opening": False,
        "spatial_outcome_used_for_row_or_feature_selection": False,
    }
    args.manifest.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

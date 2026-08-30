#!/usr/bin/env python3
"""Aggregate the frozen Chapter 1 evaluation shards without changing measurement rules.

The script is deliberately outcome-blind.  It concatenates the 36 Florence evaluation
shards, restores frozen spatial metadata by photo_id, and enforces the 6 x 120 design.
No colour threshold, class label, boundary call, or feature-tuning decision is made here.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

EXPECTED_TOTAL = 720
EXPECTED_PER_SPECIES = 120
EXPECTED_SHARDS_PER_SPECIES = 6
EXPECTED_ROWS_PER_SHARD = 20


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def scalar_id(value: Any) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if text.endswith(".0"):
        try:
            return str(int(float(text)))
        except ValueError:
            pass
    return text


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                row = json.loads(line)
            except json.JSONDecodeError as exc:
                raise ValueError(f"{path}:{line_number}: invalid JSON: {exc}") from exc
            if not isinstance(row, dict):
                raise ValueError(f"{path}:{line_number}: each JSONL row must be an object")
            row = dict(row)
            row["_source_shard_file"] = path.name
            rows.append(row)
    return rows


def load_split_metadata(path: Path) -> tuple[dict[str, dict[str, str]], list[str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"{path}: missing CSV header")
        required = {"species", "photo_id", "split", "latitude", "longitude"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"{path}: missing required columns {sorted(missing)}")
        metadata: dict[str, dict[str, str]] = {}
        species: set[str] = set()
        for row in reader:
            if str(row.get("split", "")).strip().lower() != "evaluation":
                continue
            photo_id = scalar_id(row.get("photo_id"))
            if not photo_id:
                raise ValueError(f"{path}: evaluation row with empty photo_id")
            if photo_id in metadata:
                raise ValueError(f"{path}: duplicate evaluation photo_id {photo_id}")
            metadata[photo_id] = dict(row)
            species.add(str(row["species"]).strip())
    if len(metadata) != EXPECTED_TOTAL:
        raise ValueError(
            f"{path}: expected {EXPECTED_TOTAL} frozen evaluation rows, found {len(metadata)}"
        )
    return metadata, sorted(species)


def flatten_schema(value: Any, prefix: str = "") -> Iterable[tuple[str, str]]:
    if isinstance(value, dict):
        if not value:
            yield prefix, "object(empty)"
        for key in sorted(value):
            child = f"{prefix}.{key}" if prefix else str(key)
            yield from flatten_schema(value[key], child)
    elif isinstance(value, list):
        yield prefix, "array"
        for item in value[:5]:
            yield from flatten_schema(item, f"{prefix}[]")
    elif value is None:
        yield prefix, "null"
    elif isinstance(value, bool):
        yield prefix, "boolean"
    elif isinstance(value, (int, float)) and not isinstance(value, bool):
        yield prefix, "number"
    else:
        yield prefix, "string"


def finite_coordinate(value: Any, field: str, photo_id: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"photo_id {photo_id}: invalid {field}={value!r}") from exc
    if not math.isfinite(result):
        raise ValueError(f"photo_id {photo_id}: non-finite {field}")
    return result


def fail(message: str, diagnostics: dict[str, Any], diagnostics_path: Path) -> None:
    diagnostics["validation_status"] = "failed"
    diagnostics["validation_error"] = message
    diagnostics_path.parent.mkdir(parents=True, exist_ok=True)
    diagnostics_path.write_text(
        json.dumps(diagnostics, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    raise SystemExit(message)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--split-csv", type=Path, required=True)
    parser.add_argument("--output-jsonl", type=Path, required=True)
    parser.add_argument("--manifest-json", type=Path, required=True)
    parser.add_argument("--schema-json", type=Path, required=True)
    parser.add_argument("--source-run-id", required=True)
    parser.add_argument("--source-head-sha", required=True)
    args = parser.parse_args()

    diagnostics: dict[str, Any] = {
        "protocol": "jbi_ch1_frozen_evaluation_aggregate_v1",
        "source_run_id": str(args.source_run_id),
        "source_head_sha": str(args.source_head_sha),
        "expected_total_rows": EXPECTED_TOTAL,
        "expected_rows_per_species": EXPECTED_PER_SPECIES,
        "expected_shards_per_species": EXPECTED_SHARDS_PER_SPECIES,
        "expected_rows_per_shard": EXPECTED_ROWS_PER_SHARD,
        "outcome_blind": True,
        "changes_measurement_rule": False,
    }

    try:
        split_by_photo, expected_species = load_split_metadata(args.split_csv)
    except Exception as exc:
        fail(str(exc), diagnostics, args.manifest_json)

    shard_paths = sorted(args.input_dir.rglob("*.jsonl"))
    diagnostics["shard_files"] = [str(path.relative_to(args.input_dir)) for path in shard_paths]
    diagnostics["shard_file_count"] = len(shard_paths)
    if len(shard_paths) != len(expected_species) * EXPECTED_SHARDS_PER_SPECIES:
        fail(
            "expected "
            f"{len(expected_species) * EXPECTED_SHARDS_PER_SPECIES} JSONL shard files, "
            f"found {len(shard_paths)}",
            diagnostics,
            args.manifest_json,
        )

    rows: list[dict[str, Any]] = []
    shard_row_counts: dict[str, int] = {}
    for shard_path in shard_paths:
        try:
            shard_rows = read_jsonl(shard_path)
        except Exception as exc:
            fail(str(exc), diagnostics, args.manifest_json)
        shard_row_counts[str(shard_path.relative_to(args.input_dir))] = len(shard_rows)
        if len(shard_rows) != EXPECTED_ROWS_PER_SHARD:
            fail(
                f"{shard_path}: expected {EXPECTED_ROWS_PER_SHARD} rows, found {len(shard_rows)}",
                diagnostics,
                args.manifest_json,
            )
        rows.extend(shard_rows)
    diagnostics["shard_row_counts"] = shard_row_counts

    if len(rows) != EXPECTED_TOTAL:
        fail(
            f"expected {EXPECTED_TOTAL} aggregate rows, found {len(rows)}",
            diagnostics,
            args.manifest_json,
        )

    seen_photo: set[str] = set()
    seen_blind: set[str] = set()
    enriched: list[dict[str, Any]] = []
    failures: list[dict[str, str]] = []
    by_species: Counter[str] = Counter()
    by_species_shard: defaultdict[str, Counter[int]] = defaultdict(Counter)

    for row in rows:
        photo_id = scalar_id(row.get("photo_id"))
        blind_id = scalar_id(row.get("blind_id"))
        species = str(row.get("species", "")).strip()
        if not photo_id or not blind_id or not species:
            fail(
                "every shard row must contain non-empty species, photo_id, and blind_id",
                diagnostics,
                args.manifest_json,
            )
        if photo_id in seen_photo:
            fail(f"duplicate photo_id {photo_id}", diagnostics, args.manifest_json)
        if blind_id in seen_blind:
            fail(f"duplicate blind_id {blind_id}", diagnostics, args.manifest_json)
        seen_photo.add(photo_id)
        seen_blind.add(blind_id)

        frozen = split_by_photo.get(photo_id)
        if frozen is None:
            fail(
                f"photo_id {photo_id} is absent from the frozen evaluation split",
                diagnostics,
                args.manifest_json,
            )
        frozen_species = str(frozen["species"]).strip()
        if species != frozen_species:
            fail(
                f"photo_id {photo_id}: shard species {species!r} != frozen species {frozen_species!r}",
                diagnostics,
                args.manifest_json,
            )
        if row.get("evaluation_row") is not True:
            fail(f"photo_id {photo_id}: evaluation_row is not true", diagnostics, args.manifest_json)
        if row.get("calibration_only") is not False:
            fail(f"photo_id {photo_id}: calibration_only is not false", diagnostics, args.manifest_json)
        if row.get("final_label") is not False:
            fail(f"photo_id {photo_id}: final_label is not false", diagnostics, args.manifest_json)
        if row.get("evaluation_feature_measurement") is not True:
            fail(
                f"photo_id {photo_id}: evaluation_feature_measurement is not true",
                diagnostics,
                args.manifest_json,
            )

        status = str(row.get("feature_status", "")).strip().lower()
        if status != "ok":
            failures.append(
                {
                    "species": species,
                    "photo_id": photo_id,
                    "blind_id": blind_id,
                    "feature_status": status or "missing",
                    "feature_error": str(row.get("feature_error", "")),
                }
            )

        shard_index = row.get("compute_shard_index")
        shard_count = row.get("compute_shard_count")
        try:
            shard_index_i = int(shard_index)
            shard_count_i = int(shard_count)
        except (TypeError, ValueError) as exc:
            fail(
                f"photo_id {photo_id}: invalid compute shard metadata",
                diagnostics,
                args.manifest_json,
            )
        if shard_count_i != EXPECTED_SHARDS_PER_SPECIES or not 0 <= shard_index_i < shard_count_i:
            fail(
                f"photo_id {photo_id}: invalid shard {shard_index_i}/{shard_count_i}",
                diagnostics,
                args.manifest_json,
            )

        output_row = dict(row)
        output_row["photo_id"] = photo_id
        output_row["blind_id"] = blind_id
        output_row["latitude"] = finite_coordinate(frozen["latitude"], "latitude", photo_id)
        output_row["longitude"] = finite_coordinate(frozen["longitude"], "longitude", photo_id)
        for key in (
            "observation_id",
            "photo_url",
            "photo_license",
            "attribution",
            "observed_on",
            "observed_month",
            "spatial_cell",
            "selection_hash",
            "split_rank_hash",
        ):
            if key in frozen and key not in output_row:
                output_row[key] = frozen[key]
        output_row["split"] = "evaluation"
        enriched.append(output_row)
        by_species[species] += 1
        by_species_shard[species][shard_index_i] += 1

    missing_frozen = sorted(set(split_by_photo) - seen_photo)
    if missing_frozen:
        fail(
            f"{len(missing_frozen)} frozen evaluation photo_ids are missing from artifacts",
            diagnostics,
            args.manifest_json,
        )
    if sorted(by_species) != expected_species:
        fail(
            f"aggregate species {sorted(by_species)} != frozen species {expected_species}",
            diagnostics,
            args.manifest_json,
        )
    for species in expected_species:
        if by_species[species] != EXPECTED_PER_SPECIES:
            fail(
                f"{species}: expected {EXPECTED_PER_SPECIES} rows, found {by_species[species]}",
                diagnostics,
                args.manifest_json,
            )
        expected_shards = {index: EXPECTED_ROWS_PER_SHARD for index in range(EXPECTED_SHARDS_PER_SPECIES)}
        observed_shards = dict(sorted(by_species_shard[species].items()))
        if observed_shards != expected_shards:
            fail(
                f"{species}: shard counts {observed_shards} != {expected_shards}",
                diagnostics,
                args.manifest_json,
            )

    diagnostics["feature_failures"] = failures
    diagnostics["feature_failure_count"] = len(failures)
    if failures:
        fail(
            f"{len(failures)} evaluation photographs have non-ok feature status",
            diagnostics,
            args.manifest_json,
        )

    enriched.sort(key=lambda row: (str(row["species"]), scalar_id(row["photo_id"])))
    args.output_jsonl.parent.mkdir(parents=True, exist_ok=True)
    output_text = "".join(
        json.dumps(row, sort_keys=True, ensure_ascii=False, separators=(",", ":")) + "\n"
        for row in enriched
    )
    args.output_jsonl.write_text(output_text, encoding="utf-8")

    schema_counts: defaultdict[str, Counter[str]] = defaultdict(Counter)
    key_presence: Counter[str] = Counter()
    for row in enriched:
        row_paths: set[str] = set()
        for path, value_type in flatten_schema(row):
            schema_counts[path][value_type] += 1
            row_paths.add(path)
        key_presence.update(row_paths)
    schema = {
        "protocol": "jbi_ch1_frozen_evaluation_schema_v1",
        "rows": len(enriched),
        "fields": {
            path: {
                "present_rows": key_presence[path],
                "types": dict(sorted(counts.items())),
            }
            for path, counts in sorted(schema_counts.items())
        },
    }
    args.schema_json.parent.mkdir(parents=True, exist_ok=True)
    args.schema_json.write_text(
        json.dumps(schema, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    diagnostics.update(
        {
            "validation_status": "success",
            "aggregate_rows": len(enriched),
            "species": expected_species,
            "per_species_rows": dict(sorted(by_species.items())),
            "per_species_shard_rows": {
                species: dict(sorted(counts.items()))
                for species, counts in sorted(by_species_shard.items())
            },
            "unique_photo_ids": len(seen_photo),
            "unique_blind_ids": len(seen_blind),
            "split_csv_sha256": sha256_file(args.split_csv),
            "aggregate_jsonl_sha256": sha256_file(args.output_jsonl),
            "schema_json_sha256": sha256_file(args.schema_json),
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "evaluation_opened_for_rule_tuning": False,
            "final_labels_created": False,
        }
    )
    args.manifest_json.parent.mkdir(parents=True, exist_ok=True)
    args.manifest_json.write_text(
        json.dumps(diagnostics, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(
        f"validated {len(enriched)} evaluation rows across {len(expected_species)} species; "
        f"sha256={diagnostics['aggregate_jsonl_sha256']}"
    )


if __name__ == "__main__":
    main()

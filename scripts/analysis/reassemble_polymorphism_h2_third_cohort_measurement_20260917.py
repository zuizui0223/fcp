#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import pandas as pd

from fcp_pipeline.photo_first_measurement_execution import reassemble_complete_measurement
from fcp_pipeline.third_cohort_prospective_execution_gate import (
    EXPECTED_ROWS,
    EXPECTED_SPECIES,
    EXPECTED_TERMINAL_PARTITIONS,
    MIN_MEASUREMENT_EVALUABLE_SPECIES,
    validate_stage_transition,
)

ROOT = Path(__file__).resolve().parents[2]
PROTOCOL = ROOT / "docs" / "POLYMORPHISM_H2_THIRD_COHORT_PROSPECTIVE_MEASUREMENT_PROTOCOL_20260917.md"
TARGET_PER_SPECIES = 100
MIN_CLASSIFIABLE = 40


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--results-dir", type=Path, required=True)
    p.add_argument("--worker-batch0", type=Path, required=True)
    p.add_argument("--worker-batch1", type=Path, required=True)
    p.add_argument("--metadata-join-key", type=Path, required=True)
    p.add_argument("--firewall-manifest", type=Path, required=True)
    p.add_argument("--output-csv", type=Path, required=True)
    p.add_argument("--support-csv", type=Path, required=True)
    p.add_argument("--output-json", type=Path, required=True)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    firewall = json.loads(args.firewall_manifest.read_text(encoding="utf-8"))
    if firewall.get("schema") != "third_cohort_measurement_firewall_v1":
        raise RuntimeError("unexpected third-cohort firewall schema")
    if firewall.get("status") != "frozen_before_third_cohort_pixels":
        raise RuntimeError("third-cohort firewall status mismatch")
    if int(firewall.get("frozen_rows", -1)) != EXPECTED_ROWS:
        raise RuntimeError("third-cohort frozen row denominator drifted")
    if int(firewall.get("frozen_species", -1)) != EXPECTED_SPECIES:
        raise RuntimeError("third-cohort frozen species denominator drifted")
    if int(firewall.get("total_terminal_partitions", -1)) != EXPECTED_TERMINAL_PARTITIONS:
        raise RuntimeError("third-cohort terminal partition count drifted")
    firewall_stage = firewall.get("stage_receipts", {}).get("firewall")
    start_stage = firewall.get("stage_receipts", {}).get("start")
    if not isinstance(firewall_stage, dict) or not isinstance(start_stage, dict):
        raise RuntimeError("third-cohort firewall lacks frozen stage receipts")

    workers = [
        pd.read_csv(args.worker_batch0, dtype=str).fillna(""),
        pd.read_csv(args.worker_batch1, dtype=str).fillna(""),
    ]
    worker = pd.concat(workers, ignore_index=True)
    if len(worker) != EXPECTED_ROWS or worker["measurement_id"].nunique() != EXPECTED_ROWS:
        raise RuntimeError("third-cohort worker denominator is incomplete")

    metadata_join = pd.read_csv(args.metadata_join_key, dtype={"measurement_id": str}, low_memory=False).fillna("")
    if len(metadata_join) != EXPECTED_ROWS or metadata_join["measurement_id"].nunique() != EXPECTED_ROWS:
        raise RuntimeError("third-cohort metadata join denominator is incomplete")
    if "species" not in metadata_join.columns:
        raise RuntimeError("third-cohort metadata join lacks species identity")
    if metadata_join["species"].astype(str).nunique() != EXPECTED_SPECIES:
        raise RuntimeError("third-cohort metadata join species denominator drifted")
    if not (metadata_join.groupby("species", observed=True).size().astype(int) == TARGET_PER_SPECIES).all():
        raise RuntimeError("third-cohort metadata join no longer has exactly 100 rows per species")

    result_pattern = re.compile(r"^b(\d+)_partition_s(\d{2})_p(\d{2})\.csv$")
    receipt_pattern = re.compile(r"^b(\d+)_partition_s(\d{2})_p(\d{2})\.json$")
    result_files: list[Path] = []
    result_keys: set[tuple[int, int, int]] = set()
    receipt_keys: set[tuple[int, int, int]] = set()
    receipt_rows = 0

    for path in sorted(args.results_dir.iterdir()):
        m = result_pattern.match(path.name)
        if m:
            key = tuple(int(m.group(i)) for i in (1, 2, 3))
            if key in result_keys:
                raise RuntimeError(f"duplicate third-cohort result partition {key}")
            result_keys.add(key)
            result_files.append(path)
            continue
        m = receipt_pattern.match(path.name)
        if m:
            key = tuple(int(m.group(i)) for i in (1, 2, 3))
            if key in receipt_keys:
                raise RuntimeError(f"duplicate third-cohort receipt partition {key}")
            receipt_keys.add(key)
            receipt = json.loads(path.read_text(encoding="utf-8"))
            if receipt.get("status") != "complete_random_photo_first_terminal_partition":
                raise RuntimeError(f"unexpected terminal receipt status for {key}")
            if int(receipt.get("semantic_shard", -1)) != key[1] or int(receipt.get("compute_partition", -1)) != key[2]:
                raise RuntimeError(f"terminal receipt identity mismatch for {key}")
            if receipt.get("source_urls_present") is not False or receipt.get("species_present") is not False or receipt.get("coordinates_present") is not False:
                raise RuntimeError(f"terminal partition leaked sealed metadata for {key}")
            receipt_rows += int(receipt.get("terminal_rows") or 0)

    expected_keys = {(b, s, p) for b in range(2) for s in range(32) for p in range(4)}
    if result_keys != expected_keys or receipt_keys != expected_keys:
        raise RuntimeError(
            "not_evaluable_incomplete_third_cohort_measurement_partitions: "
            f"missing_results={sorted(expected_keys-result_keys)}; missing_receipts={sorted(expected_keys-receipt_keys)}"
        )
    if receipt_rows != EXPECTED_ROWS:
        raise RuntimeError(f"third-cohort terminal receipt rows {receipt_rows} != {EXPECTED_ROWS}")

    partition_stage = {
        "stage": "PARTITION_MEASUREMENT_COMPLETE",
        "species": EXPECTED_SPECIES,
        "rows": EXPECTED_ROWS,
        "replacement_rows": 0,
        "replacement_species": 0,
        "terminal_acquisition_rows": EXPECTED_ROWS,
        "terminal_measurement_rows": EXPECTED_ROWS,
        "terminal_partitions": EXPECTED_TERMINAL_PARTITIONS,
        "persisted_image_pixels": False,
    }
    validate_stage_transition(firewall_stage, partition_stage)

    partition_results = [pd.read_csv(path, dtype={"measurement_id": str}).fillna("") for path in result_files]
    assembled = reassemble_complete_measurement(
        partition_results,
        worker,
        metadata_join,
        expected_partition_receipts=EXPECTED_TERMINAL_PARTITIONS,
    )
    joined = assembled.joined_photos.copy()
    if len(joined) != EXPECTED_ROWS or joined["measurement_id"].nunique() != EXPECTED_ROWS:
        raise RuntimeError("third-cohort joined terminal denominator is incomplete")
    if joined["species"].astype(str).nunique() != EXPECTED_SPECIES:
        raise RuntimeError("third-cohort joined species denominator drifted")
    if not (joined.groupby("species", observed=True).size().astype(int) == TARGET_PER_SPECIES).all():
        raise RuntimeError("third-cohort joined raw-photo target drifted")

    reassembly_stage = {
        "stage": "REASSEMBLY_COMPLETE",
        "species": EXPECTED_SPECIES,
        "rows": EXPECTED_ROWS,
        "replacement_rows": 0,
        "replacement_species": 0,
        "unique_measurement_ids": int(joined["measurement_id"].nunique()),
        "duplicate_measurement_ids": int(joined["measurement_id"].duplicated().sum()),
    }
    validate_stage_transition(partition_stage, reassembly_stage)

    biological = {"white", "yellow_orange", "red_pink", "blue_purple"}
    joined["global_classifiable"] = joined["morph"].isin(biological) & joined["measurement_status"].eq("classified_four_state_morph")
    fraction_cols = [
        "flower_fraction_white", "flower_fraction_yellow", "flower_fraction_orange",
        "flower_fraction_red", "flower_fraction_pink", "flower_fraction_magenta",
        "flower_fraction_purple", "flower_fraction_blue", "flower_fraction_bronze",
    ]
    missing = sorted(set(fraction_cols) - set(joined.columns))
    if missing:
        raise RuntimeError(f"third-cohort measurement lacks frozen palette fractions: {missing}")
    for col in fraction_cols:
        joined[col] = pd.to_numeric(joined[col], errors="raise")
    classifiable = joined["global_classifiable"].to_numpy(dtype=bool)
    palette_sum = joined.loc[classifiable, fraction_cols].sum(axis=1)
    if len(palette_sum) and not ((palette_sum - 1.0).abs() <= 1e-8).all():
        raise RuntimeError("classifiable third-cohort nine-colour fractions do not sum to one")

    support = (
        joined.groupby(["species", "inat_taxon_id"], observed=True)
        .agg(raw_photos=("measurement_id", "size"), classifiable_photos=("global_classifiable", "sum"))
        .reset_index()
    )
    if len(support) != EXPECTED_SPECIES:
        raise RuntimeError("third-cohort support table species denominator drifted")
    support["measurement_evaluable"] = support["classifiable_photos"].astype(int) >= MIN_CLASSIFIABLE
    evaluable_species = int(support["measurement_evaluable"].sum())
    support_decision = "PASS" if evaluable_species >= MIN_MEASUREMENT_EVALUABLE_SPECIES else "NOT_EVALUABLE"
    support_stage = {
        "stage": "SUPPORT_GATE_COMPLETE",
        "species": EXPECTED_SPECIES,
        "rows": EXPECTED_ROWS,
        "replacement_rows": 0,
        "replacement_species": 0,
        "measurement_evaluable_species": evaluable_species,
        "support_decision": support_decision,
    }
    validate_stage_transition(reassembly_stage, support_stage)

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    args.support_csv.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    joined.to_csv(args.output_csv, index=False, lineterminator="\n")
    support.to_csv(args.support_csv, index=False, lineterminator="\n")

    terminal_counts = joined["measurement_status"].astype(str).value_counts(dropna=False).sort_index()
    result = {
        "schema": "third_cohort_prospective_measurement_result_v1",
        "analysis": "polymorphism_h2_third_cohort_prospective_measurement",
        "date_jst": "2026-09-17",
        "status": "complete_third_cohort_location_blind_measurement_and_support_gate",
        "stage": "SUPPORT_GATE_COMPLETE",
        "protocol": str(PROTOCOL.relative_to(ROOT)),
        "species": EXPECTED_SPECIES,
        "rows": EXPECTED_ROWS,
        "target_rows_per_species": TARGET_PER_SPECIES,
        "partition_receipts": EXPECTED_TERMINAL_PARTITIONS,
        "terminal_result_rows": int(len(joined)),
        "classifiable_rows": int(joined["global_classifiable"].sum()),
        "nonclassifiable_rows": int((~joined["global_classifiable"]).sum()),
        "measurement_evaluable_species": evaluable_species,
        "support_decision": support_decision,
        "minimum_classifiable_photos_per_species": MIN_CLASSIFIABLE,
        "minimum_measurement_evaluable_species": MIN_MEASUREMENT_EVALUABLE_SPECIES,
        "replacement_rows": 0,
        "replacement_species": 0,
        "persisted_image_pixels": False,
        "H2_opened": False,
        "terminal_status_counts": {str(k): int(v) for k, v in terminal_counts.items()},
        "coordinate_colour_join_opened_after_complete_measurement": True,
        "stage_receipts": {
            "start": start_stage,
            "firewall": firewall_stage,
            "partition": partition_stage,
            "reassembly": reassembly_stage,
            "support": support_stage,
        },
        "lineage": {
            "firewall_manifest_sha256": sha256_file(args.firewall_manifest),
            "protocol_sha256": sha256_file(PROTOCOL),
            "measured_table_sha256": sha256_file(args.output_csv),
            "species_support_sha256": sha256_file(args.support_csv),
        },
        "files": {
            "measured_table": str(args.output_csv.relative_to(ROOT)),
            "species_support": str(args.support_csv.relative_to(ROOT)),
        },
        "hard_nonclaims": [
            "measurement support is not global polymorphism prevalence",
            "nonclassifiable rows are not a fifth biological morph",
            "measurement failure is not monomorphism",
            "this receipt does not itself test the frozen white/non-white axis",
        ],
    }
    args.output_json.write_text(
        json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    reread = json.loads(args.output_json.read_text(encoding="utf-8"))
    if reread != result:
        raise RuntimeError("third-cohort measurement result failed durable read-back")
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

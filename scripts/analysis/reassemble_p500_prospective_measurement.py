#!/usr/bin/env python3
"""Reassemble all 256 P500 blind measurement partitions and apply support gate."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
from pathlib import Path

import pandas as pd

from fcp_pipeline.p500_prospective_execution_gate import (
    EXPECTED_ROWS,
    EXPECTED_SPECIES,
    EXPECTED_TERMINAL_PARTITIONS,
    MIN_MEASUREMENT_EVALUABLE_SPECIES,
    validate_stage_transition,
)
from fcp_pipeline.photo_first_measurement_execution import reassemble_complete_measurement

ROOT = Path(__file__).resolve().parents[2]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_deterministic_gzip_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = frame.to_csv(index=False, lineterminator="\n").encode("utf-8")
    with path.open("wb") as fh:
        with gzip.GzipFile(filename="", mode="wb", compresslevel=9, mtime=0, fileobj=fh) as gz:
            gz.write(raw)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--results-dir", type=Path, required=True)
    p.add_argument("--worker-batch0", type=Path, required=True)
    p.add_argument("--worker-batch1", type=Path, required=True)
    p.add_argument("--metadata-join-key", type=Path, required=True)
    p.add_argument("--firewall-manifest", type=Path, required=True)
    p.add_argument("--output-csv-gz", type=Path, required=True)
    p.add_argument("--output-support", type=Path, required=True)
    p.add_argument("--output-manifest", type=Path, required=True)
    return p.parse_args()


def _base(stage: str) -> dict[str, object]:
    return {
        "stage": stage,
        "species": EXPECTED_SPECIES,
        "rows": EXPECTED_ROWS,
        "replacement_rows": 0,
        "replacement_species": 0,
    }


def main() -> int:
    args = parse_args()
    firewall = json.loads(args.firewall_manifest.read_text(encoding="utf-8"))
    if firewall.get("status") != "p500_measurement_firewall_frozen_before_pixels":
        raise RuntimeError("P500 firewall state mismatch")
    if firewall.get("species") != EXPECTED_SPECIES or firewall.get("rows") != EXPECTED_ROWS:
        raise RuntimeError("P500 firewall denominator mismatch")
    if firewall.get("measurement_ids") != EXPECTED_ROWS:
        raise RuntimeError("P500 firewall measurement-ID census mismatch")
    if firewall.get("terminal_partitions") != EXPECTED_TERMINAL_PARTITIONS:
        raise RuntimeError("P500 firewall partition census mismatch")
    if firewall.get("candidate_pixels_opened") is not False:
        raise RuntimeError("P500 firewall was not frozen before pixels")

    firewall_stage = _base("FIREWALL_FROZEN")
    firewall_stage.update(
        measurement_ids=EXPECTED_ROWS,
        terminal_partitions=EXPECTED_TERMINAL_PARTITIONS,
        candidate_pixels_opened=False,
    )

    workers = [
        pd.read_csv(args.worker_batch0, dtype=str).fillna(""),
        pd.read_csv(args.worker_batch1, dtype=str).fillna(""),
    ]
    worker = pd.concat(workers, ignore_index=True)
    if len(worker) != EXPECTED_ROWS or worker["measurement_id"].nunique() != EXPECTED_ROWS:
        raise RuntimeError("P500 worker denominator is incomplete")

    metadata_join = pd.read_csv(args.metadata_join_key, dtype={"measurement_id": str}).fillna("")
    if len(metadata_join) != EXPECTED_ROWS or metadata_join["measurement_id"].nunique() != EXPECTED_ROWS:
        raise RuntimeError("P500 metadata join denominator is incomplete")
    if int(metadata_join["inat_taxon_id"].nunique()) != EXPECTED_SPECIES:
        raise RuntimeError("P500 metadata join species denominator mismatch")

    result_pattern = re.compile(r"^b(\d+)_partition_s(\d{2})_p(\d{2})\.csv$")
    receipt_pattern = re.compile(r"^b(\d+)_partition_s(\d{2})_p(\d{2})\.json$")
    result_files: list[Path] = []
    result_keys: set[tuple[int, int, int]] = set()
    receipt_keys: set[tuple[int, int, int]] = set()
    terminal_rows = 0

    for path in sorted(args.results_dir.iterdir()):
        m = result_pattern.match(path.name)
        if m:
            key = tuple(int(m.group(i)) for i in (1, 2, 3))
            if key in result_keys:
                raise RuntimeError(f"duplicate P500 result partition: {key}")
            result_keys.add(key)
            result_files.append(path)
            continue
        m = receipt_pattern.match(path.name)
        if m:
            key = tuple(int(m.group(i)) for i in (1, 2, 3))
            if key in receipt_keys:
                raise RuntimeError(f"duplicate P500 receipt partition: {key}")
            receipt_keys.add(key)
            receipt = json.loads(path.read_text(encoding="utf-8"))
            if receipt.get("status") != "complete_random_photo_first_terminal_partition":
                raise RuntimeError(f"unexpected terminal receipt status: {key}")
            if int(receipt.get("semantic_shard", -1)) != key[1]:
                raise RuntimeError(f"receipt semantic shard mismatch: {key}")
            if int(receipt.get("compute_partition", -1)) != key[2]:
                raise RuntimeError(f"receipt compute partition mismatch: {key}")
            if receipt.get("source_urls_present") is not False:
                raise RuntimeError(f"source URL leaked into terminal receipt: {key}")
            if receipt.get("species_present") is not False:
                raise RuntimeError(f"species leaked into terminal receipt: {key}")
            if receipt.get("coordinates_present") is not False:
                raise RuntimeError(f"coordinates leaked into terminal receipt: {key}")
            terminal_rows += int(receipt.get("terminal_rows") or 0)

    expected_keys = {
        (b, s, p)
        for b in range(2)
        for s in range(32)
        for p in range(4)
    }
    if result_keys != expected_keys or receipt_keys != expected_keys:
        raise RuntimeError(
            "not_evaluable_p500_incomplete_measurement_partitions: "
            f"missing_results={sorted(expected_keys-result_keys)}; "
            f"missing_receipts={sorted(expected_keys-receipt_keys)}"
        )
    if terminal_rows != EXPECTED_ROWS:
        raise RuntimeError(f"P500 terminal rows {terminal_rows} != {EXPECTED_ROWS}")

    measured_stage = _base("PARTITION_MEASUREMENT_COMPLETE")
    measured_stage.update(
        terminal_acquisition_rows=terminal_rows,
        terminal_measurement_rows=terminal_rows,
        terminal_partitions=EXPECTED_TERMINAL_PARTITIONS,
        persisted_image_pixels=False,
    )
    validate_stage_transition(firewall_stage, measured_stage)

    partition_results = [
        pd.read_csv(path, dtype={"measurement_id": str}).fillna("")
        for path in result_files
    ]
    assembled = reassemble_complete_measurement(
        partition_results,
        worker,
        metadata_join,
        expected_partition_receipts=EXPECTED_TERMINAL_PARTITIONS,
    )
    joined = assembled.joined_photos.copy()
    if len(joined) != EXPECTED_ROWS or joined["measurement_id"].nunique() != EXPECTED_ROWS:
        raise RuntimeError("P500 joined measurement denominator is incomplete")
    if int(joined["inat_taxon_id"].nunique()) != EXPECTED_SPECIES:
        raise RuntimeError("P500 joined species denominator mismatch")
    if not (joined.groupby("inat_taxon_id", observed=True).size().astype(int) == 100).all():
        raise RuntimeError("P500 joined per-species raw target drifted")

    reassembly_stage = _base("REASSEMBLY_COMPLETE")
    reassembly_stage.update(unique_measurement_ids=EXPECTED_ROWS, duplicate_measurement_ids=0)
    validate_stage_transition(measured_stage, reassembly_stage)

    biological = {"white", "yellow_orange", "red_pink", "blue_purple"}
    joined["global_classifiable"] = (
        joined["morph"].isin(biological)
        & joined["measurement_status"].eq("classified_four_state_morph")
    )
    fraction_cols = [
        "flower_fraction_white", "flower_fraction_yellow", "flower_fraction_orange",
        "flower_fraction_red", "flower_fraction_pink", "flower_fraction_magenta",
        "flower_fraction_purple", "flower_fraction_blue", "flower_fraction_bronze",
    ]
    missing = sorted(set(fraction_cols) - set(joined.columns))
    if missing:
        raise RuntimeError(f"P500 measurement lacks palette fractions: {missing}")

    joined["colour_white"] = pd.to_numeric(joined["flower_fraction_white"], errors="raise")
    joined["colour_yellow_orange"] = (
        pd.to_numeric(joined["flower_fraction_yellow"], errors="raise")
        + pd.to_numeric(joined["flower_fraction_orange"], errors="raise")
        + pd.to_numeric(joined["flower_fraction_bronze"], errors="raise")
    )
    joined["colour_red_pink"] = (
        pd.to_numeric(joined["flower_fraction_red"], errors="raise")
        + pd.to_numeric(joined["flower_fraction_pink"], errors="raise")
        + pd.to_numeric(joined["flower_fraction_magenta"], errors="raise")
    )
    joined["colour_blue_purple"] = (
        pd.to_numeric(joined["flower_fraction_blue"], errors="raise")
        + pd.to_numeric(joined["flower_fraction_purple"], errors="raise")
    )

    species_col = "species" if "species" in joined.columns else "query_species"
    support = (
        joined.groupby([species_col, "inat_taxon_id"], observed=True)
        .agg(
            raw_photos=("measurement_id", "size"),
            classifiable_photos=("global_classifiable", "sum"),
        )
        .reset_index()
    )
    support["measurement_evaluable"] = support["classifiable_photos"].astype(int) >= 40
    evaluable_species = int(support["measurement_evaluable"].sum())
    support_decision = "PASS" if evaluable_species >= MIN_MEASUREMENT_EVALUABLE_SPECIES else "NOT_EVALUABLE"

    support_stage = _base("SUPPORT_GATE_COMPLETE")
    support_stage.update(
        measurement_evaluable_species=evaluable_species,
        support_decision=support_decision,
    )
    validate_stage_transition(reassembly_stage, support_stage)

    write_deterministic_gzip_csv(joined, args.output_csv_gz)
    args.output_support.parent.mkdir(parents=True, exist_ok=True)
    support.to_csv(args.output_support, index=False, lineterminator="\n")
    args.output_manifest.parent.mkdir(parents=True, exist_ok=True)

    result = {
        "schema": "p500_prospective_measurement_result_v1",
        "status": "complete_p500_location_blind_measurement_and_support_gate",
        "stage": "SUPPORT_GATE_COMPLETE",
        "species": EXPECTED_SPECIES,
        "rows": EXPECTED_ROWS,
        "partition_receipts": EXPECTED_TERMINAL_PARTITIONS,
        "terminal_rows": terminal_rows,
        "joined_rows": int(len(joined)),
        "classifiable_rows": int(joined["global_classifiable"].sum()),
        "nonclassifiable_rows": int((~joined["global_classifiable"]).sum()),
        "measurement_evaluable_species": evaluable_species,
        "support_decision": support_decision,
        "H2_opened": False,
        "replacement_rows": 0,
        "replacement_species": 0,
        "coordinate_colour_join_opened_after_complete_measurement": True,
        "persisted_image_pixels": False,
        "lineage": {
            "firewall_manifest_sha256": sha256_file(args.firewall_manifest),
            "measured_table_sha256": sha256_file(args.output_csv_gz),
            "species_support_sha256": sha256_file(args.output_support),
        },
        "stage_receipts": {
            "partition_measurement": measured_stage,
            "reassembly": reassembly_stage,
            "support": support_stage,
        },
    }
    args.output_manifest.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

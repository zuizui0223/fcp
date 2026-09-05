#!/usr/bin/env python3
"""Reassemble all 256 global blind measurement partitions and apply the frozen support gate."""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import pandas as pd

from fcp_pipeline.photo_first_measurement_execution import reassemble_complete_measurement

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_FIREWALL = Path("/tmp/firewall/measurement_firewall_manifest.json")
DEFAULT_JOIN = Path("/tmp/firewall/sealed_keys/metadata_join_key.csv")
DEFAULT_EXECUTION = ROOT / "docs/supporting/global_monte_carlo_measurement_execution_contract_v1.json"
DEFAULT_CANDIDATE = ROOT / "docs/supporting/global_monte_carlo_candidate_acquisition_manifest_v1.json"


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
    p.add_argument("--metadata-join-key", type=Path, default=DEFAULT_JOIN)
    p.add_argument("--firewall-manifest", type=Path, default=DEFAULT_FIREWALL)
    p.add_argument("--execution-contract", type=Path, default=DEFAULT_EXECUTION)
    p.add_argument("--candidate-manifest", type=Path, default=DEFAULT_CANDIDATE)
    p.add_argument("--output-csv", type=Path, required=True)
    p.add_argument("--output-manifest", type=Path, required=True)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    execution = json.loads(args.execution_contract.read_text(encoding="utf-8"))
    candidate = json.loads(args.candidate_manifest.read_text(encoding="utf-8"))
    firewall = json.loads(args.firewall_manifest.read_text(encoding="utf-8"))
    if execution.get("status") != "frozen_before_capacity_outcome_before_candidate_acquisition_outcome_and_before_global_candidate_pixels":
        raise RuntimeError("global measurement contract is not frozen pre-outcome")
    if candidate.get("premeasurement_gate", {}).get("pass") is not True:
        raise RuntimeError("global candidate acquisition gate did not pass")
    if firewall.get("status") != "global_measurement_firewall_frozen_before_pixels":
        raise RuntimeError("global measurement firewall state mismatch")

    target = int(candidate["capacity_selected_raw_photo_target"])
    full_species = int(candidate["full_target_species"])
    candidate_rows = full_species * target
    measurement_species = int(firewall["frozen_measurement_species"])
    expected_rows = int(firewall["frozen_measurement_rows"])
    if int(firewall["candidate_pool_rows"]) != candidate_rows or int(candidate["candidate_rows"]) != candidate_rows:
        raise RuntimeError("global full candidate denominator drifted")
    if int(firewall["candidate_full_target_species"]) != full_species:
        raise RuntimeError("global firewall full-target species denominator drifted")
    if expected_rows != measurement_species * target:
        raise RuntimeError("bounded measurement rows do not equal measurement species * target")
    if measurement_species != min(full_species, int(firewall["measurement_species_budget"])):
        raise RuntimeError("bounded measurement species denominator drifted")
    if int(firewall["selected_raw_photo_target"]) != target:
        raise RuntimeError("global firewall raw photo target drifted")

    workers = [
        pd.read_csv(args.worker_batch0, dtype=str).fillna(""),
        pd.read_csv(args.worker_batch1, dtype=str).fillna(""),
    ]
    worker = pd.concat(workers, ignore_index=True)
    if len(worker) != expected_rows or worker["measurement_id"].nunique() != expected_rows:
        raise RuntimeError("global worker denominator is incomplete")
    metadata_join = pd.read_csv(args.metadata_join_key, dtype={"measurement_id": str}).fillna("")
    if len(metadata_join) != expected_rows or metadata_join["measurement_id"].nunique() != expected_rows:
        raise RuntimeError("global metadata join denominator is incomplete")
    if int(metadata_join["inat_taxon_id"].nunique()) != measurement_species:
        raise RuntimeError("metadata join taxon denominator differs from bounded measurement species")

    result_pattern = re.compile(r"^b(\d+)_partition_s(\d{2})_p(\d{2})\.csv$")
    receipt_pattern = re.compile(r"^b(\d+)_partition_s(\d{2})_p(\d{2})\.json$")
    result_files: list[Path] = []
    result_keys: set[tuple[int, int, int]] = set()
    receipt_keys: set[tuple[int, int, int]] = set()
    receipt_rows = 0
    for path in sorted(args.results_dir.iterdir()):
        match = result_pattern.match(path.name)
        if match:
            key = tuple(int(match.group(i)) for i in (1, 2, 3))
            if key in result_keys:
                raise RuntimeError(f"duplicate global result partition {key}")
            result_keys.add(key)
            result_files.append(path)
            continue
        match = receipt_pattern.match(path.name)
        if match:
            key = tuple(int(match.group(i)) for i in (1, 2, 3))
            if key in receipt_keys:
                raise RuntimeError(f"duplicate global receipt partition {key}")
            receipt_keys.add(key)
            receipt = json.loads(path.read_text(encoding="utf-8"))
            if receipt.get("status") != "complete_random_photo_first_terminal_partition":
                raise RuntimeError(f"unexpected terminal receipt status for {key}")
            if int(receipt.get("semantic_shard", -1)) != key[1] or int(receipt.get("compute_partition", -1)) != key[2]:
                raise RuntimeError(f"terminal receipt identity does not match filename for {key}")
            if receipt.get("source_urls_present") is not False or receipt.get("species_present") is not False or receipt.get("coordinates_present") is not False:
                raise RuntimeError(f"terminal receipt leaked metadata for {key}")
            receipt_rows += int(receipt.get("terminal_rows") or 0)

    expected_keys = {(batch, semantic, compute) for batch in range(2) for semantic in range(32) for compute in range(4)}
    if result_keys != expected_keys or receipt_keys != expected_keys:
        raise RuntimeError(
            "not_evaluable_global_incomplete_measurement_partitions: "
            f"missing_results={sorted(expected_keys-result_keys)}; missing_receipts={sorted(expected_keys-receipt_keys)}"
        )
    if receipt_rows != expected_rows:
        raise RuntimeError(f"global terminal receipt rows {receipt_rows} != frozen bounded denominator {expected_rows}")

    partition_results = [pd.read_csv(path, dtype={"measurement_id": str}).fillna("") for path in result_files]
    assembled = reassemble_complete_measurement(
        partition_results,
        worker,
        metadata_join,
        expected_partition_receipts=int(execution["partitioning"]["total_terminal_partitions"]),
    )
    joined = assembled.joined_photos.copy()
    if len(joined) != expected_rows or joined["measurement_id"].nunique() != expected_rows:
        raise RuntimeError("global joined measurement denominator is incomplete")
    if int(joined["inat_taxon_id"].nunique()) != measurement_species:
        raise RuntimeError("global joined bounded species denominator drifted")
    if not (joined.groupby("inat_taxon_id", observed=True).size().astype(int) == target).all():
        raise RuntimeError("global joined raw photo target drifted")

    biological = {"white", "yellow_orange", "red_pink", "blue_purple"}
    joined["global_classifiable"] = joined["morph"].isin(biological) & joined["measurement_status"].eq("classified_four_state_morph")
    required_fraction_columns = {
        "flower_fraction_white", "flower_fraction_yellow", "flower_fraction_orange", "flower_fraction_bronze",
        "flower_fraction_red", "flower_fraction_pink", "flower_fraction_magenta", "flower_fraction_blue", "flower_fraction_purple",
    }
    missing_fractions = sorted(required_fraction_columns - set(joined.columns))
    if missing_fractions:
        raise RuntimeError(f"global measurement lacks frozen palette fractions: {missing_fractions}")
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
    colour_cols = ["colour_white", "colour_yellow_orange", "colour_red_pink", "colour_blue_purple"]
    classifiable = joined["global_classifiable"].to_numpy(dtype=bool)
    colour_sum = joined.loc[classifiable, colour_cols].sum(axis=1)
    if len(colour_sum) and not ((colour_sum - 1.0).abs() <= 1e-8).all():
        raise RuntimeError("classifiable global four-group colour vectors do not sum to one")

    species_support = (
        joined.groupby(["species", "inat_taxon_id"], observed=True)
        .agg(raw_photos=("measurement_id", "size"), classifiable_photos=("global_classifiable", "sum"))
        .reset_index()
    )
    minimum_classifiable = int(execution["postmeasurement_gate"]["minimum_classifiable_photos_per_species"])
    required_species = int(execution["postmeasurement_gate"]["minimum_inferential_species"])
    species_support["global_measurement_evaluable"] = species_support["classifiable_photos"].astype(int) >= minimum_classifiable
    evaluable_species = int(species_support["global_measurement_evaluable"].sum())
    passed = evaluable_species >= required_species
    decision = execution["postmeasurement_gate"]["if_pass"] if passed else execution["postmeasurement_gate"]["if_fail"]

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    args.output_manifest.parent.mkdir(parents=True, exist_ok=True)
    joined.to_csv(args.output_csv, index=False, lineterminator="\n")
    support_path = args.output_csv.with_name("global_monte_carlo_measurement_species_support_v1.csv")
    species_support.to_csv(support_path, index=False, lineterminator="\n")
    result = {
        "protocol": execution["protocol"],
        "status": "complete_global_location_blind_measurement_and_join",
        "candidate_pool_rows": candidate_rows,
        "candidate_full_target_species": full_species,
        "measurement_species_budget": int(firewall["measurement_species_budget"]),
        "measurement_species_selection_seed": int(firewall["measurement_species_selection_seed"]),
        "measurement_taxon_id_sha256": firewall["measurement_taxon_id_sha256"],
        "frozen_measurement_rows": expected_rows,
        "frozen_measurement_species": measurement_species,
        "selected_raw_photo_target": target,
        "partition_receipts": int(execution["partitioning"]["total_terminal_partitions"]),
        "terminal_result_rows": expected_rows,
        "joined_rows": expected_rows,
        "classifiable_rows": int(joined["global_classifiable"].sum()),
        "mixed_uncertain_rows": int((~joined["global_classifiable"]).sum()),
        "coordinate_colour_join_opened_after_complete_measurement": True,
        "rgfca_colour_columns": colour_cols,
        "postmeasurement_gate": {
            "minimum_classifiable_photos_per_species": minimum_classifiable,
            "evaluable_species": evaluable_species,
            "required_species": required_species,
            "pass": bool(passed),
            "decision": str(decision),
        },
        "measurement_infrastructure_source_commit": execution["validated_measurement_infrastructure"]["source_commit"],
        "g1_g3_inference_run": False,
        "external_overlay_opened": False,
        "lineage": {
            "candidate_photos_sha256": firewall["candidate_photos_sha256"],
            "measurement_budget_amendment_sha256": firewall["measurement_budget_amendment_sha256"],
            "execution_contract_sha256": sha256_file(args.execution_contract),
            "measured_table_sha256": sha256_file(args.output_csv),
            "species_support_sha256": sha256_file(support_path),
        },
        "files": {
            "measured_table": str(args.output_csv.relative_to(ROOT)),
            "species_support": str(support_path.relative_to(ROOT)),
        },
    }
    args.output_manifest.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

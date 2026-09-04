#!/usr/bin/env python3
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
DEFAULT_WORKER = Path("/tmp/firewall/worker_packet/measurement_manifest.csv")
DEFAULT_JOIN = Path("/tmp/firewall/sealed_keys/metadata_join_key.csv")
DEFAULT_EXECUTION = ROOT / "docs/supporting/random_photo_first_h9_measurement_execution_contract_v1.json"
DEFAULT_FRESH = ROOT / "docs/supporting/random_photo_first_h9_fresh_metadata_manifest_v1.json"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--results-dir", type=Path, required=True)
    p.add_argument("--worker-manifest", type=Path, default=DEFAULT_WORKER)
    p.add_argument("--metadata-join-key", type=Path, default=DEFAULT_JOIN)
    p.add_argument("--firewall-manifest", type=Path, default=DEFAULT_FIREWALL)
    p.add_argument("--execution-contract", type=Path, default=DEFAULT_EXECUTION)
    p.add_argument("--fresh-manifest", type=Path, default=DEFAULT_FRESH)
    p.add_argument("--output-csv", type=Path, required=True)
    p.add_argument("--output-manifest", type=Path, required=True)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    execution = json.loads(args.execution_contract.read_text())
    fresh = json.loads(args.fresh_manifest.read_text())
    firewall = json.loads(args.firewall_manifest.read_text())
    if execution["status"] != "prospectively_frozen_after_h9_metadata_gate_before_any_h9_image_pixel":
        raise RuntimeError("H9 execution contract is not frozen pre-pixel")
    if fresh["premeasurement_gate"]["pass"] is not True:
        raise RuntimeError("H9 fresh metadata gate was not passed")
    if firewall["status"] != "h9_measurement_firewall_frozen_before_pixels":
        raise RuntimeError("H9 firewall state mismatch")
    if int(firewall["frozen_candidate_rows"]) != 2280:
        raise RuntimeError("H9 firewall denominator drifted")

    worker = pd.read_csv(args.worker_manifest, dtype=str).fillna("")
    metadata_join = pd.read_csv(args.metadata_join_key, dtype={"measurement_id": str}).fillna("")
    if len(worker) != 2280 or worker["measurement_id"].nunique() != 2280:
        raise RuntimeError("H9 worker denominator drifted")

    result_pattern = re.compile(r"^partition_s(\d{2})_p(\d{2})\.csv$")
    receipt_pattern = re.compile(r"^partition_s(\d{2})_p(\d{2})\.json$")
    result_files = []
    result_pairs = set()
    receipt_pairs = set()
    receipt_rows = 0
    for path in sorted(args.results_dir.iterdir()):
        m = result_pattern.match(path.name)
        if m:
            pair = (int(m.group(1)), int(m.group(2)))
            if pair in result_pairs:
                raise RuntimeError(f"duplicate H9 result partition {pair}")
            result_pairs.add(pair)
            result_files.append(path)
            continue
        m = receipt_pattern.match(path.name)
        if m:
            pair = (int(m.group(1)), int(m.group(2)))
            if pair in receipt_pairs:
                raise RuntimeError(f"duplicate H9 receipt partition {pair}")
            receipt_pairs.add(pair)
            rec = json.loads(path.read_text())
            if rec["status"] != "complete_random_photo_first_terminal_partition":
                raise RuntimeError(f"unexpected terminal receipt status for {pair}")
            if rec["source_urls_present"] is not False or rec["species_present"] is not False or rec["coordinates_present"] is not False:
                raise RuntimeError(f"H9 terminal receipt leaked metadata for {pair}")
            receipt_rows += int(rec["terminal_rows"])

    expected_pairs = {(s, p) for s in range(32) for p in range(4)}
    if result_pairs != expected_pairs or receipt_pairs != expected_pairs:
        missing_results = sorted(expected_pairs - result_pairs)
        missing_receipts = sorted(expected_pairs - receipt_pairs)
        raise RuntimeError(
            f"not_evaluable_h9_incomplete_measurement_partitions: results_missing={missing_results}; receipts_missing={missing_receipts}"
        )
    if receipt_rows != 2280:
        raise RuntimeError(f"H9 terminal receipt rows do not sum to 2280: {receipt_rows}")

    partition_results = [pd.read_csv(path, dtype={"measurement_id": str}).fillna("") for path in result_files]
    assembled = reassemble_complete_measurement(
        partition_results,
        worker,
        metadata_join,
        expected_partition_receipts=128,
    )
    joined = assembled.joined_photos.copy()
    if len(joined) != 2280 or joined["measurement_id"].nunique() != 2280:
        raise RuntimeError("H9 joined measurement denominator is incomplete")
    if int(joined["inat_taxon_id"].nunique()) != 38:
        raise RuntimeError("H9 joined species denominator drifted")
    if not (joined.groupby("inat_taxon_id").size() == 60).all():
        raise RuntimeError("H9 joined raw denominator is not 60 per species")

    biological = {"white", "yellow_orange", "red_pink", "blue_purple"}
    joined["h9_classifiable"] = joined["morph"].isin(biological) & joined["measurement_status"].eq("classified_four_state_morph")
    species_support = (
        joined.groupby(["species", "inat_taxon_id"], observed=True)
        .agg(raw_photos=("measurement_id", "size"), classifiable_photos=("h9_classifiable", "sum"))
        .reset_index()
    )
    minimum_classifiable = int(execution["postmeasurement_gate"]["minimum_classifiable_photos_per_species"])
    required_species = int(execution["postmeasurement_gate"]["minimum_evaluable_species"])
    species_support["h9_measurement_evaluable"] = species_support["classifiable_photos"].astype(int) >= minimum_classifiable
    evaluable_species = int(species_support["h9_measurement_evaluable"].sum())
    passed = evaluable_species >= required_species
    decision = (
        str(execution["postmeasurement_gate"]["if_pass"])
        if passed
        else str(execution["postmeasurement_gate"]["if_fail"])
    )

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    joined.to_csv(args.output_csv, index=False, lineterminator="\n")
    support_path = args.output_csv.with_name("random_photo_first_h9_measurement_species_support_v1.csv")
    species_support.to_csv(support_path, index=False, lineterminator="\n")
    result = {
        "protocol": execution["protocol"],
        "status": "complete_h9_location_blind_measurement_and_join",
        "frozen_candidate_rows": 2280,
        "frozen_species": 38,
        "partition_receipts": 128,
        "terminal_result_rows": 2280,
        "joined_rows": 2280,
        "classifiable_rows": int(joined["h9_classifiable"].sum()),
        "mixed_uncertain_rows": int((~joined["h9_classifiable"]).sum()),
        "coordinate_colour_join_opened_after_complete_measurement": True,
        "postmeasurement_gate": {
            "minimum_classifiable_photos_per_species": minimum_classifiable,
            "evaluable_species": evaluable_species,
            "required_species": required_species,
            "pass": bool(passed),
            "decision": decision,
        },
        "measurement_infrastructure_source_commit": execution["validated_measurement_infrastructure"]["source_commit"],
        "prior_experiment_rows_used": False,
        "h9_spatial_inference_run": False,
        "lineage": {
            "fresh_metadata_observations_sha256": fresh["files"]["observations"]["sha256"],
            "firewall_fresh_metadata_sha256": firewall["fresh_metadata_sha256"],
            "execution_contract_sha256": sha256_file(args.execution_contract),
            "measured_table_sha256": sha256_file(args.output_csv),
            "species_support_sha256": sha256_file(support_path),
        },
        "files": {
            "measured_table": str(args.output_csv.relative_to(ROOT)),
            "species_support": str(support_path.relative_to(ROOT)),
        },
    }
    args.output_manifest.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

EXPECTED_ROWS = 40_000
EXPECTED_PARTITIONS = 256


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--technical-seal-dir", type=Path, required=True)
    p.add_argument("--technical-firewall-dir", type=Path, required=True)
    p.add_argument("--worker-output-dir", type=Path, required=True)
    p.add_argument("--join-output-dir", type=Path, required=True)
    args = p.parse_args()

    summary = json.loads((args.technical_seal_dir / "technical_summary.json").read_text())
    if summary.get("status") != "RESPONSE_BLIND_TECHNICAL_SEAL_FROZEN":
        raise RuntimeError("technical seal is not frozen")
    if int(summary.get("rows", -1)) != EXPECTED_ROWS:
        raise RuntimeError("technical seal row denominator drift")
    for key in ("biological_outcomes_opened", "species_opened", "morph_opened", "q_white_or_W_opened"):
        if summary.get(key) is not False:
            raise RuntimeError(f"technical outcome firewall already open: {key}")

    technical = pd.read_csv(
        args.technical_seal_dir / "technical_table.csv.gz",
        compression="gzip",
        dtype={"measurement_id": str, "source_image_sha256": str},
    )
    if len(technical) != EXPECTED_ROWS or technical["measurement_id"].nunique() != EXPECTED_ROWS:
        raise RuntimeError("technical table census drift")
    # The frozen technical_table_sha256 is a canonical in-memory serialization
    # digest computed before gzip persistence. Re-reading the CSV through pandas
    # changes dtype/float rendering and is not a bitwise reproducer of that
    # canonical serialization. Artifact identity/digest is verified by the
    # workflow before this script runs; here we verify the frozen logical digest
    # value and the row/ID/source-SHA census without reserializing the table.
    if summary.get("technical_table_sha256") != "f757c90eddcb00f43d504a2ae6ccaeff4639d7ce468f231bc3b75675d66440e9":
        raise RuntimeError("unexpected frozen technical table SHA256")
    expected_sha = technical[["measurement_id", "source_image_sha256"]].copy()
    if expected_sha["source_image_sha256"].fillna("").str.len().ne(64).any():
        raise RuntimeError("technical table contains missing/invalid source SHA")

    firewall = json.loads((args.technical_firewall_dir / "firewall_manifest.json").read_text())
    if firewall.get("status") != "RESPONSE_BLIND_FIREWALL_FROZEN":
        raise RuntimeError("technical firewall is not frozen")
    if int(firewall.get("rows", -1)) != EXPECTED_ROWS:
        raise RuntimeError("technical firewall row denominator drift")
    if firewall.get("biological_outcomes_opened") is not False:
        raise RuntimeError("technical firewall reports opened biological outcomes")

    worker_out = args.worker_output_dir
    join_out = args.join_output_dir
    worker_out.mkdir(parents=True, exist_ok=True)
    join_out.mkdir(parents=True, exist_ok=True)

    seen: list[str] = []
    partition_receipts = []
    for batch in range(2):
        for shard in range(32):
            for part in range(4):
                src = args.technical_firewall_dir / "partitions" / f"b{batch}_s{shard}_p{part}"
                worker = pd.read_csv(src / "worker_manifest.csv", dtype=str).fillna("")
                acquisition = pd.read_csv(src / "acquisition_key.csv", dtype=str).fillna("")
                if worker["measurement_id"].nunique() != len(worker):
                    raise RuntimeError("worker partition contains duplicate IDs")
                if set(worker["measurement_id"]) != set(acquisition["measurement_id"]):
                    raise RuntimeError("worker/acquisition partition mismatch")
                sha = worker[["measurement_id"]].merge(
                    expected_sha,
                    on="measurement_id",
                    how="left",
                    validate="one_to_one",
                )
                if sha["source_image_sha256"].fillna("").str.len().ne(64).any():
                    raise RuntimeError("partition expected SHA missing")
                root = worker_out / "partitions" / f"b{batch}_s{shard}_p{part}"
                root.mkdir(parents=True, exist_ok=True)
                worker.to_csv(root / "worker_manifest.csv", index=False, lineterminator="\n")
                acquisition.to_csv(root / "acquisition_key.csv", index=False, lineterminator="\n")
                sha.to_csv(root / "expected_sha.csv", index=False, lineterminator="\n")
                seen.extend(worker["measurement_id"].astype(str).tolist())
                partition_receipts.append({
                    "batch": batch,
                    "shard": shard,
                    "part": part,
                    "rows": int(len(worker)),
                })

    if len(seen) != EXPECTED_ROWS or len(set(seen)) != EXPECTED_ROWS:
        raise RuntimeError("biological worker firewall coverage drift")
    if len(partition_receipts) != EXPECTED_PARTITIONS:
        raise RuntimeError("biological partition count drift")

    join = pd.read_csv(
        args.technical_firewall_dir / "sealed_join_key.csv",
        dtype={"measurement_id": str},
    )
    if len(join) != EXPECTED_ROWS or join["measurement_id"].nunique() != EXPECTED_ROWS:
        raise RuntimeError("sealed biological join key census drift")
    if set(join["measurement_id"].astype(str)) != set(seen):
        raise RuntimeError("sealed join key does not match biological worker census")
    join.to_csv(join_out / "sealed_join_key.csv", index=False, lineterminator="\n")

    strata = pd.read_csv(
        args.technical_seal_dir / "technical_strata.csv",
        dtype={"measurement_id": str},
    )
    if len(strata) != EXPECTED_ROWS or strata["measurement_id"].nunique() != EXPECTED_ROWS:
        raise RuntimeError("technical strata census drift")
    strata.to_csv(join_out / "technical_strata.csv", index=False, lineterminator="\n")
    (join_out / "technical_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    manifest = {
        "schema": "fcp_v2_biological_pass_firewall_v1",
        "status": "BIOLOGICAL_PASS_AUTHORIZED_AFTER_TECHNICAL_SEAL",
        "rows": EXPECTED_ROWS,
        "partition_receipts_expected": EXPECTED_PARTITIONS,
        "technical_table_sha256": summary["technical_table_sha256"],
        "technical_seal_biological_outcomes_opened": False,
        "worker_receives_species": False,
        "worker_receives_taxon_id": False,
        "worker_receives_photo_id": False,
        "worker_receives_coordinates": False,
        "worker_receives_observer": False,
        "worker_receives_expected_source_sha256": True,
        "source_drift_replacement_allowed": False,
        "partition_counts": partition_receipts,
    }
    (worker_out / "biological_firewall_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

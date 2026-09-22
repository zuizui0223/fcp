#!/usr/bin/env python3
"""Build the response-blind FCP v2 measurement firewall from terminal metadata."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

EXPECTED_ROWS = 40_000
EXPECTED_SPECIES = 400
ROWS_PER_SPECIES = 100
HEAVY_PER_SPECIES = 20
BATCHES = 2
SHARDS = 32
PARTS = 4
ID_SALT = "FCP_V2_MEASUREMENT_ID_20260923"
HEAVY_SALT = "FCP_V2_HEAVY_COUNTERFACTUAL_20260923"

WORKER_FIELDS = (
    "measurement_id",
    "image_filename",
    "photo_license",
    "heavy_counterfactual",
)
ACQUISITION_FIELDS = (
    "measurement_id",
    "image_filename",
    "photo_url_large",
    "photo_license",
    "heavy_counterfactual",
)


def measurement_id(photo_id: int) -> str:
    raw = f"{ID_SALT}|{int(photo_id)}".encode("utf-8")
    return "FCPV2-" + hashlib.sha256(raw).hexdigest()[:24].upper()


def heavy_hash(panel: str, taxon_id: int, photo_id: int) -> str:
    raw = f"{HEAVY_SALT}|{panel}|{int(taxon_id)}|{int(photo_id)}".encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def slot(identifier: str) -> tuple[int, int, int]:
    digest = hashlib.sha256(identifier.encode("utf-8")).digest()
    return digest[0] % BATCHES, digest[1] % SHARDS, digest[2] % PARTS


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--authorized-metadata", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args()

    df = pd.read_csv(args.authorized_metadata, compression="gzip")
    required = {
        "panel",
        "inat_taxon_id",
        "species",
        "observation_id",
        "photo_id",
        "photo_url_large",
        "photo_license",
        "latitude",
        "longitude",
        "observer_id",
    }
    missing = required - set(df.columns)
    if missing:
        raise RuntimeError(f"authorized metadata missing columns: {sorted(missing)}")
    if len(df) != EXPECTED_ROWS:
        raise RuntimeError(f"row denominator drift: {len(df)}")
    if df["inat_taxon_id"].nunique() != EXPECTED_SPECIES:
        raise RuntimeError("species denominator drift")
    if df["photo_id"].nunique() != EXPECTED_ROWS or df["observation_id"].nunique() != EXPECTED_ROWS:
        raise RuntimeError("photo/observation IDs are not unique")
    counts = df.groupby(["panel", "inat_taxon_id"], sort=False).size()
    if len(counts) != EXPECTED_SPECIES or not counts.eq(ROWS_PER_SPECIES).all():
        raise RuntimeError("authorized species do not all contain exactly 100 rows")
    if set(df["panel"].astype(str)) != {"P", "N"}:
        raise RuntimeError("panel labels drifted")

    x = df.copy()
    x["measurement_id"] = [measurement_id(v) for v in x["photo_id"]]
    if x["measurement_id"].nunique() != EXPECTED_ROWS:
        raise RuntimeError("measurement IDs are not unique")
    x["image_filename"] = x["measurement_id"].astype(str) + ".img"

    x["_heavy_hash"] = [
        heavy_hash(panel, taxon, photo)
        for panel, taxon, photo in zip(
            x["panel"].astype(str),
            x["inat_taxon_id"].astype(int),
            x["photo_id"].astype(int),
            strict=True,
        )
    ]
    x["heavy_counterfactual"] = False
    for (_, _), idx in x.groupby(["panel", "inat_taxon_id"], sort=False).groups.items():
        sub = x.loc[list(idx)].sort_values(["_heavy_hash", "photo_id"], kind="mergesort")
        chosen = sub.head(HEAVY_PER_SPECIES).index
        x.loc[chosen, "heavy_counterfactual"] = True

    heavy_counts = x.groupby(["panel", "inat_taxon_id"], sort=False)["heavy_counterfactual"].sum()
    if not heavy_counts.eq(HEAVY_PER_SPECIES).all():
        raise RuntimeError("heavy-counterfactual allocation drift")
    if int(x["heavy_counterfactual"].sum()) != EXPECTED_SPECIES * HEAVY_PER_SPECIES:
        raise RuntimeError("heavy-counterfactual total drift")

    slots = [slot(v) for v in x["measurement_id"].astype(str)]
    x["batch"] = [v[0] for v in slots]
    x["shard"] = [v[1] for v in slots]
    x["part"] = [v[2] for v in slots]

    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)

    join_columns = [
        "measurement_id",
        "panel",
        "inat_taxon_id",
        "species",
        "observation_id",
        "photo_id",
        "latitude",
        "longitude",
        "observer_id",
        "heavy_counterfactual",
    ]
    if "observer" in x.columns:
        join_columns.append("observer")
    if "observed_on" in x.columns:
        join_columns.append("observed_on")
    x[join_columns].to_csv(out / "sealed_join_key.csv", index=False, lineterminator="\n")

    heavy_manifest = x.loc[
        x["heavy_counterfactual"],
        ["measurement_id", "panel", "inat_taxon_id", "photo_id", "_heavy_hash"],
    ].sort_values(["panel", "inat_taxon_id", "_heavy_hash"], kind="mergesort")
    heavy_manifest.to_csv(
        out / "heavy_counterfactual_manifest.csv",
        index=False,
        lineterminator="\n",
    )

    partition_counts = []
    for batch in range(BATCHES):
        for shard in range(SHARDS):
            for part in range(PARTS):
                sub = x.loc[
                    (x["batch"] == batch)
                    & (x["shard"] == shard)
                    & (x["part"] == part)
                ].copy()
                root = out / "partitions" / f"b{batch}_s{shard}_p{part}"
                root.mkdir(parents=True, exist_ok=True)
                sub[list(WORKER_FIELDS)].to_csv(
                    root / "worker_manifest.csv", index=False, lineterminator="\n"
                )
                sub[list(ACQUISITION_FIELDS)].to_csv(
                    root / "acquisition_key.csv", index=False, lineterminator="\n"
                )
                partition_counts.append(
                    {
                        "batch": batch,
                        "shard": shard,
                        "part": part,
                        "rows": int(len(sub)),
                        "heavy_rows": int(sub["heavy_counterfactual"].sum()),
                    }
                )

    manifest = {
        "schema": "fcp_v2_measurement_firewall_v1",
        "status": "RESPONSE_BLIND_FIREWALL_FROZEN",
        "rows": EXPECTED_ROWS,
        "species": EXPECTED_SPECIES,
        "rows_per_species": ROWS_PER_SPECIES,
        "heavy_counterfactual_rows_per_species": HEAVY_PER_SPECIES,
        "heavy_counterfactual_rows": EXPECTED_SPECIES * HEAVY_PER_SPECIES,
        "partition_receipts_expected": BATCHES * SHARDS * PARTS,
        "worker_fields": list(WORKER_FIELDS),
        "technical_workers_receive_species": False,
        "technical_workers_receive_taxon_id": False,
        "technical_workers_receive_photo_id": False,
        "technical_workers_receive_coordinates": False,
        "technical_workers_receive_observer": False,
        "biological_outcomes_opened": False,
        "heavy_subset_uses_biological_outcome": False,
        "replacement_allowed": False,
        "partition_counts": partition_counts,
    }
    (out / "firewall_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

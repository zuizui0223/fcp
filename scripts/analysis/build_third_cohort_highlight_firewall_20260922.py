#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

EXPECTED_ROWS = 49_900
EXPECTED_SPECIES = 499
EXPECTED_METADATA_SHA256 = "13b25d72f20ed2b09ebcf3f80e0058aede08474a7e9051f7fb6ce1e521a16290"
SALT = "FCP_H2_THIRD_COHORT_HIGHLIGHT_VALIDITY_20260922_V1"
BATCHES = 2
SHARDS = 32
PARTS = 4


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def hid(photo_id: object) -> str:
    payload = f"{SALT}|{int(photo_id)}".encode()
    return "FCPH2HL-" + hashlib.sha256(payload).hexdigest()[:24].upper()


def slot(identifier: str) -> tuple[int, int, int]:
    digest = hashlib.sha256(identifier.encode()).digest()
    return digest[0] % BATCHES, digest[1] % SHARDS, digest[2] % PARTS


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--authorized-metadata", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args()

    observed_sha = sha256(args.authorized_metadata)
    if observed_sha != EXPECTED_METADATA_SHA256:
        raise RuntimeError(f"authorized metadata SHA mismatch: {observed_sha}")

    df = pd.read_csv(args.authorized_metadata, compression="gzip")
    required = {"species", "photo_id", "photo_url_large", "photo_license"}
    if not required.issubset(df.columns):
        raise RuntimeError(f"authorized metadata missing columns: {sorted(required-set(df.columns))}")
    if len(df) != EXPECTED_ROWS:
        raise RuntimeError(f"row denominator drift: {len(df)}")
    if df["species"].nunique() != EXPECTED_SPECIES:
        raise RuntimeError("species denominator drift")
    if df["photo_id"].nunique() != EXPECTED_ROWS:
        raise RuntimeError("photo_id denominator is not unique")
    if df["photo_url_large"].isna().any() or df["photo_url_large"].astype(str).eq("").any():
        raise RuntimeError("authorized metadata contains missing photo URL")

    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for row in df.itertuples(index=False):
        identifier = hid(row.photo_id)
        batch, shard, part = slot(identifier)
        rows.append({
            "measurement_id": identifier,
            "image_filename": identifier + ".img",
            "photo_license": str(row.photo_license),
            "photo_url_large": str(row.photo_url_large),
            "photo_id": int(row.photo_id),
            "batch": batch,
            "shard": shard,
            "part": part,
        })
    x = pd.DataFrame(rows)
    if x["measurement_id"].nunique() != EXPECTED_ROWS:
        raise RuntimeError("highlight measurement IDs are not unique")

    join_key = x[["measurement_id", "photo_id"]].copy()
    join_key.to_csv(out / "sealed_join_key.csv", index=False, lineterminator="\n")

    partition_counts = []
    for b in range(BATCHES):
        for s in range(SHARDS):
            for q in range(PARTS):
                sub = x.loc[(x.batch == b) & (x.shard == s) & (x.part == q)].copy()
                root = out / "partitions" / f"b{b}_s{s}_p{q}"
                root.mkdir(parents=True, exist_ok=True)
                sub[["measurement_id", "image_filename", "photo_license"]].to_csv(
                    root / "worker_manifest.csv", index=False, lineterminator="\n"
                )
                sub[["measurement_id", "image_filename", "photo_url_large", "photo_license"]].to_csv(
                    root / "acquisition_key.csv", index=False, lineterminator="\n"
                )
                partition_counts.append({
                    "batch": b, "shard": s, "part": q, "rows": int(len(sub))
                })

    manifest = {
        "schema": "third_cohort_highlight_firewall_v1",
        "status": "response_blind_firewall_frozen",
        "authorized_metadata_sha256": observed_sha,
        "rows": EXPECTED_ROWS,
        "species": EXPECTED_SPECIES,
        "photo_ids_unique": EXPECTED_ROWS,
        "measurement_ids_unique": EXPECTED_ROWS,
        "partition_receipts_expected": BATCHES * SHARDS * PARTS,
        "salt": SALT,
        "biological_outcomes_opened": False,
        "selection_uses_biological_outcome": False,
        "replacement_allowed": False,
        "partition_counts": partition_counts,
    }
    (out / "firewall_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()

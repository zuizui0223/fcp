#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "docs/supporting/random_photo_first_h9_individual_distance_contract_v1.json"
H7_EXCLUSION = ROOT / "data/frozen/random_photo_first_h7_exclusion_ledger_v1.csv"
H7_EXCLUSION_MANIFEST = ROOT / "docs/supporting/random_photo_first_h7_exclusion_manifest_v1.json"
H7_FRESH = ROOT / "data/frozen/random_photo_first_h7_fresh_metadata_v1.csv"
H7_FRESH_MANIFEST = ROOT / "docs/supporting/random_photo_first_h7_fresh_metadata_manifest_v1.json"
OUT = ROOT / "data/frozen/random_photo_first_h9_exclusion_ledger_v1.csv"
MANIFEST = ROOT / "docs/supporting/random_photo_first_h9_exclusion_manifest_v1.json"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    contract = json.loads(CONTRACT.read_text())
    h7e = json.loads(H7_EXCLUSION_MANIFEST.read_text())
    h7m = json.loads(H7_FRESH_MANIFEST.read_text())
    if sha256_file(H7_EXCLUSION) != h7e["exclusion_table_sha256"]:
        raise RuntimeError("H7 exclusion ledger lineage mismatch")
    if sha256_file(H7_FRESH) != h7m["files"]["observations"]["sha256"]:
        raise RuntimeError("H7 fresh metadata lineage mismatch")
    if bool(h7m["candidate_image_pixels_opened"]):
        raise RuntimeError("H7 pixels unexpectedly opened")

    prior = pd.read_csv(H7_EXCLUSION, usecols=["observation_id", "photo_id"])
    h7 = pd.read_csv(H7_FRESH, usecols=["observation_id", "photo_id"])
    prior["source"] = "pre_h7_experiments"
    h7["source"] = "h7_fresh_metadata"
    combined = pd.concat([prior, h7], ignore_index=True)
    combined["observation_id"] = pd.to_numeric(combined["observation_id"], errors="raise").astype("int64")
    combined["photo_id"] = pd.to_numeric(combined["photo_id"], errors="raise").astype("int64")
    pairs = (
        combined.groupby(["observation_id", "photo_id"], observed=True)["source"]
        .agg(lambda s: ";".join(sorted(set(map(str, s)))))
        .reset_index()
        .sort_values(["observation_id", "photo_id"], kind="mergesort")
    )
    expected_obs = int(contract["freshness_firewall"]["expected_unique_prior_observations"])
    expected_photo = int(contract["freshness_firewall"]["expected_unique_prior_photos"])
    if int(pairs["observation_id"].nunique()) != expected_obs:
        raise RuntimeError(f"unexpected H9 exclusion observation count: {pairs['observation_id'].nunique()}")
    if int(pairs["photo_id"].nunique()) != expected_photo:
        raise RuntimeError(f"unexpected H9 exclusion photo count: {pairs['photo_id'].nunique()}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    pairs.to_csv(OUT, index=False, lineterminator="\n")
    manifest = {
        "protocol": contract["protocol"],
        "status": "h9_prior_id_exclusion_frozen_before_h9_query",
        "source_rows": {
            "pre_h7_unique_pairs": int(len(prior)),
            "h7_fresh_rows": int(len(h7)),
        },
        "unique_exclusion_pairs": int(len(pairs)),
        "unique_exclusion_observation_ids": int(pairs["observation_id"].nunique()),
        "unique_exclusion_photo_ids": int(pairs["photo_id"].nunique()),
        "h7_exclusion_sha256": sha256_file(H7_EXCLUSION),
        "h7_fresh_sha256": sha256_file(H7_FRESH),
        "h9_exclusion_sha256": sha256_file(OUT),
        "h9_api_queries_opened": False,
        "h9_pixels_opened": False,
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

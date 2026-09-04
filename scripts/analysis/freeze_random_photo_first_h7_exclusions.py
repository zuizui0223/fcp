#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "docs/supporting/random_photo_first_h7_balanced_itv_contract_v1.json"
OLD_RANDOM = ROOT / "data/frozen/random_photo_first_candidate_pool_v1.csv"
CH1 = ROOT / "data/frozen/jbi_ch1_photo_source_manifest.csv"
OUT = ROOT / "data/frozen/random_photo_first_h7_exclusion_ledger_v1.csv"
MANIFEST = ROOT / "docs/supporting/random_photo_first_h7_exclusion_manifest_v1.json"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git_blob_sha1(path: Path) -> str:
    data = path.read_bytes()
    header = f"blob {len(data)}\0".encode()
    return hashlib.sha1(header + data).hexdigest()


def load_pairs(path: Path, source: str) -> pd.DataFrame:
    df = pd.read_csv(path, usecols=["observation_id", "photo_id"])
    df["observation_id"] = pd.to_numeric(df["observation_id"], errors="raise").astype("int64")
    df["photo_id"] = pd.to_numeric(df["photo_id"], errors="raise").astype("int64")
    df["source"] = source
    return df


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pr21-manifest", required=True)
    args = parser.parse_args()
    pr21 = Path(args.pr21_manifest)

    contract = json.loads(CONTRACT.read_text())
    pr21_spec = contract["freshness_firewall"]["exclude_pr21_terminal_manifest"]
    if git_blob_sha1(pr21) != pr21_spec["git_blob_sha"]:
        raise RuntimeError("PR21 manifest git blob mismatch")

    random_df = load_pairs(OLD_RANDOM, "random_photo_first_v1")
    ch1_df = load_pairs(CH1, "six_species_ch1")
    pr21_df = load_pairs(pr21, "pr21_terminal_60000")
    if len(pr21_df) != int(pr21_spec["expected_rows"]):
        raise RuntimeError(f"expected 60000 PR21 rows, found {len(pr21_df)}")

    combined = pd.concat([random_df, ch1_df, pr21_df], ignore_index=True)
    pair_sources = (
        combined.groupby(["observation_id", "photo_id"], observed=True)["source"]
        .agg(lambda s: ";".join(sorted(set(map(str, s)))))
        .reset_index()
        .sort_values(["observation_id", "photo_id"], kind="mergesort")
    )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    pair_sources.to_csv(OUT, index=False)

    manifest = {
        "protocol": contract["protocol"],
        "status": "h7_prior_photo_exclusion_ledger_frozen_before_fresh_query",
        "source_rows": {
            "random_photo_first_v1": int(len(random_df)),
            "six_species_ch1": int(len(ch1_df)),
            "pr21_terminal_60000": int(len(pr21_df)),
        },
        "unique_exclusion_pairs": int(len(pair_sources)),
        "unique_exclusion_observation_ids": int(pair_sources["observation_id"].nunique()),
        "unique_exclusion_photo_ids": int(pair_sources["photo_id"].nunique()),
        "old_random_sha256": sha256_file(OLD_RANDOM),
        "six_species_git_blob_sha": git_blob_sha1(CH1),
        "pr21_git_blob_sha": git_blob_sha1(pr21),
        "exclusion_table_sha256": sha256_file(OUT),
        "fresh_h7_api_queries_opened": False,
        "fresh_h7_pixels_opened": False,
        "enforcement": "H7 acquisition rejects a returned row if either its observation_id or its photo_id appears in this ledger.",
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Reconstruct pre-H7 historical observer identities from frozen metadata only."""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import pandas as pd

PROTOCOL_COMMIT = "6424b32c580f1f7c232d8876fdb43a86dd9af1a4"
RANDOM_ROWS = 20845
SIX_ROWS = 1200
PR21_ROWS = 60000
EXPECTED_PAIRS = 81803
RANDOM_SHA = "c124af7229d8ccb26d8ab20886772fbf61f12b8dabbcb529f83012bafcd5413d"
SIX_BLOB = "d8442bed4973556110a0678bee6ac580ddf90242"
PR21_SHA = "25835fb30510297d14b0863b204367e972c266fef3736cbed8f289116cf422cc"
H7_SHA = "b877a5fe5d36e0a5f84e73b9eae14b3c33b3f0fddc0ebc1a5f69331cc88b811a"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def norm(series: pd.Series) -> pd.Series:
    return series.astype("string").str.strip().str.replace(r"\.0$", "", regex=True).fillna("")


def load_source(path: Path, source: str, expected_rows: int) -> pd.DataFrame:
    frame = pd.read_csv(path, usecols=["observation_id", "photo_id", "observer_id"])
    if len(frame) != expected_rows:
        raise RuntimeError(f"{source} row count {len(frame)} != {expected_rows}")
    for col in ("observation_id", "photo_id", "observer_id"):
        frame[col] = norm(frame[col])
        if frame[col].eq("").any():
            raise RuntimeError(f"{source} has missing {col}")
    if frame[["observation_id", "photo_id"]].duplicated().any():
        raise RuntimeError(f"{source} contains duplicate observation/photo pairs")
    frame["source"] = source
    return frame


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--random", type=Path, required=True)
    p.add_argument("--six", type=Path, required=True)
    p.add_argument("--pr21", type=Path, required=True)
    p.add_argument("--h7-ledger", type=Path, required=True)
    p.add_argument("--output-ledger", type=Path, required=True)
    p.add_argument("--output-result", type=Path, required=True)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    for p in (args.random, args.six, args.pr21, args.h7_ledger):
        if not p.exists():
            raise RuntimeError(f"missing source {p}")

    if sha256_file(args.random) != RANDOM_SHA:
        raise RuntimeError("random-photo-first candidate pool SHA drifted")
    six_blob = subprocess.check_output(["git", "hash-object", str(args.six)], text=True).strip()
    if six_blob != SIX_BLOB:
        raise RuntimeError(f"six-species Git blob drifted: {six_blob}")
    if sha256_file(args.pr21) != PR21_SHA:
        raise RuntimeError("PR21 terminal manifest SHA drifted")
    if sha256_file(args.h7_ledger) != H7_SHA:
        raise RuntimeError("H7 exclusion ledger SHA drifted")

    sources = [
        load_source(args.random, "random_photo_first_candidate_pool_v1", RANDOM_ROWS),
        load_source(args.six, "six_species_ch1", SIX_ROWS),
        load_source(args.pr21, "pr21_terminal_60000", PR21_ROWS),
    ]
    union = pd.concat(sources, ignore_index=True)

    conflicts = (
        union.groupby(["observation_id", "photo_id"])["observer_id"]
        .nunique()
        .rename("n_observer_ids")
        .reset_index()
    )
    conflict_rows = conflicts.loc[conflicts["n_observer_ids"] != 1]
    if len(conflict_rows):
        raise RuntimeError(f"historical observer conflicts found for {len(conflict_rows)} pairs")

    collapsed = (
        union.sort_values(["observation_id", "photo_id", "source"], kind="mergesort")
        .groupby(["observation_id", "photo_id"], as_index=False)
        .agg(
            observer_id=("observer_id", "first"),
            source=("source", lambda x: ";".join(sorted(set(map(str, x))))),
        )
    )
    if len(collapsed) != EXPECTED_PAIRS:
        raise RuntimeError(f"source union unique pairs {len(collapsed)} != {EXPECTED_PAIRS}")

    frozen = pd.read_csv(args.h7_ledger, usecols=["observation_id", "photo_id"])
    for col in ("observation_id", "photo_id"):
        frozen[col] = norm(frozen[col])
    if len(frozen) != EXPECTED_PAIRS or frozen.drop_duplicates().shape[0] != EXPECTED_PAIRS:
        raise RuntimeError("frozen H7 exclusion ledger pair denominator drifted")

    joined = frozen.merge(
        collapsed,
        on=["observation_id", "photo_id"],
        how="left",
        validate="one_to_one",
        indicator=True,
    )
    missing = int((joined["_merge"] != "both").sum())
    if missing:
        raise RuntimeError(f"{missing} frozen H7 pairs lack historical observer reconstruction")
    joined = joined.drop(columns="_merge")
    if joined["observer_id"].eq("").any() or joined["observer_id"].isna().any():
        raise RuntimeError("reconstructed observer ledger contains missing observer IDs")

    reverse = collapsed.merge(
        frozen,
        on=["observation_id", "photo_id"],
        how="left",
        validate="one_to_one",
        indicator=True,
    )
    extra = int((reverse["_merge"] != "both").sum())
    if extra:
        raise RuntimeError(f"historical source union contains {extra} extra pairs outside frozen H7 ledger")

    joined = joined.sort_values(["observation_id", "photo_id"], kind="mergesort").reset_index(drop=True)
    args.output_ledger.parent.mkdir(parents=True, exist_ok=True)
    args.output_result.parent.mkdir(parents=True, exist_ok=True)
    joined.to_csv(args.output_ledger, index=False, lineterminator="\n")

    overlap_pairs = int(sum(len(x) - 1 for _, x in union.groupby(["observation_id", "photo_id"])))
    result = {
        "protocol": "rgfca-pre-h7-observer-identity-reconstruction-v1",
        "protocol_commit": PROTOCOL_COMMIT,
        "status": "complete_exact_historical_observer_reconstruction",
        "pre_h7_observer_identity_complete": True,
        "current_api_identity_queries_used": False,
        "image_pixels_opened": False,
        "flower_colour_used": False,
        "M1_M3_scores_used": False,
        "source_rows": {
            "random_photo_first_candidate_pool_v1": RANDOM_ROWS,
            "six_species_ch1": SIX_ROWS,
            "pr21_terminal_60000": PR21_ROWS,
            "total_before_source_overlap_collapse": int(len(union)),
        },
        "source_overlap_duplicate_pair_rows": overlap_pairs,
        "reconstructed_pairs": int(len(joined)),
        "frozen_h7_pairs": EXPECTED_PAIRS,
        "pair_coverage_missing": missing,
        "source_union_extra_pairs": extra,
        "observer_conflicts": int(len(conflict_rows)),
        "unique_historical_observers": int(joined["observer_id"].nunique()),
        "lineage": {
            "random_sha256": sha256_file(args.random),
            "six_git_blob_sha": six_blob,
            "pr21_sha256": sha256_file(args.pr21),
            "h7_exclusion_ledger_sha256": sha256_file(args.h7_ledger),
            "reconstructed_ledger_sha256": sha256_file(args.output_ledger),
        },
    }
    args.output_result.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

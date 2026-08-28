#!/usr/bin/env python3
"""Freeze the Chapter 1 480/720 photograph split from a canonical acquisition CSV."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess

import pandas as pd

from fcp_pipeline.photo_split import (
    SPLIT_VERSION,
    SplitSpec,
    assignment_hash,
    canonical_id_hash,
    freeze_photo_split,
)


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def resolve_column(frame: pd.DataFrame, explicit: str | None, candidates: tuple[str, ...], kind: str) -> str:
    if explicit:
        if explicit not in frame.columns:
            raise ValueError(f"requested {kind} column not found: {explicit}")
        return explicit
    for candidate in candidates:
        if candidate in frame.columns:
            return candidate
    raise ValueError(f"could not auto-detect {kind} column; tried {candidates}")


def current_commit() -> str | None:
    if os.environ.get("GITHUB_SHA"):
        return os.environ["GITHUB_SHA"]
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True, stderr=subprocess.DEVNULL
        ).strip()
    except Exception:
        return None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    parser.add_argument("manifest_json", type=Path)
    parser.add_argument("--species-col")
    parser.add_argument("--photo-id-col")
    parser.add_argument("--salt", default="fcp-jbi-ch1-photo-split-v1")
    args = parser.parse_args()

    raw = args.input_csv.read_bytes()
    frame = pd.read_csv(args.input_csv)
    species_col = resolve_column(
        frame,
        args.species_col,
        ("species", "taxon_name", "scientific_name", "accepted_species"),
        "species",
    )
    photo_id_col = resolve_column(
        frame,
        args.photo_id_col,
        ("photo_id", "image_id", "media_id"),
        "photo ID",
    )

    spec = SplitSpec(salt=args.salt)
    frozen = freeze_photo_split(
        frame,
        species_col=species_col,
        photo_id_col=photo_id_col,
        spec=spec,
    )

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    args.manifest_json.parent.mkdir(parents=True, exist_ok=True)
    frozen.to_csv(args.output_csv, index=False, lineterminator="\n")
    output_bytes = args.output_csv.read_bytes()

    per_species = (
        frozen.groupby([species_col, "split"], sort=True)
        .size()
        .unstack(fill_value=0)
        .reset_index()
        .to_dict(orient="records")
    )
    payload = {
        "status": "frozen",
        "protocol": SPLIT_VERSION,
        "source_file": args.input_csv.name,
        "source_raw_sha256": sha256_bytes(raw),
        "source_species_photo_id_sha256": canonical_id_hash(
            frame, species_col=species_col, photo_id_col=photo_id_col
        ),
        "assignment_sha256": assignment_hash(
            frozen, species_col=species_col, photo_id_col=photo_id_col
        ),
        "output_csv_sha256": sha256_bytes(output_bytes),
        "software_commit_sha": current_commit(),
        "species_column": species_col,
        "photo_id_column": photo_id_col,
        "split_basis": "SHA256(salt, species, photo_id) only; all other metadata ignored",
        "salt": spec.salt,
        "expected_species": spec.expected_species,
        "photographs_per_species": spec.photographs_per_species,
        "calibration_per_species": spec.calibration_per_species,
        "evaluation_per_species": spec.evaluation_per_species,
        "total_rows": int(len(frozen)),
        "per_species_split_counts": per_species,
        "outcome_blind": True,
        "evaluation_opened_for_rule_tuning": False,
    }
    args.manifest_json.write_text(
        json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

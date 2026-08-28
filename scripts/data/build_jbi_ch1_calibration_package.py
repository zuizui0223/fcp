#!/usr/bin/env python3
"""Build the blinded Chapter 1 calibration sheet and draft species codebook."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import os

import pandas as pd

from fcp_pipeline.photo_calibration import (
    CALIBRATION_PACKAGE_VERSION,
    build_calibration_sheet,
    calibration_summary,
    validate_calibration_sheet,
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


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
    parser.add_argument(
        "split_csv",
        nargs="?",
        type=Path,
        default=Path("data/frozen/jbi_ch1_photo_split_v1.csv"),
    )
    parser.add_argument(
        "--review-csv",
        type=Path,
        default=Path("data/calibration/jbi_ch1_calibration_review_v1.csv"),
    )
    parser.add_argument(
        "--codebook-json",
        type=Path,
        default=Path("data/calibration/jbi_ch1_species_colour_codebook_v1.json"),
    )
    parser.add_argument(
        "--manifest-json",
        type=Path,
        default=Path("data/calibration/jbi_ch1_calibration_package_manifest_v1.json"),
    )
    args = parser.parse_args()

    split = pd.read_csv(args.split_csv)
    review = build_calibration_sheet(split)
    validate_calibration_sheet(review)

    args.review_csv.parent.mkdir(parents=True, exist_ok=True)
    review.to_csv(args.review_csv, index=False, lineterminator="\n")

    species = sorted(review["species"].astype(str).unique())
    codebook = {
        "protocol": CALIBRATION_PACKAGE_VERSION,
        "status": "draft_calibration",
        "species": {
            sp: {
                "measurement_status": "draft",
                "visibility_rule_version": "draft-v1",
                "segmentation_rule_version": "draft-v1",
                "colour_code_version": "draft-v1",
                "allowed_colour_states": [],
                "unresolved_rule": "assign unresolved when no frozen state can be supported reliably",
                "not_evaluable_rule": "retain as measurement failure; never coerce to a colour state",
                "notes": "",
            }
            for sp in species
        },
    }
    args.codebook_json.write_text(
        json.dumps(codebook, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    summary = calibration_summary(review).to_dict(orient="records")
    manifest = {
        "protocol": CALIBRATION_PACKAGE_VERSION,
        "status": "generated_unmeasured",
        "software_commit_sha": current_commit(),
        "source_split": str(args.split_csv),
        "source_split_sha256": sha256(args.split_csv),
        "review_csv": str(args.review_csv),
        "review_csv_sha256": sha256(args.review_csv),
        "codebook_json": str(args.codebook_json),
        "codebook_json_sha256": sha256(args.codebook_json),
        "n_rows": int(len(review)),
        "species_count": len(species),
        "per_species": summary,
        "contains_evaluation_rows": False,
        "geography_observer_date_hidden_by_default": True,
        "evaluation_opened_for_rule_tuning": False,
    }
    args.manifest_json.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Build the exact 21,424-photo matched background-control reacquisition firewall."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from fcp_pipeline.photo_first_measurement import REFERENCE_RGB
from fcp_pipeline.photo_first_measurement_execution import (
    ACQUISITION_FIELDS,
    WORKER_FIELDS,
    semantic_shard,
    compute_partition,
)

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONTRACT = ROOT / "docs/supporting/global_rgfca_within_species_spatial_background_control_contract_v2.json"
DEFAULT_AMENDMENT = ROOT / "docs/supporting/global_rgfca_background_control_postoutcome_status_amendment_v2a.json"
DEFAULT_PREFLIGHT = ROOT / "docs/supporting/global_rgfca_background_control_preflight_v2.json"
DEFAULT_MEASURED = ROOT / "data/derived/global_monte_carlo_measured_photos_v1.csv"
DEFAULT_CANDIDATE = ROOT / "data/frozen/global_monte_carlo_candidate_photos_v1.csv"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def bool_series(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.casefold().isin({"true", "1"})


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    ap.add_argument("--amendment", type=Path, default=DEFAULT_AMENDMENT)
    ap.add_argument("--preflight", type=Path, default=DEFAULT_PREFLIGHT)
    ap.add_argument("--measured", type=Path, default=DEFAULT_MEASURED)
    ap.add_argument("--candidate", type=Path, default=DEFAULT_CANDIDATE)
    ap.add_argument("--output-dir", type=Path, required=True)
    args = ap.parse_args()

    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    amendment = json.loads(args.amendment.read_text(encoding="utf-8"))
    preflight = json.loads(args.preflight.read_text(encoding="utf-8"))
    if amendment.get("status") != "frozen_after_flower_only_omnibus_outcome_but_before_any_background_colour_recovery":
        raise RuntimeError("background timing amendment is not valid")
    if preflight.get("status") != "pass_background_control_preflight_without_pixels":
        raise RuntimeError("background preflight did not authorize reacquisition")
    if preflight.get("background_colour_opened") is not False or preflight.get("pixels_opened") is not False:
        raise RuntimeError("preflight unexpectedly opened background outcome")
    if sha256_file(args.measured) != contract["matched_frame"]["measured_table_sha256"]:
        raise RuntimeError("measured table differs from frozen matched frame")

    measured = pd.read_csv(args.measured, dtype={"measurement_id": str, "photo_id": str}).fillna("")
    classifiable = bool_series(measured["global_classifiable"])
    cf = measured.loc[classifiable].copy()
    counts = cf.groupby("species", observed=True).size()
    eligible = counts[counts >= 40].index
    frame = cf.loc[cf["species"].isin(eligible)].copy()
    if len(frame) != 21424 or frame["species"].nunique() != 369:
        raise RuntimeError("matched diagnostic frame is not exactly 21,424 photos / 369 species")
    if frame["measurement_id"].nunique() != len(frame) or frame["photo_id"].nunique() != len(frame):
        raise RuntimeError("matched diagnostic frame contains duplicate IDs")

    candidate = pd.read_csv(args.candidate, dtype={"photo_id": str}).fillna("")
    if candidate["photo_id"].nunique() != len(candidate):
        raise RuntimeError("candidate photo IDs are not unique")
    joined = frame.merge(
        candidate[["photo_id", "photo_url_large", "photo_license"]],
        on="photo_id",
        how="left",
        validate="one_to_one",
    )
    if len(joined) != 21424 or joined["photo_url_large"].astype(str).str.len().eq(0).any():
        raise RuntimeError("matched diagnostic URLs are incomplete")

    worker = pd.DataFrame({
        "measurement_id": joined["measurement_id"].astype(str),
        "image_filename": joined["measurement_id"].astype(str) + ".jpg",
        "photo_license": joined["photo_license"].astype(str),
    })
    acquisition = pd.DataFrame({
        "measurement_id": worker["measurement_id"],
        "image_filename": worker["image_filename"],
        "photo_url_large": joined["photo_url_large"].astype(str),
        "photo_license": worker["photo_license"],
    })
    if tuple(worker.columns) != WORKER_FIELDS or tuple(acquisition.columns) != ACQUISITION_FIELDS:
        raise RuntimeError("matched-control blinded/acquisition interfaces drifted")

    expected_cols = [
        "measurement_id",
        "image_sha256",
        "background_effective_pixels",
        "mask_pixels",
        *[f"palette_count_{name}" for name in REFERENCE_RGB],
    ]
    expected = joined[expected_cols].copy()
    expected = expected.rename(columns={
        "image_sha256": "expected_image_sha256",
        "background_effective_pixels": "expected_background_effective_pixels",
        "mask_pixels": "expected_flower_mask_pixels",
        **{f"palette_count_{name}": f"expected_flower_palette_count_{name}" for name in REFERENCE_RGB},
    })
    forbidden_expected = {"species", "latitude", "longitude", "photo_id", "photo_url_large"}
    if forbidden_expected.intersection(expected.columns):
        raise RuntimeError("safe expected key leaks source/species/geography")

    assignments = pd.DataFrame({
        "measurement_id": worker["measurement_id"].astype(str),
        "semantic_shard": [semantic_shard(x, 32) for x in worker["measurement_id"].astype(str)],
        "compute_partition": [compute_partition(x, 4) for x in worker["measurement_id"].astype(str)],
    })
    partition_counts = assignments.groupby(["semantic_shard", "compute_partition"], observed=True).size()
    if len(partition_counts) != 128 or int(partition_counts.sum()) != 21424:
        raise RuntimeError("matched-control partition census is not complete")

    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    (out / "sealed").mkdir(parents=True, exist_ok=True)
    worker.to_csv(out / "worker_manifest.csv", index=False, lineterminator="\n")
    acquisition.to_csv(out / "sealed/acquisition_key.csv", index=False, lineterminator="\n")
    expected.to_csv(out / "safe_expected_reproduction_key.csv", index=False, lineterminator="\n")
    assignments.to_csv(out / "partition_assignments.csv", index=False, lineterminator="\n")

    manifest = {
        "protocol": contract["protocol"],
        "status": "background_control_reacquisition_firewall_ready_before_pixels",
        "rows": 21424,
        "species_hidden_from_workers": True,
        "coordinates_hidden_from_workers": True,
        "source_urls_only_in_sealed_acquisition_key": True,
        "expected_reproduction_key_contains_species_or_coordinates": False,
        "semantic_shards": 32,
        "compute_partitions_per_shard": 4,
        "nonempty_partitions": int(len(partition_counts)),
        "minimum_partition_rows": int(partition_counts.min()),
        "maximum_partition_rows": int(partition_counts.max()),
        "pixels_opened": False,
        "background_colour_opened": False,
        "lineage": {
            "contract_sha256": sha256_file(args.contract),
            "amendment_sha256": sha256_file(args.amendment),
            "preflight_sha256": sha256_file(args.preflight),
            "measured_sha256": sha256_file(args.measured),
            "candidate_sha256": sha256_file(args.candidate),
            "worker_sha256": sha256_file(out / "worker_manifest.csv"),
            "acquisition_key_sha256": sha256_file(out / "sealed/acquisition_key.csv"),
            "expected_key_sha256": sha256_file(out / "safe_expected_reproduction_key.csv"),
        }
    }
    (out / "firewall_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

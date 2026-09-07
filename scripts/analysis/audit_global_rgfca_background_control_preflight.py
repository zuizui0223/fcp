#!/usr/bin/env python3
"""Preflight the frozen RGFCA matched-background recovery without opening pixels."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from fcp_pipeline.photo_first_measurement import REFERENCE_RGB

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONTRACT = ROOT / "docs/supporting/global_rgfca_within_species_spatial_background_control_contract_v2.json"
DEFAULT_MEASURED = ROOT / "data/derived/global_monte_carlo_measured_photos_v1.csv"
DEFAULT_CANDIDATE = ROOT / "data/frozen/global_monte_carlo_candidate_photos_v1.csv"
DEFAULT_MEASUREMENT_RESULT = ROOT / "docs/supporting/global_monte_carlo_measurement_result_v1.json"
DEFAULT_CANDIDATE_MANIFEST = ROOT / "docs/supporting/global_monte_carlo_candidate_acquisition_manifest_v1.json"


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
    ap.add_argument("--measured", type=Path, default=DEFAULT_MEASURED)
    ap.add_argument("--candidate", type=Path, default=DEFAULT_CANDIDATE)
    ap.add_argument("--measurement-result", type=Path, default=DEFAULT_MEASUREMENT_RESULT)
    ap.add_argument("--candidate-manifest", type=Path, default=DEFAULT_CANDIDATE_MANIFEST)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    measurement = json.loads(args.measurement_result.read_text(encoding="utf-8"))
    candidate_manifest = json.loads(args.candidate_manifest.read_text(encoding="utf-8"))
    if contract.get("status") != "frozen_before_any_within_species_omnibus_permutation_outcome_and_before_any_background_colour_recovery":
        raise RuntimeError("background-control contract is not frozen pre-recovery")
    expected_measured_sha = contract["matched_frame"]["measured_table_sha256"]
    observed_measured_sha = sha256_file(args.measured)
    if observed_measured_sha != expected_measured_sha:
        raise RuntimeError("measured table SHA differs from frozen background-control contract")
    if measurement.get("lineage", {}).get("measured_table_sha256") != expected_measured_sha:
        raise RuntimeError("measurement manifest and background-control contract disagree on measured table")
    expected_candidate_sha = str(measurement.get("lineage", {}).get("candidate_photos_sha256") or "")
    if not expected_candidate_sha or sha256_file(args.candidate) != expected_candidate_sha:
        raise RuntimeError("candidate table SHA differs from frozen measurement lineage")
    if candidate_manifest.get("lineage", {}).get("candidate_photos_sha256") not in {None, expected_candidate_sha}:
        raise RuntimeError("candidate manifest lineage disagrees with measurement lineage")

    measured = pd.read_csv(args.measured, dtype={"measurement_id": str, "photo_id": str}).fillna("")
    required_measured = {
        "measurement_id", "photo_id", "species", "latitude", "longitude", "image_sha256",
        "global_classifiable", "background_effective_pixels", "mask_pixels",
        *[f"palette_count_{name}" for name in REFERENCE_RGB],
    }
    missing = sorted(required_measured - set(measured.columns))
    if missing:
        raise RuntimeError(f"measured table lacks matched-control fields: {missing}")
    if measured["measurement_id"].nunique() != len(measured) or measured["photo_id"].nunique() != len(measured):
        raise RuntimeError("measured table IDs are not unique")

    classifiable = bool_series(measured["global_classifiable"])
    class_frame = measured.loc[classifiable].copy()
    minimum = int(contract["matched_frame"]["selection_rule"].split("at least ", 1)[1].split(" classifiable", 1)[0])
    counts = class_frame.groupby("species", observed=True).size()
    eligible = counts[counts >= minimum].index
    frame = class_frame.loc[class_frame["species"].isin(eligible)].copy()
    expected_rows = int(contract["matched_frame"]["photos"])
    expected_species = int(contract["matched_frame"]["species"])
    if len(frame) != expected_rows or frame["species"].nunique() != expected_species:
        raise RuntimeError(
            f"matched frame drifted: rows={len(frame)} species={frame['species'].nunique()}"
        )
    if int(frame.groupby("species", observed=True).size().min()) < minimum:
        raise RuntimeError("matched frame contains a species below its frozen old-G3 gate")

    candidate = pd.read_csv(args.candidate, dtype={"photo_id": str}).fillna("")
    required_candidate = {"photo_id", "photo_url_large", "photo_license"}
    missing_candidate = sorted(required_candidate - set(candidate.columns))
    if missing_candidate:
        raise RuntimeError(f"candidate table lacks recovery fields: {missing_candidate}")
    if candidate["photo_id"].nunique() != len(candidate):
        raise RuntimeError("candidate photo IDs are not unique")
    recovery = frame.merge(
        candidate[["photo_id", "photo_url_large", "photo_license"]],
        on="photo_id",
        how="left",
        validate="one_to_one",
    )
    if len(recovery) != expected_rows:
        raise RuntimeError("candidate join changed matched-control denominator")
    if recovery["photo_url_large"].astype(str).str.len().eq(0).any():
        raise RuntimeError("one or more matched-control photos lack the frozen source URL")
    if recovery["image_sha256"].astype(str).str.fullmatch(r"[0-9a-f]{64}").fillna(False).sum() != expected_rows:
        raise RuntimeError("one or more matched-control photos lack a frozen image SHA256")

    palette_cols = [f"palette_count_{name}" for name in REFERENCE_RGB]
    palette = recovery[palette_cols].apply(pd.to_numeric, errors="raise").astype(int)
    if (palette < 0).any().any():
        raise RuntimeError("stored flower palette counts contain a negative value")
    mask_pixels = pd.to_numeric(recovery["mask_pixels"], errors="raise").astype(int)
    if not palette.sum(axis=1).eq(mask_pixels).all():
        raise RuntimeError("stored twelve-anchor flower counts do not sum to stored mask pixels")
    background_pixels = pd.to_numeric(recovery["background_effective_pixels"], errors="raise").astype(int)
    if (background_pixels <= 0).any():
        raise RuntimeError("matched old-G3 frame contains a photo with no stored background pixels")

    report = {
        "protocol": contract["protocol"],
        "status": "pass_background_control_preflight_without_pixels",
        "pixels_opened": False,
        "background_colour_opened": False,
        "source_urls_persisted": False,
        "matched_rows": int(len(recovery)),
        "matched_species": int(recovery["species"].nunique()),
        "minimum_classifiable_photos_per_species_for_this_diagnostic_only": minimum,
        "candidate_url_join_complete": True,
        "original_image_sha_available_for_every_row": True,
        "stored_twelve_anchor_flower_counts_complete": True,
        "stored_flower_palette_counts_equal_mask_pixels": True,
        "stored_background_pixel_count_positive_for_every_row": True,
        "background_effective_pixels": {
            "minimum": int(background_pixels.min()),
            "median": float(background_pixels.median()),
            "maximum": int(background_pixels.max()),
        },
        "lineage": {
            "contract_sha256": sha256_file(args.contract),
            "measured_table_sha256": observed_measured_sha,
            "candidate_table_sha256": expected_candidate_sha,
            "measurement_result_sha256": sha256_file(args.measurement_result),
        },
        "next_gate": "Only exact-image-SHA matched re-acquisition using frozen ROI-v4 may open background pixels; no replacement photos or denominator adaptation are permitted."
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

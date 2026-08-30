#!/usr/bin/env python3
"""Freeze the prospective 12 x 200 Chapter 1 scale-up source and split.

The program accepts only a completed metadata-only feasibility report in which all 12
species passed. Assignment uses species and stable photo ID only. No image pixels or
measurement outcomes are read.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil

import pandas as pd

from fcp_pipeline.photo_split import (
    SplitSpec,
    assignment_hash,
    canonical_id_hash,
    freeze_photo_split,
)


PROTOCOL = "jbi-ch1-scaleup-photo-split-v1"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_preconditions(
    contract: dict,
    report: dict,
    candidate_path: Path,
) -> None:
    if contract.get("protocol") != PROTOCOL:
        raise ValueError(f"unexpected split protocol: {contract.get('protocol')!r}")
    if contract.get("status") != "frozen_before_scaleup_colour_measurement":
        raise ValueError("scale-up split contract is not frozen before colour measurement")

    expected = contract.get("preconditions", {})
    checks = {
        "feasibility_status": report.get("status"),
        "candidate_manifest_valid_for_final_freeze": report.get(
            "candidate_manifest_valid_for_final_freeze"
        ),
        "candidate_manifest_rows": report.get("candidate_manifest_rows"),
        "species_passed": report.get("species_passed"),
        "failed_species": report.get("failed_species"),
    }
    for key, observed in checks.items():
        required = expected.get(key)
        if observed != required:
            raise ValueError(
                f"feasibility precondition failed for {key}: observed={observed!r}, required={required!r}"
            )

    for key in (
        "candidate_images_downloaded",
        "flower_colour_pixels_inspected",
        "stage_a_effects_used",
        "stage_b_surfaces_used",
        "environmental_layers_used",
    ):
        if report.get(key) is not False:
            raise ValueError(f"feasibility report permits forbidden input/action: {key}")

    observed_candidate_hash = sha256(candidate_path)
    if report.get("candidate_manifest_sha256") != observed_candidate_hash:
        raise ValueError(
            "candidate manifest SHA256 does not match the completed feasibility report"
        )


def build_spec(contract: dict) -> SplitSpec:
    assignment = contract.get("assignment_rule", {})
    return SplitSpec(
        expected_species=int(contract["expected_species"]),
        photographs_per_species=int(contract["photographs_per_species"]),
        calibration_per_species=int(contract["calibration_per_species"]),
        evaluation_per_species=int(contract["evaluation_per_species"]),
        salt=str(assignment["salt"]),
    )


def freeze_scaleup(
    *,
    contract_path: Path,
    report_path: Path,
    candidate_path: Path,
    source_output: Path,
    split_output: Path,
    manifest_output: Path,
) -> dict:
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    report = json.loads(report_path.read_text(encoding="utf-8"))
    verify_preconditions(contract, report, candidate_path)

    source = pd.read_csv(candidate_path, dtype={"photo_id": str})
    required = {"species", "photo_id"}
    missing = sorted(required - set(source.columns))
    if missing:
        raise ValueError(f"candidate manifest missing required columns: {missing}")

    spec = build_spec(contract)
    frozen = freeze_photo_split(
        source,
        species_col="species",
        photo_id_col="photo_id",
        spec=spec,
    )

    if len(source) != int(contract["total_photographs"]):
        raise ValueError("candidate manifest row count differs from frozen contract")
    if int((frozen["split"] == "calibration").sum()) != int(
        contract["total_calibration"]
    ):
        raise ValueError("frozen calibration total differs from contract")
    if int((frozen["split"] == "evaluation").sum()) != int(
        contract["total_evaluation"]
    ):
        raise ValueError("frozen evaluation total differs from contract")

    source_output.parent.mkdir(parents=True, exist_ok=True)
    split_output.parent.mkdir(parents=True, exist_ok=True)
    manifest_output.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(candidate_path, source_output)
    frozen.to_csv(split_output, index=False, lineterminator="\n")

    per_species = (
        frozen.groupby(["species", "split"], sort=True)
        .size()
        .unstack(fill_value=0)
        .reset_index()
        .to_dict(orient="records")
    )
    manifest = {
        "protocol": PROTOCOL,
        "status": "scaleup_photo_source_and_split_frozen",
        "contract_sha256": sha256(contract_path),
        "feasibility_report_sha256": sha256(report_path),
        "candidate_manifest_sha256": sha256(candidate_path),
        "frozen_source_manifest_sha256": sha256(source_output),
        "source_species_photo_id_sha256": canonical_id_hash(
            source, species_col="species", photo_id_col="photo_id"
        ),
        "assignment_sha256": assignment_hash(
            frozen, species_col="species", photo_id_col="photo_id"
        ),
        "frozen_split_sha256": sha256(split_output),
        "expected_species": spec.expected_species,
        "photographs_per_species": spec.photographs_per_species,
        "calibration_per_species": spec.calibration_per_species,
        "evaluation_per_species": spec.evaluation_per_species,
        "total_rows": int(len(frozen)),
        "total_calibration": int((frozen["split"] == "calibration").sum()),
        "total_evaluation": int((frozen["split"] == "evaluation").sum()),
        "per_species_split_counts": per_species,
        "split_basis": "SHA256(salt, species, photo_id) only; all other metadata ignored",
        "salt": spec.salt,
        "outcome_blind": True,
        "image_pixels_read": False,
        "evaluation_opened_for_rule_tuning": False,
        "next_gate": "run 960-image calibration-only visibility, localization and continuous-colour measurement development",
    }
    manifest_output.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--contract",
        type=Path,
        default=Path("docs/supporting/jbi_ch1_scaleup_photo_split_contract_v1.json"),
    )
    parser.add_argument(
        "--report",
        type=Path,
        default=Path("docs/supporting/jbi_ch1_scaleup_inat_feasibility_v1.json"),
    )
    parser.add_argument(
        "--candidate",
        type=Path,
        default=Path("data/scaleup/jbi_ch1_scaleup_inat_candidate_manifest_v1.csv"),
    )
    parser.add_argument(
        "--source-output",
        type=Path,
        default=Path("data/scaleup/frozen/jbi_ch1_scaleup_photo_source_manifest_v1.csv"),
    )
    parser.add_argument(
        "--split-output",
        type=Path,
        default=Path("data/scaleup/frozen/jbi_ch1_scaleup_photo_split_v1.csv"),
    )
    parser.add_argument(
        "--manifest-output",
        type=Path,
        default=Path("data/scaleup/frozen/jbi_ch1_scaleup_photo_split_manifest_v1.json"),
    )
    args = parser.parse_args()

    manifest = freeze_scaleup(
        contract_path=args.contract,
        report_path=args.report,
        candidate_path=args.candidate,
        source_output=args.source_output,
        split_output=args.split_output,
        manifest_output=args.manifest_output,
    )
    print(json.dumps(manifest, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

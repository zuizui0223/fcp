#!/usr/bin/env python3
"""Audit iNaturalist metadata feasibility for the frozen 12-species scale-up cohort.

This program queries taxon, observation and photo metadata only. It does not download or
inspect image pixels. A failed species is reported, not silently replaced.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
from typing import Any


PROTOCOL = "jbi-ch1-scaleup-inat-feasibility-v1"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_acquisition_core():
    path = Path(__file__).with_name("acquire_jbi_ch1_inat_photos.py")
    spec = importlib.util.spec_from_file_location("jbi_ch1_inat_acquisition_core", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load acquisition core from {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_cohort(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"cohort_order", "canonical_name", "family", "rank"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"cohort CSV missing columns: {sorted(missing)}")
        rows = [
            {key: (value or "").strip() for key, value in row.items()}
            for row in reader
        ]
    if len(rows) != 12:
        raise ValueError(f"expected 12 frozen scale-up species, found {len(rows)}")
    orders = [int(row["cohort_order"]) for row in rows]
    if orders != list(range(1, 13)):
        raise ValueError(f"cohort order is not 1..12: {orders}")
    names = [row["canonical_name"] for row in rows]
    if len(set(names)) != 12 or any(not name for name in names):
        raise ValueError("cohort names must be 12 unique non-empty species")
    return rows


def audit_cohort(
    core,
    cohort: list[dict[str, str]],
    config: dict[str, Any],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    if config.get("protocol") != PROTOCOL:
        raise ValueError(f"unexpected feasibility protocol: {config.get('protocol')!r}")
    if config.get("status") != "frozen_before_candidate_api_queries":
        raise ValueError("feasibility contract is not frozen before API queries")
    audit_rules = config.get("audit_rules", {})
    forbidden_true = [
        key
        for key in (
            "candidate_images_downloaded",
            "flower_colour_pixels_inspected",
            "stage_a_effects_used",
            "stage_b_surfaces_used",
            "environmental_layers_used",
            "failed_species_replaced_automatically",
        )
        if audit_rules.get(key) is not False
    ]
    if forbidden_true:
        raise ValueError(f"forbidden feasibility inputs/actions enabled: {forbidden_true}")
    if audit_rules.get("all_species_must_pass_before_final_source_freeze") is not True:
        raise ValueError("contract must require all species to pass before final source freeze")

    reports: list[dict[str, Any]] = []
    passing_rows: list[dict[str, Any]] = []

    for cohort_row in cohort:
        species = cohort_row["canonical_name"]
        base = {
            "cohort_order": int(cohort_row["cohort_order"]),
            "rank": int(cohort_row["rank"]),
            "species": species,
            "family": cohort_row["family"],
            "candidate_images_downloaded": False,
            "flower_colour_pixels_inspected": False,
        }
        candidates: list[dict[str, Any]] = []
        try:
            taxon, candidates = core.fetch_candidates(config, species)
        except Exception as exc:
            reports.append(
                {
                    **base,
                    "status": "taxon_or_api_failed",
                    "gate_pass": False,
                    "error_type": type(exc).__name__,
                    "error": str(exc),
                    "candidate_count": len(candidates),
                }
            )
            continue

        taxon_summary = {
            "inat_taxon_id": int(taxon["id"]),
            "inat_taxon_name": str(taxon.get("name") or ""),
            "inat_taxon_rank": str(taxon.get("rank") or ""),
            "candidate_count": len(candidates),
        }
        try:
            selected = core.select_rows(config, species, candidates)
        except Exception as exc:
            reports.append(
                {
                    **base,
                    **taxon_summary,
                    "status": "insufficient_balanced_metadata",
                    "gate_pass": False,
                    "error_type": type(exc).__name__,
                    "error": str(exc),
                }
            )
            continue

        qc = core.species_qc(config, species, taxon, candidates, selected)
        gate_pass = bool(qc.get("gate_pass")) and len(selected) == int(
            config["selection"]["target_photographs_per_species"]
        )
        report = {
            **base,
            **taxon_summary,
            "status": "pass" if gate_pass else "metadata_qc_failed",
            "gate_pass": gate_pass,
            "selected_count": len(selected),
            "unique_observers": int(qc.get("unique_observers", 0)),
            "unique_spatial_cells": int(qc.get("unique_spatial_cells", 0)),
            "unique_calendar_months": int(qc.get("unique_calendar_months", 0)),
            "maximum_observer_fraction": qc.get("maximum_observer_fraction"),
            "maximum_spatial_cell_fraction": qc.get("maximum_spatial_cell_fraction"),
            "maximum_month_fraction": qc.get("maximum_month_fraction"),
            "gate_failures": list(qc.get("gate_failures", [])),
            "observer_counts_top10": qc.get("observer_counts_top10", []),
            "month_counts": qc.get("month_counts", []),
        }
        reports.append(report)
        if gate_pass:
            for row in selected:
                passing_rows.append(
                    {
                        "scaleup_cohort_order": int(cohort_row["cohort_order"]),
                        "scaleup_literature_rank": int(cohort_row["rank"]),
                        "scaleup_family": cohort_row["family"],
                        **row,
                    }
                )

    reports.sort(key=lambda row: int(row["cohort_order"]))
    passing_rows.sort(
        key=lambda row: (
            int(row["scaleup_cohort_order"]),
            str(row.get("selection_hash", "")),
            str(row.get("photo_id", "")),
        )
    )
    return reports, passing_rows


def write_candidate_manifest(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    preferred = [
        "scaleup_cohort_order",
        "scaleup_literature_rank",
        "scaleup_family",
        "species",
        "inat_taxon_id",
        "inat_taxon_name",
        "observation_id",
        "photo_id",
        "photo_url",
        "photo_url_api",
        "photo_license",
        "attribution",
        "observation_license",
        "quality_grade",
        "latitude",
        "longitude",
        "positional_accuracy_m",
        "geoprivacy",
        "obscured",
        "observed_on",
        "observed_month",
        "time_observed_at",
        "created_at",
        "observer_id",
        "observer",
        "place_guess",
        "spatial_cell",
        "selection_hash",
    ]
    available = {key for row in rows for key in row}
    fields = [key for key in preferred if key in available]
    fields.extend(sorted(available - set(fields)))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--contract",
        type=Path,
        default=Path("docs/supporting/jbi_ch1_scaleup_inat_feasibility_contract_v1.json"),
    )
    parser.add_argument(
        "--cohort",
        type=Path,
        default=Path("docs/supporting/jbi_ch1_scaleup_cohort_v1.csv"),
    )
    parser.add_argument(
        "--cohort-manifest",
        type=Path,
        default=Path("docs/supporting/jbi_ch1_scaleup_cohort_manifest_v1.json"),
    )
    parser.add_argument(
        "--candidate-manifest",
        type=Path,
        default=Path("data/scaleup/jbi_ch1_scaleup_inat_candidate_manifest_v1.csv"),
    )
    parser.add_argument(
        "--report",
        type=Path,
        default=Path("docs/supporting/jbi_ch1_scaleup_inat_feasibility_v1.json"),
    )
    args = parser.parse_args()

    config = json.loads(args.contract.read_text(encoding="utf-8"))
    cohort_manifest = json.loads(args.cohort_manifest.read_text(encoding="utf-8"))
    if cohort_manifest.get("status") != "scaleup_cohort_selected_before_photo_acquisition":
        raise ValueError("cohort manifest is not frozen before photo acquisition")
    cohort = read_cohort(args.cohort)
    if [row["canonical_name"] for row in cohort] != cohort_manifest.get("selected_species"):
        raise ValueError("cohort CSV and cohort manifest species order differ")

    core = load_acquisition_core()
    reports, passing_rows = audit_cohort(core, cohort, config)
    write_candidate_manifest(args.candidate_manifest, passing_rows)

    passing_species = [row["species"] for row in reports if row["gate_pass"]]
    failed_species = [row["species"] for row in reports if not row["gate_pass"]]
    all_pass = len(passing_species) == 12
    per_species_rows = Counter(row.get("species") for row in passing_rows)
    candidate_manifest_valid_for_final_freeze = (
        all_pass
        and len(passing_rows) == 2400
        and len(per_species_rows) == 12
        and all(count == 200 for count in per_species_rows.values())
        and len({str(row.get("photo_id")) for row in passing_rows}) == 2400
    )

    report = {
        "protocol": PROTOCOL,
        "status": "pass_all_12" if candidate_manifest_valid_for_final_freeze else "audit_complete_replacement_required",
        "contract_sha256": sha256(args.contract),
        "cohort_csv_sha256": sha256(args.cohort),
        "cohort_manifest_sha256": sha256(args.cohort_manifest),
        "candidate_manifest_sha256": sha256(args.candidate_manifest),
        "candidate_images_downloaded": False,
        "flower_colour_pixels_inspected": False,
        "stage_a_effects_used": False,
        "stage_b_surfaces_used": False,
        "environmental_layers_used": False,
        "species_audited": len(reports),
        "species_passed": len(passing_species),
        "passing_species": passing_species,
        "failed_species": failed_species,
        "candidate_manifest_rows": len(passing_rows),
        "candidate_manifest_valid_for_final_freeze": candidate_manifest_valid_for_final_freeze,
        "species_results": reports,
        "next_gate": (
            "freeze the 2400-row source and deterministic 960/1440 split"
            if candidate_manifest_valid_for_final_freeze
            else "version a literature-led replacement amendment for failed species before any colour measurement"
        ),
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

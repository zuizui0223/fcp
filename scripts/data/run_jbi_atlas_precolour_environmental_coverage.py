#!/usr/bin/env python3
"""Evaluate environmental coverage on a frozen atlas geometry before pixels.

The live-API feasibility result is deliberately unable to authorize image
acquisition.  The same runner must be applied to the final dated-source cohort,
where it emits the status consumed by the protected colour-join runner.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
from typing import Any, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.atlas_environment import evaluate_environmental_coverage_gate


SCALES = (100, 250, 500)
SOURCE_STAGES = ("live-feasibility", "final-dated-source")


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def parse_false(value: object, *, field: str) -> None:
    if value is False or str(value).strip().casefold() == "false":
        return
    raise ValueError(f"{field} must remain false before atlas image acquisition")


def selected_geometry(
    feasibility: Mapping[str, Any], panels: Sequence[Mapping[str, Any]]
) -> tuple[list[str], list[dict[str, Any]]]:
    """Validate the 8 x 25 panel denominator and return its geometry rows."""

    if feasibility.get("status") != "pass_live_api_scaleout_feasibility":
        raise ValueError("metadata scale-out did not pass its frozen source gate")
    parse_false(
        feasibility.get("candidate_image_pixels_opened"),
        field="candidate_image_pixels_opened",
    )
    parse_false(feasibility.get("continuous_colour_used"), field="continuous_colour_used")
    if len(panels) != 200:
        raise ValueError("atlas panels must contain exactly 200 species")

    taxa: list[str] = []
    counts: dict[str, int] = {}
    for row in panels:
        taxon = str(row.get("taxon_id", "")).strip()
        cohort = str(row.get("cohort_id", "")).strip()
        if not taxon or not cohort:
            raise ValueError("panel row lacks taxon_id or cohort_id")
        parse_false(
            row.get("candidate_image_pixels_opened", False),
            field="panel candidate_image_pixels_opened",
        )
        taxa.append(taxon)
        counts[cohort] = counts.get(cohort, 0) + 1
    expected = {f"C{index:02d}": 25 for index in range(1, 9)}
    if len(set(taxa)) != 200 or counts != expected:
        raise ValueError("atlas panels changed their disjoint 8 x 25 denominator")

    by_taxon: dict[str, Mapping[str, Any]] = {}
    for row in feasibility.get("species_results", ()):
        taxon = str(row.get("taxon_id", "")).strip()
        if not taxon or taxon in by_taxon:
            raise ValueError("metadata audit has a missing or duplicate taxon ID")
        by_taxon[taxon] = row
    missing = set(taxa) - set(by_taxon)
    if missing:
        raise ValueError(f"selected taxa lack metadata audit evidence: {sorted(missing)}")

    geometry: list[dict[str, Any]] = []
    for taxon in taxa:
        row = by_taxon[taxon]
        if row.get("status") != "geometry_eligible":
            raise ValueError(f"selected taxon {taxon} was not geometry eligible")
        scale_results = row.get("geometry_scale_results")
        if not isinstance(scale_results, list):
            raise ValueError(f"selected taxon {taxon} lacks geometry scale results")
        geometry.append({"taxon_id": taxon, "scale_results": scale_results})
    return taxa, geometry


def load_frozen_boundaries(
    environment_dir: Path, manifest: Mapping[str, Any]
) -> tuple[dict[int, list[dict[str, str]]], dict[str, str]]:
    if (
        manifest.get("status") != "pass_precolour_environmental_boundary_freeze"
        or manifest.get("scaleout_colour_opened") is not False
        or manifest.get("environment_colour_join_performed") is not False
    ):
        raise ValueError("environmental boundary freeze did not pass before colour")
    expected_hashes = manifest.get("output_sha256", {})
    rows_by_scale: dict[int, list[dict[str, str]]] = {}
    hashes: dict[str, str] = {}
    for scale in SCALES:
        name = f"environmental_boundary_cells_{scale}km.csv"
        path = environment_dir / name
        digest = file_sha256(path)
        if digest != expected_hashes.get(name):
            raise ValueError(f"environmental boundary hash changed: {name}")
        rows = read_csv(path)
        if not rows or any(int(row["scale_km"]) != scale for row in rows):
            raise ValueError(f"environmental boundary scale changed: {name}")
        rows_by_scale[scale] = rows
        hashes[name] = digest
    return rows_by_scale, hashes


def build_evidence(
    *,
    feasibility: Mapping[str, Any],
    panels: Sequence[Mapping[str, Any]],
    boundary_rows_by_scale: Mapping[int, Sequence[Mapping[str, Any]]],
    environment_contract: Mapping[str, Any],
    source_stage: str,
    dated_source_reconciliation: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    if source_stage not in SOURCE_STAGES:
        raise ValueError("unknown source stage")
    if source_stage == "live-feasibility":
        if dated_source_reconciliation is not None:
            raise ValueError("live feasibility must not consume dated-source outcomes")
    else:
        dated = dated_source_reconciliation
        if not isinstance(dated, Mapping):
            raise ValueError("final coverage requires dated-source reconciliation")
        if (
            dated.get("status") != "pass_dated_source_m2m_scaleout_freeze"
            or dated.get("candidate_image_pixels_opened") is not False
            or dated.get("continuous_colour_used") is not False
            or dated.get("selected_species") != 200
            or dated.get("selected_photo_assets") != 60000
            or dated.get("frozen_observations") != 60000
            or dated.get("replacement_permitted") is not False
            or dated.get("image_acquisition_authorized") is not False
        ):
            raise ValueError("dated-source reconciliation did not pass unchanged")
    taxa, geometry = selected_geometry(feasibility, panels)
    gate = evaluate_environmental_coverage_gate(
        geometry,
        taxa,
        boundary_rows_by_scale,
        environment_contract,
    )
    passed = gate["status"] == "pass_precolour_environmental_coverage"
    if source_stage == "live-feasibility":
        status = (
            "pass_live_api_precolour_environmental_coverage_feasibility"
            if passed
            else "not_evaluable_live_api_precolour_environmental_coverage"
        )
    else:
        status = gate["status"]
    return {
        **gate,
        "status": status,
        "coverage_gate_status": gate["status"],
        "source_stage": source_stage,
        "final_dated_source_required": source_stage != "final-dated-source",
        "image_acquisition_authorized": False,
        "selected_taxon_ids": taxa,
        "claim_ceiling": (
            "Pre-colour opportunity-cell coverage only; no image pixels, flower colour, "
            "environment-colour concordance, or biological result."
        ),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata-feasibility", type=Path, required=True)
    parser.add_argument("--species-panels", type=Path, required=True)
    parser.add_argument(
        "--environment-dir", type=Path, default=ROOT / "data/atlas/environment"
    )
    parser.add_argument(
        "--environment-contract",
        type=Path,
        default=ROOT / "docs/supporting/jbi_atlas_environmental_overlay_contract_v1.json",
    )
    parser.add_argument(
        "--boundary-manifest",
        type=Path,
        default=ROOT / "data/atlas/environment/environmental_boundary_freeze_manifest.json",
    )
    parser.add_argument("--source-stage", choices=SOURCE_STAGES, required=True)
    parser.add_argument("--dated-source-reconciliation", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    feasibility = json.loads(args.metadata_feasibility.read_text(encoding="utf-8"))
    panels = read_csv(args.species_panels)
    contract = json.loads(args.environment_contract.read_text(encoding="utf-8"))
    boundary_manifest = json.loads(args.boundary_manifest.read_text(encoding="utf-8"))
    dated_source_reconciliation = None
    if args.dated_source_reconciliation is not None:
        dated_source_reconciliation = json.loads(
            args.dated_source_reconciliation.read_text(encoding="utf-8")
        )
    rows, boundary_hashes = load_frozen_boundaries(
        args.environment_dir, boundary_manifest
    )
    result = build_evidence(
        feasibility=feasibility,
        panels=panels,
        boundary_rows_by_scale=rows,
        environment_contract=contract,
        source_stage=args.source_stage,
        dated_source_reconciliation=dated_source_reconciliation,
    )
    result["parents"] = {
        "metadata_feasibility_sha256": file_sha256(args.metadata_feasibility),
        "species_panels_sha256": file_sha256(args.species_panels),
        "environment_contract_sha256": file_sha256(args.environment_contract),
        "boundary_manifest_sha256": file_sha256(args.boundary_manifest),
        "boundary_files_sha256": boundary_hashes,
    }
    if args.dated_source_reconciliation is not None:
        result["parents"]["dated_source_reconciliation_sha256"] = file_sha256(
            args.dated_source_reconciliation
        )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["coverage_gate_status"] == "pass_precolour_environmental_coverage" else 2


if __name__ == "__main__":
    raise SystemExit(main())

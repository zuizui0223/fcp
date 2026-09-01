#!/usr/bin/env python3
"""Resolve the exact 60k atlas through iNaturalist's photo-observation M:M table."""

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

from fcp_pipeline.atlas_dated_source_m2m import (
    reconcile_rows_m2m,
    scan_snapshot_m2m,
    validate_m2m_amendment,
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        return
    fields = list(rows[0])
    union = {key for row in rows for key in row}
    fields.extend(sorted(union - set(fields)))
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--snapshot-archive", type=Path, required=True)
    parser.add_argument("--metadata-feasibility", type=Path, required=True)
    parser.add_argument("--species-panels", type=Path, required=True)
    parser.add_argument("--selected-observations", type=Path, required=True)
    parser.add_argument(
        "--amendment",
        type=Path,
        default=ROOT
        / "docs/supporting/jbi_atlas_dated_source_m2m_amendment_v2.json",
    )
    parser.add_argument(
        "--snapshot-receipt",
        type=Path,
        default=ROOT
        / "docs/supporting/jbi_atlas_inaturalist_snapshot_receipt_v1.json",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    amendment = json.loads(args.amendment.read_text(encoding="utf-8"))
    validate_m2m_amendment(amendment)
    receipt = json.loads(args.snapshot_receipt.read_text(encoding="utf-8"))
    snapshot = amendment["immutable_parents"]["snapshot"]
    if (
        receipt.get("status") != "pass_exact_archive_identity_and_tar_integrity"
        or receipt.get("source", {}).get("sha256") != snapshot["sha256"]
        or args.snapshot_archive.stat().st_size != snapshot["content_length_bytes"]
    ):
        raise RuntimeError("dated snapshot receipt or archive identity changed")
    archive_sha = sha256(args.snapshot_archive)
    if archive_sha != snapshot["sha256"]:
        raise RuntimeError("dated snapshot archive SHA-256 changed")
    feasibility = json.loads(args.metadata_feasibility.read_text(encoding="utf-8"))
    if (
        feasibility.get("status") != "pass_live_api_scaleout_feasibility"
        or feasibility.get("candidate_image_pixels_opened") is not False
        or feasibility.get("continuous_colour_used") is not False
    ):
        raise RuntimeError("live metadata feasibility did not pass before images")
    panels = read_csv(args.species_panels)
    observations = read_csv(args.selected_observations)
    taxon_ids = {str(row["taxon_id"]) for row in panels}
    photo_ids = {str(row["photo_id"]) for row in observations}
    observer_ids = {str(row["observer_id"]) for row in observations}
    genus_ids = {str(row["inat_genus_id"]) for row in observations}
    technical_failure = None
    scanned: dict[str, Any] | None = None
    try:
        scanned = scan_snapshot_m2m(
            args.snapshot_archive,
            taxon_ids=taxon_ids,
            photo_ids=photo_ids,
            observer_ids=observer_ids,
            genus_ids=genus_ids,
        )
        audit, frozen_rows = reconcile_rows_m2m(
            panels, observations, scanned, amendment
        )
    except (OSError, UnicodeError, ValueError) as exc:
        technical_failure = f"{type(exc).__name__}: {exc}"
        frozen_rows = []
        audit = {
            "protocol": amendment["protocol"],
            "status": "not_evaluable_dated_source_m2m_reconciliation",
            "candidate_image_pixels_opened": False,
            "continuous_colour_used": False,
            "selected_species": len(taxon_ids),
            "selected_photo_assets": len(photo_ids),
            "frozen_observations": 0,
            "technical_failure": technical_failure,
            "replacement_permitted": False,
            "image_acquisition_authorized": False,
            "claim_ceiling": amendment["claim_ceiling"],
        }
    audit["snapshot"] = {
        **snapshot,
        "computed_sha256": archive_sha,
        "archive_members": scanned["members"] if scanned is not None else [],
        "second_pass_for_observations": (
            scanned["second_pass_for_observations"] if scanned is not None else None
        ),
    }
    audit["parents"] = {
        "amendment_sha256": sha256(args.amendment),
        "snapshot_receipt_sha256": sha256(args.snapshot_receipt),
        "metadata_feasibility_sha256": sha256(args.metadata_feasibility),
        "species_panels_sha256": sha256(args.species_panels),
        "selected_observations_sha256": sha256(args.selected_observations),
    }
    reconciliation_path = args.output_dir / "dated_source_m2m_reconciliation.json"
    write_json(reconciliation_path, audit)
    frozen_path = args.output_dir / "dated_source_m2m_observation_manifest.csv"
    if frozen_rows:
        write_csv(frozen_path, frozen_rows)
    manifest = {
        "protocol": amendment["protocol"],
        "status": audit["status"],
        "candidate_image_pixels_opened": False,
        "files": {reconciliation_path.name: sha256(reconciliation_path)},
    }
    if frozen_path.exists():
        manifest["files"][frozen_path.name] = sha256(frozen_path)
    manifest_path = args.output_dir / "dated_source_m2m_manifest.json"
    write_json(manifest_path, manifest)
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0 if audit["status"] == "pass_dated_source_m2m_scaleout_freeze" else 2


if __name__ == "__main__":
    raise SystemExit(main())

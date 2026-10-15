#!/usr/bin/env python3
"""Select ONLY literature-documented local flower-colour polymorphism taxa.

Canonical source: the immutable JBI v2.2 34-species strict single-pass
C/S evidence freeze. C=1 explicitly requires multiple natural *floral display*
colour morphs in the same population/site. S=1 alone is geographic colour
differentiation, NOT enough to establish local FCP.

This script NEVER scores FCP from an iNaturalist photo, a 42,111-taxon
opportunity frame, a climate association, or an old model output.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

SOURCE = "data/frozen/jbi_cs_evidence_freeze_v22.csv"
SOURCE_SHA = "fa4cd099b875531256964b101a8db7fbbe8a2598ad3dfdb1a83145af98c4a2bd"
SOURCE_MANIFEST = "data/frozen/jbi_cs_evidence_freeze_v22_manifest.json"
STRICT_STATES = {
    "local_coexistence_only": (1, 0),
    "spatial_segregation_only": (0, 1),
    "coexistence_and_segregation": (1, 1),
}
SELECTED_FIELDS = [
    "canonical_name", "family", "organization_state",
    "C_local_coexistence_documented", "S_spatial_segregation_documented",
    "n_resolved_sources", "n_C_positive_sources", "n_S_positive_sources",
    "n_conflict_sources", "C_source_ids", "S_source_ids",
    "selection_reason", "evidence_review_status",
]
EXCLUDED_FIELDS = [
    "canonical_name", "family", "organization_state",
    "C_local_coexistence_documented", "S_spatial_segregation_documented",
    "S_source_ids", "exclusion_reason",
]
STATUS = "DOCUMENTED_LOCAL_FCP_C_POSITIVE_ONLY_NOT_IMAGE_VALIDATED"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build(source_file: Path, source_manifest: Path) -> tuple[list[dict], list[dict], dict]:
    if sha256(source_file) != SOURCE_SHA:
        raise ValueError("Immutable original strict FCP source SHA256 mismatch")
    freeze = json.loads(source_manifest.read_text(encoding="utf-8"))
    if (freeze.get("status") != "frozen"
        or freeze.get("dataset") != SOURCE
        or freeze.get("dataset_sha256") != SOURCE_SHA
        or freeze.get("n_species") != 34
        or freeze.get("C_positive") != 26
        or freeze.get("state_counts") != {
            "local_coexistence_only": 11,
            "spatial_segregation_only": 8,
            "coexistence_and_segregation": 15,
        }
        or freeze.get("classification") !=
            "single-pass conservative documented-evidence audit"):
        raise ValueError("Original literature evidence freeze changed")
    selected: list[dict] = []
    excluded: list[dict] = []
    seen: set[str] = set()
    counts = {label: 0 for label in STRICT_STATES}
    with source_file.open("r", encoding="utf-8", newline="") as handle:
        rows = csv.DictReader(handle)
        required = {
            "canonical_name", "family", "organization_state",
            "C_local_coexistence_documented", "S_spatial_segregation_documented",
            "n_resolved_sources", "n_C_positive_sources", "n_S_positive_sources",
            "n_conflict_sources", "C_source_ids", "S_source_ids",
        }
        if set(rows.fieldnames or ()) != required:
            raise ValueError("Strict source evidence columns changed or were lost")
        for row in rows:
            name = row["canonical_name"].strip()
            state = row["organization_state"]
            if name in seen or not name or not row["family"].strip() or state not in STRICT_STATES:
                raise ValueError("Invalid duplicate/name/state in original FCP literature freeze")
            seen.add(name)
            counts[state] += 1
            c, s = STRICT_STATES[state]
            try:
                numeric = {
                    col: int(row[col]) for col in (
                        "C_local_coexistence_documented", "S_spatial_segregation_documented",
                        "n_resolved_sources", "n_C_positive_sources", "n_S_positive_sources",
                        "n_conflict_sources",
                    )
                }
            except (ValueError, TypeError) as err:
                raise ValueError("Noninteger original evidence flags/counts") from err
            if (numeric["C_local_coexistence_documented"] != c
                or numeric["S_spatial_segregation_documented"] != s
                or numeric["n_resolved_sources"] < 1
                or numeric["n_C_positive_sources"] < c
                or numeric["n_S_positive_sources"] < s
                or any(v < 0 for v in numeric.values())
                or (c and not row["C_source_ids"].strip())
                or (s and not row["S_source_ids"].strip())):
                raise ValueError("Original explicit source evidence is internally inconsistent")
            if c:
                out = {key: row[key] for key in SELECTED_FIELDS if key in row}
                out["selection_reason"] = "explicit_natural_within_population_discrete_floral_colour_morphs"
                out["evidence_review_status"] = "strict_single_pass_source_positive_no_independent_human_duplicate_review"
                selected.append(out)
            else:
                out = {key: row[key] for key in EXCLUDED_FIELDS if key in row}
                out["exclusion_reason"] = "spatial_difference_only_local_coexistence_not_documented"
                excluded.append(out)
    if (len(seen) != 34 or counts != freeze["state_counts"]
        or len(selected) != 26 or len(excluded) != 8
        or sum(int(x["S_spatial_segregation_documented"]) for x in selected) != 15):
        raise ValueError("FCP inclusion would change the original 26 confirmed-by-source C-positive species")
    return selected, excluded, {
        "schema": "fcp_strict_documented_local_coexistence_species_selection_v1",
        "status": STATUS,
        "source": SOURCE,
        "original_evidence_sha256": SOURCE_SHA,
        "source_workflow_run_id": freeze["source_workflow_run_id"],
        "selection_rule": "C_local_coexistence_documented == 1",
        "selected_species": len(selected),
        "local_coexistence_only": counts["local_coexistence_only"],
        "coexistence_and_segregation": counts["coexistence_and_segregation"],
        "excluded_S_only": len(excluded),
        "source_supported_not_independently_reviewer_verified": True,
        "not_source_photo_validation": True,
        "all_42111_species_are_not_polymorphic_cases": True,
        "previous_93209351_photo_ID_count_excluded_from_strict_FCP_measurement": True,
        "warning": ("Study estimates structure/heterogeneity conditional on documented FCP; "
                    "it cannot estimate prevalence among all plants, or compare FCP to "
                    "true biological absence without an independently validated negative group."),
    }


def write_csv(path: Path, fields: list[str], rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, default=Path(SOURCE))
    p.add_argument("--source-manifest", type=Path, default=Path(SOURCE_MANIFEST))
    p.add_argument("--output-dir", required=True, type=Path)
    args = p.parse_args(argv)
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise ValueError("Never overwrite an existing polymorphism-evidence cohort")
    keep, drop, result = build(args.source, args.source_manifest)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    selected = args.output_dir / "strict_local_fcp_26_source_documented.csv"
    excluded = args.output_dir / "geographic_only_8_excluded_from_strict_fcp.csv"
    write_csv(selected, SELECTED_FIELDS, keep)
    write_csv(excluded, EXCLUDED_FIELDS, drop)
    result["selected_csv_sha256"] = sha256(selected)
    result["geographic_only_exclusion_csv_sha256"] = sha256(excluded)
    (args.output_dir / "cohort_receipt.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")
    print("FCP_STRICT_26_SOURCE_DOCUMENTED", json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

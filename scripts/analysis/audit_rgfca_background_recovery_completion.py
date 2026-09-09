"""Audit the completed discovery-background artifacts without colour inference.

Consumes only run 34094784607's closed artifacts and literal discovery Git
objects. Never imports an image model, downloads images, or reads reserve data.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import io
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[2]
EXECUTION = "3ba0c11a66b1b0c6aef09483ebc7eb0666f244a3"
RESULT_COMMIT = "2f00847ccc4d15be8637a1e3589239ce2ab5c114"
RESULT_PATH = "docs/supporting/global_rgfca_background_falsification_result_v2.json"
FLAGS = (
    "exact_image_sha_match", "roi_admission_reproduced",
    "flower_mask_pixels_exact", "flower_palette_exact", "background_pixels_exact",
)
EXACT = "exact_matched_background_recovered"
FAILED = "original_roi_or_flower_palette_not_exactly_reproduced"


def git_bytes(revision: str, path: str) -> bytes:
    return subprocess.check_output(["git", "show", f"{revision}:{path}"], cwd=ROOT)


def csv_rows(raw: bytes) -> list[dict[str, str]]:
    return list(csv.DictReader(io.StringIO(raw.decode("utf-8"))))


def audit(full_dir: Path, partitions_dir: Path) -> dict:
    summary_bytes = (full_dir / "summary.json").read_bytes()
    assert summary_bytes == git_bytes(RESULT_COMMIT, RESULT_PATH)
    summary = json.loads(summary_bytes)
    assert summary["status"] == "not_evaluable_incomplete_exact_background_recovery"
    for key in ("primary_statistic_computed", "replacement_photos_used",
                "denominator_adapted_after_recovery"):
        assert summary[key] is False, key
    assert {p.name for p in full_dir.iterdir()} == {
        "summary.json", "background_recovery_rows_v2.csv"
    }, "Unexpected inference output in a non-evaluable artifact"
    lineage = summary["lineage"]
    for key, path in (
        ("contract_sha256", "docs/supporting/global_rgfca_within_species_spatial_background_control_contract_v2.json"),
        ("amendment_sha256", "docs/supporting/global_rgfca_background_control_postoutcome_status_amendment_v2a.json"),
        ("measured_sha256", "data/derived/global_monte_carlo_measured_photos_v1.csv"),
    ):
        assert hashlib.sha256(git_bytes(EXECUTION, path)).hexdigest() == lineage[key]
    recovery_raw = (full_dir / "background_recovery_rows_v2.csv").read_bytes()
    assert hashlib.sha256(recovery_raw).hexdigest() == lineage["recovery_rows_sha256"]
    rows = csv_rows(recovery_raw)
    files = sorted(partitions_dir.rglob("recovery_s*_p*.csv"), key=lambda p: p.name)
    expected_names = {f"recovery_s{s:02d}_p{p:02d}.csv" for s in range(32) for p in range(4)}
    assert len(files) == 128 and {p.name for p in files} == expected_names
    concatenated = []
    for path in files:
        part = csv_rows(path.read_bytes())
        manifest = json.loads(path.with_suffix(".json").read_bytes())
        counts = Counter(row["recovery_status"] for row in part)
        assert manifest["status"] == "complete_background_control_partition_recovery"
        assert manifest["selected_rows"] == len(part)
        for key in (EXACT, FAILED, "acquisition_failed_no_replacement", "image_sha_mismatch_no_replacement"):
            assert manifest[key] == counts[key], (path.name, key)
        assert manifest["roi_runtime_or_admission_failure"] == 0
        for key in ("species_opened_to_worker", "coordinates_opened_to_worker", "replacement_photos_used"):
            assert manifest[key] is False, (path.name, key)
        concatenated.extend(part)
    assert concatenated == rows, "Partition rows differ from the finalized census"
    assert len(rows) == len({row["measurement_id"] for row in rows}) == 21424
    measured = csv_rows(git_bytes(EXECUTION, "data/derived/global_monte_carlo_measured_photos_v1.csv"))
    classifiable = [r for r in measured if r["global_classifiable"].strip().casefold() in {"true", "1"}]
    species_counts = Counter(r["species"] for r in classifiable)
    frame = [r for r in classifiable if species_counts[r["species"]] >= 40]
    assert len(frame) == 21424 and len({r["species"] for r in frame}) == 369
    assert {r["measurement_id"] for r in frame} == {r["measurement_id"] for r in rows}
    assert dict(Counter(r["recovery_status"] for r in rows)) == summary["status_counts"]
    failures = Counter()
    combinations = Counter()
    for row in rows:
        assert all(row[key] in {"True", "False"} for key in FLAGS)
        bad = tuple(key for key in FLAGS if row[key] == "False")
        assert row["recovery_status"] == (FAILED if bad else EXACT)
        failures.update(bad)
        if bad:
            combinations[" + ".join(bad)] += 1
    assert summary["exact_rows"] == 21339 and summary["failed_rows"] == 85
    assert sum(combinations.values()) == 85
    return {
        "status": "verified_completed_non_evaluable_discovery_background_recovery",
        "run_id": 34094784607, "execution_head": EXECUTION,
        "result_commit": RESULT_COMMIT,
        "partitions_verified": 128, "rows_verified": len(rows), "species_in_fixed_frame": 369,
        "exact_rows": summary["exact_rows"], "failed_rows": summary["failed_rows"],
        "failed_flags": {key: failures[key] for key in FLAGS},
        "failure_combinations": dict(combinations),
        "lineage": lineage,
        "primary_statistic_computed": False, "reserve_outcomes_read": False,
        "raw_images_reacquired_by_audit": False,
        "scope": "Saved worker flags and exact artifact/ID census; not an independent rerun of image measurement or a diagnosis of the root cause.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full-dir", type=Path, required=True)
    parser.add_argument("--partitions-dir", type=Path, required=True)
    parser.add_argument("--output-json", type=Path)
    args = parser.parse_args()
    result = json.dumps(audit(args.full_dir, args.partitions_dir), indent=2) + "\n"
    if args.output_json:
        args.output_json.write_text(result, encoding="utf-8", newline="\n")
    print(result, end="")


if __name__ == "__main__":
    main()

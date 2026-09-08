"""Verify the closed metadata-only run from saved bytes; no provider requests."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

from scripts.analysis.audit_rgfca_flowermask_pair_metadata import (
    ROOT, require, validate_dispatch_history, verify_saved_result,
)


RUN = 34201981494
HEAD = "49d24a1b3ca060608d629577a01ce337b2bb45fe"
PLAN = "docs/supporting/rgfca_flowermask_pair_metadata_plan_v1.json"
PLAN_SHA = "2ce2a578db6c470a9021189241e2e5c9fae561c2435896f64e61d50e0ccd4d83"
INNER = "rgfca-flowermask-pair-metadata-v1"
RECEIPT_HASHES = {
    f"{INNER}/result.json": "27c79f60b4088c4c2942778c43c9089509ac3dcce5ec825166d97804166359b3",
    f"{INNER}/execution_start.json": "b133ba62b29a0ed5539b28ed54218b5e74da0069d669aec8598b7f101aa08793",
    "flowermask-pair-dispatch-history.json": "cc1ba4ef73ec7f7a57feffe4907cc591ab0782f57e94040436993e70999f3da0",
}


def verify(bundle: Path, root: Path = ROOT):
    plan_raw = (root / PLAN).read_bytes()
    require(hashlib.sha256(plan_raw).hexdigest() == PLAN_SHA, "frozen plan bytes")
    plan = json.loads(plan_raw)
    expected_paths = set(RECEIPT_HASHES) | {
        f"{INNER}/responses/{row['annotation_id']}.json" for row in plan["rows"]
    }
    actual_paths = {p.relative_to(bundle).as_posix() for p in bundle.rglob("*") if p.is_file()}
    require(actual_paths == expected_paths and len(actual_paths) == 303, "exact 303-file census")
    for relative, expected in RECEIPT_HASHES.items():
        require(hashlib.sha256((bundle / relative).read_bytes()).hexdigest() == expected,
                f"immutable receipt drift: {relative}")
    start = json.loads((bundle / INNER / "execution_start.json").read_bytes())
    require(start == {
        "head_sha": HEAD, "run_id": str(RUN), "started_utc": "2026-09-08T07:58:31Z",
        "canonical_plan_sha256": PLAN_SHA, "all_rows_checkpointed_before_requests": True,
    }, "execution identity")
    history = json.loads((bundle / "flowermask-pair-dispatch-history.json").read_bytes())
    validate_dispatch_history(history, RUN, HEAD, 1)
    result = verify_saved_result(plan, bundle / INNER, root)
    require(result["status_counts"] == {"metadata_inspected": 300}, "terminal census")
    require(result["new_requests_started"] == 299, "request census")
    unresolved = []
    group_counts = Counter()
    for row in result["terminal"]:
        group_counts[row["group"]] += 1
        if not row["summary"]["required_pair_metadata_present"]:
            names = {item["filename"] for item in row["summary"]["files"]}
            unresolved.append({
                "annotation_id": row["annotation_id"], "group": row["group"],
                "annotation_filename": row["annotation_filename"],
                "candidate_folder_id": row["candidate_folder_id"],
                "candidate_folder_name": row["candidate_folder_name"],
                "missing_required_filenames": sorted(set(plan["rules"]["required_pair_filenames"]) - names),
                "observed_filenames": sorted(names),
            })
    require(result["pair_metadata_present_count"] == 299 and len(unresolved) == 1,
            "closed pairing outcome")
    require(unresolved[0]["annotation_id"] == "bf1604ee-fe76-4ec2-b2b7-51c50a732090"
            and unresolved[0]["missing_required_filenames"] == ["img.png"], "unresolved identity")
    return {
        "schema": "rgfca-flowermask-pair-metadata-retained-summary-v1",
        "run_id": RUN, "execution_head_sha": HEAD, "plan_sha256": PLAN_SHA,
        "result_sha256": RECEIPT_HASHES[f"{INNER}/result.json"],
        "artifact_id": 10046355184,
        "provider_artifact_zip_sha256": "a26dabf38c4a51274303dad0191fb2284df65b0806832f9c988210018df950e6",
        "provider_zip_digest_locally_verified": False,
        "retained_files_verified": 303, "metadata_rows_inspected": 300,
        "new_requests_started": 299, "cached_responses_reused": 1,
        "prior_schema_opened_rows_retained": sum(r["prior_schema_opened"] for r in result["terminal"]),
        "source_group_row_counts": dict(sorted(group_counts.items())),
        "pair_metadata_present_count": 299, "all_300_pair_metadata_present": False,
        "unresolved": unresolved, "status": "completed_metadata_only_pairing_incomplete",
        "acquisition_source_failures": 0, "retained_verification_passed": True,
        "image_pixels_decoded": False, "model_run": False, "coordinates_joined": False,
        "benchmark_execution_authorized": False,
        "claim_ceiling": "Documentary metadata only; one returned listing lacks img.png; no subset benchmark.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bundle", type=Path)
    parser.add_argument("--verify-summary", type=Path)
    args = parser.parse_args()
    out = verify(args.bundle)
    if args.verify_summary:
        require(json.loads(args.verify_summary.read_bytes()) == out, "retained summary drift")
    print(json.dumps(out, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

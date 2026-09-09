"""One bounded, metadata-only FlowerMask paired-manifest inspection.

Never follows any download/view URL in returned metadata. No image, mask,
annotation body, model or geographic coordinate is requested or interpreted.
"""
from __future__ import annotations

import argparse
from collections import Counter
from copy import deepcopy
import hashlib
import json
import os
from pathlib import Path
import re
import time
from urllib.request import HTTPRedirectHandler, Request, build_opener

from scripts.analysis.audit_rgfca_reference_schema import validate_plan as validate_frame


ROOT = Path(__file__).resolve().parents[2]
FRAME = "docs/supporting/rgfca_reference_schema_plan_v1.json"
INPUTS = "data/validation/flowermask_pair_metadata_inputs_v1"
SOURCE_HASHES = {
    FRAME: "48436a74160b4700e2d9b2bbc803aeac7fcfe8a1f09885140a1eade034c5aade",
    f"{INPUTS}/folders.json": "fd30ad1ba248c43eca04c1a53cec31e87617dbd815b5a98b752eb461dec2b892",
    f"{INPUTS}/cached_folder_files.json": "a227790321584118998d390850cca1b66ab4f7ac54276ebc0c0262f566c77485",
}
CACHED_FOLDER = "e677f918-a695-4472-9a03-39e6827d007a"
BASE = "https://data.mendeley.com/public-api/datasets/3pw57gdcj2/files?folder_id="
UUID = r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}"
MAX_RESPONSE_BYTES = 250_000
RULES = {
    "maximum_new_requests": 299, "maximum_response_bytes": MAX_RESPONSE_BYTES,
    "request_timeout_seconds": 30, "minimum_request_start_interval_seconds": 1,
    "retries": 0, "redirects": 0, "pagination_followups": 0,
    "stop_on_transport_or_schema_failure": True,
    "required_pair_filenames": ["img.png", "label.png", "label_names.txt"],
    "benchmark_execution_authorized": False,
}


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_sources(root=ROOT):
    raw = {name: (root / name).read_bytes() for name in SOURCE_HASHES}
    for name, expected in SOURCE_HASHES.items():
        require(sha(raw[name]) == expected, f"immutable source drift: {name}")
    return raw


def derive_plan(frame, folders):
    validate_frame(frame)
    by_id = {folder["id"]: folder for folder in folders}
    require(len(by_id) == len(folders) == 319, "folder census or duplicate IDs")
    rows = []
    for annotation in frame["public_annotation_files"]:
        parent = by_id[annotation["folder_id"]]["parent_id"]
        require(by_id[parent]["name"] == annotation["group"], "provider group mismatch")
        containers = [folder for folder in folders if folder.get("parent_id") == parent
                      and re.fullmatch(r".+_masking(?:_files)?", folder["name"], re.I)]
        require(len(containers) == 1, "ambiguous masking container")
        stem = Path(annotation["filename"]).stem
        matches = [folder for folder in folders if folder.get("parent_id") == containers[0]["id"]
                   and folder["name"] == stem + "_dataset"]
        require(len(matches) == 1, "ambiguous filename-based correspondence")
        folder = matches[0]
        require(re.fullmatch(UUID, folder["id"]), "invalid folder ID")
        rows.append({
            "annotation_id": annotation["id"], "group": annotation["group"],
            "annotation_filename": annotation["filename"],
            "candidate_folder_id": folder["id"], "candidate_folder_name": folder["name"],
            "prior_schema_opened": annotation["id"] in frame["selected_ids"],
            "source": "cached" if folder["id"] == CACHED_FOLDER else "request",
            "metadata_url": BASE + folder["id"] + "&version=1",
        })
    require(len({r["candidate_folder_id"] for r in rows}) == 300, "nonunique folder mapping")
    require(Counter(r["source"] for r in rows) == {"cached": 1, "request": 299}, "request census")
    return {"schema": "rgfca-flowermask-pair-metadata-plan-v1", "source_hashes": dict(SOURCE_HASHES),
            "rules": deepcopy(RULES), "rows": rows,
            "scope": "300_previously_exposed_annotation_ids_not_complete_3600_dataset"}


def build_plan(root=ROOT):
    raw = read_sources(root)
    return derive_plan(json.loads(raw[FRAME]), json.loads(raw[f"{INPUTS}/folders.json"]))


def validate_plan(plan, root=ROOT):
    require(plan == build_plan(root), "plan differs from exact source-derived frame/rules")
    return plan["rows"]


def summarize_metadata(raw, row):
    require(0 < len(raw) <= MAX_RESPONSE_BYTES, "metadata response size")
    entries = json.loads(raw)
    require(isinstance(entries, list) and len(entries) <= 50, "metadata list schema/cap")
    files = []
    for entry in entries:
        require(isinstance(entry, dict), "metadata entry type")
        # A closed metadata grammar rejects an annotation/imageData payload.
        require(set(entry) <= {"filename", "id", "folder_id", "content_details", "size",
                               "last_modified_date", "status"}, "unexpected metadata fields")
        detail = entry["content_details"]
        require(isinstance(detail, dict), "content details type")
        require(set(detail) <= {"id", "sha256_hash", "content_type", "size", "created_date",
                               "download_url", "view_url", "download_expiry_time"},
                "unexpected content-detail fields")
        require(entry["folder_id"] == row["candidate_folder_id"], "wrong response folder")
        require(re.fullmatch(UUID, entry["id"]), "file ID")
        require(entry["status"] == "COMPLETED", "incomplete provider file")
        require(isinstance(entry["filename"], str) and entry["filename"]
                and len(entry["filename"]) <= 250, "filename")
        require(type(entry["size"]) is int and entry["size"] > 0
                and type(detail["size"]) is int and detail["size"] == entry["size"], "size disagreement")
        require(re.fullmatch(r"[0-9a-f]{64}", detail["sha256_hash"]), "provider hash")
        require(isinstance(detail["content_type"], str), "content type")
        files.append({"filename": entry["filename"], "provider_file_id": entry["id"],
                      "advertised_bytes": entry["size"], "advertised_sha256": detail["sha256_hash"],
                      "content_type": detail["content_type"]})
    require(len({f["provider_file_id"] for f in files}) == len(files), "duplicate file ID")
    require(len({f["filename"] for f in files}) == len(files), "duplicate exact filename")
    require(len({f["filename"].casefold() for f in files}) == len(files), "casefold filename collision")
    by_name = {f["filename"]: f for f in files}
    expected_types = {"img.png": "image/png", "label.png": "image/png", "label_names.txt": "text/plain"}
    pair_present = all(name in by_name and by_name[name]["content_type"] == kind
                       for name, kind in expected_types.items())
    return {"files": sorted(files, key=lambda f: f["filename"]),
            "file_count": len(files), "required_pair_metadata_present": pair_present,
            "pagination_completeness_verified": False,
            "advertised_payload_hashes_verified": False}


class NoRedirect(HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        raise ValueError("metadata redirects prohibited; no alternate endpoint")


def fetch_metadata(row):
    expected = BASE + row["candidate_folder_id"] + "&version=1"
    require(row["metadata_url"] == expected
            and re.fullmatch(UUID, row["candidate_folder_id"]), "non-metadata request")
    request = Request(expected, headers={"User-Agent": "FCP-reference-metadata-audit/1",
                                         "Accept": "application/json"})
    with build_opener(NoRedirect()).open(request, timeout=RULES["request_timeout_seconds"]) as response:
        require(response.status == 200 and response.url == expected, "unexpected response/redirect")
        require(response.headers.get_content_type() == "application/json", "non-JSON metadata response")
        return response.read(MAX_RESPONSE_BYTES + 1)


def write_json(path, value, *, exclusive=False):
    text = json.dumps(value, indent=2, ensure_ascii=True) + "\n"
    target = path if exclusive else path.with_suffix(path.suffix + ".tmp")
    with target.open("x" if exclusive else "w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)
        handle.flush()
        os.fsync(handle.fileno())
    if not exclusive:
        os.replace(target, path)


def run_inspection(plan, output, root=ROOT, fetch=fetch_metadata, pause=time.sleep, clock=time.monotonic):
    rows = validate_plan(plan, root)
    cached = read_sources(root)[f"{INPUTS}/cached_folder_files.json"]
    # Validate the existing response before any new request or output claim.
    summarize_metadata(cached, next(row for row in rows if row["source"] == "cached"))
    output.mkdir(parents=True, exist_ok=False)
    (output / "responses").mkdir()
    result = {"schema": "rgfca-flowermask-pair-metadata-result-v1", "status": "started",
              "terminal": [{**row, "status": "not_attempted"} for row in rows],
              "new_requests_started": 0, "image_pixels_decoded": False, "model_run": False,
              "coordinates_joined": False, "benchmark_execution_authorized": False}
    result_path = output / "result.json"
    write_json(result_path, result, exclusive=True)
    write_json(output / "execution_start.json", {"head_sha": os.environ.get("GITHUB_SHA"),
               "run_id": os.environ.get("GITHUB_RUN_ID"), "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
               "canonical_plan_sha256": sha((json.dumps(plan, indent=2, ensure_ascii=True) + "\n").encode()),
               "all_rows_checkpointed_before_requests": True}, exclusive=True)
    last_start, stopped = None, False
    for target in result["terminal"]:
        if stopped:
            target["status"] = "not_attempted_after_source_failure"
            continue
        raw = None
        try:
            if target["source"] == "cached":
                raw = cached
            else:
                if last_start is not None:
                    pause(max(0.0, RULES["minimum_request_start_interval_seconds"] - (clock() - last_start)))
                result["new_requests_started"] += 1
                require(result["new_requests_started"] <= RULES["maximum_new_requests"], "request cap")
                target["status"] = "request_started"
                write_json(result_path, result)
                last_start = clock()
                raw = fetch(target)
            summary = summarize_metadata(raw, target)
            # Save valid metadata bytes only. Malformed/opaque payloads are not retained.
            with (output / "responses" / (target["annotation_id"] + ".json")).open("xb") as handle:
                handle.write(raw)
                handle.flush()
                os.fsync(handle.fileno())
            target.update(status="metadata_inspected", summary=summary,
                          response_sha256=sha(raw), response_bytes=len(raw))
        except Exception as exc:
            target.update(status="source_inspection_failed", error_type=type(exc).__name__)
            if raw is not None:
                target.update(response_sha256=sha(raw), response_bytes=len(raw))
            stopped = True
        write_json(result_path, result)
    result["status"] = "stopped_source_failure" if stopped else "completed_metadata_only"
    result["status_counts"] = dict(Counter(row["status"] for row in result["terminal"]))
    result["pair_metadata_present_count"] = sum(row.get("summary", {}).get("required_pair_metadata_present", False)
                                                 for row in result["terminal"])
    result["all_300_pair_metadata_present"] = result["pair_metadata_present_count"] == 300
    write_json(result_path, result)
    return result


def verify_saved_result(plan, output, root=ROOT):
    rows = validate_plan(plan, root)
    result = json.loads((output / "result.json").read_bytes())
    require(len(result["terminal"]) == 300, "incomplete row census")
    stopped = False
    requests = 0
    for expected, target in zip(rows, result["terminal"]):
        require(all(target.get(k) == v for k, v in expected.items()), "row identity drift")
        status = target["status"]
        if stopped:
            require(status == "not_attempted_after_source_failure", "work after source failure")
        else:
            require(status in {"metadata_inspected", "source_inspection_failed"}, "unfinished/unknown terminal")
            requests += target["source"] == "request"
        if status != "metadata_inspected":
            require("summary" not in target, "uninspected row has summary")
        if status == "source_inspection_failed":
            require(not stopped and isinstance(target.get("error_type"), str), "failure census")
            stopped = True
        if target["status"] == "metadata_inspected":
            raw = (output / "responses" / (target["annotation_id"] + ".json")).read_bytes()
            require(sha(raw) == target["response_sha256"] and len(raw) == target["response_bytes"], "response drift")
            require(summarize_metadata(raw, target) == target["summary"], "summary drift")
            if target["source"] == "cached":
                require(sha(raw) == SOURCE_HASHES[f"{INPUTS}/cached_folder_files.json"], "cached source drift")
    require(result["status"] == ("stopped_source_failure" if stopped else "completed_metadata_only"), "completion status")
    require(type(result["new_requests_started"]) is int and result["new_requests_started"] == requests, "request census")
    counts = dict(Counter(row["status"] for row in result["terminal"]))
    require(counts == result["status_counts"], "status counts")
    pair_count = sum(row.get("summary", {}).get("required_pair_metadata_present", False) for row in result["terminal"])
    require(pair_count == result["pair_metadata_present_count"], "pair count")
    require((pair_count == 300) is result["all_300_pair_metadata_present"], "pair gate")
    for key in ("image_pixels_decoded", "model_run", "coordinates_joined", "benchmark_execution_authorized"):
        require(result[key] is False, "unsupported claim escalation")
    return result


def validate_dispatch_history(history, run_id, head_sha, run_attempt):
    """Single dispatch across this v1 workflow; fail closed on retries/history gaps."""
    runs = history["workflow_runs"]
    require(run_attempt == 1, "workflow rerun prohibited")
    require(history["total_count"] == len(runs) == 1, "prior/concurrent dispatch or incomplete history")
    run = runs[0]
    require(run["id"] == run_id and run["head_sha"] == head_sha, "dispatch identity mismatch")
    require(run["event"] == "workflow_dispatch" and run["run_attempt"] == 1, "dispatch attempt mismatch")
    require(run["path"].split("@")[0] == ".github/workflows/rgfca-flowermask-pair-metadata.yml", "wrong workflow history")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["plan", "preflight", "inspect", "verify", "history"])
    parser.add_argument("--plan", type=Path)
    parser.add_argument("--expected-plan-sha256")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--history", type=Path)
    args = parser.parse_args()
    if args.command == "history":
        validate_dispatch_history(json.loads(args.history.read_bytes()), int(os.environ["GITHUB_RUN_ID"]),
                                  os.environ["GITHUB_SHA"], int(os.environ["GITHUB_RUN_ATTEMPT"]))
        print("Verified one first-attempt metadata-only dispatch; no prior/concurrent run.")
        return
    require(args.plan is not None, "plan required")
    if args.command == "plan":
        write_json(args.plan, build_plan(), exclusive=True)
        print(sha(args.plan.read_bytes()))
        return
    raw = args.plan.read_bytes()
    require(args.expected_plan_sha256 and sha(raw) == args.expected_plan_sha256, "plan byte hash")
    plan = json.loads(raw)
    validate_plan(plan)
    if args.command == "preflight":
        print(json.dumps({"status": "qualified_metadata_only_plan", "rows": 300, "cached": 1, "new_requests": 299,
                          "image_pixels_decoded": False, "requests_performed": 0}))
        return
    require(args.output is not None, "output directory required")
    result = run_inspection(plan, args.output) if args.command == "inspect" else verify_saved_result(plan, args.output)
    print(json.dumps({key: value for key, value in result.items() if key != "terminal"}, sort_keys=True))
    if result["status"] != "completed_metadata_only":
        raise SystemExit("Source inspection stopped; preserve every row; no retry or replacement.")


if __name__ == "__main__":
    main()

"""Bounded FlowerMask annotation-schema inspection, never image evaluation.

Build a plan from public metadata receipts before acquiring annotation bodies.
Choose the smallest provider file ID in each of six groups (not by size/colour).
Documents may contain opaque imageData: transfer/hash it, but never decode it.
No image URL/path is followed, no model is loaded, no coordinates are joined.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import time
from urllib.request import Request, urlopen


GROUPS = (
    "Butterfly Pea", "Caesalpinia Pulcherrima", "Jatropha Integerimma",
    "Plumeria", "Rose", "Tecoma Stans",
)
BASE = "https://data.mendeley.com/public-files/datasets/3pw57gdcj2/files/"
MAX_BYTES = 10_000_000


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def validate_plan(plan):
    require(plan["schema"] == "rgfca-reference-schema-plan-v1", "plan schema")
    require(plan["dataset_version"] == "10.17632/3pw57gdcj2.1", "dataset version")
    rows = plan["public_annotation_files"]
    require(len(rows) == 300, "public inventory must contain all 300 entries")
    require(Counter(r["group"] for r in rows) == Counter({g: 50 for g in GROUPS}),
            "six provider groups, 50 files each; preserve names verbatim")
    require(len({r["id"] for r in rows}) == 300, "duplicate file ID")
    require(len({(r["group"], r["filename"]) for r in rows}) == 300,
            "duplicate group/filename")
    for row in rows:
        require(re.fullmatch(r"[0-9a-f-]{36}", row["id"]) is not None, "file ID")
        require(re.fullmatch(r"[0-9a-f]{64}", row["sha256"]) is not None, "hash")
        require(type(row["bytes"]) is int and 0 < row["bytes"] <= MAX_BYTES,
                "declared size exceeds inspection bound")
        require(row["filename"].endswith(".json"), "not an annotation JSON")
        require(row["url"] == BASE + row["id"] + "/file_downloaded", "source URL")
    chosen = [min(r["id"] for r in rows if r["group"] == g) for g in GROUPS]
    require(plan["selected_ids"] == chosen, "fixed smallest-ID selection changed")
    require(plan["selection"] == "smallest_provider_file_id_per_group", "selection")
    return [next(r for r in rows if r["id"] == key) for key in chosen]


def build_plan(receipts):
    folder_path = receipts / "flowermask-v1-public-folders.json"
    folders = {r["id"]: r for r in json.loads(folder_path.read_bytes())}
    sources = {folder_path.name: digest(folder_path.read_bytes())}
    rows = []
    for path in sorted(receipts.glob("flowermask-v1-json-files-*.json")):
        raw = path.read_bytes()
        sources[path.name] = digest(raw)
        for row in json.loads(raw):
            folder = folders[row["folder_id"]]
            detail = row["content_details"]
            require("json" in folder["name"].lower(), "wrong folder")
            require(row["status"] == "COMPLETED", "incomplete provider resource")
            require(detail["content_type"] == "application/json", "content type")
            require(row["size"] == detail["size"], "provider size disagreement")
            rows.append({"group": folders[folder["parent_id"]]["name"],
                         "folder_id": folder["id"], "id": row["id"],
                         "filename": row["filename"], "bytes": detail["size"],
                         "sha256": detail["sha256_hash"], "url": detail["download_url"]})
    rows.sort(key=lambda r: (r["group"], r["id"]))
    plan = {"schema": "rgfca-reference-schema-plan-v1",
            "dataset_version": "10.17632/3pw57gdcj2.1",
            "source_receipt_sha256": sources,
            "selection": "smallest_provider_file_id_per_group",
            "public_annotation_files": rows,
            "selected_ids": [min(r["id"] for r in rows if r["group"] == g) for g in GROUPS]}
    validate_plan(plan)
    return plan


def summarize_document(raw, row):
    require(len(raw) == row["bytes"], "download size mismatch")
    require(digest(raw) == row["sha256"], "download SHA-256 mismatch")
    obj = json.loads(raw)
    require(isinstance(obj, dict), "annotation root is not an object")
    shapes = obj.get("shapes")
    require(isinstance(shapes, list), "shapes is not a list")
    labels, kinds, extra_fields, points = Counter(), Counter(), set(), []
    nonempty_groups = 0
    for shape in shapes:
        require(isinstance(shape, dict), "shape is not an object")
        label, kind = shape.get("label"), shape.get("shape_type")
        require(isinstance(label, str) and 0 < len(label) <= 256, "label format")
        require(isinstance(kind, str) and 0 < len(kind) <= 64, "shape type format")
        require(isinstance(shape.get("points"), list), "points is not a list")
        labels[label] += 1
        kinds[kind] += 1
        points.append(len(shape["points"]))
        nonempty_groups += shape.get("group_id") is not None
        extra_fields.update(set(shape) - {"label", "shape_type", "points", "group_id", "flags"})
    embedded = obj.get("imageData")
    image_path = obj.get("imagePath")
    require(image_path is None or isinstance(image_path, str), "imagePath type")
    require(embedded is None or isinstance(embedded, str), "imageData type")
    dims = [obj.get("imageHeight"), obj.get("imageWidth")]
    require(all(type(v) is int and v > 0 for v in dims), "image dimensions")
    # No imageData decoding, point coordinates, flag values or path resolution.
    return {"top_level_keys": sorted(obj), "shape_count": len(shapes),
            "labels": dict(sorted(labels.items())), "shape_types": dict(sorted(kinds.items())),
            "point_counts": points, "nonempty_group_id_count": nonempty_groups,
            "extra_shape_keys": sorted(extra_fields), "image_height_width": dims,
            "embedded_image_string_characters": len(embedded) if embedded else 0,
            "image_path_sha256": digest(image_path.encode()) if image_path else None,
            "image_path_suffix": Path(image_path).suffix.lower() if image_path else None}


def fetch_document(row):
    # Only the six fixed provider URLs are requested; no linked image requests.
    req = Request(row["url"], headers={"User-Agent": "FCP-reference-schema-audit/1"})
    with urlopen(req, timeout=45) as response:
        require(response.url.startswith("https://"), "non-HTTPS redirect")
        return response.read(row["bytes"] + 1)


def run_audit(plan, fetch=fetch_document):
    chosen = validate_plan(plan)
    terminal = []
    for row in chosen:
        result = {"group": row["group"], "id": row["id"],
                  "expected_sha256": row["sha256"], "expected_bytes": row["bytes"]}
        try:
            raw = fetch(row)
            result.update({"observed_sha256": digest(raw), "observed_bytes": len(raw)})
            result["summary"] = summarize_document(raw, row)
            result["status"] = "schema_inspected"
        except Exception as exc:
            # Retain all six terminal entries; no retry, replacement or subset promotion.
            result.update(status="inspection_failed", error_type=type(exc).__name__)
        terminal.append(result)
    return {"schema": "rgfca-reference-schema-audit-v1", "terminal": terminal,
            "all_six_schema_inspected": all(r["status"] == "schema_inspected" for r in terminal),
            "image_pixels_decoded": False, "coordinates_joined": False,
            "model_run": False, "measurement_validation_passed": False,
            "ecological_inference_performed": False,
            "scope": "six_preselected_schema_examples_not_annotation_accuracy_or_complete_payload_census"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    build = sub.add_parser("plan")
    build.add_argument("--receipts", type=Path, required=True)
    build.add_argument("--output", type=Path, required=True)
    inspect = sub.add_parser("inspect")
    inspect.add_argument("--plan", type=Path, required=True)
    inspect.add_argument("--expected-plan-sha256", required=True)
    inspect.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(not args.output.exists(), "output already exists; no reacquisition or overwrite")
    if args.command == "plan":
        result = build_plan(args.receipts)
    else:
        raw = args.plan.read_bytes()
        require(digest(raw) == args.expected_plan_sha256, "plan byte hash mismatch")
        def paced_fetch(row):
            time.sleep(1)
            return fetch_document(row)
        result = run_audit(json.loads(raw), paced_fetch)
        result["plan_sha256"] = digest(raw)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    # Exclusive creation prevents silently replacing a previous plan or outcome.
    with args.output.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(result, handle, indent=2, ensure_ascii=True)
        handle.write("\n")
    if args.command == "inspect" and not result["all_six_schema_inspected"]:
        raise SystemExit("One or more reference-schema inspections failed; see full terminal ledger.")


if __name__ == "__main__":
    main()

"""Read-only ZIP/COCO intake. No pixel decoding, rasterization or model calls."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path, PurePosixPath
import re
import stat
import zipfile


EXPECTED_SHA = "7c9d213514e319b11cc15623747277f88eb37c453722f665b8570b40b9f0dd0d"


def digest(data):
    return hashlib.sha256(data).hexdigest()


def inspect_archive(path: Path):
    before = path.stat()
    archive_digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            archive_digest.update(chunk)
    archive_sha = archive_digest.hexdigest()
    if archive_sha != EXPECTED_SHA:
        raise ValueError("Not the received, pinned Monarda export")
    checks = Counter()
    inventory, members, splits = [], [], {}
    hashes, original_names, paths_used = defaultdict(list), defaultdict(list), set()
    with zipfile.ZipFile(path) as archive:
        entries = archive.infolist()
        names = [e.filename for e in entries]
        if len(names) != len(set(names)) or len(names) != len({n.casefold() for n in names}):
            raise ValueError("Duplicate or case-colliding archive member")
        if sum(e.file_size for e in entries) > 500_000_000:
            raise ValueError("Archive exceeds bounded intake size")
        for entry in entries:
            member = PurePosixPath(entry.filename)
            # ZipInfo normalizes backslashes on Windows and truncates NULs.
            # Reject transformed source names rather than trusting that repair.
            if entry.orig_filename != entry.filename or member.is_absolute() or ".." in member.parts or "\\" in entry.filename or ":" in entry.filename:
                raise ValueError("Unsafe archive path")
            if stat.S_ISLNK(entry.external_attr >> 16) or entry.flag_bits & 1:
                raise ValueError("Linked or encrypted member")
            if entry.file_size > 20_000_000:
                raise ValueError("Member exceeds bounded intake size")
            # Reading every complete member verifies its ZIP CRC; JPEG bytes
            # are hashed, not decoded. No member is extracted or executed.
            raw = archive.read(entry)
            members.append({"path": entry.filename, "bytes": len(raw), "sha256": digest(raw)})
        member_index = {m["path"]: m for m in members}
        readmes = {n: archive.read(n).decode("utf-8-sig") for n in ("README.dataset.txt", "README.roboflow.txt")}
        for split, expected in (("train", 77), ("valid", 22), ("test", 11)):
            coco = json.loads(archive.read(f"{split}/_annotations.coco.json"))
            image_ids = [r["id"] for r in coco["images"]]
            annotation_ids = [r["id"] for r in coco["annotations"]]
            category_ids = [r["id"] for r in coco["categories"]]
            if len(image_ids) != len(set(image_ids)) or len(annotation_ids) != len(set(annotation_ids)):
                raise ValueError("Duplicate COCO key within split")
            if len(category_ids) != len(set(category_ids)):
                raise ValueError("Duplicate category key")
            if len(image_ids) != expected:
                raise ValueError("Split count differs from the pinned v1 page")
            images = {r["id"]: r for r in coco["images"]}
            for row in images.values():
                if any(type(row[k]) is not int or row[k] <= 0 for k in ("width", "height")):
                    raise ValueError("Invalid declared image dimensions")
            per_image, used_categories, segmentation_types = Counter(), Counter(), Counter()
            split_checks = Counter()
            for annotation in coco["annotations"]:
                if annotation["image_id"] not in images or annotation["category_id"] not in category_ids:
                    raise ValueError("Orphan annotation foreign key")
                image = images[annotation["image_id"]]
                per_image[annotation["image_id"]] += 1
                used_categories[str(annotation["category_id"])] += 1
                segmentation = annotation["segmentation"]
                segmentation_types[type(segmentation).__name__] += 1
                if isinstance(segmentation, list):
                    if not segmentation:
                        split_checks["empty_polygon_segmentations"] += 1
                    split_checks["polygon_components"] += len(segmentation)
                    for polygon in segmentation:
                        valid = isinstance(polygon, list) and len(polygon) >= 6 and len(polygon) % 2 == 0
                        valid = valid and all(isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v) for v in polygon)
                        if not valid:
                            split_checks["invalid_polygon_arrays"] += 1
                        elif any(x < 0 or x > image["width"] or y < 0 or y > image["height"] for x, y in zip(polygon[::2], polygon[1::2])):
                            split_checks["polygons_outside_declared_dimensions"] += 1
                else:
                    split_checks["non_polygon_segmentations_not_assessed"] += 1
                box = annotation["bbox"]
                if len(box) != 4 or any(not isinstance(v, (int, float)) or isinstance(v, bool) or not math.isfinite(v) for v in box):
                    split_checks["invalid_bbox"] += 1
                else:
                    x, y, width, height = box
                    if x < 0 or y < 0 or width < 0 or height < 0 or x + width > image["width"] or y + height > image["height"]:
                        split_checks["bbox_outside_declared_dimensions"] += 1
            for image in coco["images"]:
                filename = image["file_name"]
                if PurePosixPath(filename).name != filename:
                    raise ValueError("COCO filename is not a basename")
                member_path = f"{split}/{filename}"
                if member_path not in member_index or member_path in paths_used:
                    raise ValueError("Missing or multiply referenced image member")
                paths_used.add(member_path)
                original = image.get("extra", {}).get("name")
                hashes[member_index[member_path]["sha256"]].append(member_path)
                original_names[original].append(member_path)
                row = {"split": split, "coco_image_id": image["id"], "path": member_path,
                       "original_name": original, "width_declared": image["width"], "height_declared": image["height"],
                       "bytes": member_index[member_path]["bytes"], "sha256": member_index[member_path]["sha256"],
                       "annotation_count": per_image[image["id"]], "license_id": image["license"],
                       "export_date_captured": image.get("date_captured"),
                       "date_captured_equals_export_date_created": image.get("date_captured") == coco["info"].get("date_created")}
                inventory.append(row)
                if per_image[image["id"]] == 0:
                    split_checks["images_without_annotations"] += 1
                if not original or re.fullmatch(r"[0-9]+\.jpg", original) is None:
                    split_checks["unmapped_original_name_format"] += 1
                if row["date_captured_equals_export_date_created"]:
                    split_checks["date_captured_equals_export_date_created"] += 1
            checks.update(split_checks)
            splits[split] = {"images": len(images), "annotations": len(coco["annotations"]),
                             "categories": coco["categories"], "used_category_ids": dict(used_categories),
                             "segmentation_types": dict(segmentation_types), "checks": dict(split_checks),
                             "info": coco["info"], "licenses": coco["licenses"]}
        archive_images = {n for n in names if PurePosixPath(n).suffix.lower() in {".jpg", ".jpeg", ".png"}}
        if archive_images != paths_used:
            raise ValueError("Image members and COCO image census differ")
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError("Input changed during audit")
    return {
        "schema_version": "fcp-monarda-archive-intake-v1", "archive_sha256": archive_sha,
        "archive_bytes": before.st_size, "members": sorted(members, key=lambda r: r["path"]),
        "archive_crc_all_members_verified": True, "readmes_verbatim_untrusted_data": readmes,
        "splits": splits, "total_images": len(inventory),
        "total_annotations": sum(r["annotations"] for r in splits.values()), "checks": dict(checks),
        "duplicate_image_byte_groups": [v for v in hashes.values() if len(v) > 1],
        "duplicate_original_name_groups": [v for v in original_names.values() if len(v) > 1],
        "image_inventory": sorted(inventory, key=lambda r: r["path"]),
        "read_depth": {"zip_bytes_and_coco_json_read": True, "pixel_decode": False, "rasterized_masks": False,
                       "model_execution": False, "colour_measurement": False, "geographic_join": False},
        "limits": ["COCO dimensions are declarations, not JPEG header or decoded-pixel validation.",
                   "Exact-byte and original-name uniqueness are not perceptual, event or observer independence.",
                   "Export dates are not verified observation or phenology dates.",
                   "COCO image IDs are split-local; use split plus id.",
                   "Numeric original names require author-source mapping before use as original photo identifiers.",
                   "Generic flowers annotations are not verified focal-taxon petal truth.",
                   "No external model-training, FCP frame or JRC overlap audit has yet been performed.",
                   "Provider CC BY 4.0 declarations do not verify every original photo's rights."],
        "decision": "complete_archive_and_coco_pair_census_received_not_admitted_for_model_evaluation"
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    result = inspect_archive(args.archive)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({k: result[k] for k in ("archive_sha256", "archive_bytes", "total_images", "total_annotations", "checks", "duplicate_image_byte_groups", "duplicate_original_name_groups", "decision")}, indent=2))

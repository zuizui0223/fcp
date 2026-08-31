#!/usr/bin/env python3
"""Hash-freeze JRC and Oxford ROI sources without decoding image pixels."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path
import sys
from typing import Any, Mapping, Sequence

from scipy.io import loadmat


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.segformer_roi import validate_roi_v3_contract


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_sha256(path: Path) -> str:
    payload = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(payload).hexdigest()


def write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError("cannot write an empty source inventory")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def _verified_json(path: Path, expected_sha256: str) -> dict[str, Any]:
    if sha256(path) != expected_sha256:
        raise RuntimeError(f"source hash mismatch: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def freeze_jrc(
    source_root: Path, contract: Mapping[str, Any]
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    gate = contract["jrc_field_gate"]
    rows: list[dict[str, Any]] = []
    all_hashes: dict[str, str] = {}
    for split, expected_images, expected_boxes in (
        ("train", int(gate["train_images"]), int(gate["train_boxes"])),
        ("test", int(gate["test_images"]), int(gate["test_boxes"])),
    ):
        annotation_path = source_root / "annotations" / f"{split}.json"
        data = _verified_json(annotation_path, gate[f"{split}_annotation_sha256"])
        if data.get("categories") != [{"id": 1, "name": "flower", "supercategory": ""}]:
            raise RuntimeError("JRC flower category identity changed")
        if len(data.get("images", ())) != expected_images or len(
            data.get("annotations", ())
        ) != expected_boxes:
            raise RuntimeError("JRC official split denominator changed")
        annotations = Counter(int(row["image_id"]) for row in data["annotations"])
        for image in sorted(data["images"], key=lambda row: int(row["id"])):
            file_name = str(image["file_name"])
            if Path(file_name).name != file_name:
                raise RuntimeError("JRC annotation contains an unsafe image path")
            image_path = source_root / "images" / split / file_name
            if not image_path.is_file():
                raise RuntimeError(f"missing JRC source image: {split}/{file_name}")
            digest = sha256(image_path)
            if digest in all_hashes:
                raise RuntimeError(
                    f"JRC exact image duplicate: {split}/{file_name} and {all_hashes[digest]}"
                )
            all_hashes[digest] = f"{split}/{file_name}"
            rows.append(
                {
                    "split": split,
                    "image_id": int(image["id"]),
                    "file_name": file_name,
                    "width": int(image["width"]),
                    "height": int(image["height"]),
                    "annotation_boxes": annotations[int(image["id"])],
                    "image_bytes": image_path.stat().st_size,
                    "image_sha256": digest,
                    "image_pixels_decoded": False,
                    "estimator_prediction_run": False,
                }
            )
    for name, expected in (
        ("readme.txt", gate["readme_sha256"]),
        ("copyright.txt", gate["copyright_sha256"]),
    ):
        if sha256(source_root / name) != expected:
            raise RuntimeError(f"JRC {name} hash mismatch")
    return rows, {
        "images": len(rows),
        "boxes": sum(int(row["annotation_boxes"]) for row in rows),
        "split_images": dict(Counter(str(row["split"]) for row in rows)),
        "image_bytes": sum(int(row["image_bytes"]) for row in rows),
        "exact_internal_duplicates": 0,
    }


def freeze_oxford102(
    source_root: Path,
    oxford17_images_dir: Path,
    contract: Mapping[str, Any],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    proxy = contract["oxford102_proxy"]
    for name, expected in (
        ("102flowers.tgz", proxy["images_sha256"]),
        ("102segmentations.tgz", proxy["segmentations_sha256"]),
        ("imagelabels.mat", proxy["labels_sha256"]),
        ("setid.mat", proxy["setid_sha256"]),
    ):
        if sha256(source_root / name) != expected:
            raise RuntimeError(f"Oxford-102 source hash mismatch: {name}")
    labels = loadmat(source_root / "imagelabels.mat")["labels"][0]
    split_data = loadmat(source_root / "setid.mat")
    split_by_id = {
        int(image_id): split
        for split in ("trnid", "valid", "tstid")
        for image_id in split_data[split][0]
    }
    if len(labels) != 8189 or len(split_by_id) != 8189:
        raise RuntimeError("Oxford-102 labels or split denominator changed")
    old_hashes = {
        sha256(path): path.stem for path in sorted(oxford17_images_dir.glob("image_*.jpg"))
    }
    if len(old_hashes) != 1360:
        raise RuntimeError("complete Oxford-17 archive is required for deduplication")

    raw_rows: list[dict[str, Any]] = []
    groups: dict[str, list[int]] = defaultdict(list)
    for image_id in range(1, 8190):
        image_path = source_root / "jpg" / f"image_{image_id:05d}.jpg"
        mask_path = source_root / "segmim" / f"segmim_{image_id:05d}.jpg"
        if not image_path.is_file() or not mask_path.is_file():
            raise RuntimeError(f"Oxford-102 image-mask pair missing: {image_id}")
        image_hash = sha256(image_path)
        groups[image_hash].append(image_id)
        raw_rows.append(
            {
                "image_id": image_id,
                "class_id": int(labels[image_id - 1]),
                "official_split": split_by_id[image_id],
                "image_bytes": image_path.stat().st_size,
                "image_sha256": image_hash,
                "proxy_mask_bytes": mask_path.stat().st_size,
                "proxy_mask_sha256": sha256(mask_path),
            }
        )
    internal_duplicate_ids = {
        image_id for ids in groups.values() if len(ids) > 1 for image_id in ids
    }
    old_overlap_ids = {
        int(row["image_id"]) for row in raw_rows if row["image_sha256"] in old_hashes
    }
    if sorted(old_overlap_ids) != proxy["drop_exact_duplicates_against_full_oxford17_ids"]:
        raise RuntimeError("Oxford-17/102 exact-overlap identity changed")
    excluded = internal_duplicate_ids | old_overlap_ids
    salt = str(proxy["locked_selection_salt"])
    for row in raw_rows:
        payload = (
            f"{salt}\x1f{row['class_id']}\x1f{row['image_sha256']}\x1f{row['image_id']}"
        ).encode("utf-8")
        row["selection_sha256"] = hashlib.sha256(payload).hexdigest()
    eligible_test_by_class: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for row in raw_rows:
        if int(row["image_id"]) not in excluded and row["official_split"] == "tstid":
            eligible_test_by_class[int(row["class_id"])].append(row)
    locked: set[int] = set()
    for class_id in range(1, 103):
        candidates = sorted(
            eligible_test_by_class[class_id],
            key=lambda row: (str(row["selection_sha256"]), int(row["image_id"])),
        )
        if len(candidates) < 20:
            raise RuntimeError(f"Oxford-102 class {class_id} lacks 20 unique test images")
        locked.update(int(row["image_id"]) for row in candidates[:20])
    if len(locked) != int(proxy["locked_images"]):
        raise RuntimeError("Oxford-102 locked proxy denominator changed")

    output: list[dict[str, Any]] = []
    for row in raw_rows:
        image_id = int(row["image_id"])
        if image_id in excluded:
            role = "excluded_exact_duplicate"
        elif row["official_split"] in {"trnid", "valid"}:
            role = "development"
        elif image_id in locked:
            role = "locked_proxy"
        else:
            role = "quarantine"
        output.append(
            {
                **row,
                "full_oxford17_exact_overlap": row["image_sha256"] in old_hashes,
                "full_oxford17_overlap_image_id": old_hashes.get(row["image_sha256"], ""),
                "internal_exact_duplicate_group_size": len(groups[row["image_sha256"]]),
                "role": role,
                "image_pixels_decoded_for_estimator": False,
                "estimator_prediction_run": False,
            }
        )
    return output, {
        "images": len(output),
        "classes": len({int(row["class_id"]) for row in output}),
        "roles": dict(Counter(str(row["role"]) for row in output)),
        "full_oxford17_exact_overlap_ids": sorted(old_overlap_ids),
        "internal_exact_duplicate_ids": sorted(internal_duplicate_ids),
        "unique_image_hashes": len(groups),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--contract",
        type=Path,
        default=Path("docs/supporting/jbi_atlas_roi_estimator_contract_v3.json"),
    )
    parser.add_argument("--jrc-root", type=Path, required=True)
    parser.add_argument("--oxford102-root", type=Path, required=True)
    parser.add_argument("--oxford17-images-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    validate_roi_v3_contract(contract)
    jrc_rows, jrc_summary = freeze_jrc(args.jrc_root, contract)
    oxford_rows, oxford_summary = freeze_oxford102(
        args.oxford102_root, args.oxford17_images_dir, contract
    )
    jrc_hashes = {str(row["image_sha256"]) for row in jrc_rows}
    oxford_hashes = {str(row["image_sha256"]) for row in oxford_rows}
    if jrc_hashes & oxford_hashes:
        raise RuntimeError("JRC contains an exact Oxford-102 image duplicate")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    jrc_path = args.output_dir / "jrc_flower_detection_source_inventory_v1.csv"
    oxford_path = args.output_dir / "oxford102_roi_proxy_source_inventory_v1.csv"
    write_csv(jrc_path, jrc_rows)
    write_csv(oxford_path, oxford_rows)
    manifest = {
        "protocol": contract["protocol"],
        "status": "pass_roi_v3_source_inventory_freeze",
        "contract_path": args.contract.as_posix(),
        "contract_sha256_lf_canonical_v1": canonical_sha256(args.contract),
        "scaleout_candidate_pixels_opened": False,
        "estimator_predictions_run": False,
        "locked_images_decoded": False,
        "jrc": jrc_summary,
        "oxford102": oxford_summary,
        "cross_source_exact_duplicates": 0,
        "files": {
            jrc_path.name: sha256(jrc_path),
            oxford_path.name: sha256(oxford_path),
        },
    }
    manifest_path = args.output_dir / "roi_v3_source_inventory_manifest.json"
    write_json(manifest_path, manifest)
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

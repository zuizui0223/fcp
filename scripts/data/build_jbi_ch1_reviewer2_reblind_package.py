#!/usr/bin/env python3
"""Build a globally re-shuffled reviewer-2 ROI/condition package for all 480 calibration rows.

The review sheets and blank review CSV expose only a fresh reviewer-2 identifier and the
Florence-selected crop. Species, original blind ID, reviewer-1 decisions, colour scores,
geography, observer and date are withheld. A separate mapping file is written for later
post-review reconciliation but is never uploaded in the reviewer-facing artifact.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import time
from urllib.request import Request, urlopen

from PIL import Image, ImageDraw, ImageFont, ImageOps

PROTOCOL = "jbi-ch1-reviewer2-reblind-package-v1"
SALT = "jbi-ch1-reviewer2-reblind-v1-20260830"
USER_AGENT = "zuizui0223-fcp-jbi-ch1-r2-reblind/1.0 (research reproducibility)"
ALLOWED_ROI = ["usable", "rescue_segment", "invalid", "ambiguous"]
ALLOWED_CONDITION = ["fresh", "senescent", "damaged", "mixed_or_ambiguous", "not_evaluable"]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_features(path: Path) -> list[dict]:
    rows = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(rows) != 480:
        raise RuntimeError(f"expected 480 feature rows, found {len(rows)}")
    if len({str(row["blind_id"]) for row in rows}) != 480 or len({str(row["photo_id"]) for row in rows}) != 480:
        raise RuntimeError("feature IDs are not unique")
    if any(row.get("calibration_only") is not True for row in rows):
        raise RuntimeError("non-calibration feature row present")
    if any(row.get("evaluation_row") is not False or row.get("final_label") is not False for row in rows):
        raise RuntimeError("evaluation/final-label firewall violation")
    return rows


def make_r2_id(blind_id: str) -> str:
    return hashlib.sha256(f"{SALT}|{blind_id}".encode()).hexdigest()[:16]


def shuffle_key(blind_id: str) -> str:
    return hashlib.sha256(f"{SALT}|order|{blind_id}".encode()).hexdigest()


def build_assignments(rows: list[dict]) -> list[dict]:
    assignments = []
    for row in rows:
        blind_id = str(row["blind_id"])
        assignments.append({
            "r2_id": make_r2_id(blind_id),
            "blind_id": blind_id,
            "photo_id": str(row["photo_id"]),
            "species": str(row["species"]),
            "downloaded_from": str(row.get("downloaded_from", "")),
            "selected_bbox": list(row.get("selected_bbox", [])),
            "order_key": shuffle_key(blind_id),
        })
    assignments.sort(key=lambda row: row["order_key"])
    for index, row in enumerate(assignments, start=1):
        row["review_order"] = index
    if len(assignments) != 480:
        raise RuntimeError("assignment count changed")
    if len({row["r2_id"] for row in assignments}) != 480:
        raise RuntimeError("r2_id collision")
    return assignments


def download_image(url: str) -> Image.Image:
    if not url:
        raise RuntimeError("missing frozen calibration image URL")
    errors = []
    for attempt in range(1, 4):
        try:
            req = Request(url, headers={"User-Agent": USER_AGENT, "Accept": "image/*"})
            with urlopen(req, timeout=60) as response:
                payload = response.read()
            if len(payload) < 1024:
                raise RuntimeError(f"response too small: {len(payload)}")
            with Image.open(io.BytesIO(payload)) as image:
                return ImageOps.exif_transpose(image).convert("RGB").copy()
        except Exception as exc:
            errors.append(f"attempt={attempt}: {type(exc).__name__}: {exc}")
            time.sleep(0.5 * attempt)
    raise RuntimeError("; ".join(errors[-3:]))


def crop_box_with_context(image: Image.Image, box: list[float], pad_fraction: float = 0.12) -> Image.Image:
    if len(box) != 4:
        raise RuntimeError("selected_bbox is missing")
    width, height = image.size
    x0, y0, x1, y1 = map(float, box)
    x0, x1 = sorted((max(0.0, min(width, x0)), max(0.0, min(width, x1))))
    y0, y1 = sorted((max(0.0, min(height, y0)), max(0.0, min(height, y1))))
    if x1 <= x0 or y1 <= y0:
        raise RuntimeError("invalid selected_bbox")
    pad = max(x1 - x0, y1 - y0) * pad_fraction
    left = max(0, int(x0 - pad))
    top = max(0, int(y0 - pad))
    right = min(width, int(x1 + pad + 0.999))
    bottom = min(height, int(y1 + pad + 0.999))
    return image.crop((left, top, right, bottom)).convert("RGB")


def tile(crop: Image.Image, r2_id: str, review_order: int) -> Image.Image:
    tile_w, tile_h, header = 260, 220, 34
    out = Image.new("RGB", (tile_w, tile_h + header), "white")
    draw = ImageDraw.Draw(out)
    font = ImageFont.load_default()
    draw.text((6, 5), f"{review_order:03d} | {r2_id}", fill="black", font=font)
    fit = ImageOps.contain(crop, (tile_w, tile_h), Image.Resampling.LANCZOS)
    x = (tile_w - fit.width) // 2
    y = header + (tile_h - fit.height) // 2
    out.paste(fit, (x, y))
    return out


def make_pages(tiles: list[Image.Image], output_dir: Path) -> list[str]:
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    per_page = 20
    cols = 4
    rows_per_page = 5
    cell_w, cell_h = 260, 254
    for start in range(0, len(tiles), per_page):
        subset = tiles[start:start + per_page]
        page = Image.new("RGB", (cols * cell_w, rows_per_page * cell_h), "white")
        for offset, image in enumerate(subset):
            x = (offset % cols) * cell_w
            y = (offset // cols) * cell_h
            page.paste(image, (x, y))
        path = output_dir / f"reviewer2_reblind_page_{start // per_page + 1:02d}.jpg"
        page.save(path, quality=92)
        paths.append(str(path))
    if len(paths) != 24:
        raise RuntimeError(f"expected 24 pages, found {len(paths)}")
    return paths


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--features", type=Path, default=Path("data/calibration/jbi_ch1_florence_calibration_features_v1.jsonl"))
    parser.add_argument("--review-csv", type=Path, default=Path("data/calibration/jbi_ch1_reviewer2_reblind_queue_v1.csv"))
    parser.add_argument("--mapping-csv", type=Path, default=Path("data/calibration/jbi_ch1_reviewer2_reblind_mapping_v1.csv"))
    parser.add_argument("--sheets-dir", type=Path, default=Path("artifacts/jbi_ch1_reviewer2_reblind_v1/contact_sheets"))
    parser.add_argument("--manifest", type=Path, default=Path("docs/supporting/jbi_ch1_reviewer2_reblind_manifest_v1.json"))
    args = parser.parse_args()

    features = load_features(args.features)
    assignments = build_assignments(features)
    tiles = []
    failures = []
    for row in assignments:
        try:
            image = download_image(row["downloaded_from"])
            crop = crop_box_with_context(image, row["selected_bbox"])
            tiles.append(tile(crop, row["r2_id"], int(row["review_order"])))
        except Exception as exc:
            failures.append({"r2_id": row["r2_id"], "error": f"{type(exc).__name__}: {exc}"})
            blank = Image.new("RGB", (260, 220), "white")
            ImageDraw.Draw(blank).text((8, 8), "image/crop failure", fill="black")
            tiles.append(tile(blank, row["r2_id"], int(row["review_order"])))
        if int(row["review_order"]) % 40 == 0:
            print(f"built {row['review_order']}/480 reviewer2 crops", flush=True)

    pages = make_pages(tiles, args.sheets_dir)
    review_rows = [
        {
            "review_order": row["review_order"],
            "r2_id": row["r2_id"],
            "target_roi_validity": "",
            "condition_review": "",
            "reviewer2_notes": "",
        }
        for row in assignments
    ]
    mapping_rows = [
        {
            "review_order": row["review_order"],
            "r2_id": row["r2_id"],
            "blind_id": row["blind_id"],
            "photo_id": row["photo_id"],
            "species": row["species"],
            "evaluation_row": False,
            "final_label": False,
        }
        for row in assignments
    ]
    write_csv(args.review_csv, review_rows, ["review_order", "r2_id", "target_roi_validity", "condition_review", "reviewer2_notes"])
    write_csv(args.mapping_csv, mapping_rows, ["review_order", "r2_id", "blind_id", "photo_id", "species", "evaluation_row", "final_label"])

    manifest = {
        "protocol": PROTOCOL,
        "status": "reviewer2_reblind_package_generated_pending_independent_review",
        "calibration_only": True,
        "evaluation_rows_opened": False,
        "final_label": False,
        "n_rows": 480,
        "n_sheets": len(pages),
        "n_failures": len(failures),
        "failures": failures,
        "reblind_salt": SALT,
        "source_feature_sha256": sha256(args.features),
        "review_csv_sha256": sha256(args.review_csv),
        "mapping_csv_sha256": sha256(args.mapping_csv),
        "reviewer_facing_artifact_includes_mapping": False,
        "reviewer_facing_columns": ["review_order", "r2_id", "target_roi_validity", "condition_review", "reviewer2_notes"],
        "reviewer_facing_sheet_exposes_species": False,
        "reviewer_facing_sheet_exposes_original_blind_id": False,
        "reviewer_facing_sheet_exposes_reviewer1_decision": False,
        "reviewer_facing_sheet_exposes_colour_scores": False,
        "geography_exposed": False,
        "observer_exposed": False,
        "date_exposed": False,
        "allowed_target_roi_validity": ALLOWED_ROI,
        "allowed_condition_review": ALLOWED_CONDITION,
        "reviewer2_labels_created_by_script": False,
        "independent_second_review_completed": False,
        "contact_sheets": pages,
        "next_gate": "complete reviewer2 decisions without opening mapping or reviewer1 decisions, then reconcile agreement and disagreements before any colour-rule freeze"
    }
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(manifest, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

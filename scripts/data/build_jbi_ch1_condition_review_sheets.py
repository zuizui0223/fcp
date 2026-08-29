#!/usr/bin/env python3
"""Build blinded flower-condition review sheets from frozen calibration ROI features.

Only 480 calibration IDs are allowed. Images are re-downloaded from their frozen URLs,
cropped to the Florence-selected flower ROI, and shown with blind IDs only. Geography,
observer, date, colour candidate scores and environmental metadata are not included in
review sheets. No condition label is created by this script.
"""
from __future__ import annotations

import argparse
import io
import json
from pathlib import Path
import re
import time
from urllib.request import Request, urlopen

import pandas as pd
from PIL import Image, ImageDraw, ImageFont, ImageOps

PROTOCOL = "jbi-ch1-condition-review-sheets-v1"
USER_AGENT = "zuizui0223-fcp-condition-review/1.0 (research reproducibility)"
TILES_PER_PAGE = 20
TILE_W = 300
TILE_H = 300
HEADER_H = 36


def slug(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_")


def load_features(path: Path) -> list[dict]:
    rows = [json.loads(x) for x in path.read_text(encoding="utf-8").splitlines() if x.strip()]
    if len(rows) != 480:
        raise RuntimeError(f"expected 480 calibration feature rows, found {len(rows)}")
    if len({str(r["photo_id"]) for r in rows}) != 480:
        raise RuntimeError("duplicate photo IDs")
    if any(r.get("evaluation_row") is not False or r.get("final_label") is not False for r in rows):
        raise RuntimeError("feature firewall violation")
    return rows


def download(url: str) -> Image.Image:
    last = None
    for attempt in range(1, 4):
        try:
            req = Request(url, headers={"User-Agent": USER_AGENT, "Accept": "image/*"})
            with urlopen(req, timeout=60) as response:
                payload = response.read()
            if len(payload) < 1024:
                raise RuntimeError(f"response too small: {len(payload)} bytes")
            with Image.open(io.BytesIO(payload)) as image:
                return ImageOps.exif_transpose(image).convert("RGB").copy()
        except Exception as exc:
            last = exc
            time.sleep(0.5 * attempt)
    raise RuntimeError(f"download failed: {last}")


def crop_with_context(image: Image.Image, bbox: list[float], pad_fraction: float = 0.12) -> Image.Image:
    if len(bbox) != 4:
        raise ValueError("bbox must have four coordinates")
    w, h = image.size
    x0, y0, x1, y1 = map(float, bbox)
    bw, bh = x1 - x0, y1 - y0
    if bw <= 0 or bh <= 0:
        raise ValueError("invalid bbox")
    pad = max(bw, bh) * pad_fraction
    x0 = max(0, int(x0 - pad))
    y0 = max(0, int(y0 - pad))
    x1 = min(w, int(x1 + pad + 0.999))
    y1 = min(h, int(y1 + pad + 0.999))
    return image.crop((x0, y0, x1, y1)).convert("RGB")


def tile_for_crop(crop: Image.Image, blind_id: str) -> Image.Image:
    tile = Image.new("RGB", (TILE_W, TILE_H + HEADER_H), "white")
    canvas = ImageOps.contain(crop, (TILE_W, TILE_H), Image.Resampling.LANCZOS)
    x = (TILE_W - canvas.width) // 2
    y = HEADER_H + (TILE_H - canvas.height) // 2
    tile.paste(canvas, (x, y))
    draw = ImageDraw.Draw(tile)
    draw.text((8, 10), blind_id, fill="black", font=ImageFont.load_default())
    return tile


def save_pages(species: str, tiles: list[Image.Image], out_dir: Path) -> list[str]:
    cols, rows_per_page = 5, 4
    page_paths = []
    for start in range(0, len(tiles), TILES_PER_PAGE):
        subset = tiles[start:start + TILES_PER_PAGE]
        page = Image.new("RGB", (cols * TILE_W, rows_per_page * (TILE_H + HEADER_H)), "white")
        for i, tile in enumerate(subset):
            col = i % cols
            row = i // cols
            page.paste(tile, (col * TILE_W, row * (TILE_H + HEADER_H)))
        page_no = start // TILES_PER_PAGE + 1
        path = out_dir / f"{slug(species)}_condition_page_{page_no:02d}.jpg"
        page.save(path, quality=92)
        page_paths.append(str(path))
    return page_paths


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--features",
        type=Path,
        default=Path("data/calibration/jbi_ch1_florence_calibration_features_v1.jsonl"),
    )
    parser.add_argument(
        "--split",
        type=Path,
        default=Path("data/frozen/jbi_ch1_photo_split_v1.csv"),
    )
    parser.add_argument(
        "--sheets-dir",
        type=Path,
        default=Path("artifacts/jbi_ch1_condition_review_v1/contact_sheets"),
    )
    parser.add_argument(
        "--review-csv",
        type=Path,
        default=Path("data/calibration/jbi_ch1_condition_review_queue_v1.csv"),
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("docs/supporting/jbi_ch1_condition_review_manifest_v1.json"),
    )
    args = parser.parse_args()

    features = load_features(args.features)
    split = pd.read_csv(args.split)
    calibration = split.loc[split["split"].astype(str).eq("calibration")].copy()
    evaluation_ids = set(split.loc[split["split"].astype(str).eq("evaluation"), "photo_id"].astype(str))
    if len(calibration) != 480:
        raise RuntimeError("frozen split does not contain 480 calibration rows")
    by_id = {str(r["photo_id"]): r for _, r in calibration.iterrows()}
    if set(by_id) != {str(r["photo_id"]) for r in features}:
        raise RuntimeError("feature IDs differ from frozen calibration IDs")
    if set(by_id) & evaluation_ids:
        raise RuntimeError("evaluation leakage")

    args.sheets_dir.mkdir(parents=True, exist_ok=True)
    review_rows = []
    sheet_paths = []
    failures = []

    species_set = sorted({str(r["species"]) for r in features})
    for species in species_set:
        group = sorted([r for r in features if r["species"] == species], key=lambda r: r["blind_id"])
        if len(group) != 80:
            raise RuntimeError(f"{species}: expected 80 rows")
        tiles = []
        for index, row in enumerate(group, start=1):
            pid = str(row["photo_id"])
            blind_id = str(row["blind_id"])
            split_row = by_id[pid]
            status = row.get("feature_status")
            bbox = row.get("selected_bbox")
            crop_status = "ok"
            if status != "ok" or not isinstance(bbox, list):
                crop_status = "not_evaluable_localization"
                crop = Image.new("RGB", (TILE_W, TILE_H), "white")
                ImageDraw.Draw(crop).text((10, 10), "ROI unavailable", fill="black")
            else:
                try:
                    urls = [
                        str(split_row[c]).strip()
                        for c in ("photo_url", "photo_url_api")
                        if c in split_row.index and not pd.isna(split_row[c]) and str(split_row[c]).strip()
                    ]
                    if not urls:
                        raise RuntimeError("no image URL")
                    image = download(urls[0])
                    crop = crop_with_context(image, bbox)
                except Exception as exc:
                    crop_status = "not_evaluable_download_or_crop"
                    failures.append({"blind_id": blind_id, "species": species, "error": f"{type(exc).__name__}: {exc}"})
                    crop = Image.new("RGB", (TILE_W, TILE_H), "white")
                    ImageDraw.Draw(crop).text((10, 10), "Crop unavailable", fill="black")
            tiles.append(tile_for_crop(crop, blind_id))
            review_rows.append({
                "species": species,
                "blind_id": blind_id,
                "photo_id": pid,
                "review_order_within_species": index,
                "crop_status": crop_status,
                "condition_review": "",
                "condition_failure_reason": "",
                "review_notes": "",
                "final_condition_label": False,
                "evaluation_row": False,
            })
        sheet_paths.extend(save_pages(species, tiles, args.sheets_dir))

    review = pd.DataFrame(review_rows)
    if len(review) != 480 or review["blind_id"].nunique() != 480:
        raise RuntimeError("review queue is not exactly 480 unique calibration rows")
    if review["evaluation_row"].astype(bool).any():
        raise RuntimeError("evaluation leakage in review queue")
    args.review_csv.parent.mkdir(parents=True, exist_ok=True)
    review.to_csv(args.review_csv, index=False, lineterminator="\n")

    manifest = {
        "protocol": PROTOCOL,
        "status": "blinded_condition_review_package_generated",
        "n_rows": 480,
        "species_count": 6,
        "sheets_per_species": 4,
        "n_sheets": len(sheet_paths),
        "contact_sheets": sheet_paths,
        "review_csv": str(args.review_csv),
        "calibration_only": True,
        "evaluation_rows_opened": False,
        "geography_exposed_to_reviewer": False,
        "observer_exposed_to_reviewer": False,
        "date_exposed_to_reviewer": False,
        "colour_candidate_scores_exposed_to_reviewer": False,
        "condition_labels_created_by_script": False,
        "allowed_review_values": ["fresh", "senescent", "damaged", "mixed_or_ambiguous", "not_evaluable"],
        "failures": failures,
    }
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(manifest, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

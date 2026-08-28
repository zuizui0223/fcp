#!/usr/bin/env python3
"""Materialize and technically audit only the frozen 480-photo calibration set.

This script never selects or opens evaluation rows. It downloads the already-frozen
calibration photo URLs, records byte hashes and technical image metrics, and builds
blind contact sheets for measurement-rule calibration. Technical metrics are not used
to assign flower colour and do not alter the frozen split.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import io
import json
from pathlib import Path
import re
import time
from typing import Iterable
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont, ImageOps


PROTOCOL = "jbi-ch1-calibration-images-v1"
USER_AGENT = "zuizui0223-fcp-jbi-ch1-calibration/1.0 (research reproducibility)"


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_path(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def blind_id(species: str, photo_id: str) -> str:
    payload = f"jbi-ch1-calibration-v1\x1f{species}\x1f{photo_id}".encode("utf-8")
    return hashlib.sha256(payload).hexdigest()[:16]


def slugify(value: str) -> str:
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    return value.strip("_") or "species"


def calibration_rows(split: pd.DataFrame) -> pd.DataFrame:
    required = {"species", "photo_id", "split"}
    missing = required - set(split.columns)
    if missing:
        raise ValueError(f"split missing required columns: {sorted(missing)}")
    rows = split.loc[split["split"].astype(str) == "calibration"].copy()
    if len(rows) != 480:
        raise ValueError(f"expected exactly 480 calibration rows, found {len(rows)}")
    counts = rows.groupby("species").size()
    if len(counts) != 6 or not (counts == 80).all():
        raise ValueError("expected six species with exactly 80 calibration rows each")
    if rows["photo_id"].astype(str).duplicated().any():
        raise ValueError("calibration rows contain duplicate photo IDs")
    if not any(c in rows.columns for c in ("photo_url", "photo_url_api")):
        raise ValueError("calibration rows contain no photo URL column")
    return rows


def candidate_urls(row: pd.Series) -> list[str]:
    urls: list[str] = []
    for column in ("photo_url", "photo_url_api"):
        if column not in row.index:
            continue
        value = row.get(column)
        if pd.isna(value):
            continue
        text = str(value).strip()
        if text and text not in urls:
            urls.append(text)
    return urls


def download_bytes(
    urls: Iterable[str], *, retries: int = 3, pause_seconds: float = 0.5
) -> tuple[bytes, str]:
    errors: list[str] = []
    for url in urls:
        for attempt in range(1, retries + 1):
            try:
                req = Request(url, headers={"User-Agent": USER_AGENT, "Accept": "image/*"})
                with urlopen(req, timeout=60) as response:
                    payload = response.read()
                if len(payload) < 1024:
                    raise RuntimeError(f"response too small: {len(payload)} bytes")
                return payload, url
            except Exception as exc:
                errors.append(f"{url} attempt={attempt}: {type(exc).__name__}: {exc}")
                if pause_seconds > 0:
                    time.sleep(pause_seconds * attempt)
    raise RuntimeError("; ".join(errors[-8:]))


def decode_rgb(payload: bytes) -> Image.Image:
    with Image.open(io.BytesIO(payload)) as image:
        image = ImageOps.exif_transpose(image)
        image.load()
        return image.convert("RGB")


def technical_metrics(image: Image.Image) -> dict[str, float | int]:
    width, height = image.size
    thumb = image.copy()
    thumb.thumbnail((384, 384), Image.Resampling.LANCZOS)
    arr = np.asarray(thumb, dtype=np.float32) / 255.0
    if arr.ndim != 3 or arr.shape[2] != 3:
        raise ValueError("expected RGB image")

    vmax = arr.max(axis=2)
    vmin = arr.min(axis=2)
    saturation = np.divide(
        vmax - vmin,
        np.maximum(vmax, 1e-6),
        out=np.zeros_like(vmax),
        where=vmax > 1e-6,
    )
    luminance = 0.2126 * arr[:, :, 0] + 0.7152 * arr[:, :, 1] + 0.0722 * arr[:, :, 2]
    dx = np.diff(luminance, axis=1)
    dy = np.diff(luminance, axis=0)
    edge_energy = float((np.mean(dx * dx) + np.mean(dy * dy)) / 2.0)

    h, w = luminance.shape
    y0, y1 = int(h * 0.25), max(int(h * 0.75), int(h * 0.25) + 1)
    x0, x1 = int(w * 0.25), max(int(w * 0.75), int(w * 0.25) + 1)
    central = luminance[y0:y1, x0:x1]

    return {
        "width_px": int(width),
        "height_px": int(height),
        "megapixels": float(width * height / 1_000_000.0),
        "min_dimension_px": int(min(width, height)),
        "mean_luminance": float(luminance.mean()),
        "sd_luminance": float(luminance.std()),
        "central_sd_luminance": float(central.std()) if central.size else 0.0,
        "dark_clip_fraction": float((luminance <= 0.02).mean()),
        "bright_clip_fraction": float((luminance >= 0.98).mean()),
        "mean_saturation": float(saturation.mean()),
        "p90_saturation": float(np.quantile(saturation, 0.90)),
        "edge_energy": edge_energy,
    }


def technical_flags(metrics: dict[str, float | int]) -> list[str]:
    """Conservative technical flags only; these are not biological exclusions."""
    flags: list[str] = []
    if int(metrics["min_dimension_px"]) < 256:
        flags.append("low_resolution")
    if float(metrics["dark_clip_fraction"]) > 0.50:
        flags.append("severe_dark_clipping")
    if float(metrics["bright_clip_fraction"]) > 0.50:
        flags.append("severe_bright_clipping")
    if float(metrics["edge_energy"]) < 0.00005:
        flags.append("very_low_detail")
    return flags


def fit_thumbnail(image: Image.Image, size: tuple[int, int]) -> Image.Image:
    tile = Image.new("RGB", size, "white")
    work = image.copy()
    work.thumbnail(size, Image.Resampling.LANCZOS)
    x = (size[0] - work.width) // 2
    y = (size[1] - work.height) // 2
    tile.paste(work, (x, y))
    return tile


def make_contact_sheet(
    entries: list[tuple[str, Path]],
    output: Path,
    *,
    columns: int = 5,
    rows: int = 4,
    tile_size: tuple[int, int] = (240, 210),
    label_height: int = 28,
) -> None:
    capacity = columns * rows
    if len(entries) > capacity:
        raise ValueError("too many entries for one contact sheet")
    canvas = Image.new(
        "RGB",
        (columns * tile_size[0], rows * (tile_size[1] + label_height)),
        "white",
    )
    draw = ImageDraw.Draw(canvas)
    font = ImageFont.load_default()
    for index, (label, image_path) in enumerate(entries):
        r, c = divmod(index, columns)
        x = c * tile_size[0]
        y = r * (tile_size[1] + label_height)
        with Image.open(image_path) as source:
            source = ImageOps.exif_transpose(source).convert("RGB")
            tile = fit_thumbnail(source, tile_size)
        canvas.paste(tile, (x, y))
        draw.rectangle(
            [x, y + tile_size[1], x + tile_size[0] - 1, y + tile_size[1] + label_height - 1],
            outline="black",
            width=1,
        )
        draw.text((x + 5, y + tile_size[1] + 7), label, fill="black", font=font)
    output.parent.mkdir(parents=True, exist_ok=True)
    canvas.save(output, format="JPEG", quality=88, optimize=True)


def materialize(
    split_csv: Path,
    image_root: Path,
    qc_csv: Path,
    manifest_json: Path,
    contact_root: Path,
    *,
    retries: int = 3,
    pause_seconds: float = 0.5,
) -> dict[str, object]:
    split = pd.read_csv(split_csv)
    rows = calibration_rows(split)
    image_root.mkdir(parents=True, exist_ok=True)
    contact_root.mkdir(parents=True, exist_ok=True)

    qc_rows: list[dict[str, object]] = []
    contact_entries: dict[str, list[tuple[str, Path]]] = defaultdict(list)
    failures: list[dict[str, str]] = []

    for index, (_, row) in enumerate(rows.iterrows(), start=1):
        species = str(row["species"]).strip()
        photo_id = str(row["photo_id"]).strip()
        bid = blind_id(species, photo_id)
        out_path = image_root / slugify(species) / f"{bid}.jpg"
        out_path.parent.mkdir(parents=True, exist_ok=True)
        record: dict[str, object] = {
            "species": species,
            "blind_id": bid,
            "photo_id": photo_id,
            "materialization_status": "failed",
            "technical_flags": "",
        }
        try:
            payload, used_url = download_bytes(
                candidate_urls(row), retries=retries, pause_seconds=pause_seconds
            )
            image = decode_rgb(payload)
            metrics = technical_metrics(image)
            flags = technical_flags(metrics)
            out_path.write_bytes(payload)
            record.update(metrics)
            record.update(
                {
                    "materialization_status": "ok",
                    "downloaded_from": used_url,
                    "raw_image_sha256": sha256_bytes(payload),
                    "raw_image_bytes": len(payload),
                    "technical_flags": ";".join(flags),
                    "artifact_path": str(out_path),
                }
            )
            contact_entries[species].append((bid, out_path))
            image.close()
        except Exception as exc:
            record["error"] = f"{type(exc).__name__}: {exc}"
            failures.append({"species": species, "blind_id": bid, "error": record["error"]})
        qc_rows.append(record)
        if index % 40 == 0:
            print(f"materialized {index}/480", flush=True)

    qc = pd.DataFrame(qc_rows).sort_values(["species", "blind_id"], kind="mergesort")
    qc_csv.parent.mkdir(parents=True, exist_ok=True)
    qc.to_csv(qc_csv, index=False, lineterminator="\n")

    sheet_paths: list[str] = []
    for species, entries in sorted(contact_entries.items()):
        entries = sorted(entries, key=lambda item: item[0])
        for page_index, start in enumerate(range(0, len(entries), 20), start=1):
            page = entries[start : start + 20]
            path = contact_root / f"{slugify(species)}_page_{page_index:02d}.jpg"
            make_contact_sheet(page, path)
            sheet_paths.append(str(path))

    per_species = []
    for species, group in qc.groupby("species", sort=True):
        ok = group["materialization_status"].astype(str).eq("ok")
        flagged = group["technical_flags"].fillna("").astype(str).ne("")
        per_species.append(
            {
                "species": species,
                "n_expected": int(len(group)),
                "n_materialized": int(ok.sum()),
                "n_failed": int((~ok).sum()),
                "n_technical_flagged": int((ok & flagged).sum()),
            }
        )

    manifest: dict[str, object] = {
        "protocol": PROTOCOL,
        "status": "pass" if not failures else "failed",
        "source_split": str(split_csv),
        "source_split_sha256": sha256_path(split_csv),
        "qc_csv": str(qc_csv),
        "qc_csv_sha256": sha256_path(qc_csv),
        "image_root": str(image_root),
        "contact_sheet_root": str(contact_root),
        "contact_sheets": sheet_paths,
        "n_expected": 480,
        "n_materialized": int((qc["materialization_status"] == "ok").sum()),
        "n_failed": len(failures),
        "per_species": per_species,
        "evaluation_rows_opened": False,
        "technical_metrics_used_for_colour_assignment": False,
        "technical_flags_are_exclusion_rules": False,
        "failures": failures,
    }
    manifest_json.parent.mkdir(parents=True, exist_ok=True)
    manifest_json.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--split", type=Path, default=Path("data/frozen/jbi_ch1_photo_split_v1.csv"))
    parser.add_argument(
        "--image-root",
        type=Path,
        default=Path("artifacts/jbi_ch1_calibration_images_v1/images"),
    )
    parser.add_argument(
        "--qc-csv",
        type=Path,
        default=Path("data/calibration/jbi_ch1_calibration_image_qc_v1.csv"),
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("docs/supporting/jbi_ch1_calibration_image_manifest_v1.json"),
    )
    parser.add_argument(
        "--contact-root",
        type=Path,
        default=Path("artifacts/jbi_ch1_calibration_images_v1/contact_sheets"),
    )
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument("--pause-seconds", type=float, default=0.5)
    args = parser.parse_args()

    manifest = materialize(
        args.split,
        args.image_root,
        args.qc_csv,
        args.manifest,
        args.contact_root,
        retries=args.retries,
        pause_seconds=args.pause_seconds,
    )
    print(json.dumps(manifest, indent=2, ensure_ascii=False))
    if manifest["n_failed"]:
        raise SystemExit(f"{manifest['n_failed']} frozen calibration images failed to materialize")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

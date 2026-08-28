#!/usr/bin/env python3
"""Materialize one deterministic calibration image per species for a vision pilot.

The pilot is strictly calibration-only and exists to validate semantic image-review
automation before any 480-image run is considered.
"""

from __future__ import annotations

import argparse
import hashlib
import io
from pathlib import Path
import re
import time
from urllib.request import Request, urlopen

import pandas as pd
from PIL import Image, ImageOps


USER_AGENT = "zuizui0223-fcp-jbi-ch1-vision-pilot/1.0"


def blind_id(species: str, photo_id: str) -> str:
    payload = f"jbi-ch1-calibration-v1\x1f{species}\x1f{photo_id}".encode("utf-8")
    return hashlib.sha256(payload).hexdigest()[:16]


def slugify(value: str) -> str:
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    return value.strip("_") or "species"


def pick_pilot(split: pd.DataFrame) -> pd.DataFrame:
    required = {"species", "photo_id", "split"}
    missing = required - set(split.columns)
    if missing:
        raise ValueError(f"split missing columns: {sorted(missing)}")
    calibration = split.loc[split["split"].astype(str).eq("calibration")].copy()
    evaluation = split.loc[split["split"].astype(str).eq("evaluation")].copy()
    if len(calibration) != 480 or len(evaluation) != 720:
        raise ValueError("expected frozen 480/720 split")
    if calibration.groupby("species").size().ne(80).any() or calibration["species"].nunique() != 6:
        raise ValueError("expected 80 calibration rows for each of six species")
    calibration["blind_id"] = [
        blind_id(str(s).strip(), str(p).strip())
        for s, p in calibration[["species", "photo_id"]].itertuples(index=False, name=None)
    ]
    pilot = (
        calibration.sort_values(["species", "blind_id"], kind="mergesort")
        .groupby("species", sort=True, as_index=False)
        .head(1)
        .copy()
    )
    if len(pilot) != 6:
        raise ValueError("pilot must contain exactly one image per species")
    if set(pilot["photo_id"].astype(str)) & set(evaluation["photo_id"].astype(str)):
        raise ValueError("evaluation leakage detected")
    return pilot


def urls_for(row: pd.Series) -> list[str]:
    urls = []
    for column in ("photo_url", "photo_url_api"):
        if column in row.index and not pd.isna(row[column]):
            url = str(row[column]).strip()
            if url and url not in urls:
                urls.append(url)
    if not urls:
        raise ValueError("pilot row has no image URL")
    return urls


def download_decodable(urls: list[str], *, retries: int = 3) -> tuple[bytes, str, tuple[int, int]]:
    errors = []
    for url in urls:
        for attempt in range(1, retries + 1):
            try:
                req = Request(url, headers={"User-Agent": USER_AGENT, "Accept": "image/*"})
                with urlopen(req, timeout=60) as response:
                    payload = response.read()
                if len(payload) < 1024:
                    raise RuntimeError(f"response too small: {len(payload)}")
                with Image.open(io.BytesIO(payload)) as image:
                    image = ImageOps.exif_transpose(image)
                    image.load()
                    size = image.size
                return payload, url, size
            except Exception as exc:
                errors.append(f"{url} attempt={attempt}: {type(exc).__name__}: {exc}")
                time.sleep(0.5 * attempt)
    raise RuntimeError("; ".join(errors[-8:]))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--split", type=Path, default=Path("data/frozen/jbi_ch1_photo_split_v1.csv"))
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("artifacts/jbi_ch1_vision_pilot_v1/images"),
    )
    parser.add_argument(
        "--manifest-csv",
        type=Path,
        default=Path("artifacts/jbi_ch1_vision_pilot_v1/pilot_manifest.csv"),
    )
    args = parser.parse_args()

    split = pd.read_csv(args.split)
    pilot = pick_pilot(split)
    records = []
    for _, row in pilot.iterrows():
        species = str(row["species"]).strip()
        photo_id = str(row["photo_id"]).strip()
        bid = str(row["blind_id"])
        payload, used_url, size = download_decodable(urls_for(row))
        path = args.output_root / slugify(species) / f"{bid}.jpg"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(payload)
        records.append(
            {
                "species": species,
                "blind_id": bid,
                "photo_id": photo_id,
                "image_path": str(path),
                "image_sha256": hashlib.sha256(payload).hexdigest(),
                "width_px": int(size[0]),
                "height_px": int(size[1]),
                "downloaded_from": used_url,
                "pilot_only": True,
                "evaluation_row": False,
            }
        )
        print(f"pilot materialized: {species} {bid}", flush=True)

    out = pd.DataFrame(records).sort_values("species", kind="mergesort")
    args.manifest_csv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.manifest_csv, index=False, lineterminator="\n")
    if len(out) != 6 or out["blind_id"].nunique() != 6:
        raise RuntimeError("vision pilot manifest failed final cardinality check")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

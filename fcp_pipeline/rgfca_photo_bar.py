"""Deterministic, licensed discovery-only photographic illustration.

The estimator is imported by the execution entry point only; selection and
display tests require no model, network, reserve data or source image.
"""
from __future__ import annotations

import hashlib
import io
import json
from pathlib import Path
import subprocess
from urllib.parse import urlsplit

import numpy as np
import pandas as pd
from PIL import Image, ImageOps

ROOT = Path(__file__).resolve().parents[1]
SOURCE_COMMIT = "29584f3ad7ae0cd99a1d8f43459252af38f615da"
SALT = "fcp-rgfca-photo-bar-v1"
SLOTS = 24
CC0_URL = "https://creativecommons.org/publicdomain/zero/1.0/"
METADATA = ["measurement_id", "photo_id", "observation_id", "species", "observer_id",
            "photo_license", "attribution", "longitude", "latitude"]
PALETTE = ("white", "yellow", "orange", "red", "pink", "magenta", "purple", "blue",
           "bronze", "green", "brown", "black")


def sha256(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def canonical(value) -> bytes:
    return (json.dumps(value, ensure_ascii=False, sort_keys=True, indent=2, allow_nan=False) + "\n").encode("utf-8")


def committed(path: str) -> bytes:
    return subprocess.check_output(["git", "show", f"{SOURCE_COMMIT}:{path}"], cwd=ROOT)


def discovery_sources():
    manifest_path = "docs/supporting/global_monte_carlo_measurement_result_v1.json"
    manifest = json.loads(committed(manifest_path))
    paths = {
        "measured": "data/derived/global_monte_carlo_measured_photos_v1.csv",
        "candidate": "data/frozen/global_monte_carlo_candidate_photos_v1.csv",
    }
    raw = {key: committed(path) for key, path in paths.items()}
    expected = {"measured": manifest["lineage"]["measured_table_sha256"],
                "candidate": manifest["lineage"]["candidate_photos_sha256"]}
    if {k: sha256(v) for k, v in raw.items()} != expected:
        raise ValueError("Exact discovery source hash differs")
    frames = {k: pd.read_csv(io.BytesIO(v), dtype=str, keep_default_na=False) for k, v in raw.items()}
    measured, candidate = frames["measured"], frames["candidate"]
    if len(measured) != 50000 or measured.photo_id.nunique() != 50000 or measured.species.nunique() != 500:
        raise ValueError("Discovery measurement census differs")
    if not measured.groupby("species").size().eq(100).all():
        raise ValueError("Discovery species budget differs")
    if len(candidate) != 100000 or candidate.photo_id.nunique() != 100000:
        raise ValueError("Candidate metadata census differs")
    if not measured.global_classifiable.isin(["True", "False"]).all():
        raise ValueError("Unknown classifiable flag")
    classified = measured.loc[measured.global_classifiable.eq("True")]
    counts = classified.groupby("species").size()
    eligible = counts[counts >= 40].index
    pool = classified.loc[classified.species.isin(eligible)].copy()
    if len(classified) != 25377 or len(eligible) != 369 or len(pool) != 21424:
        raise ValueError("Eligible discovery census differs")
    return pool, candidate, {paths[k]: expected[k] for k in paths}


def select_metadata(metadata: pd.DataFrame, slots: int = SLOTS) -> list[dict]:
    if set(metadata.columns) != set(METADATA):
        raise ValueError("Selection must receive metadata only")
    if metadata.photo_id.duplicated().any() or metadata.measurement_id.duplicated().any():
        raise ValueError("Duplicate photo or measurement ID")
    pool = metadata.loc[metadata.photo_license.eq("cc0") & metadata.attribution.str.strip().ne("")].copy()
    if len(pool) < slots or slots <= 0:
        raise ValueError("Insufficient CC0 display pool")
    for field in ("photo_id", "observation_id", "observer_id"):
        if not pool[field].str.fullmatch(r"[1-9][0-9]*").all():
            raise ValueError(f"Invalid {field}")
    pool["longitude"] = pd.to_numeric(pool.longitude, errors="raise")
    pool["latitude"] = pd.to_numeric(pool.latitude, errors="raise")
    if not np.isfinite(pool[["longitude", "latitude"]].to_numpy()).all() or not pool.longitude.between(-180, 180).all() or not pool.latitude.between(-90, 90).all():
        raise ValueError("Invalid display geography")
    pool["numeric_photo_id"] = pool.photo_id.map(int)
    pool = pool.sort_values(["longitude", "numeric_photo_id"]).reset_index(drop=True)
    pool["rank_bin"] = np.arange(len(pool)) * slots // len(pool)
    pool["selection_rank"] = pool.photo_id.map(lambda value: sha256(f"{SALT}|{value}".encode()))
    used_species, used_observers, selected = set(), set(), []
    for index in range(slots):
        group = pool.loc[pool.rank_bin.eq(index)].sort_values(["selection_rank", "numeric_photo_id"])
        admissible = group.loc[~group.species.isin(used_species) & ~group.observer_id.isin(used_observers)]
        if admissible.empty:
            raise ValueError(f"No metadata-admissible photo in rank bin {index + 1}")
        row = admissible.iloc[0]
        used_species.add(row.species)
        used_observers.add(row.observer_id)
        selected.append({**row[METADATA].to_dict(), "slot": index + 1,
                         "selection_rank": row.selection_rank, "rank_bin_rows": len(group)})
    return selected


def validate_photo_url(url: str, photo_id: str):
    parsed = urlsplit(url)
    hosts = {"inaturalist-open-data.s3.amazonaws.com", "static.inaturalist.org"}
    if parsed.scheme != "https" or parsed.hostname not in hosts or parsed.username or parsed.password or parsed.port or parsed.query or parsed.fragment:
        raise ValueError("Unexpected photo source URL")
    if parsed.path not in {f"/photos/{photo_id}/large.jpg", f"/photos/{photo_id}/large.jpeg", f"/photos/{photo_id}/large.png"}:
        raise ValueError("Photo URL does not identify the frozen large image")


def make_plan() -> dict:
    pool, candidate, hashes = discovery_sources()
    selected = select_metadata(pool[METADATA])
    measured_by_id = pool.set_index("photo_id")
    candidates_by_id = candidate.set_index("photo_id")
    for row in selected:
        photo_id = row["photo_id"]
        original, metadata = measured_by_id.loc[photo_id], candidates_by_id.loc[photo_id]
        for field in ("observation_id", "photo_license", "species", "observer_id"):
            if str(original[field]) != str(metadata[field]):
                raise ValueError("Candidate and measurement metadata differ")
        url = str(metadata.photo_url_large)
        validate_photo_url(url, photo_id)
        row.update({"photo_url_large": url, "image_sha256": str(original.image_sha256),
                    "expected_mask_pixels": str(original.mask_pixels),
                    "expected_background_pixels": str(original.background_effective_pixels),
                    "expected_palette": {name: str(original[f"palette_count_{name}"]) for name in PALETTE},
                    "photo_page": f"https://www.inaturalist.org/photos/{photo_id}",
                    "observation_page": f"https://www.inaturalist.org/observations/{row['observation_id']}",
                    "license_url": CC0_URL})
    return {"protocol": "rgfca-discovery-photo-bar-v1", "source_commit": SOURCE_COMMIT,
            "status": "metadata_only_display_plan_before_reacquisition", "slots": SLOTS,
            "source_sha256": hashes, "eligible_discovery_photos": len(pool),
            "eligible_cc0_photos_with_credit": int((pool.photo_license.eq("cc0") & pool.attribution.str.strip().ne("")).sum()),
            "selection_salt": SALT, "selected": selected,
            "new_inference": False, "reserve_outcomes_read": False, "replacement_allowed": False}


def current_photo_rights(observation: dict, selected: dict) -> dict:
    results = observation.get("results", [])
    if len(results) != 1 or str(results[0].get("id")) != selected["observation_id"]:
        raise ValueError("Observation response identity differs")
    matches = [p for p in results[0].get("photos", []) if str(p.get("id")) == selected["photo_id"]]
    if len(matches) != 1 or matches[0].get("license_code") != "cc0" or not str(matches[0].get("attribution") or "").strip():
        raise ValueError("Selected photo-level CC0 declaration or credit not verified")
    return {"photo_id": selected["photo_id"], "license_code": "cc0", "license_url": CC0_URL,
            "attribution": matches[0]["attribution"], "photo_page": selected["photo_page"]}


def exact_integer(value) -> int:
    number = float(value)
    if not np.isfinite(number) or number < 0 or not number.is_integer():
        raise ValueError("Expected pixel/palette count must be an exact nonnegative integer")
    return int(number)


def make_crop(image: Image.Image, measured: dict, expected: dict):
    from fcp_pipeline.photo_first_measurement import nearest_palette_counts

    if measured.get("automated_colour_state_status") != "automated_colour_state_admitted":
        raise ValueError("Original ROI admission not reproduced")
    oriented = ImageOps.exif_transpose(image).convert("RGB")
    rgb = np.asarray(oriented, dtype=np.uint8)
    flower, background = np.asarray(measured["flower_mask"]), np.asarray(measured["background_mask"])
    if flower.dtype != bool or background.dtype != bool or flower.shape != rgb.shape[:2] or background.shape != flower.shape:
        raise ValueError("Expected original-resolution boolean masks")
    if np.any(flower & background):
        raise ValueError("Flower/background masks overlap")
    if np.count_nonzero(flower) != exact_integer(expected["expected_mask_pixels"]) or np.count_nonzero(background) != exact_integer(expected["expected_background_pixels"]):
        raise ValueError("Original mask pixel counts not reproduced")
    counts = nearest_palette_counts(rgb[flower])
    if set(expected["expected_palette"]) != set(PALETTE) or any(int(counts[k]) != exact_integer(expected["expected_palette"][k]) for k in PALETTE):
        raise ValueError("Original flower palette not reproduced")
    ys, xs = np.where(flower)
    if not len(xs):
        raise ValueError("Empty flower mask")
    bbox = [int(xs.min()), int(ys.min()), int(xs.max()) + 1, int(ys.max()) + 1]
    rgba = np.dstack([rgb, flower.astype(np.uint8) * 255])
    # Remove even hidden non-flower RGB: publication crops cannot retain background pixels.
    rgba[~flower] = 0
    crop = Image.fromarray(rgba).crop(tuple(bbox))
    return crop, {"oriented_size": list(oriented.size), "bbox_xyxy": bbox,
                  "full_flower_mask_sha256": sha256(flower.astype(np.uint8).tobytes()),
                  "flower_mask_pixels": int(len(xs)), "palette_counts_reproduced": True}

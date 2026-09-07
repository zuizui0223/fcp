#!/usr/bin/env python3
"""Plan, replay, or render a fixed discovery-only CC0 flower photo strip."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import io
import importlib.metadata
import json
from pathlib import Path
import time
import urllib.request

import numpy as np
from PIL import Image

from fcp_pipeline.rgfca_photo_bar import (
    PALETTE, SLOTS, canonical, current_photo_rights, exact_integer, make_crop,
    make_plan, sha256, validate_photo_url,
)


def checked_plan(path):
    raw = path.read_bytes()
    plan = json.loads(raw)
    if raw != canonical(plan) or plan != make_plan():
        raise ValueError("Display plan differs from deterministic fixed-source selection")
    return plan, sha256(raw)


def fetch_bytes(url, limit):
    class NoRedirect(urllib.request.HTTPRedirectHandler):
        def redirect_request(self, req, fp, code, msg, headers, newurl):
            raise ValueError("Unexpected redirect from the frozen publication source")

    request = urllib.request.Request(url, headers={"User-Agent": "FCP-RGFCA-publication-photo-audit/1.0"})
    with urllib.request.build_opener(NoRedirect()).open(request, timeout=60) as response:
        if response.status != 200:
            raise ValueError("Unexpected HTTP response")
        data = response.read(limit + 1)
    if len(data) > limit:
        raise ValueError("Response exceeds bounded publication acquisition")
    return data


def safe_credit(text):
    for char in "\\`*_{}[]<>#":
        text = text.replace(char, "\\" + char)
    return text.replace("\n", " ").replace("\r", " ")


def acquire(plan_path, output, detector_weight, efficient_sam_dir):
    from fcp_pipeline.flower_roi_v4 import validate_roi_v4_contract
    from fcp_pipeline.flower_roi_v4_runtime import (
        FrozenFlowerColourEstimator, file_sha256, validate_scaleout_authorization,
    )

    plan, plan_sha = checked_plan(plan_path)
    root = Path(__file__).resolve().parents[2]
    execution = json.loads((root / "docs/supporting/random_photo_first_measurement_execution_v1.json").read_text())
    roi = json.loads((root / "docs/supporting/jbi_atlas_roi_estimator_contract_v4.json").read_text())
    locked = json.loads((root / "data/atlas/qualification/roi_v4_locked_test/jrc_roi_v4_locked_test_result.json").read_text())
    validate_roi_v4_contract(roi)
    detector_sha = file_sha256(detector_weight)
    if detector_sha != execution["runtime"]["detector_sha256"]:
        raise ValueError("Detector weight differs from fixed measurement")
    validate_scaleout_authorization(locked, trained_weight_sha256=detector_sha)
    if output.exists() and any(output.iterdir()):
        raise ValueError("Publication output must be new; preserve previous attempts")
    output.mkdir(parents=True, exist_ok=True)
    result = {"protocol": plan["protocol"], "source_commit": plan["source_commit"],
              "plan_sha256": plan_sha, "started_at_utc": datetime.now(timezone.utc).isoformat(),
              "status": "publication_rights_verification", "rights": [], "photos": [],
              "new_inference": False, "reserve_outcomes_read": False,
              "source_images_persisted": False, "replacement_used": False,
              "detector_sha256": detector_sha}
    result["runtime_versions"] = {name: importlib.metadata.version(name) for name in
                                  ("numpy", "pandas", "Pillow", "torch", "ultralytics", "onnxruntime", "scikit-image")}
    result_path = output / "photo_bar_execution.json"

    # Complete the metadata-only rights gate before opening any source image.
    for row in plan["selected"]:
        endpoint = f"https://api.inaturalist.org/v1/observations/{row['observation_id']}"
        try:
            payload = fetch_bytes(endpoint, 5 * 1024 * 1024)
            rights = current_photo_rights(json.loads(payload), row)
            rights.update({"status": "cc0_photo_verified", "endpoint": endpoint,
                           "response_sha256": sha256(payload),
                           "verified_at_utc": datetime.now(timezone.utc).isoformat()})
        except Exception as error:
            rights = {"photo_id": row["photo_id"], "status": "rights_not_verified",
                      "error_type": type(error).__name__}
        result["rights"].append(rights)
        result_path.write_bytes(canonical(result))
        time.sleep(1.1)
    if any(r["status"] != "cc0_photo_verified" for r in result["rights"]):
        result["status"] = "photo_bar_not_evaluable_rights_gate_failed_before_pixels"
        result_path.write_bytes(canonical(result))
        return 2

    estimator = FrozenFlowerColourEstimator(detector_weight, efficient_sam_dir, roi,
                                           torch_threads=int(execution["runtime"]["torch_threads_per_worker"]))
    crops_dir = output / "crops"
    crops_dir.mkdir()
    for row in plan["selected"]:
        receipt = {"slot": row["slot"], "photo_id": row["photo_id"], "measurement_id": row["measurement_id"]}
        try:
            validate_photo_url(row["photo_url_large"], row["photo_id"])
            raw = fetch_bytes(row["photo_url_large"], 20 * 1024 * 1024)
            if sha256(raw) != row["image_sha256"]:
                raise ValueError("Original image SHA-256 not reproduced")
            with Image.open(io.BytesIO(raw)) as image:
                measured = estimator.measure(image)
                crop, geometry = make_crop(image, measured, row)
            name = f"slot_{row['slot']:02d}.png"
            path = crops_dir / name
            # Reconstructed, masked RGB only. No source EXIF or hidden background RGB.
            crop.save(path, format="PNG")
            receipt.update(geometry)
            receipt.update({"status": "original_summaries_reproduced_crop_created",
                            "source_image_sha256": sha256(raw), "crop_name": name,
                            "crop_sha256": sha256(path.read_bytes())})
        except Exception as error:
            receipt.update({"status": "display_reconstruction_failed_no_replacement",
                            "error_type": type(error).__name__, "error": str(error)[:300]})
        result["photos"].append(receipt)
        result_path.write_bytes(canonical(result))
        print(f"photo_bar_slot={row['slot']}/{SLOTS} {receipt['status']}", flush=True)
        time.sleep(1.1)
    complete = len(result["photos"]) == SLOTS and all(r["status"] == "original_summaries_reproduced_crop_created" for r in result["photos"])
    result["status"] = "complete_reconstructed_cc0_discovery_photo_bar_inputs" if complete else "photo_bar_not_evaluable_incomplete_reproduction"
    result["finished_at_utc"] = datetime.now(timezone.utc).isoformat()
    result_path.write_bytes(canonical(result))
    if complete:
        credits = ["# RGFCA discovery photo bar: source credits", "",
                   "Illustrative CC0 discovery sample; not independent validation or a global flower-colour frequency estimate.", "",
                   "Changes: reconstructed flower mask, transparent non-flower pixels, bounding-rectangle crop and display resizing. Original RGB is retained within the mask. No photographer endorsement is implied.", ""]
        for row, rights in zip(plan["selected"], result["rights"], strict=True):
            credits.append(f"- {row['slot']:02d}: {safe_credit(rights['attribution'])}. [Photo {row['photo_id']}]({row['photo_page']}); [observation]({row['observation_page']}); [CC0]({row['license_url']}).")
        (output / "RGFCA_PHOTO_BAR_CREDITS.md").write_text("\n".join(credits) + "\n", encoding="utf-8")
    return 0 if complete else 2


def load_render_inputs(plan_path, inputs):
    """Verify the complete receipt and every crop before producing a figure."""
    from fcp_pipeline.photo_first_measurement import nearest_palette_counts

    plan, plan_sha = checked_plan(plan_path)
    result = json.loads((inputs / "photo_bar_execution.json").read_bytes())
    if result["status"] != "complete_reconstructed_cc0_discovery_photo_bar_inputs" or result["plan_sha256"] != plan_sha:
        raise ValueError("Complete exact-plan photo bar inputs required")
    if len(result["photos"]) != SLOTS or len(result["rights"]) != SLOTS:
        raise ValueError("Photo bar input census differs")
    if result["source_commit"] != plan["source_commit"] or result["protocol"] != plan["protocol"] or any(
        result.get(key) is not False for key in
        ("new_inference", "reserve_outcomes_read", "source_images_persisted", "replacement_used")
    ):
        raise ValueError("Unexpected publication scope or source")
    crops = []
    for row, receipt, rights in zip(plan["selected"], result["photos"], result["rights"], strict=True):
        if receipt["slot"] != row["slot"] or receipt["photo_id"] != row["photo_id"] or rights["photo_id"] != row["photo_id"] or rights["status"] != "cc0_photo_verified":
            raise ValueError("Slot/photo/rights identity mismatch")
        if rights["license_code"] != "cc0" or rights["license_url"] != row["license_url"] or rights["photo_page"] != row["photo_page"] or not rights["attribution"].strip():
            raise ValueError("Photo-level CC0 receipt differs")
        if receipt["status"] != "original_summaries_reproduced_crop_created" or receipt["measurement_id"] != row["measurement_id"] or receipt["source_image_sha256"] != row["image_sha256"] or receipt["palette_counts_reproduced"] is not True:
            raise ValueError("Original measurement reproduction receipt differs")
        name = f"slot_{row['slot']:02d}.png"
        if receipt["crop_name"] != name:
            raise ValueError("Unexpected crop path")
        raw = (inputs / "crops" / name).read_bytes()
        if sha256(raw) != receipt["crop_sha256"]:
            raise ValueError("Crop SHA-256 differs")
        with Image.open(io.BytesIO(raw)) as image:
            if image.mode != "RGBA":
                raise ValueError("Expected lossless RGBA crop")
            pixels = np.asarray(image)
        opaque = pixels[:, :, 3] == 255
        if not np.isin(pixels[:, :, 3], [0, 255]).all() or np.any(pixels[~opaque]):
            raise ValueError("Crop contains hidden RGB or non-boolean alpha")
        if np.count_nonzero(opaque) != exact_integer(row["expected_mask_pixels"]):
            raise ValueError("Crop mask census differs")
        if receipt["flower_mask_pixels"] != exact_integer(row["expected_mask_pixels"]):
            raise ValueError("Receipt mask census differs")
        width, height = receipt["oriented_size"]
        x0, y0, x1, y1 = receipt["bbox_xyxy"]
        if any(type(v) is not int for v in (width, height, x0, y0, x1, y1)) or not (0 <= x0 < x1 <= width <= 8192 and 0 <= y0 < y1 <= height <= 8192):
            raise ValueError("Invalid original image/crop geometry")
        if pixels.shape[:2] != (y1 - y0, x1 - x0):
            raise ValueError("Crop size and original bounding rectangle differ")
        full_mask = np.zeros((height, width), dtype=np.uint8)
        full_mask[y0:y1, x0:x1] = opaque
        if sha256(full_mask.tobytes()) != receipt["full_flower_mask_sha256"]:
            raise ValueError("Reconstructed crop mask digest differs")
        counts = nearest_palette_counts(pixels[:, :, :3][opaque])
        if any(counts[k] != exact_integer(row["expected_palette"][k]) for k in PALETTE):
            raise ValueError("Crop does not retain original measured RGB palette")
        crops.append(pixels)
    return plan, plan_sha, crops


def render(plan_path, inputs, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plan, plan_sha, crops = load_render_inputs(plan_path, inputs)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "text.color": "#30343A", "pdf.fonttype": 42})
    fig, axes = plt.subplots(2, 12, figsize=(12, 5), facecolor="white")
    for ax, row, pixels in zip(axes.flat, plan["selected"], crops, strict=True):
        ax.set_facecolor("#EEEEEE")
        ax.imshow(pixels, interpolation="nearest")
        height, width = pixels.shape[:2]
        side = max(height, width)
        ax.set_xlim((width - side) / 2 - .5, (width + side) / 2 - .5)
        ax.set_ylim((height + side) / 2 - .5, (height - side) / 2 - .5)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        lon = row["longitude"]
        label = ax.set_title(f"{row['slot']:02d}  {abs(lon):.0f}°{'W' if lon < 0 else 'E'}", fontsize=8, loc="left", pad=5)
        label.set_url(row["photo_page"])
    fig.suptitle("Flower regions from discovery photographs", x=.025, y=.97, ha="left", fontsize=15)
    fig.text(.025, .89, "24 CC0 photos · west-to-east order across rows · unequal geographic spacing · species names suppressed", fontsize=9)
    fig.text(.025, .065, "Original RGB within reconstructed flower masks; non-flower pixels transparent. Cropped/resized for display only.", fontsize=8)
    fig.text(.025, .025, "Illustrative sample, not global colour frequencies or independent validation. Numbered sources: RGFCA_PHOTO_BAR_CREDITS.md.", fontsize=8)
    fig.subplots_adjust(left=.025, right=.99, top=.82, bottom=.15, wspace=.20, hspace=.35)
    output.mkdir(parents=True, exist_ok=True)
    outputs = []
    for suffix in ("png", "pdf"):
        path = output / f"rgfca_figure3_discovery_photo_bar.{suffix}"
        meta = {"Software": "FCP RGFCA fixed photo bar"} if suffix == "png" else {"CreationDate": None, "ModDate": None, "Creator": "FCP RGFCA fixed photo bar"}
        fig.savefig(path, dpi=240, metadata=meta)
        outputs.append({"name": path.name, "sha256": sha256(path.read_bytes()), "bytes": path.stat().st_size})
    plt.close(fig)
    return {"plan_sha256": plan_sha, "execution_sha256": sha256((inputs / "photo_bar_execution.json").read_bytes()), "outputs": outputs,
            "status": "rendered_complete_discovery_photo_bar", "new_inference": False, "reserve_outcomes_read": False}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["plan", "acquire", "render"])
    ap.add_argument("--plan", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path)
    ap.add_argument("--inputs-dir", type=Path)
    ap.add_argument("--detector-weight", type=Path)
    ap.add_argument("--efficient-sam-dir", type=Path)
    args = ap.parse_args()
    if args.stage == "plan":
        plan = make_plan()
        args.plan.parent.mkdir(parents=True, exist_ok=True)
        args.plan.write_bytes(canonical(plan))
        print(json.dumps({"plan_sha256": sha256(canonical(plan)), "slots": SLOTS,
                          "eligible_cc0_photos": plan["eligible_cc0_photos_with_credit"], "pixels_opened": False}))
        return 0
    if args.stage == "acquire":
        return acquire(args.plan, args.output_dir, args.detector_weight, args.efficient_sam_dir)
    manifest = render(args.plan, args.inputs_dir, args.output_dir)
    (args.output_dir / "photo_bar_figure_manifest.json").write_bytes(canonical(manifest))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

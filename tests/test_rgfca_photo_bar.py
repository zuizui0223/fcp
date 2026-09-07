"""Photo-bar guards use metadata and synthetic owned pixels, never live photos."""
import copy
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

import numpy as np
import pytest
from PIL import Image

from fcp_pipeline import rgfca_photo_bar as bar
from fcp_pipeline.photo_first_measurement import nearest_palette_counts
from scripts.analysis import make_rgfca_photo_bar as cli


@pytest.fixture(scope="module")
def sources():
    return bar.discovery_sources()


def test_frozen_plan_metadata_census_and_reproduction(sources):
    pool, _, _ = sources
    plan = bar.make_plan()
    saved = json.loads((bar.ROOT / "docs/supporting/rgfca_photo_bar_plan_v1.json").read_bytes())
    assert plan == saved
    assert plan["eligible_cc0_photos_with_credit"] == 659
    selected = plan["selected"]
    assert len(selected) == len({r["photo_id"] for r in selected}) == 24
    assert len({r["species"] for r in selected}) == len({r["observer_id"] for r in selected}) == 24
    assert all(r["photo_license"] == "cc0" for r in selected)
    assert all(r["photo_id"] in set(pool.photo_id) for r in selected)
    assert not plan["new_inference"] and not plan["reserve_outcomes_read"] and not plan["replacement_allowed"]


def test_selection_is_row_order_independent_and_rejects_colour_columns(sources):
    pool, _, _ = sources
    meta = pool[bar.METADATA]
    assert bar.select_metadata(meta) == bar.select_metadata(meta.iloc[::-1])
    with pytest.raises(ValueError, match="metadata only"):
        bar.select_metadata(pool)


def test_selection_does_not_replace_an_empty_metadata_bin(sources):
    pool, _, _ = sources
    meta = pool[bar.METADATA].copy()
    meta.loc[:, "observer_id"] = "1"
    with pytest.raises(ValueError, match="No metadata-admissible"):
        bar.select_metadata(meta)


@pytest.mark.parametrize("url", [
    "http://static.inaturalist.org/photos/123/large.jpg",
    "https://example.org/photos/123/large.jpg",
    "https://static.inaturalist.org/photos/999/large.jpg",
    "https://static.inaturalist.org/photos/123/original.jpg",
    "https://static.inaturalist.org/photos/123/large.jpg?override=1",
    "https://user:secret@static.inaturalist.org/photos/123/large.jpg",
])
def test_unfrozen_image_urls_rejected(url):
    with pytest.raises(ValueError):
        bar.validate_photo_url(url, "123")


def test_photo_level_rights_not_observation_rights():
    selected = {"observation_id": "10", "photo_id": "20", "photo_page": "https://www.inaturalist.org/photos/20"}
    response = {"results": [{"id": 10, "license_code": "cc0", "photos": [{"id": 20, "license_code": "cc-by-nc", "attribution": "Example"}]}]}
    with pytest.raises(ValueError, match="photo-level CC0"):
        bar.current_photo_rights(response, selected)
    response["results"][0]["photos"][0]["license_code"] = "cc0"
    assert bar.current_photo_rights(response, selected)["attribution"] == "Example"
    response["results"][0]["photos"][0]["id"] = 21
    with pytest.raises(ValueError):
        bar.current_photo_rights(response, selected)


@pytest.mark.parametrize("value", ["1.1", "-1", "nan", "inf", ""])
def test_numeric_reproduction_never_rounds(value):
    with pytest.raises(ValueError):
        bar.exact_integer(value)


@pytest.fixture
def synthetic():
    rgb = np.zeros((6, 7, 3), dtype=np.uint8)
    rgb[:] = (200, 10, 20)
    mask = np.zeros((6, 7), dtype=bool)
    mask[1:5, 2:6] = True
    mask[2, 3] = False
    background = ~mask
    expected = {"expected_mask_pixels": "15.0", "expected_background_pixels": "27.0",
                "expected_palette": nearest_palette_counts(rgb[mask])}
    measured = {"flower_mask": mask, "background_mask": background,
                "automated_colour_state_status": "automated_colour_state_admitted"}
    return Image.fromarray(rgb), measured, expected


def test_crop_keeps_original_flower_rgb_and_no_hidden_background(synthetic):
    image, measured, expected = synthetic
    crop, receipt = bar.make_crop(image, measured, expected)
    assert crop.mode == "RGBA" and crop.size == (4, 4)
    assert receipt["bbox_xyxy"] == [2, 1, 6, 5]
    assert receipt["flower_mask_pixels"] == 15
    pixels = np.asarray(crop)
    assert not np.any(pixels[pixels[:, :, 3] == 0])
    assert np.all(pixels[pixels[:, :, 3] == 255][:, :3] == [200, 10, 20])


@pytest.mark.parametrize("fault", ["not_admitted", "overlap", "pixels", "palette", "not_boolean"])
def test_crop_reproduction_failures_cannot_be_presented(synthetic, fault):
    image, measured, expected = copy.deepcopy(synthetic)
    if fault == "not_admitted":
        measured["automated_colour_state_status"] = "uncertain"
    elif fault == "overlap":
        measured["background_mask"][:] = True
    elif fault == "pixels":
        expected["expected_mask_pixels"] = "14.0"
    elif fault == "palette":
        expected["expected_palette"]["white"] += 1
    else:
        measured["flower_mask"] = measured["flower_mask"].astype(int)
    with pytest.raises(ValueError):
        bar.make_crop(image, measured, expected)


@pytest.fixture
def synthetic_render_inputs(synthetic, tmp_path, monkeypatch):
    """Only owned synthetic pixels; mock the plan boundary, not crop validation."""
    image, measured, expected = synthetic
    crop, geometry = bar.make_crop(image, measured, expected)
    inputs = tmp_path / "synthetic-inputs"
    (inputs / "crops").mkdir(parents=True)
    plan = {"source_commit": "synthetic-only", "protocol": "synthetic-only", "selected": []}
    receipt = {"source_commit": "synthetic-only", "protocol": "synthetic-only",
               "plan_sha256": "synthetic-only", "photos": [], "rights": [],
               "status": "complete_reconstructed_cc0_discovery_photo_bar_inputs",
               "new_inference": False, "reserve_outcomes_read": False,
               "source_images_persisted": False, "replacement_used": False}
    for slot in range(1, 25):
        row = {**expected, "slot": slot, "photo_id": str(slot), "measurement_id": str(slot),
               "image_sha256": bar.sha256(image.tobytes()), "longitude": slot - 12.5,
               "license_url": bar.CC0_URL, "photo_page": f"https://example.invalid/{slot}"}
        plan["selected"].append(row)
        name = f"slot_{slot:02d}.png"
        crop.save(inputs / "crops" / name)
        receipt["photos"].append({**geometry, "slot": slot, "photo_id": str(slot),
                                  "measurement_id": str(slot), "source_image_sha256": row["image_sha256"],
                                  "status": "original_summaries_reproduced_crop_created",
                                  "crop_name": name, "crop_sha256": bar.sha256((inputs / "crops" / name).read_bytes())})
        receipt["rights"].append({"photo_id": str(slot), "status": "cc0_photo_verified",
                                  "license_code": "cc0", "license_url": bar.CC0_URL,
                                  "photo_page": row["photo_page"], "attribution": "Synthetic test only"})
    (inputs / "photo_bar_execution.json").write_bytes(bar.canonical(receipt))
    monkeypatch.setattr(cli, "checked_plan", lambda path: (plan, "synthetic-only"))
    return inputs, receipt


def test_complete_synthetic_strip_renders_identically_twice(synthetic_render_inputs, tmp_path):
    inputs, _ = synthetic_render_inputs
    first = cli.render(Path("synthetic"), inputs, tmp_path / "first")
    second = cli.render(Path("synthetic"), inputs, tmp_path / "second")
    assert first == second
    assert [r["name"] for r in first["outputs"]] == [
        "rgfca_figure3_discovery_photo_bar.png", "rgfca_figure3_discovery_photo_bar.pdf"]
    assert not first["new_inference"] and not first["reserve_outcomes_read"]


@pytest.mark.parametrize("fault", ["partial", "license", "failed_photo", "source_sha", "mask_sha", "bbox", "hidden_rgb", "scope"])
def test_incomplete_or_tampered_inputs_never_render(synthetic_render_inputs, tmp_path, fault):
    inputs, receipt = synthetic_render_inputs
    if fault == "partial":
        receipt["photos"].pop()
    elif fault == "license":
        receipt["rights"][0]["license_code"] = "cc-by-nc"
    elif fault == "failed_photo":
        receipt["photos"][0]["status"] = "display_reconstruction_failed_no_replacement"
    elif fault == "source_sha":
        receipt["photos"][0]["source_image_sha256"] = "0" * 64
    elif fault == "mask_sha":
        receipt["photos"][0]["full_flower_mask_sha256"] = "0" * 64
    elif fault == "bbox":
        receipt["photos"][0]["bbox_xyxy"][0] = -1
    elif fault == "scope":
        receipt["reserve_outcomes_read"] = True
    else:
        path = inputs / "crops/slot_01.png"
        with Image.open(path) as image:
            pixels = np.array(image)
        pixels[pixels[:, :, 3] == 0] = (255, 0, 0, 0)
        Image.fromarray(pixels).save(path)
        # Even updating the file digest must not let hidden RGB pass.
        receipt["photos"][0]["crop_sha256"] = bar.sha256(path.read_bytes())
    (inputs / "photo_bar_execution.json").write_bytes(bar.canonical(receipt))
    with pytest.raises(ValueError):
        cli.render(Path("synthetic"), inputs, tmp_path / "blocked-output")
    assert not (tmp_path / "blocked-output").exists()


def test_committed_display_has_complete_source_crop_and_figure_provenance():
    release = json.loads((bar.ROOT / "docs/supporting/rgfca_photo_bar_release_v1.json").read_bytes())
    for field in ("execution", "figure_manifest", "credits"):
        assert bar.sha256((bar.ROOT / release[field + "_path"]).read_bytes()) == release[field + "_sha256"]
    inputs = (bar.ROOT / release["execution_path"]).parent
    plan, plan_sha, crops = cli.load_render_inputs(bar.ROOT / "docs/supporting/rgfca_photo_bar_plan_v1.json", inputs)
    assert plan_sha == release["plan_sha256"]
    assert len(crops) == release["verified_photo_count"] == release["current_photo_level_cc0_count"] == 24
    manifest = json.loads((bar.ROOT / release["figure_manifest_path"]).read_bytes())
    assert manifest["execution_sha256"] == release["execution_sha256"]
    assert manifest["plan_sha256"] == plan_sha
    for row in manifest["outputs"]:
        raw = (bar.ROOT / "docs/figures" / row["name"]).read_bytes()
        assert len(raw) == row["bytes"] and bar.sha256(raw) == row["sha256"]
    credits = (bar.ROOT / release["credits_path"]).read_text(encoding="utf-8")
    assert len([line for line in credits.splitlines() if line.startswith("- ")]) == 24
    assert all(row["photo_page"] in credits and row["observation_page"] in credits for row in plan["selected"])
    for flag in ("replacement_used", "reserve_outcomes_read", "new_inference",
                 "independent_validation_claim_allowed", "segmentation_accuracy_validated_by_display"):
        assert release[flag] is False


def test_real_committed_crops_render_twice_without_reacquisition(tmp_path, monkeypatch):
    monkeypatch.setattr(cli, "fetch_bytes", lambda *args: pytest.fail("Rendering must not reacquire photographs"))
    plan = bar.ROOT / "docs/supporting/rgfca_photo_bar_plan_v1.json"
    inputs = bar.ROOT / "docs/figures/rgfca_photo_bar_v1"
    assert cli.render(plan, inputs, tmp_path / "a") == cli.render(plan, inputs, tmp_path / "b")


def test_actual_publication_export_command_without_inherited_pythonpath(tmp_path):
    workflow = (bar.ROOT / ".github/workflows/rgfca-publication-figures.yml").read_text(encoding="utf-8")
    commands = [line.strip() for line in workflow.splitlines()
                if line.strip().startswith("python ") and "make_rgfca_photo_bar" in line and " render " in line]
    assert len(commands) == 1
    args = shlex.split(commands[0])
    args[0] = sys.executable
    output = tmp_path / "actual-cli-export"
    args[args.index("--output-dir") + 1] = str(output)
    env = dict(os.environ)
    env.pop("PYTHONPATH", None)
    env["PYTHONUTF8"] = "1"
    result = subprocess.run(args, cwd=bar.ROOT, env=env, capture_output=True, text=True, encoding="utf-8", timeout=60)
    assert result.returncode == 0, result.stderr
    manifest = json.loads((output / "photo_bar_figure_manifest.json").read_bytes())
    assert manifest["status"] == "rendered_complete_discovery_photo_bar"
    assert len(manifest["outputs"]) == 2 and not manifest["reserve_outcomes_read"]

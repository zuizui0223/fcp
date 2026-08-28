import importlib.util
import io
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "data" / "materialize_jbi_ch1_calibration_images.py"
spec = importlib.util.spec_from_file_location("calibration_image_materializer", SCRIPT)
assert spec is not None and spec.loader is not None
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def make_split() -> pd.DataFrame:
    rows = []
    for s in range(6):
        for i in range(200):
            rows.append(
                {
                    "species": f"Species {s+1}",
                    "photo_id": f"p{s+1:02d}_{i+1:03d}",
                    "photo_url": f"https://example.invalid/{s+1}/{i+1}.jpg",
                    "photo_url_api": f"https://fallback.invalid/{s+1}/{i+1}.jpg",
                    "latitude": s,
                    "longitude": i,
                    "observer": "hidden",
                    "split": "calibration" if i < 80 else "evaluation",
                }
            )
    return pd.DataFrame(rows)


def jpeg_bytes() -> bytes:
    y, x = np.mgrid[0:320, 0:320]
    arr = np.stack(
        [
            (x % 256).astype(np.uint8),
            (y % 256).astype(np.uint8),
            ((x + y) % 256).astype(np.uint8),
        ],
        axis=2,
    )
    image = Image.fromarray(arr, mode="RGB")
    buf = io.BytesIO()
    image.save(buf, format="JPEG", quality=90)
    return buf.getvalue()


def test_calibration_rows_never_include_evaluation():
    rows = mod.calibration_rows(make_split())
    assert len(rows) == 480
    assert rows.groupby("species").size().eq(80).all()
    assert set(rows["split"]) == {"calibration"}


def test_blind_id_is_deterministic_and_species_conditioned():
    a = mod.blind_id("Species A", "123")
    b = mod.blind_id("Species A", "123")
    c = mod.blind_id("Species B", "123")
    assert a == b
    assert a != c
    assert len(a) == 16


def test_technical_metrics_are_finite():
    image = Image.open(io.BytesIO(jpeg_bytes())).convert("RGB")
    metrics = mod.technical_metrics(image)
    assert metrics["width_px"] == 320
    assert metrics["height_px"] == 320
    for key in (
        "mean_luminance",
        "sd_luminance",
        "central_sd_luminance",
        "dark_clip_fraction",
        "bright_clip_fraction",
        "mean_saturation",
        "p90_saturation",
        "edge_energy",
    ):
        assert np.isfinite(metrics[key])


def test_technical_flags_do_not_claim_biological_evaluability():
    flags = mod.technical_flags(
        {
            "min_dimension_px": 100,
            "dark_clip_fraction": 0.8,
            "bright_clip_fraction": 0.0,
            "edge_energy": 0.0,
        }
    )
    assert "low_resolution" in flags
    assert "severe_dark_clipping" in flags
    assert "very_low_detail" in flags
    assert all("flower" not in flag and "evaluable" not in flag for flag in flags)


def test_full_materialization_uses_exactly_480_frozen_rows(tmp_path, monkeypatch):
    split_path = tmp_path / "split.csv"
    make_split().to_csv(split_path, index=False)
    payload = jpeg_bytes()

    def fake_download(urls, *, retries=3, pause_seconds=0.5):
        urls = list(urls)
        assert urls
        return payload, urls[0]

    monkeypatch.setattr(mod, "download_bytes", fake_download)
    manifest = mod.materialize(
        split_path,
        tmp_path / "images",
        tmp_path / "qc.csv",
        tmp_path / "manifest.json",
        tmp_path / "contact",
        retries=1,
        pause_seconds=0,
    )
    assert manifest["status"] == "pass"
    assert manifest["n_expected"] == 480
    assert manifest["n_materialized"] == 480
    assert manifest["n_failed"] == 0
    assert manifest["evaluation_rows_opened"] is False
    assert len(manifest["contact_sheets"]) == 24

    qc = pd.read_csv(tmp_path / "qc.csv")
    assert len(qc) == 480
    assert qc["blind_id"].is_unique
    assert qc.groupby("species").size().eq(80).all()
    assert not any("latitude" in c or "longitude" in c or "observer" in c for c in qc.columns)

"""Purely synthetic source gates for exact 149900-method fourth cohort."""
import hashlib
import json

import pandas as pd
import pytest

from scripts.analysis import run_fcp_wave4_original_roi_and_statistics_20261017 as m


def source(tmp_path, *, authorized=True):
    p = tmp_path / "source"
    (p / "selection").mkdir(parents=True)
    (p / "metadata").mkdir()
    rows = []
    for s in range(300):
        for k in range(4):
            pid = s * 4 + k + 1
            rows.append({
                "inat_taxon_id": 20000+s, "query_species": f"Testus {s}", "species": f"Testus {s}",
                "photo_id": pid, "observation_id": 100000 + pid, "photo_url_large": f"https://example.org/{pid}.jpg",
                "photo_license": "cc-by", "latitude": float(s % 90), "longitude": float((s % 300)-150),
                "observer_id": f"observer-{s}-{k}", "observed_on": "2026-05-12",
            })
    data = p / "metadata/wave4_authorized_metadata.csv.gz"
    pd.DataFrame(rows).to_csv(data, index=False, compression={"method":"gzip","mtime":0}, lineterminator="\n")
    (p / "selection/selection_result.json").write_text(json.dumps({
        "wave4_selected_species": 500, "remaining_candidate_species": 2730,
        "new_selection_manifest_sha256": "synthetic-source-sha"
    }))
    (p / "metadata/metadata_result.json").write_text(json.dumps({
        "selection_sha256": "synthetic-source-sha", "metadata_gate_pass": authorized,
        "image_pixels_opened": False, "new_selected_species": 500, "original_photos": 149900,
        "new_exact_100_species": 300, "new_exact_100_photo_rows": 1200,
        "outputs": {"wave4_authorized_metadata.csv.gz": m.sha(data)}
    }))
    return p


def test_new_photos_are_blind_and_exactly_partitioned(tmp_path, monkeypatch):
    monkeypatch.setattr(m, "TARGET", 4)
    src = source(tmp_path)
    result = m.build_firewall(src, tmp_path / "fw")
    assert result["new_photo_rows"] == 1200
    assert result["partitions"] == 256
    assert not result["photo_pixels_opened"]
    assert result["blind_worker_columns"] == ["measurement_id","image_filename","photo_license"]
    file = tmp_path / "fw/batch_0/worker_packet/measurement_manifest.csv"
    assert list(pd.read_csv(file).columns) == result["blind_worker_columns"]
    assert "species" not in file.read_text() and "latitude" not in file.read_text()
    expected = set(pd.read_csv(tmp_path / "fw/batch_0/sealed_keys/acquisition_key.csv").measurement_id.astype(str))
    expected |= set(pd.read_csv(tmp_path / "fw/batch_1/sealed_keys/acquisition_key.csv").measurement_id.astype(str))
    assert len(expected) == 1200
    assert all(len(pd.read_csv(tmp_path / f"fw/batch_{b}/worker_packet/measurement_manifest.csv")) for b in (0,1))


def test_no_photo_opening_without_original_source_gate(tmp_path, monkeypatch):
    monkeypatch.setattr(m, "TARGET", 4)
    src = source(tmp_path, authorized=False)
    with pytest.raises(RuntimeError, match="unqualified"):
        m.build_firewall(src, tmp_path / "fw")


def test_detect_tampered_frozen_new_original_photo_csv(tmp_path, monkeypatch):
    monkeypatch.setattr(m, "TARGET", 4)
    src = source(tmp_path)
    f = src / "metadata/wave4_authorized_metadata.csv.gz"
    f.write_bytes(f.read_bytes() + b"tampered")
    with pytest.raises(RuntimeError, match="SHA drift"):
        m.authenticated_source(src)


def test_original_partition_geometry_stays_fixed():
    assert (m.N_BATCHES, m.N_SHARDS, m.N_PARTS, m.TOTAL_PARTS, m.SOURCE_RUN) == (2,32,4,256,38989666486)


def test_real_original_roi_colour_module_imports_before_pixel_download():
    """Catch missing original photo_first_atlas import chain before 256 ROI jobs."""
    import importlib
    from pathlib import Path
    modules = (
        "fcp_pipeline.photo_first_atlas",
        "fcp_pipeline.photo_first_atlas_v2",
        "fcp_pipeline.photo_first_measurement",
    )
    loaded = [importlib.import_module(mod) for mod in modules]
    assert all(Path(module.__file__).is_file() for module in loaded)
    assert callable(loaded[-1].coarse_morph_from_palette)


def test_old_50km_summary_string_key_is_reported_not_all_nclass40_species():
    # Historical cohort_summary() uses "50", NOT the float 50.0.
    summary = {
        "50": {"evaluable": True, "n_species": 17},
        "100": {"evaluable": True, "n_species": 23},
    }
    assert m.original_spatial_primary_50km_n(summary) == 17
    assert m.original_spatial_primary_50km_n({"50": {"evaluable": False, "n_species": 0}}) == 0
    assert m.original_spatial_primary_50km_n({}) == 0

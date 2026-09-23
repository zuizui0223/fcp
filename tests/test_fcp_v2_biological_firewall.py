from __future__ import annotations

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "analysis" / "build_fcp_v2_biological_firewall_20260923.py"


def load_module():
    spec = importlib.util.spec_from_file_location("fcp_v2_bio_firewall", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_pass_b_measurement_id_and_slot_are_deterministic() -> None:
    m = load_module()
    a = m.measurement_id(123456)
    b = m.measurement_id(123456)
    assert a == b
    assert a.startswith("FCPV2-")
    assert m.slot(a) == m.slot(b)
    batch, shard, part = m.slot(a)
    assert 0 <= batch < 2
    assert 0 <= shard < 32
    assert 0 <= part < 4


def test_pass_b_worker_schema_contains_source_sha_but_no_identity_leak() -> None:
    m = load_module()
    assert m.WORKER_FIELDS == (
        "measurement_id",
        "image_filename",
        "photo_license",
        "heavy_counterfactual",
        "expected_source_sha256",
    )
    forbidden = {
        "species",
        "inat_taxon_id",
        "photo_id",
        "observation_id",
        "latitude",
        "longitude",
        "observer_id",
        "observer",
        "D",
        "W",
        "q_white",
    }
    assert not (set(m.WORKER_FIELDS) & forbidden)

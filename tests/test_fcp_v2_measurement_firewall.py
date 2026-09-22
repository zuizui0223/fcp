from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "analysis" / "build_fcp_v2_measurement_firewall_20260923.py"


def load_module():
    spec = importlib.util.spec_from_file_location("fcp_v2_firewall", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_measurement_id_and_partition_are_deterministic() -> None:
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


def test_heavy_hash_uses_only_frozen_identity_fields() -> None:
    m = load_module()
    a = m.heavy_hash("P", 1001, 90001)
    b = m.heavy_hash("P", 1001, 90001)
    c = m.heavy_hash("P", 1001, 90002)
    assert a == b
    assert a != c


def test_worker_schema_has_no_biological_or_location_identity() -> None:
    m = load_module()
    assert m.WORKER_FIELDS == (
        "measurement_id",
        "image_filename",
        "photo_license",
        "heavy_counterfactual",
    )
    forbidden = ("species", "taxon", "photo_id", "latitude", "longitude", "observer", "morph", "q_white", "spatial_outcome")
    for field in m.WORKER_FIELDS:
        assert not any(token.casefold() in field.casefold() for token in forbidden)

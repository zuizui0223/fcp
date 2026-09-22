from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "analysis" / "select_fcp_v2_measurement_validity_species_20260923.py"


def load_module():
    spec = importlib.util.spec_from_file_location("fcp_v2_select", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_canonical_identity_reads_only_allowed_columns() -> None:
    m = load_module()
    frame = pd.DataFrame(
        {
            "inat_taxon_id": [2, 1],
            "species": ["Beta beta", "Alpha alpha"],
            "after_observer_cap": [220, 180],
            "morph": ["white", "blue_purple"],
            "D": [0.2, 0.8],
            "family": ["X", "Y"],
        }
    )
    out = m.canonical_identity(frame)
    assert list(out.columns) == list(m.ALLOWED)
    assert "morph" not in out.columns
    assert "D" not in out.columns
    assert "family" not in out.columns


def test_hash_rank_is_deterministic_and_input_order_invariant() -> None:
    m = load_module()
    frame = pd.DataFrame(
        {
            "inat_taxon_id": [4, 1, 3, 2],
            "species": ["D d", "A a", "C c", "B b"],
            "after_observer_cap": [250, 250, 250, 250],
        }
    )
    a = m.hash_rank(frame, salt=m.P_SALT)
    b = m.hash_rank(frame.iloc[::-1].reset_index(drop=True), salt=m.P_SALT)
    assert a["inat_taxon_id"].tolist() == b["inat_taxon_id"].tolist()
    assert a["selection_hash"].tolist() == b["selection_hash"].tolist()
    assert a["selection_rank"].tolist() == [1, 2, 3, 4]


def test_panel_salts_are_distinct() -> None:
    m = load_module()
    assert m.P_SALT != m.N_SALT
    frame = pd.DataFrame(
        {
            "inat_taxon_id": [1, 2],
            "species": ["A a", "B b"],
            "after_observer_cap": [300, 300],
        }
    )
    p = m.hash_rank(frame, salt=m.P_SALT)
    n = m.hash_rank(frame, salt=m.N_SALT)
    assert p["selection_hash"].tolist() != n["selection_hash"].tolist()


def test_canonical_identity_rejects_duplicate_species_or_taxa() -> None:
    m = load_module()
    dup_taxon = pd.DataFrame(
        {
            "inat_taxon_id": [1, 1],
            "species": ["A a", "B b"],
            "after_observer_cap": [200, 200],
        }
    )
    with pytest.raises(ValueError):
        m.canonical_identity(dup_taxon)

    dup_species = pd.DataFrame(
        {
            "inat_taxon_id": [1, 2],
            "species": ["A a", "A a"],
            "after_observer_cap": [200, 200],
        }
    )
    with pytest.raises(ValueError):
        m.canonical_identity(dup_species)

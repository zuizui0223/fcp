from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "analysis" / "run_fcp_v2_measurement_validity_summary_20260923.py"


def load_module():
    spec = importlib.util.spec_from_file_location("fcp_v2_mv", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_qwhite_is_unit_zero_sum() -> None:
    m = load_module()
    q = m._qwhite()
    assert np.isclose(q.sum(), 0.0, atol=1e-12)
    assert np.isclose(np.linalg.norm(q), 1.0, atol=1e-12)


def test_technical_displacement_exactly_on_qwhite_has_T_one() -> None:
    m = load_module()
    q = m._qwhite()
    p0 = np.full(9, 1.0 / 9.0)
    p1 = p0 + 0.05 * q
    assert np.all(p1 > 0)
    p1 = p1 / p1.sum()
    row = {
        "panel": "P",
        "species": "Synthetic species",
    }
    for j, c in enumerate(m.FRACTIONS):
        row[c] = p0[j]
        row[f"cf_{c}"] = p1[j]
    pairs = pd.DataFrame([row])
    summary, per = m._mv3a_for_scope(pairs, "all")
    assert summary["nonzero_technical_displacements"] == 1
    assert summary["T_white_image"]["mean"] == pytest.approx(1.0)
    assert per.loc[0, "mean_T_white"] == pytest.approx(1.0)


def test_species_D_half_for_equal_two_morphs() -> None:
    m = load_module()
    rows = []
    for i in range(40):
        rows.append(
            {
                "condition_id": "baseline",
                "panel": "P",
                "species": "Synthetic species",
                "morph": "white" if i < 20 else "red_pink",
                "measurement_status": "classified_four_state_morph",
            }
        )
    out = m._species_D(pd.DataFrame(rows), "baseline")
    assert len(out) == 1
    assert out.loc[0, "n_classifiable"] == 40
    assert out.loc[0, "D"] == pytest.approx(0.5)


def test_primary_vector_table_uses_fixed_white_axis_without_refit() -> None:
    m = load_module()
    rows = []
    for i in range(40):
        white = i < 20
        row = {
            "condition_id": "baseline",
            "panel": "P",
            "species": "Synthetic species",
            "morph": "white" if white else "red_pink",
            "measurement_status": "classified_four_state_morph",
        }
        for c in m.FRACTIONS:
            row[c] = 0.0
        if white:
            row["flower_fraction_white"] = 1.0
        else:
            row["flower_fraction_red"] = 1.0
        rows.append(row)
    out = m._vector_table(pd.DataFrame(rows), "baseline", "all")
    assert len(out) == 1
    assert out.loc[0, "n_classifiable"] == 40
    assert out.loc[0, "minor_cluster_fraction"] == pytest.approx(0.5)
    assert 0.0 <= out.loc[0, "projection_sq"] <= 1.0


def test_hellinger_identical_palette_is_zero() -> None:
    m = load_module()
    p = np.array([[0.5, 0.5] + [0.0] * 7], dtype=float)
    assert m._hellinger_rows(p, p)[0] == pytest.approx(0.0)

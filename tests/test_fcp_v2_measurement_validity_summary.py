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


def test_four_colour_collapse_matches_frozen_groups() -> None:
    m = load_module()
    row = {c: 0.0 for c in m.FRACTIONS}
    row["flower_fraction_white"] = 0.1
    row["flower_fraction_yellow"] = 0.1
    row["flower_fraction_orange"] = 0.1
    row["flower_fraction_bronze"] = 0.1
    row["flower_fraction_red"] = 0.1
    row["flower_fraction_pink"] = 0.1
    row["flower_fraction_magenta"] = 0.1
    row["flower_fraction_blue"] = 0.15
    row["flower_fraction_purple"] = 0.15
    x = m._four_colour_matrix(pd.DataFrame([row]))
    assert x.shape == (1, 4)
    assert x[0].tolist() == pytest.approx([0.1, 0.3, 0.3, 0.3])


def test_spatial_rho_detects_ordered_colour_change() -> None:
    m = load_module()
    latitude = np.array([0.0, 1.0, 2.0, 3.0])
    longitude = np.zeros(4)
    t = np.array([0.0, 0.25, 0.75, 1.0])
    traits = np.column_stack([1.0 - t, t, np.zeros(4), np.zeros(4)])
    rho = m._spatial_rho(latitude, longitude, traits)
    assert rho > 0.9


def test_spatial_pair_identity_has_zero_change() -> None:
    m = load_module()
    baseline = pd.DataFrame(
        {
            "condition_id": ["baseline", "baseline"],
            "panel": ["P", "P"],
            "species": ["A a", "B b"],
            "n_classifiable": [40, 40],
            "spatial_rho": [0.1, 0.4],
        }
    )
    other = baseline.copy()
    other["condition_id"] = "fixed_ev_p0_5"
    result = m._spatial_pair_summary(
        baseline,
        other,
        scope="P",
        condition_id="fixed_ev_p0_5",
    )
    assert result["species"] == 2
    assert result["spearman_rho"] == pytest.approx(1.0)
    assert result["lin_ccc"] == pytest.approx(1.0)
    assert result["signed_spatial_rho_change"]["mean"] == pytest.approx(0.0)

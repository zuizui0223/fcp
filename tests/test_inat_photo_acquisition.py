from collections import Counter

from scripts.data.acquire_jbi_ch1_inat_photos import select_rows, species_qc


def config():
    return {
        "selection": {
            "target_photographs_per_species": 200,
            "stable_hash_salt": "test-salt",
            "spatial_cell_degrees": 1.0,
            "hard_caps": {
                "maximum_per_observer": 20,
                "maximum_per_spatial_cell": 50,
                "maximum_per_calendar_month": 100,
            },
            "qc_gates": {
                "minimum_unique_observers": 10,
                "minimum_unique_spatial_cells": 3,
                "minimum_unique_calendar_months": 2,
                "maximum_observer_fraction": 0.25,
            },
        }
    }


def candidates(n=1200):
    rows = []
    for i in range(n):
        rows.append(
            {
                "species": "Species alpha",
                "observation_id": 100000 + i,
                "photo_id": str(500000 + i),
                "latitude": -30 + (i % 60) * 1.1,
                "longitude": -150 + (i % 120) * 2.1,
                "observer": f"observer_{i % 40}",
                "observer_id": i % 40,
                "observed_month": (i % 6) + 1,
            }
        )
    return rows


def assignment(rows):
    return [row["photo_id"] for row in rows]


def test_selection_is_exact_and_deterministic():
    source = candidates()
    a = select_rows(config(), "Species alpha", source)
    b = select_rows(config(), "Species alpha", list(reversed(source)))
    assert len(a) == 200
    assert assignment(a) == assignment(b)
    assert len(set(assignment(a))) == 200


def test_selection_respects_hard_caps():
    selected = select_rows(config(), "Species alpha", candidates())
    observers = Counter(row["observer"] for row in selected)
    cells = Counter(row["spatial_cell"] for row in selected)
    months = Counter(row["observed_month"] for row in selected)
    assert max(observers.values()) <= 20
    assert max(cells.values()) <= 50
    assert max(months.values()) <= 100


def test_qc_gate_passes_balanced_selection():
    selected = select_rows(config(), "Species alpha", candidates())
    qc = species_qc(
        config(),
        "Species alpha",
        {"id": 123},
        candidates(),
        selected,
    )
    assert qc["gate_pass"] is True
    assert qc["selected_count"] == 200
    assert qc["unique_observers"] >= 10
    assert qc["unique_spatial_cells"] >= 3
    assert qc["unique_calendar_months"] >= 2


def test_selection_fails_when_observer_cap_makes_200_impossible():
    bad = candidates(300)
    for row in bad:
        row["observer"] = "one_observer"
        row["observer_id"] = 1
    try:
        select_rows(config(), "Species alpha", bad)
    except RuntimeError as exc:
        assert "only" in str(exc)
    else:
        raise AssertionError("expected selection to fail closed")

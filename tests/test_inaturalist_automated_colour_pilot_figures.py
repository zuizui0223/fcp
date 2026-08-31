from __future__ import annotations

import numpy as np

from scripts.analysis.make_jbi_inaturalist_automated_colour_pilot_figures import (
    ADMITTED,
    display_rgb,
    mass_box,
    select_photo_bar_encounters,
)


def test_species_free_photo_bar_selection_is_deterministic_and_spans_longitude() -> None:
    rows = [
        {
            "encounter_status": ADMITTED,
            "longitude": str(value),
            "encounter_blind_id": f"id-{value:03d}",
        }
        for value in range(30)
    ]
    first = select_photo_bar_encounters(rows, count=6)
    second = select_photo_bar_encounters(list(reversed(rows)), count=6)
    assert first == second
    assert [float(row["longitude"]) for row in first] == [0, 6, 12, 17, 23, 29]


def test_mass_box_tracks_soft_weight_and_stays_in_bounds() -> None:
    weights = np.zeros((20, 30), dtype=float)
    weights[5:15, 10:20] = 1.0
    x0, y0, x1, y1 = mass_box(weights)
    assert 0 <= x0 < 10 < 20 < x1 <= 30
    assert 0 <= y0 < 5 < 15 < y1 <= 20


def test_display_rgb_is_bounded() -> None:
    colours = display_rgb(
        [{"flower_L_mean": "50", "flower_a_mean": "25", "flower_b_mean": "-30"}]
    )
    assert colours.shape == (1, 3)
    assert np.all((colours >= 0) & (colours <= 1))

from __future__ import annotations

import numpy as np

from scripts.analysis.make_jbi_atlas_species_free_figure import (
    central_mask_box,
    select_photo_bar_rows,
)


def test_photo_bar_selection_is_geometry_only_deterministic_and_complete() -> None:
    rows = [
        {
            "measurement_id": f"M-{index:03d}",
            "longitude": str(-179.0 + index * 3.0),
            "latitude": str((index % 17) - 8),
            "flower_L_mean": str(index),
            "flower_a_mean": "0",
            "flower_b_mean": "0",
        }
        for index in range(100)
    ]
    selected = select_photo_bar_rows(rows)
    reversed_selected = select_photo_bar_rows(list(reversed(rows)))
    assert len(selected) == 48
    assert [row["measurement_id"] for row in selected] == [
        row["measurement_id"] for row in reversed_selected
    ]
    assert selected[0]["measurement_id"] == "M-000"
    assert selected[-1]["measurement_id"] == "M-099"


def test_photo_bar_crop_uses_central_mask_mass_with_padding() -> None:
    mask = np.zeros((100, 120), dtype=bool)
    mask[20:80, 30:90] = True
    x0, y0, x1, y1 = central_mask_box(mask)
    assert x0 < 30 and y0 < 20
    assert x1 > 90 and y1 > 80
    assert 0 <= x0 < x1 <= 120
    assert 0 <= y0 < y1 <= 100

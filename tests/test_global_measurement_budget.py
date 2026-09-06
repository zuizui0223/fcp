from __future__ import annotations

import pandas as pd

from fcp_pipeline.global_measurement_budget import (
    select_hashed_taxa,
    select_measurement_rows,
    select_measurement_taxa,
)


def _candidate(n_species: int, target: int) -> pd.DataFrame:
    rows = []
    for taxon in range(1, n_species + 1):
        for p in range(target):
            rows.append({
                "inat_taxon_id": taxon,
                "photo_id": taxon * 100000 + p,
                "species": f"species_{taxon}",
            })
    return pd.DataFrame(rows)


def test_measure_all_when_candidate_species_below_budget():
    taxa = select_measurement_taxa(range(1, 321), maximum_species=500, seed=7)
    assert len(taxa) == 320
    assert set(taxa) == set(range(1, 321))


def test_candidate_budget_caps_at_exactly_1000_and_is_reproducible():
    a = select_hashed_taxa(range(1, 5001), maximum_species=1000, seed=20260916)
    b = select_hashed_taxa(reversed(range(1, 5001)), maximum_species=1000, seed=20260916)
    c = select_hashed_taxa(range(1, 5001), maximum_species=1000, seed=20260917)
    assert len(a) == len(set(a)) == 1000
    assert a == b
    assert a != c


def test_measurement_budget_caps_at_exactly_500_and_is_reproducible():
    a = select_measurement_taxa(range(1, 4001), maximum_species=500, seed=20260918)
    b = select_measurement_taxa(reversed(range(1, 4001)), maximum_species=500, seed=20260918)
    c = select_measurement_taxa(range(1, 4001), maximum_species=500, seed=20260919)
    assert len(a) == len(set(a)) == 500
    assert a == b
    assert a != c


def test_selected_rows_keep_every_photo_for_selected_taxa():
    frame = _candidate(800, 5)
    selected = select_measurement_rows(
        frame,
        target_photos_per_species=5,
        maximum_species=500,
        seed=20260918,
    )
    assert selected["inat_taxon_id"].nunique() == 500
    assert len(selected) == 2500
    assert (selected.groupby("inat_taxon_id").size() == 5).all()
    assert selected["photo_id"].is_unique


def test_input_row_order_does_not_change_measurement_subset():
    frame = _candidate(650, 4)
    left = select_measurement_rows(frame, target_photos_per_species=4, maximum_species=500, seed=44)
    right = select_measurement_rows(
        frame.sample(frac=1.0, random_state=3),
        target_photos_per_species=4,
        maximum_species=500,
        seed=44,
    )
    assert left[["inat_taxon_id", "photo_id"]].equals(right[["inat_taxon_id", "photo_id"]])

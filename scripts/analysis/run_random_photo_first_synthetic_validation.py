#!/usr/bin/env python3
"""Run a non-biological synthetic validation of the quality-safe photo-first H1 machinery."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from fcp_pipeline.photo_first_atlas_v2 import persistence_null_test
from fcp_pipeline.shared_transition_surface import EqualAreaGrid, equal_area_cell_centers


OUTPUT = Path("data/derived/random_photo_first_synthetic_validation_v1.json")


def synthetic_photos() -> pd.DataFrame:
    grid = EqualAreaGrid(n_lon=8, n_sinlat=4)
    cell_id, latitude, longitude = equal_area_cell_centers(grid)
    rng = np.random.default_rng(20260903)
    rows: list[dict[str, object]] = []
    for species_index in range(12):
        species = f"synthetic_species_{species_index:02d}"
        for cid, lat, lon in zip(cell_id, latitude, longitude, strict=True):
            col = int(cid) % grid.n_lon
            base_morph = "red_pink" if col < grid.n_lon // 2 else "blue_purple"
            for photo_index in range(5):
                morph = base_morph
                draw = rng.random()
                if draw < 0.03:
                    morph = "blue_purple" if base_morph == "red_pink" else "red_pink"
                elif draw < 0.08:
                    # Planted spatially unstructured measurement uncertainty must not
                    # become a fifth biological morph or be shuffled by the null.
                    morph = "mixed_uncertain"
                rows.append(
                    {
                        "species": species,
                        "latitude": float(lat + rng.normal(0.0, 0.25)),
                        "longitude": float(lon + rng.normal(0.0, 0.25)),
                        "morph": morph,
                        "synthetic_cell": int(cid),
                        "synthetic_photo": int(photo_index),
                    }
                )
    return pd.DataFrame(rows)


def main() -> None:
    photos = synthetic_photos()
    grid = EqualAreaGrid(n_lon=8, n_sinlat=4)
    observed, null, p_upper = persistence_null_test(
        photos,
        grid=grid,
        target_n=600,
        n_replicates=30,
        species_cap_per_cell=2,
        min_photos_per_cell=8,
        transition_quantile=0.90,
        n_permutations=39,
        sampling_seed=20260903,
        permutation_seed=20260904,
        morph_levels=("red_pink", "blue_purple"),
    )

    supported = observed.edge_table[observed.edge_table["opportunities"] > 0].copy()
    supported = supported.sort_values(
        ["persistence", "opportunities", "edge_id"], ascending=[False, False, True]
    )
    top_edges = supported.head(12).to_dict(orient="records")
    null_median = float(np.median(null))
    if not observed.concentration > null_median:
        raise AssertionError(
            "synthetic planted boundary did not exceed the median species-conditioned null"
        )

    payload = {
        "protocol": "random-photo-first-boundary-persistence-v1",
        "implementation": "photo_first_atlas_v2_quality_safe",
        "role": "software_validation_only_not_biological_evidence",
        "synthetic_design": {
            "species": 12,
            "grid": {"n_lon": 8, "n_sinlat": 4},
            "photos": int(len(photos)),
            "planted_transition": "red_pink in western half versus blue_purple in eastern half",
            "label_noise_fraction": 0.03,
            "measurement_uncertainty_fraction": 0.05,
        },
        "observed": {
            "persistence_concentration": observed.concentration,
            "realized_transition_rate": observed.transition_rate,
            "mean_sampled_photos": observed.mean_sampled_photos,
            "n_replicates": observed.n_replicates,
            "top_edges": top_edges,
        },
        "null": {
            "permutations": int(len(null)),
            "mean": float(np.mean(null)),
            "median": null_median,
            "maximum": float(np.max(null)),
            "p_upper": p_upper,
        },
        "validation_pass": True,
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()

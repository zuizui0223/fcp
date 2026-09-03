#!/usr/bin/env python3
"""Non-biological equivalence and timing check for the cached H1 null engine."""

from __future__ import annotations

import json
from pathlib import Path
import time

import numpy as np
import pandas as pd

from fcp_pipeline.photo_first_atlas_v2 import (
    persistence_null_test as persistence_null_test_slow,
)
from fcp_pipeline.photo_first_h1_fast import persistence_null_test_cached
from fcp_pipeline.shared_transition_surface import EqualAreaGrid, equal_area_cell_centers


OUTPUT = Path("data/derived/random_photo_first_h1_cache_benchmark_v1.json")


def synthetic_photos() -> pd.DataFrame:
    grid = EqualAreaGrid(n_lon=8, n_sinlat=4)
    cell_id, latitude, longitude = equal_area_cell_centers(grid)
    rows = []
    for species_index in range(16):
        species = f"benchmark_species_{species_index:02d}"
        for cid, lat, lon in zip(cell_id, latitude, longitude, strict=True):
            col = int(cid) % grid.n_lon
            for repeat in range(5):
                if repeat == 0 and (int(cid) + species_index) % 7 == 0:
                    morph = "mixed_uncertain"
                elif col < grid.n_lon // 2:
                    morph = "red_pink"
                else:
                    morph = "blue_purple"
                rows.append(
                    {
                        "species": species,
                        "latitude": float(lat),
                        "longitude": float(lon),
                        "morph": morph,
                    }
                )
    return pd.DataFrame(rows)


def main() -> None:
    photos = synthetic_photos()
    grid = EqualAreaGrid(n_lon=8, n_sinlat=4)
    kwargs = dict(
        grid=grid,
        target_n=900,
        n_replicates=24,
        species_cap_per_cell=2,
        min_photos_per_cell=8,
        transition_quantile=0.90,
        n_permutations=9,
        sampling_seed=20260903,
        permutation_seed=20260904,
    )

    started = time.perf_counter()
    slow_observed, slow_null, slow_p = persistence_null_test_slow(photos, **kwargs)
    slow_seconds = time.perf_counter() - started

    started = time.perf_counter()
    fast_observed, fast_null, fast_p = persistence_null_test_cached(photos, **kwargs)
    fast_seconds = time.perf_counter() - started

    np.testing.assert_array_equal(
        slow_observed.edge_table["opportunities"].to_numpy(),
        fast_observed.edge_table["opportunities"].to_numpy(),
    )
    np.testing.assert_array_equal(
        slow_observed.edge_table["transition_count"].to_numpy(),
        fast_observed.edge_table["transition_count"].to_numpy(),
    )
    if not np.isclose(
        slow_observed.concentration,
        fast_observed.concentration,
        rtol=0,
        atol=1e-15,
    ):
        raise AssertionError("cached engine changed the observed H1 statistic")

    # The two null arrays need not be elementwise identical: the cached engine
    # uses an equivalent vectorized random-key permutation within each species.
    # Both preserve the same null constraints and frozen seeds, but the RNG draw
    # consumption differs from the prototype implementation.
    payload = {
        "role": "software_equivalence_and_timing_only_not_biological_evidence",
        "photos": int(len(photos)),
        "replicates": kwargs["n_replicates"],
        "permutations": kwargs["n_permutations"],
        "observed_equivalence": {
            "slow_concentration": slow_observed.concentration,
            "fast_concentration": fast_observed.concentration,
            "same_transition_counts": True,
            "same_opportunities": True,
        },
        "timing_seconds": {
            "slow": slow_seconds,
            "cached": fast_seconds,
            "speedup": slow_seconds / fast_seconds if fast_seconds > 0 else None,
        },
        "null_sanity": {
            "slow_p_upper": slow_p,
            "cached_p_upper": fast_p,
            "slow_null_mean": float(np.mean(slow_null)),
            "cached_null_mean": float(np.mean(fast_null)),
        },
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()

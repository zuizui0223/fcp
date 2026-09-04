#!/usr/bin/env python3
"""Audit capacity for a post-H1/H2 exploratory within-species climate analysis.

No threshold is selected as a favourable outcome. The script reports a complete
small grid of photo-count / occupied-cell / morph-variation eligibility thresholds
for all classifiable fresh measurements. It does not join climate, fit a model, or
make a biological H4 claim.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from fcp_pipeline.photo_first_atlas import prepare_photo_grid
from fcp_pipeline.shared_transition_surface import EqualAreaGrid


ROOT = Path(__file__).resolve().parents[2]
MEASURED = ROOT / "data/derived/random_photo_first_measured_photos_v1.csv"
MEASUREMENT_RESULT = ROOT / "docs/supporting/random_photo_first_measurement_result_v1.json"
H1_RESULT = ROOT / "docs/supporting/random_photo_first_h1_result_v1.json"
OUT_SPECIES = ROOT / "data/derived/random_photo_first_h4_within_species_capacity_v1.csv"
OUT_THRESHOLDS = ROOT / "data/derived/random_photo_first_h4_within_species_capacity_thresholds_v1.csv"
OUT_JSON = ROOT / "docs/supporting/random_photo_first_h4_within_species_capacity_v1.json"

BIOLOGICAL_MORPHS = ("white", "yellow_orange", "red_pink", "blue_purple")
PHOTO_THRESHOLDS = (3, 5, 10, 20)
CELL_THRESHOLDS = (2, 3, 5)
MORPH_THRESHOLDS = (1, 2)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def shannon_bits(values: pd.Series) -> float:
    counts = values.astype(str).value_counts().to_numpy(dtype=float)
    p = counts / counts.sum()
    return float(-np.sum(p * np.log2(p)))


def main() -> int:
    measurement = load_json(MEASUREMENT_RESULT)
    h1 = load_json(H1_RESULT)
    if measurement.get("status") != "complete_fresh_location_blind_measurement_and_join":
        raise RuntimeError("complete fresh measurement required")
    if measurement.get("measurement_table_sha256") != sha256(MEASURED):
        raise RuntimeError("measurement table SHA256 drifted")
    if h1.get("status") != "complete_h1_evaluable":
        raise RuntimeError("completed H1 required so this capacity audit is unambiguously post-H1")

    photos = pd.read_csv(MEASURED)
    required = {"measurement_id", "species", "latitude", "longitude", "morph"}
    missing = sorted(required.difference(photos.columns))
    if missing:
        raise RuntimeError(f"measured table lacks required fields: {missing}")
    if len(photos) != int(measurement["joined_rows"]):
        raise RuntimeError("measurement denominator drifted")

    classifiable = photos.loc[photos["morph"].astype(str).isin(BIOLOGICAL_MORPHS)].copy()
    if len(classifiable) != int(measurement["classified_rows"]):
        raise RuntimeError("classifiable denominator drifted from measurement manifest")
    if len(classifiable) != int(measurement["terminal_status_counts"]["classified_four_state_morph"]):
        raise RuntimeError("classified_rows disagrees with terminal status counts")

    grid = EqualAreaGrid(n_lon=18, n_sinlat=9)
    classifiable = prepare_photo_grid(classifiable, grid=grid)

    rows: list[dict[str, Any]] = []
    for species, group in classifiable.groupby("species", sort=True, observed=True):
        morph_counts = group["morph"].astype(str).value_counts()
        rows.append(
            {
                "species": str(species),
                "classifiable_photos": int(len(group)),
                "occupied_h1_cells": int(group["cell_id"].nunique()),
                "morph_levels": int(morph_counts.size),
                "morph_entropy_bits": shannon_bits(group["morph"]),
                "dominant_morph_fraction": float(morph_counts.iloc[0] / len(group)),
                "white_n": int(morph_counts.get("white", 0)),
                "yellow_orange_n": int(morph_counts.get("yellow_orange", 0)),
                "red_pink_n": int(morph_counts.get("red_pink", 0)),
                "blue_purple_n": int(morph_counts.get("blue_purple", 0)),
            }
        )
    species_table = pd.DataFrame(rows).sort_values(
        ["classifiable_photos", "occupied_h1_cells", "morph_levels", "species"],
        ascending=[False, False, False, True],
        kind="mergesort",
    ).reset_index(drop=True)

    threshold_rows: list[dict[str, Any]] = []
    for min_photos in PHOTO_THRESHOLDS:
        for min_cells in CELL_THRESHOLDS:
            for min_morphs in MORPH_THRESHOLDS:
                eligible = (
                    (species_table["classifiable_photos"] >= min_photos)
                    & (species_table["occupied_h1_cells"] >= min_cells)
                    & (species_table["morph_levels"] >= min_morphs)
                )
                sub = species_table.loc[eligible]
                threshold_rows.append(
                    {
                        "min_classifiable_photos": min_photos,
                        "min_h1_cells": min_cells,
                        "min_morph_levels": min_morphs,
                        "eligible_species": int(len(sub)),
                        "eligible_photos": int(sub["classifiable_photos"].sum()),
                        "median_photos_per_eligible_species": float(sub["classifiable_photos"].median()) if len(sub) else np.nan,
                        "median_cells_per_eligible_species": float(sub["occupied_h1_cells"].median()) if len(sub) else np.nan,
                    }
                )
    thresholds = pd.DataFrame(threshold_rows)

    OUT_SPECIES.parent.mkdir(parents=True, exist_ok=True)
    species_table.to_csv(OUT_SPECIES, index=False, lineterminator="\n")
    thresholds.to_csv(OUT_THRESHOLDS, index=False, lineterminator="\n")

    summary = {
        "protocol": "random-photo-first-h4-within-species-capacity-audit-v1",
        "status": "complete_exploratory_capacity_audit_no_h4_model_fit",
        "introduced_after_h1_h2_outcomes": True,
        "claim_role": "capacity_only_no_biological_h4_claim",
        "measurement": {
            "all_measured_rows": int(len(photos)),
            "classifiable_rows": int(len(classifiable)),
            "classifiable_species": int(species_table["species"].nunique()),
            "occupied_h1_cells": int(classifiable["cell_id"].nunique()),
        },
        "species_repetition": {
            "n_ge_2_photos": int((species_table["classifiable_photos"] >= 2).sum()),
            "n_ge_3_photos": int((species_table["classifiable_photos"] >= 3).sum()),
            "n_ge_5_photos": int((species_table["classifiable_photos"] >= 5).sum()),
            "n_ge_10_photos": int((species_table["classifiable_photos"] >= 10).sum()),
            "n_ge_20_photos": int((species_table["classifiable_photos"] >= 20).sum()),
            "n_ge_2_cells": int((species_table["occupied_h1_cells"] >= 2).sum()),
            "n_ge_3_cells": int((species_table["occupied_h1_cells"] >= 3).sum()),
            "n_ge_5_cells": int((species_table["occupied_h1_cells"] >= 5).sum()),
            "n_ge_2_morphs": int((species_table["morph_levels"] >= 2).sum()),
            "n_ge_3_morphs": int((species_table["morph_levels"] >= 3).sum()),
        },
        "joint_capacity_examples": {},
        "threshold_grid_is_complete_not_outcome_selected": True,
        "lineage": {
            "measurement_table_sha256": sha256(MEASURED),
            "h1_decision_unchanged": h1["decision"],
            "h1_p_upper": float(h1["primary"]["p_upper"]),
        },
        "files": {
            "species_capacity": str(OUT_SPECIES.relative_to(ROOT)),
            "threshold_grid": str(OUT_THRESHOLDS.relative_to(ROOT)),
        },
    }
    for min_photos, min_cells, min_morphs in ((3, 2, 2), (5, 2, 2), (10, 3, 2), (20, 3, 2), (10, 5, 2)):
        row = thresholds.loc[
            (thresholds["min_classifiable_photos"] == min_photos)
            & (thresholds["min_h1_cells"] == min_cells)
            & (thresholds["min_morph_levels"] == min_morphs)
        ].iloc[0]
        summary["joint_capacity_examples"][f"n{min_photos}_cells{min_cells}_morphs{min_morphs}"] = {
            "eligible_species": int(row["eligible_species"]),
            "eligible_photos": int(row["eligible_photos"]),
        }

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

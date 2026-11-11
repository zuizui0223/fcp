"""Replay the frozen 50-km observed spatial colour-depletion estimand.

Computes *observed* species metrics (not 199 permutation worlds). The trusted
2026-10-07 result JSON provides the independent n/mean numerical identity gate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from fcp_pipeline.u100_comparison import MORPHS, opportunity, parse_cohorts, photo_tables

EARTH_RADIUS_KM = 6371.0088
PRIMARY_RADIUS_KM = 50.0
MIN_PAIRS = 30
COARSE_CODE = {x: i for i, x in enumerate(MORPHS)}
EXPECTED_COHORTS = {"discovery": 369, "validation": 363, "prospective": 377}
REFERENCE_KEY = {"discovery": "discovery", "validation": "validation", "prospective": "third"}


def spatial_row(g: pd.DataFrame) -> dict:
    """Match the frozen original geographic unit-vector and upper-triangle logic."""
    g = g.loc[g["usable_colour"]].dropna(subset=["latitude", "longitude"]).copy()
    n = len(g)
    if n < 40:
        raise ValueError("D-evaluable species lost below 40 classifiable geocoded rows")
    lat = pd.to_numeric(g["latitude"], errors="raise").to_numpy(float)
    lon = pd.to_numeric(g["longitude"], errors="raise").to_numpy(float)
    if (not np.isfinite(lat).all() or not np.isfinite(lon).all()
            or (np.abs(lat) > 90).any() or (np.abs(lon) > 180).any()):
        raise ValueError("invalid photo coordinates")
    latr, lonr = np.deg2rad(lat), np.deg2rad(lon)
    xyz = np.column_stack((np.cos(latr) * np.cos(lonr),
                           np.cos(latr) * np.sin(lonr), np.sin(latr)))
    dots = np.clip(xyz @ xyz.T, -1.0, 1.0)
    geo = np.arccos(dots) * EARTH_RADIUS_KM
    u, v = np.triu_indices(n, k=1)
    label = np.array([COARSE_CODE[x] for x in g["morph"]], dtype=np.int8)
    different = label[u] != label[v]
    d_pair = float(np.mean(different))
    near = geo[u, v] <= PRIMARY_RADIUS_KM
    local_pairs = int(np.sum(near))
    depletion = d_pair - float(np.mean(different[near])) if local_pairs >= MIN_PAIRS else np.nan
    return {"inat_taxon_id": int(g["inat_taxon_id"].iloc[0]),
            "cohort": str(g["cohort"].iloc[0]), "D_pair": d_pair,
            "local_pairs_50km": local_pairs, "depletion_50km": depletion}


def observed_spatial_table(photos: pd.DataFrame) -> pd.DataFrame:
    good = photos.loc[photos["usable_colour"]].copy()
    counts = good.groupby(["cohort", "inat_taxon_id"]).size()
    valid = set(counts[counts >= 40].index)
    rows = []
    for key, g in good.groupby(["cohort", "inat_taxon_id"], sort=True):
        if key in valid:
            rows.append(spatial_row(g))
    return pd.DataFrame(rows, columns=["inat_taxon_id", "cohort", "D_pair", "local_pairs_50km", "depletion_50km"])


def frozen_check(rows: pd.DataFrame, ref: dict) -> dict:
    if ref.get("primary_radius_km") != PRIMARY_RADIUS_KM:
        raise ValueError("wrong frozen primary radius")
    if ref.get("confirmatory_decisions_changed") is not False:
        raise ValueError("frozen decision boundary violated")
    observed = {}
    for cohort, n_expected in EXPECTED_COHORTS.items():
        part = rows.loc[rows["cohort"].eq(cohort)]
        if len(part) != n_expected:
            raise ValueError(f"{cohort}: incorrect D-evaluable species count {len(part)}")
        ev = part.loc[part["depletion_50km"].notna()]
        target = ref["primary"][REFERENCE_KEY[cohort]]
        obs = {"n_species": len(ev),
               "mean_local_depletion": float(ev["depletion_50km"].mean()),
               "mean_specieswide_pair_diversity": float(ev["D_pair"].mean()),
               "mean_local_pair_diversity": float((ev["D_pair"] - ev["depletion_50km"]).mean())}
        for field, value in obs.items():
            expected = target[field]
            if not np.isclose(value, expected, rtol=0, atol=1e-11):
                raise ValueError(f"{cohort}: frozen observed spatial {field} mismatch: {value} vs {expected}")
        observed[cohort] = obs
    return observed


def verified_export(frame_path: Path, paths: dict[str, Path], reference: Path, out: Path) -> dict:
    if set(paths) != set(EXPECTED_COHORTS):
        raise ValueError("exact discovery/validation/prospective frozen cohorts required")
    frame = opportunity(frame_path, expected=4730)
    photos, hashes = photo_tables(paths, frame)  # hard-fails if original measured bytes differ
    rows = observed_spatial_table(photos)
    result = frozen_check(rows, json.loads(reference.read_text(encoding="utf-8")))
    out.parent.mkdir(parents=True, exist_ok=True)
    rows.sort_values(["cohort", "inat_taxon_id"], kind="stable").to_csv(out, index=False, float_format="%.17g")
    receipt = {"status": "frozen_observed_50km_spatial_means_exactly_reproduced",
               "cohort_summary": result, "n_species_rows": len(rows),
               "source_photo_sha256": hashes,
               "reference_sha256": hashlib.sha256(reference.read_bytes()).hexdigest(),
               "output_sha256": hashlib.sha256(out.read_bytes()).hexdigest(),
               "permutation_null_recomputed": False,
               "original_confirmatory_decisions_changed": False}
    out.with_suffix(".audit.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    return receipt


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--u100", type=Path, required=True)
    p.add_argument("--measured", action="append", required=True, help="cohort=path, repeat three times")
    p.add_argument("--reference", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args(argv)
    res = verified_export(a.u100, parse_cohorts(a.measured), a.reference, a.out)
    print(json.dumps({"spatial_species_rows": res["n_species_rows"],
                      "observed_50km_reproduced": True, "null_recomputed": False}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

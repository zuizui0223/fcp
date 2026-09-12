#!/usr/bin/env python3
"""Freeze H3a sampling-opportunity covariates without reading colour outcomes.

Only species, observer_id, latitude, and longitude are read. No morph,
global_classifiable, fine_state, or D column is accessed.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "polymorphism_h3a_covariate_preflight_20260912"
OUT.mkdir(parents=True, exist_ok=True)
FILES = {
    "discovery": ROOT / "data" / "derived" / "global_monte_carlo_measured_photos_v1.csv",
    "reserve": ROOT / "data" / "derived" / "rgfca_reserve_replication_measured_photos_v1.csv",
}
EARTH_RADIUS_KM = 6371.0088


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def max_haversine_span_km(lat_deg: np.ndarray, lon_deg: np.ndarray) -> float:
    keep = np.isfinite(lat_deg) & np.isfinite(lon_deg)
    lat = np.deg2rad(lat_deg[keep].astype(float))
    lon = np.deg2rad(lon_deg[keep].astype(float))
    n = len(lat)
    if n < 2:
        return float("nan")
    maxd = 0.0
    # Typical cohort size is small enough that exact pairwise distance is preferred.
    for i in range(n - 1):
        dlat = lat[i + 1:] - lat[i]
        dlon = lon[i + 1:] - lon[i]
        a = np.sin(dlat / 2.0) ** 2 + np.cos(lat[i]) * np.cos(lat[i + 1:]) * np.sin(dlon / 2.0) ** 2
        a = np.clip(a, 0.0, 1.0)
        d = 2.0 * EARTH_RADIUS_KM * np.arcsin(np.sqrt(a))
        if len(d):
            maxd = max(maxd, float(np.nanmax(d)))
    return maxd


def build_panel(cohort: str, path: Path) -> pd.DataFrame:
    cols = ["species", "observer_id", "latitude", "longitude"]
    df = pd.read_csv(path, usecols=cols)
    df["species"] = df["species"].astype(str).str.strip()
    df["latitude"] = pd.to_numeric(df["latitude"], errors="coerce")
    df["longitude"] = pd.to_numeric(df["longitude"], errors="coerce")
    rows = []
    for species, g in df.groupby("species", sort=True):
        obs = g["observer_id"].dropna().astype(str).str.strip()
        obs = obs[obs.ne("")]
        rows.append({
            "cohort": cohort,
            "species": species,
            "n_images_all_measured": int(len(g)),
            "n_observers_all_measured": int(obs.nunique()),
            "n_georeferenced_all_measured": int((g["latitude"].notna() & g["longitude"].notna()).sum()),
            "maximum_span_km_all_measured": max_haversine_span_km(g["latitude"].to_numpy(float), g["longitude"].to_numpy(float)),
        })
    return pd.DataFrame(rows)


def main() -> None:
    panels = [build_panel(cohort, path) for cohort, path in FILES.items()]
    panel = pd.concat(panels, ignore_index=True).sort_values(["cohort", "species"]).reset_index(drop=True)
    panel.to_csv(OUT / "sampling_opportunity_preoutcome.csv", index=False)
    summary = {}
    for cohort in FILES:
        x = panel.loc[panel["cohort"].eq(cohort)].copy()
        summary[cohort] = {
            "species_rows": int(len(x)),
            "n_images_nonmissing": int(x["n_images_all_measured"].notna().sum()),
            "n_observers_nonmissing": int(x["n_observers_all_measured"].notna().sum()),
            "span_nonmissing": int(x["maximum_span_km_all_measured"].notna().sum()),
            "zero_observer_species": int(x["n_observers_all_measured"].eq(0).sum()),
            "fewer_than_two_georeferenced_species": int(x["n_georeferenced_all_measured"].lt(2).sum()),
        }
    result = {
        "analysis": "polymorphism_h3a_covariate_preflight",
        "date_jst": "2026-09-12",
        "outcome_firewall": {
            "columns_read": ["species", "observer_id", "latitude", "longitude"],
            "morph_read": False,
            "fine_state_read": False,
            "global_classifiable_read": False,
            "D_computed": False,
            "association_computed": False,
        },
        "definitions": {
            "n_images": "all measured rows per species; no outcome/classifiability filtering",
            "n_observers": "distinct non-empty observer_id among all measured rows per species",
            "sampled_geographic_span": "exact maximum pairwise haversine great-circle distance (km) among all rows with finite latitude and longitude; Earth radius 6371.0088 km",
        },
        "source_sha256": {cohort: sha256(path) for cohort, path in FILES.items()},
        "panel_sha256": sha256(OUT / "sampling_opportunity_preoutcome.csv"),
        "summary": summary,
    }
    (OUT / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

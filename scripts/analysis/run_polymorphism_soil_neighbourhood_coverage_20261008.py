#!/usr/bin/env python3
"""Outcome-blind SoilGrids 5-km nearest-valid coverage feasibility audit.

Post hoc sensitivity ONLY. Original complete-case FCP results and >=200 species
gate are immutable. For an invalid soil pixel, search fixed +/-8 pixels and
accept the nearest common-valid SoilGrids prediction centre within 10 km
geodesic distance. This does not substitute an observed soil measurement.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import rowcol, xy
from rasterio.warp import transform as rio_transform

from run_polymorphism_full_gradient_partition_20261007 import (
    CLIMATE, LEGACY_COLS, THIRD_COLS, MORPHS, SOIL, SOIL_RAW, SOIL_DEPTHS,
    SHA, attach_climate_elevation, attach_soil, file_sha256, load_measured,
)

MAX_KM = 10.0
WINDOW = 8
MIN_ROWS = 40
FROZEN_GATE = 200
FROZEN_OLD = {"discovery": 172, "validation": 173, "third": 180}
R_EARTH_KM = 6371.0088


def valid_grid(soil_dir: Path):
    """Intersection of all 30 official soil layers and physical AWC condition."""
    common = None
    reference = None
    water_sum = {}
    input_grid = {}
    for prop in SOIL_RAW:
        for depth, weight in SOIL_DEPTHS:
            path = soil_dir / prop / f"{prop}_{depth}_mean_5000.tif"
            with rasterio.open(path) as ds:
                meta = (ds.width, ds.height, str(ds.crs), tuple(ds.transform)[:6])
                if reference is None:
                    reference = (ds.transform, ds.crs, ds.width, ds.height)
                else:
                    affine, crs, width, height = reference
                    if (ds.width != width or ds.height != height or ds.crs != crs
                            or not ds.transform.almost_equals(affine)):
                        raise RuntimeError(f"SoilGrids grid alignment mismatch: {path}")
                a = ds.read(1, masked=False)
                good = np.isfinite(a) & (a >= 0)
                if ds.nodata is not None:
                    good &= a != ds.nodata
                if common is None:
                    common = good.copy()
                else:
                    common &= good
                if prop in {"wv0033", "wv1500"}:
                    # Sum layers separately, only physical / non-nodata values.
                    if prop not in water_sum:
                        water_sum[prop] = np.zeros(a.shape, dtype=np.float32)
                    water_sum[prop] += np.where(good, a, 0).astype(np.float32) * np.float32(weight)
                input_grid[f"{prop}_{depth}"] = {"width":meta[0], "height":meta[1]}
    if common is None or reference is None:
        raise RuntimeError("No SoilGrids rasters found")
    common &= water_sum["wv0033"] >= water_sum["wv1500"]
    return common, reference


def haversine_km(lon: float, lat: float, lons: np.ndarray, lats: np.ndarray) -> np.ndarray:
    p1 = np.deg2rad(lat)
    p2 = np.deg2rad(np.asarray(lats, float))
    dp = p2 - p1
    dl = np.deg2rad(np.asarray(lons, float) - lon)
    a = np.sin(dp / 2)**2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2)**2
    return 2 * R_EARTH_KM * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def get_cell_rows_cols(lon: np.ndarray, lat: np.ndarray, affine, crs):
    rows = np.full(len(lon), -1, dtype=np.int64)
    cols = np.full(len(lon), -1, dtype=np.int64)
    ok = np.isfinite(lon) & np.isfinite(lat)
    ids = np.flatnonzero(ok)
    if not len(ids):
        return rows, cols
    if crs is not None and str(crs).upper() not in {"EPSG:4326", "OGC:CRS84"}:
        xx, yy = rio_transform("EPSG:4326", crs, lon[ids].tolist(), lat[ids].tolist())
    else:
        xx, yy = lon[ids], lat[ids]
    rr, cc = rowcol(affine, xx, yy)
    rows[ids] = np.asarray(rr, dtype=np.int64)
    cols[ids] = np.asarray(cc, dtype=np.int64)
    return rows, cols


def nearest_common_valid(
    lon: np.ndarray,
    lat: np.ndarray,
    common: np.ndarray,
    reference,
    *,
    max_km: float = MAX_KM,
    window: int = WINDOW,
):
    """Fixed raster-neighbourhood search; never accesses biological outcomes."""
    affine, crs, width, height = reference
    rows, cols = get_cell_rows_cols(lon, lat, affine, crs)
    selection_r = np.full(len(lon), -1, dtype=np.int64)
    selection_c = np.full(len(lon), -1, dtype=np.int64)
    dist_km = np.full(len(lon), np.nan, dtype=float)
    statuses = np.full(len(lon), "not_represented", dtype=object)
    substitutions = 0

    for i in range(len(lon)):
        if not np.isfinite(lon[i]) or not np.isfinite(lat[i]):
            statuses[i] = "invalid_coordinate"
            continue
        r, c = int(rows[i]), int(cols[i])
        if r < 0 or r >= height or c < 0 or c >= width:
            statuses[i] = "out_of_grid"
            continue
        if common[r, c]:
            selection_r[i], selection_c[i] = r, c
            dist_km[i] = 0.0  # no spatial substitution
            statuses[i] = "original_valid"
            continue

        r0, r1 = max(0, r - window), min(height, r + window + 1)
        c0, c1 = max(0, c - window), min(width, c + window + 1)
        candidate_local = np.argwhere(common[r0:r1, c0:c1])
        if len(candidate_local) == 0:
            statuses[i] = "no_valid_pixel_in_window"
            continue
        rr = candidate_local[:, 0].astype(np.int64) + r0
        cc = candidate_local[:, 1].astype(np.int64) + c0
        xs, ys = xy(affine, rr, cc, offset="center")
        if crs is not None and str(crs).upper() not in {"EPSG:4326", "OGC:CRS84"}:
            lons, lats = rio_transform(crs, "EPSG:4326", list(xs), list(ys))
        else:
            lons, lats = xs, ys
        km = haversine_km(float(lon[i]), float(lat[i]), np.asarray(lons), np.asarray(lats))
        best = int(np.argmin(km))
        if float(km[best]) > max_km:
            statuses[i] = "no_valid_pixel_within_10km"
            continue
        selection_r[i], selection_c[i] = int(rr[best]), int(cc[best])
        dist_km[i] = float(km[best])
        statuses[i] = "proxy_within_10km"
        substitutions += 1

    return {
        "row": selection_r,
        "col": selection_c,
        "distance_km": dist_km,
        "status": statuses,
        "original_valid": statuses == "original_valid",
        "proxy_valid": (statuses == "original_valid") | (statuses == "proxy_within_10km"),
        "substitutions": substitutions,
    }


def species_metadata(d: pd.DataFrame, old: np.ndarray, new: np.ndarray) -> tuple[dict, pd.DataFrame]:
    orig = pd.Series(old, index=d.index).groupby(d["inat_taxon_id"]).sum()
    prox = pd.Series(new, index=d.index).groupby(d["inat_taxon_id"]).sum()
    species_rows = []
    for taxon, g in d.groupby("inat_taxon_id", sort=True):
        counts = g["morph"].value_counts()
        p = np.array([counts.get(m, 0) / len(g) for m in MORPHS], float)
        old_n = int(orig.loc[taxon])
        new_n = int(prox.loc[taxon])
        species_rows.append({
            "inat_taxon_id": int(taxon),
            "species": str(g["species"].iloc[0]),
            "D": float(1 - np.sum(p*p)),
            "absolute_latitude": float(np.nanmedian(np.abs(pd.to_numeric(g["latitude"], errors="coerce")))),
            "old_complete_rows": old_n,
            "proxy_complete_rows": new_n,
            "old_eligible": bool(old_n >= MIN_ROWS),
            "new_eligible": bool(new_n >= MIN_ROWS),
        })
    out = pd.DataFrame(species_rows).sort_values("inat_taxon_id").reset_index(drop=True)
    def comp(mask):
        x = out.loc[mask]
        return {
            "species": int(len(x)),
            "median_D": float(x["D"].median()) if len(x) else None,
            "median_absolute_latitude": float(x["absolute_latitude"].median()) if len(x) else None,
        }
    old_ok = out["old_eligible"]
    new_ok = out["new_eligible"]
    return {
        "n_old_species": int(old_ok.sum()),
        "n_proxy_species": int(new_ok.sum()),
        "n_recovered_species": int((~old_ok & new_ok).sum()),
        "n_still_excluded_species": int((~new_ok).sum()),
        "representativeness": {
            "old_retained": comp(old_ok),
            "old_excluded": comp(~old_ok),
            "proxy_retained": comp(new_ok),
            "proxy_recovered": comp(~old_ok & new_ok),
            "proxy_excluded": comp(~new_ok),
        },
    }, out


def analyse_cohort(cohort, measured, selection, bio_dir, srad_dir, elevation, soil_dir):
    d = load_measured(measured, cohort)
    d = attach_climate_elevation(d, bio_dir, srad_dir, elevation)
    # Original soil values are recomputed for an exact old-coverage cross-check.
    original = attach_soil(d, soil_dir)
    col = THIRD_COLS if cohort == "third" else LEGACY_COLS
    shared = np.isfinite(d[CLIMATE + ["elevation_m", "latitude", "longitude"] + col].to_numpy(float)).all(axis=1)
    soil_cols = [f"soil_{x}" for x in SOIL]
    old = shared & np.isfinite(original[soil_cols].to_numpy(float)).all(axis=1)
    new = shared & selection["proxy_valid"]
    summary, per_species = species_metadata(d, old, new)
    if summary["n_old_species"] != FROZEN_OLD[cohort]:
        raise RuntimeError(f"{cohort}: frozen old coverage not reproduced: {summary['n_old_species']} vs {FROZEN_OLD[cohort]}")
    if not np.all(~old | new):
        raise RuntimeError(f"{cohort}: proxy coverage lost rows from the original complete case")
    status = pd.Series(selection["status"]).value_counts().to_dict()
    proxy = selection["distance_km"][selection["status"] == "proxy_within_10km"]
    return {
        "cohort": cohort,
        "original_evaluable_species": int(d["inat_taxon_id"].nunique()),
        "original_evaluable_rows": int(len(d)),
        "original_complete_rows": int(old.sum()),
        "proxy_complete_rows": int(new.sum()),
        "original_complete_row_fraction": float(old.mean()),
        "proxy_complete_row_fraction": float(new.mean()),
        **summary,
        "row_status_counts": {str(k): int(v) for k, v in status.items()},
        "proxy_substitution_count": int(selection["substitutions"]),
        "proxy_substitution_fraction_of_source_rows": float(selection["substitutions"] / len(d)),
        "proxy_distance_km": {
            "median": float(np.median(proxy)) if len(proxy) else None,
            "q95": float(np.quantile(proxy, .95)) if len(proxy) else None,
            "maximum": float(np.max(proxy)) if len(proxy) else None,
        },
        "coverage_gate_200_species_pass": bool(summary["n_proxy_species"] >= FROZEN_GATE),
    }, per_species


def main():
    ap = argparse.ArgumentParser()
    for cohort in ("discovery", "validation", "third"):
        ap.add_argument(f"--{cohort}", type=Path, required=True)
    ap.add_argument("--soil-dir", type=Path, required=True)
    ap.add_argument("--bio-dir", type=Path, required=True)
    ap.add_argument("--srad-dir", type=Path, required=True)
    ap.add_argument("--elevation", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    a = ap.parse_args()

    common, reference = valid_grid(a.soil_dir)
    root = a.outdir
    root.mkdir(parents=True, exist_ok=True)
    result = {
        "schema": "fcp_soil_neighbourhood_coverage_posthoc_v1",
        "date_jst": "2026-10-08",
        "role": "post_outcome_coverage_feasibility_only",
        "confirmatory_decisions_changed": False,
        "original_complete_case_hold_unchanged": True,
        "parameters": {
            "maximum_proxy_distance_km": MAX_KM,
            "window_radius_pixels": WINDOW,
            "minimum_complete_rows_per_species": MIN_ROWS,
            "minimum_analysed_species_per_cohort": FROZEN_GATE,
            "original_retained_species": FROZEN_OLD,
        },
        "soil_grid": {
            "crs": str(reference[1]),
            "width": int(reference[2]),
            "height": int(reference[3]),
            "common_valid_fraction_of_pixels": float(np.mean(common)),
            "common_valid_pixel_count": int(common.sum()),
        },
        "source_sha256": SHA,
    }
    all_rows = []
    for cohort in ("discovery", "validation", "third"):
        measured = getattr(a, cohort)
        d = load_measured(measured, cohort)
        lon = pd.to_numeric(d["longitude"], errors="coerce").to_numpy(float)
        lat = pd.to_numeric(d["latitude"], errors="coerce").to_numpy(float)
        selection = nearest_common_valid(lon, lat, common, reference)
        summary, species_rows = analyse_cohort(
            cohort, measured, selection, a.bio_dir, a.srad_dir, a.elevation, a.soil_dir
        )
        result[cohort] = summary
        species_rows.insert(0, "cohort", cohort)
        all_rows.append(species_rows)
    result["all_three_coverage_gates_pass"] = bool(
        all(result[c]["coverage_gate_200_species_pass"] for c in ("discovery", "validation", "third"))
    )
    result["status"] = (
        "COVERAGE_FEASIBLE_FOR_SEPARATE_PROXY_MODEL"
        if result["all_three_coverage_gates_pass"]
        else "HOLD_SOIL_NEIGHBOURHOOD_COVERAGE"
    )
    result["hard_nonclaims"] = [
        "does not estimate any colour-environment association",
        "a nearby modelled prediction cell is not soil measured at the flower",
        "cannot change the original complete-case HOLD",
        "no evidence for local adaptation, drift, gene flow or selection",
    ]
    pd.concat(all_rows, ignore_index=True).to_csv(root / "species_proxy_coverage.csv", index=False)
    (root / "result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

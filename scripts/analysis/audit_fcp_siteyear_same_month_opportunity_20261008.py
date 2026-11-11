#!/usr/bin/env python3
"""Pre-climate-download opportunity audit for repeated FCP photo sites by year.

Entirely retrospective, phenotype-LABEL-blind after the fixed source's
classifiability threshold. No environmental download, phenotypic association,
individual re-identification, or local adaptation estimate occurs here.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

SHA256 = {
    "discovery": "ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4",
    "validation": "0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6",
    "third": "57630fc9f281bce94a0c40a70aaf7bce879dde93d6154175adcd021e8f5c1186",
}
SOURCE_ELIGIBLE_SPECIES = {"discovery": 369, "validation": 363, "third": 377}
VALID_COLOURS = ("white", "yellow_orange", "red_pink", "blue_purple")
RADII_KM = (10.0, 25.0, 50.0)
PRIMARY_RADIUS_KM = 10.0
YEAR_MIN, YEAR_MAX = 1990, 2026
MIN_SOURCE_COLOUR_PHOTOS = 40
MIN_PHOTOS_PER_SITE_YEAR_MONTH = 2
MIN_OBSERVERS_PER_SITE_YEAR_MONTH = 2
MIN_YEARS_PER_SITE_MONTH = 2
FEASIBILITY_SPECIES_PER_COHORT = 30
EARTH_KM = 6371.0088


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def boolean_labels(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False)
    return s.fillna("").astype(str).str.strip().str.lower().isin(("true", "1", "yes", "y"))


def load_photo_opportunity(path: Path, cohort: str) -> tuple[pd.DataFrame, dict]:
    if sha256_file(path) != SHA256[cohort]:
        raise RuntimeError(f"{cohort}: frozen source SHA256 mismatch")
    fields = (
        "inat_taxon_id", "species", "photo_id", "observer_id",
        "observed_on", "latitude", "longitude", "morph", "global_classifiable",
    )
    d = pd.read_csv(path, low_memory=False)
    if not set(fields).issubset(d.columns):
        raise ValueError(f"Missing photo source fields: {sorted(set(fields)-set(d.columns))}")
    d = d.loc[
        boolean_labels(d.global_classifiable)
        & d.morph.astype(str).isin(VALID_COLOURS),
        list(fields),
    ].copy()
    d["latitude"] = pd.to_numeric(d.latitude, errors="coerce")
    d["longitude"] = pd.to_numeric(d.longitude, errors="coerce")
    d = d.loc[d.latitude.between(-90, 90) & d.longitude.between(-180, 180)].copy()
    counts = d.groupby("inat_taxon_id").size()
    d = d.loc[d.inat_taxon_id.isin(counts[counts >= MIN_SOURCE_COLOUR_PHOTOS].index)].copy()
    if d.inat_taxon_id.nunique() != SOURCE_ELIGIBLE_SPECIES[cohort]:
        raise RuntimeError(f"{cohort}: unexpectedly changed fixed >=40-photo species denominator")
    if d.photo_id.duplicated().any():
        raise RuntimeError(f"{cohort}: duplicated photo identities")
    dates = pd.to_datetime(d.observed_on, errors="coerce")
    good = dates.notna() & dates.dt.year.between(YEAR_MIN, YEAR_MAX)
    d["year"] = dates.dt.year.where(good)
    d["month"] = dates.dt.month.where(good)
    d["observer"] = d.observer_id.fillna("").astype(str).str.strip()
    d.loc[d.observer.isin(("nan", "None", "<NA>")), "observer"] = ""
    n_classified = len(d)
    d = d.loc[good].copy()
    d["year"] = d["year"].astype(int)
    d["month"] = d["month"].astype(int)
    d = d.sort_values(["inat_taxon_id", "photo_id"], kind="stable").reset_index(drop=True)
    # Morph labels are used for the pre-existing classifiability/sample gate ONLY.
    # Remove them before selection, preventing outcome-dependent locality choice.
    d = d.drop(columns=["morph", "global_classifiable", "observed_on", "observer_id"])
    return d, {
        "n_species_source_classifiable_40": SOURCE_ELIGIBLE_SPECIES[cohort],
        "n_rows_source_classified": n_classified,
        "n_rows_valid_date": len(d),
        "n_species_any_valid_date": int(d.inat_taxon_id.nunique()),
    }


def gc_matrix_km(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    lat, lon = np.deg2rad(lat.astype(float)), np.deg2rad(lon.astype(float))
    xyz = np.column_stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])
    return EARTH_KM * np.arccos(np.clip(xyz @ xyz.T, -1, 1))


def best_site_month(g: pd.DataFrame, radius_km: float) -> dict | None:
    """Choose a calendar-month/site anchor using ONLY coordinates, dates and observers.

    Each included photo is <=radius/2 from the source-photo anchor; consequently
    every within-site pair is <=radius by the triangle inequality.
    Month is EXACTLY matched between years; local year×month blocks must each
    have at least two independently observed classifiable photographs.
    """
    if radius_km <= 0:
        raise ValueError("Radius must be positive")
    required = {"latitude", "longitude", "year", "month", "observer", "photo_id"}
    if not required.issubset(g):
        raise ValueError(f"Missing fields: {sorted(required-set(g))}")
    if len(g) < 4:
        return None
    g = g.sort_values("photo_id", kind="stable").reset_index(drop=True)
    lat = g.latitude.to_numpy(float)
    lon = g.longitude.to_numpy(float)
    year = g.year.to_numpy(int)
    month = g.month.to_numpy(int)
    obs = g.observer.fillna("").astype(str).to_numpy()
    d = gc_matrix_km(lat, lon)
    best = None
    rank = None
    for anchor in range(len(g)):
        nearby = d[anchor] <= (radius_km / 2.0 + 1e-8)
        if int(nearby.sum()) < MIN_PHOTOS_PER_SITE_YEAR_MONTH * MIN_YEARS_PER_SITE_MONTH:
            continue
        for mo in sorted(set(month[nearby])):
            same_month = nearby & (month == mo) & (obs != "")
            if int(same_month.sum()) < 4:
                continue
            retained = []
            yrs = []
            for y in sorted(set(year[same_month])):
                ix = np.flatnonzero(same_month & (year == y))
                if (len(ix) >= MIN_PHOTOS_PER_SITE_YEAR_MONTH
                        and len(set(obs[ix])) >= MIN_OBSERVERS_PER_SITE_YEAR_MONTH):
                    retained.extend(ix.tolist())
                    yrs.append(int(y))
            if len(yrs) < MIN_YEARS_PER_SITE_MONTH:
                continue
            ix = np.asarray(retained, int)
            max_diameter = float(d[np.ix_(ix, ix)].max())
            if max_diameter > radius_km + 1e-6:
                raise RuntimeError("Site geometry violated its diameter limit")
            candidate_rank = (len(yrs), len(ix), len(set(obs[ix])), -anchor, -int(mo))
            if rank is None or candidate_rank > rank:
                rank = candidate_rank
                best = {
                    "eligible": True,
                    "radius_km": radius_km,
                    "month": int(mo),
                    "n_distinct_years": len(yrs),
                    "years": yrs,
                    "n_photos_qualified_site_month_years": len(ix),
                    "n_observers_qualified_site_month_years": len(set(obs[ix])),
                    "maximum_within_site_pair_distance_km": max_diameter,
                    "anchor_photo_id": str(g.photo_id.iloc[anchor]),
                    "anchor_latitude": float(lat[anchor]),
                    "anchor_longitude": float(lon[anchor]),
                }
    return best


def summarize_cohort(df: pd.DataFrame, cohort: str) -> tuple[dict, list[dict]]:
    species_rows = []
    for taxon, g in df.groupby("inat_taxon_id", sort=True):
        meta = {"cohort": cohort, "inat_taxon_id": str(taxon),
                "species": str(g.species.iloc[0]), "n_dated_classified_photos": int(len(g))}
        for radius in RADII_KM:
            candidate = best_site_month(g, radius)
            key = "r" + str(int(radius))
            meta[key + "_eligible"] = candidate is not None
            meta[key + "_best_month"] = candidate["month"] if candidate else None
            meta[key + "_qualified_years"] = candidate["n_distinct_years"] if candidate else 0
            meta[key + "_qualified_photos"] = candidate["n_photos_qualified_site_month_years"] if candidate else 0
            meta[key + "_qualified_observers"] = candidate["n_observers_qualified_site_month_years"] if candidate else 0
            meta[key + "_anchor_photo_id"] = candidate["anchor_photo_id"] if candidate else None
            meta[key + "_max_diameter_km"] = candidate["maximum_within_site_pair_distance_km"] if candidate else None
        species_rows.append(meta)
    total = SOURCE_ELIGIBLE_SPECIES[cohort]
    summary = {
        "source_species": total,
        "n_species_with_valid_dates": len(species_rows),
        "radii": {},
    }
    for radius in RADII_KM:
        key = "r" + str(int(radius))
        n = sum(bool(s[key + "_eligible"]) for s in species_rows)
        summary["radii"][key] = {
            "radius_km": radius,
            "n_species_eligible": n,
            "fraction_of_original_highdepth_species": n / total,
            "n_species_not_eligible_or_unobserved": total - n,
            "minimum_30_species_feasibility_gate_pass": n >= FEASIBILITY_SPECIES_PER_COHORT,
            "interpretation": "photo calendar/site opportunity only; no flower-colour environment association",
        }
    return summary, species_rows


def audit(paths: dict[str, Path]) -> tuple[dict, pd.DataFrame]:
    report = {
        "schema": "fcp_siteyear_same_month_photo_opportunity_v1",
        "date_jst": "2026-10-08",
        "status": "coverage_only_no_climate_outcomes",
        "role": "posthoc_same_source_observation_feasibility_gate",
        "source_sha256": SHA256,
        "confirmatory_decisions_changed": False,
        "primary_radius_km": PRIMARY_RADIUS_KM,
        "sensitivity_radii_km": [25, 50],
        "minimum_source_classifiable_photos": MIN_SOURCE_COLOUR_PHOTOS,
        "minimum_year_month_photos": MIN_PHOTOS_PER_SITE_YEAR_MONTH,
        "minimum_year_month_distinct_observers": MIN_OBSERVERS_PER_SITE_YEAR_MONTH,
        "minimum_matched_calendar_years": MIN_YEARS_PER_SITE_MONTH,
        "maximum_anchor_to_site_photo_km": "one-half of nominal site diameter",
        "minimum_species_per_cohort_feasibility_gate": FEASIBILITY_SPECIES_PER_COHORT,
        "cohorts": {},
        "hard_nonclaims": [
            "classifiability selection itself may depend on flower phenotype and imaging",
            "dates are dates of photographs, not within-plant phenotype change",
            "each selected site is photographic proximity, not a genetically connected population",
            "observer-disjoint photographs within each year are not confirmed distinct plants",
            "enlarging radius after a 10-km HOLD cannot rescue the primary estimand",
            "no year-specific climate was attached or downloaded",
            "no causal temporal selection, phenotypic plasticity, genotype or fitness inferred",
        ],
    }
    combined = []
    taxon_sets = []
    photo_sets = []
    for name in ("discovery", "validation", "third"):
        d, meta = load_photo_opportunity(paths[name], name)
        taxon_sets.append(set(d.inat_taxon_id.astype(str)))
        photo_sets.append(set(d.photo_id.astype(str)))
        result, rows = summarize_cohort(d, name)
        report["cohorts"][name] = {**meta, **result}
        combined.extend(rows)
    for i in range(3):
        for j in range(i):
            if taxon_sets[i] & taxon_sets[j]:
                raise RuntimeError("Historical taxon cohorts are not independent")
            if photo_sets[i] & photo_sets[j]:
                raise RuntimeError("Historical photo IDs overlap cohorts")
    report["primary_10km_all_cohort_feasibility_pass"] = all(
        report["cohorts"][c]["radii"]["r10"]["minimum_30_species_feasibility_gate_pass"]
        for c in ("discovery", "validation", "third")
    )
    report["next_action"] = (
        "qualify_independent_anomaly_input_and_hold_photo_exposure_checks"
        if report["primary_10km_all_cohort_feasibility_pass"] else
        "HOLD_10KM_REPEATED_SITE_YEAR_MONTH_PHOTO_OPPORTUNITY"
    )
    return report, pd.DataFrame(combined)


def main() -> None:
    p = argparse.ArgumentParser()
    for c in ("discovery", "validation", "third"):
        p.add_argument("--" + c, type=Path, required=True)
    p.add_argument("--outdir", type=Path, required=True)
    a = p.parse_args()
    report, table = audit({c: getattr(a, c) for c in ("discovery", "validation", "third")})
    a.outdir.mkdir(parents=True, exist_ok=True)
    table.to_csv(a.outdir / "species_same_month_siteyear_opportunity.csv", index=False)
    (a.outdir / "result.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

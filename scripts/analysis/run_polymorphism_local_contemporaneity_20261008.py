#!/usr/bin/env python3
"""Bounded contemporaneous photo evidence for visible white/nonwhite flowers.

Tests a necessary discriminator of a strictly season-separated, synchronised
species-wide colour switch. It CANNOT classify genetic polymorphism, fitness,
colour chemistry or adaptive maintenance.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

EXPECTED_SHA = {
    "discovery": "ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4",
    "validation": "0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6",
    "third": "57630fc9f281bce94a0c40a70aaf7bce879dde93d6154175adcd021e8f5c1186",
}
STATES = ("white", "yellow_orange", "red_pink", "blue_purple")
NONWHITE = STATES[1:]
DIAMETERS_KM = (10.0, 25.0, 50.0)
TIME_WINDOWS_DAYS = (14, 30)
PRIMARY_DIAMETER_KM = 10
PRIMARY_TIME_DAYS = 14
MIN_CLASSIFIABLE = 40
MIN_GLOBAL_EACH = 5
MIN_STRICT_EACH = 2
EARTH_KM = 6371.0088


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fp:
        for chunk in iter(lambda: fp.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def bool_series(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False).astype(bool)
    return s.fillna("").astype(str).str.strip().str.lower().isin(("true", "yes", "y", "1"))


def load(path: Path, cohort: str) -> pd.DataFrame:
    if sha(path) != EXPECTED_SHA[cohort]:
        raise ValueError(f"{cohort}: frozen input SHA256 mismatch")
    d = pd.read_csv(path, low_memory=False)
    required = {
        "inat_taxon_id", "species", "photo_id", "latitude", "longitude",
        "observer_id", "observed_on", "morph", "global_classifiable",
    }
    if not required.issubset(d.columns):
        raise ValueError(f"{cohort}: missing fields: {sorted(required-set(d.columns))}")
    d = d.loc[
        bool_series(d.global_classifiable) & d.morph.astype(str).isin(STATES),
        sorted(required),
    ].copy()
    for name in ("latitude", "longitude"):
        d[name] = pd.to_numeric(d[name], errors="coerce")
    d = d.loc[d.latitude.between(-90, 90) & d.longitude.between(-180, 180)].copy()
    counts = d.groupby("inat_taxon_id").size()
    d = d.loc[d.inat_taxon_id.isin(counts[counts >= MIN_CLASSIFIABLE].index)].copy()
    if d.photo_id.duplicated().any():
        raise ValueError(f"{cohort}: duplicated classifiable photo IDs")
    parsed = pd.to_datetime(d.observed_on, errors="coerce")
    d["year"] = parsed.dt.year
    d["day"] = parsed.dt.dayofyear
    d.loc[~d.year.between(1990, 2026), ["year", "day"]] = np.nan
    d["observer"] = d.observer_id.fillna("").astype(str).str.strip()
    d.loc[d.observer.isin(("nan", "None", "<NA>")), "observer"] = ""
    d["cohort"] = cohort
    return d.sort_values(["inat_taxon_id", "photo_id"], kind="stable").reset_index(drop=True)


def geographic_distance(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    lat = np.deg2rad(np.asarray(lat, dtype=float))
    lon = np.deg2rad(np.asarray(lon, dtype=float))
    c = np.cos(lat)
    xyz = np.column_stack((c*np.cos(lon), c*np.sin(lon), np.sin(lat)))
    return np.arccos(np.clip(xyz @ xyz.T, -1.0, 1.0))*EARTH_KM


def two_independent_observers(photos: np.ndarray, observers: np.ndarray) -> bool:
    if int(photos.sum()) < MIN_STRICT_EACH:
        return False
    values = observers[photos]
    values = values[values != ""]
    return int(len(np.unique(values))) >= MIN_STRICT_EACH


def pair_has_independent_observers(
    white_mask: np.ndarray, other_mask: np.ndarray, observers: np.ndarray
) -> bool:
    w = np.unique(observers[white_mask])
    q = np.unique(observers[other_mask])
    w, q = w[w != ""], q[q != ""]
    return bool(len(w) and len(q) and (len(w) > 1 or len(q) > 1 or w[0] != q[0]))


def one_species_evidence(g: pd.DataFrame) -> dict:
    """Anchor-based strict phenological concurrence; each anchor stays one species."""
    g = g.reset_index(drop=True)
    colours = g.morph.astype(str).to_numpy()
    y = g.year.to_numpy(dtype=float)
    day = g.day.to_numpy(dtype=float)
    obs = g.observer.astype(str).to_numpy()
    dist = geographic_distance(g.latitude.to_numpy(float), g.longitude.to_numpy(float))
    w_all = colours == "white"
    nwhite, ncolour = int(w_all.sum()), int((~w_all).sum())
    nonwhite_frequent = int(sum((colours == state).sum() >= MIN_GLOBAL_EACH
                                for state in NONWHITE))
    wc_eligible = nwhite >= MIN_GLOBAL_EACH and ncolour >= MIN_GLOBAL_EACH
    hue_eligible = nonwhite_frequent >= 2
    base = {
        "cohort": str(g.cohort.iloc[0]),
        "inat_taxon_id": int(g.inat_taxon_id.iloc[0]),
        "species": str(g.species.iloc[0]),
        "n_classifiable": int(len(g)),
        "n_valid_date": int(np.isfinite(y).sum()),
        "n_valid_date_observer": int(np.count_nonzero(np.isfinite(y) & (obs != ""))),
        "eligible_white_nonwhite": bool(wc_eligible),
        "eligible_two_nonwhite_hues": bool(hue_eligible),
    }
    all_stats = {}
    valid_year = np.isfinite(y)
    for diam in DIAMETERS_KM:
        for days in TIME_WINDOWS_DAYS:
            tag = f"d{int(diam)}_t{days}"
            opp = False
            simple_wc = False
            simple_hue = False
            strict_wc = False
            strict_hue = False
            strict_wc_years = set()
            strict_hue_years = set()
            nc_window = 0
            nc_opportunity = 0
            nstrict_wc = 0
            nstrict_hue = 0
            # Taking half the diameter and half the temporal window
            # guarantees <= stated pairwise spans, rather than merely
            # comparing each photo to a focal anchor within full span.
            for a in np.flatnonzero(valid_year & (obs != "")):
                nearby = (
                    (dist[a] <= diam/2)
                    & valid_year
                    & (y == y[a])
                    & (np.abs(day - day[a]) <= days/2)
                    & (obs != "")
                )
                if int(nearby.sum()) < 2:
                    continue
                nc_window += 1
                if int(nearby.sum()) >= 4 and len(np.unique(obs[nearby])) >= 2:
                    opp = True
                    nc_opportunity += 1
                white = nearby & w_all
                chromatic = nearby & ~w_all
                if wc_eligible and pair_has_independent_observers(white, chromatic, obs):
                    simple_wc = True
                if wc_eligible and (two_independent_observers(white, obs)
                                    and two_independent_observers(chromatic, obs)):
                    strict_wc = True
                    nstrict_wc += 1
                    strict_wc_years.add(int(y[a]))
                if hue_eligible:
                    hue_groups = [nearby & (colours == state) for state in NONWHITE]
                    present_hue = [grp for grp in hue_groups if grp.sum() > 0]
                    if any(pair_has_independent_observers(x, z, obs)
                           for i, x in enumerate(present_hue) for z in present_hue[i+1:]):
                        simple_hue = True
                    well = sum(two_independent_observers(grp, obs)
                               for grp in hue_groups)
                    if well >= 2:
                        strict_hue = True
                        nstrict_hue += 1
                        strict_hue_years.add(int(y[a]))
            all_stats[tag] = {
                "has_colour_blind_four_photo_local_year_window": bool(opp),
                "number_four_photo_windows": int(nc_opportunity),
                "number_two_photo_windows": int(nc_window),
                "has_independent_observer_white_colour_pair": bool(simple_wc),
                "has_strict_white_colour_contemporaneity": bool(strict_wc),
                "n_strict_white_colour_windows": int(nstrict_wc),
                "n_years_strict_white_colour": int(len(strict_wc_years)),
                "has_independent_observer_two_nonwhite_pair": bool(simple_hue),
                "has_strict_two_nonwhite_contemporaneity": bool(strict_hue),
                "n_strict_two_nonwhite_windows": int(nstrict_hue),
                "n_years_strict_two_nonwhite": int(len(strict_hue_years)),
            }
    for tag, stats in all_stats.items():
        for name, value in stats.items():
            base[f"{tag}_{name}"] = value
    return base


def read_prior(path: Path) -> pd.DataFrame:
    r = pd.read_csv(path, low_memory=False)
    fields = {"cohort", "inat_taxon_id", "species"}
    fields.update(f"diameter_{int(d)}km_white_nonwhite_multi_year_two_observers"
                  for d in DIAMETERS_KM)
    if not fields.issubset(r.columns):
        raise ValueError("Missing expected prior source columns")
    if r.duplicated(["cohort", "inat_taxon_id"]).any():
        raise ValueError("Duplicate cohort/taxon source in prior local recurrence data")
    return r


def summary_one_cohort(c: pd.DataFrame, prior: pd.DataFrame | None) -> dict:
    result = {
        "n_evaluable_species": len(c),
        "n_white_colour_global_at_least_five_each": int(c.eligible_white_nonwhite.sum()),
        "n_two_nonwhite_hues_global_at_least_five_each": int(c.eligible_two_nonwhite_hues.sum()),
        "n_both_global_state_comparisons_eligible": int(
            (c.eligible_white_nonwhite & c.eligible_two_nonwhite_hues).sum()),
        "scales": {},
    }
    for diam in DIAMETERS_KM:
        for days in TIME_WINDOWS_DAYS:
            tag = f"d{int(diam)}_t{days}"
            opp = c[f"{tag}_has_colour_blind_four_photo_local_year_window"].astype(bool)
            eligible_wc = c.eligible_white_nonwhite.astype(bool)
            eligible_hue = c.eligible_two_nonwhite_hues.astype(bool)
            strict_wc = eligible_wc & c[f"{tag}_has_strict_white_colour_contemporaneity"].astype(bool)
            strict_hue = eligible_hue & c[f"{tag}_has_strict_two_nonwhite_contemporaneity"].astype(bool)
            simple_wc = eligible_wc & c[f"{tag}_has_independent_observer_white_colour_pair"].astype(bool)
            simple_hue = eligible_hue & c[f"{tag}_has_independent_observer_two_nonwhite_pair"].astype(bool)
            both = eligible_wc & eligible_hue & opp
            n_wc_opp = int((eligible_wc & opp).sum())
            n_hue_opp = int((eligible_hue & opp).sum())
            pa = {
                "n_species_wc_globally_eligible_and_four_photo_opportunity": n_wc_opp,
                "n_species_hue_globally_eligible_and_four_photo_opportunity": n_hue_opp,
                "n_independent_observer_wc_pair": int(simple_wc.sum()),
                "n_independent_observer_hue_pair": int(simple_hue.sum()),
                "n_strict_wc_two_photos_two_observers_each": int(strict_wc.sum()),
                "n_strict_hue_two_photos_two_observers_each": int(strict_hue.sum()),
                "n_wc_strict_in_multiple_years_somewhere": int(
                    (strict_wc & (c[f"{tag}_n_years_strict_white_colour"] >= 2)).sum()),
                "n_hue_strict_in_multiple_years_somewhere": int(
                    (strict_hue & (c[f"{tag}_n_years_strict_two_nonwhite"] >= 2)).sum()),
                "n_species_both_comparisons_and_opportunity": int(both.sum()),
                "paired_eligible_both_contemporaneous": int((both & strict_wc & strict_hue).sum()),
                "paired_eligible_only_wc_contemporaneous": int((both & strict_wc & ~strict_hue).sum()),
                "paired_eligible_only_hue_contemporaneous": int((both & ~strict_wc & strict_hue).sum()),
                "paired_eligible_neither_contemporaneous": int((both & ~strict_wc & ~strict_hue).sum()),
                "conditional_wc_given_minimal_photo_opportunity":
                    float(strict_wc.sum()/n_wc_opp) if n_wc_opp else None,
                "comparison_role": "observational opportunity-gated counts; NOT adaptation, population polymorphism, or genotype",
                "example_wc_species_unvalidated": c.loc[
                    strict_wc, ["species", "inat_taxon_id"]
                ].sort_values("species").head(15).to_dict("records"),
            }
            if prior is not None:
                key = f"diameter_{int(diam)}km_white_nonwhite_multi_year_two_observers"
                pflags = prior.set_index("inat_taxon_id")[key].reindex(c.inat_taxon_id)
                if pflags.isna().any():
                    raise ValueError("Prior local-recurrence ledger missing a species")
                source_recurrence = bool_series(pflags.reset_index(drop=True))
                n_source = int(source_recurrence.sum())
                # Same sampled species frame, so opportunity and sign tests
                # have noninterchangeable denominators.
                pa["prior_strict_multiyear_recurrence_species"] = n_source
                pa["prior_multiyear_with_contemporaneous_photographic_opportunity"] = int(
                    (source_recurrence & opp).sum())
                pa["prior_multiyear_with_strict_contemporaneous_wc"] = int(
                    (source_recurrence & strict_wc).sum())
                pa["prior_multiyear_without_contemporaneous_opportunity"] = int(
                    (source_recurrence & ~opp).sum())
                pa["prior_multiyear_with_opportunity_but_no_strict_wc"] = int(
                    (source_recurrence & opp & ~strict_wc).sum())
            result["scales"][tag] = pa
    return result


def main() -> None:
    p = argparse.ArgumentParser()
    for cohort in EXPECTED_SHA:
        p.add_argument("--"+cohort, required=True, type=Path)
    p.add_argument("--prior-species", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    args = p.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    old = read_prior(args.prior_species)
    records = []
    all_results = {
        "schema": "fcp_local_within_year_temporal_overlap_posthoc_v1",
        "date_jst": "2026-10-08",
        "status": "complete_source_verified_descriptive_screen",
        "role": "retrospective_local_phenological_contemporaneity_not_fitness",
        "source_sha256": EXPECTED_SHA,
        "confirmatory_decisions_changed": False,
        "pre_outcome_thresholds": False,
        "primary_scale": {"maximum_pair_distance_km": PRIMARY_DIAMETER_KM,
                          "maximum_pair_date_gap_days": PRIMARY_TIME_DAYS},
        "scales_diameter_km": DIAMETERS_KM,
        "scales_gap_days": TIME_WINDOWS_DAYS,
        "minimums": {
            "global_classifiable": MIN_CLASSIFIABLE,
            "global_per_comparison_arm": MIN_GLOBAL_EACH,
            "strict_near_contemporaneous_per_arm": MIN_STRICT_EACH,
            "distinct_observers_per_arm": MIN_STRICT_EACH,
            "necessary_colour_blind_same_window_photos": 4,
        },
        "hard_nonclaims": [
            "within 10/25/50 km photographs are not a verified natural mating population",
            "two visible states near same dates and places do not identify genetic segregation",
            "seasonal plasticity can overlap among individuals or microhabitats",
            "ageing, flower developmental stage, image exposure and observer/site errors remain",
            "seasonal segregation can also occur between genetically distinct colour morphs",
            "selected spatial maximin photo sample has poor local site-year density",
            "no direction of adaptation, fitness, pollinator selection, or pigment chemistry",
            "all colours are photographic classes, not anthocyanin or UV measurements",
            "results post-date main discovery; not independent-source replication",
        ],
        "cohorts": {},
    }
    for cohort in EXPECTED_SHA:
        d = load(getattr(args, cohort), cohort)
        rows = [one_species_evidence(group) for _, group in d.groupby("inat_taxon_id", sort=True)]
        frame = pd.DataFrame(rows)
        pri = old.loc[old.cohort == cohort].copy()
        all_results["cohorts"][cohort] = summary_one_cohort(frame, pri)
        records.append(frame)
    pd.concat(records, ignore_index=True).to_csv(
        args.outdir / "species_local_contemporaneity.csv", index=False)
    (args.outdir / "result.json").write_text(
        json.dumps(all_results, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
    print(json.dumps(all_results, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()

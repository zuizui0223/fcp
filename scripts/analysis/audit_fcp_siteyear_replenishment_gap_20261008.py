#!/usr/bin/env python3
"""Outcome-label-blind gap audit for supplementation of already-sampled FCP taxa.

Diagnose how many *distinct new observer-photo slots* would be needed to
complete a repeated same-month, two-year, <=10km photo-centered locality.
All candidate years/months/sites must already occur in the original source;
no photo API requests or colour outcomes are used. Lower bounds, not promises.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from audit_fcp_siteyear_same_month_opportunity_20261008 import (
    load_photo_opportunity, gc_matrix_km, best_site_month,
    SOURCE_ELIGIBLE_SPECIES
)

DIAMETER_KM = 10.0
ANCHOR_RADIUS_KM = DIAMETER_KM / 2
MIN_DISTINCT_OBSERVERS_PER_YEAR = 2
MIN_YEARS = 2
COHORTS = ("discovery", "validation", "third")


def min_observer_gap(g: pd.DataFrame) -> dict | None:
    """Minimum missing observer-photo slots across observed years and months.

    A 'potential' site must have at least one photo IN EACH of two years for
    one exact calendar month, and within 5km of a real source photo anchor.
    New photos are hypothetical; availability in iNaturalist is NOT checked.
    """
    required = {"latitude", "longitude", "year", "month", "observer", "photo_id"}
    if not required.issubset(g):
        raise ValueError(f"Missing fields {sorted(required-set(g))}")
    if len(g) < 2:
        return None
    g = g.sort_values("photo_id", kind="stable").reset_index(drop=True)
    lat = g.latitude.to_numpy(float)
    lon = g.longitude.to_numpy(float)
    year = g.year.to_numpy(int)
    month = g.month.to_numpy(int)
    obs = g.observer.fillna("").astype(str).to_numpy()
    dm = gc_matrix_km(lat, lon)
    best_rank = None
    best = None
    for anchor in range(len(g)):
        mask = dm[anchor] <= ANCHOR_RADIUS_KM + 1e-8
        if mask.sum() < MIN_YEARS:
            continue
        for mo in np.unique(month[mask]):
            month_ix = np.flatnonzero(mask & (month == mo))
            yr = sorted(np.unique(year[month_ix]).tolist())
            if len(yr) < MIN_YEARS:
                continue
            counts = {}
            photo_counts = {}
            for y in yr:
                indices = month_ix[year[month_ix] == y]
                photo_counts[int(y)] = len(indices)
                counts[int(y)] = len(set(obs[indices]) - {""})
            # Smallest possible number of *new distinct-observer photographs*
            # among two previously observed years, not a search for new dates.
            for ai,ya in enumerate(yr):
                for yb in yr[ai+1:]:
                    have_a, have_b = counts[ya], counts[yb]
                    gap_a = max(0, MIN_DISTINCT_OBSERVERS_PER_YEAR-have_a)
                    gap_b = max(0, MIN_DISTINCT_OBSERVERS_PER_YEAR-have_b)
                    gap = gap_a+gap_b
                    rank = (gap, -(min(have_a,2)+min(have_b,2)),
                            -(photo_counts[ya]+photo_counts[yb]), int(mo),
                            str(g.photo_id.iloc[anchor]), int(ya), int(yb))
                    if best_rank is None or rank < best_rank:
                        best_rank = rank
                        best = {
                            "gap_min_additional_distinct_observer_photos":int(gap),
                            "month":int(mo), "years":[int(ya),int(yb)],
                            "existing_distinct_observers_per_year":[int(have_a),int(have_b)],
                            "existing_photos_per_year":[int(photo_counts[ya]),int(photo_counts[yb])],
                            "anchor_photo_id":str(g.photo_id.iloc[anchor]),
                            "sampled_site_half_radius_km":ANCHOR_RADIUS_KM
                        }
    return best


def audit(paths: dict[str,Path]) -> tuple[dict,pd.DataFrame]:
    report = {
        "schema":"fcp_same_siteyear_same_month_replenishment_gap_v1",
        "role":"post_outcome_label_blind_observation_metadata_sufficiency_only",
        "source_group":"historical_highdepth_three_cohorts_only",
        "status":"NO_NEW_PHOTOS_NO_ENVIRONMENTAL_RESPONSE",
        "date_jst":"2026-10-08",
        "original_10km_eligibility_hold_unchanged":True,
        "confirmatory_decisions_changed":False,
        "original_minimal_four_photo_yearmonth_condition":"two photographs from distinct identified observers in each of two different years, same month and anchored <=10km diameter",
        "gap_definition":"min across pre-existing year pairs and source-photo-anchored 5km half-radius localities of sum(max(0,2-n_distinct_observers))",
        "cohorts":{},
        "hard_nonclaims":[
            "a zero gap is not a biological genetic polymorphism estimate",
            "year/month presence in the selected source does not imply unused photos are available",
            "a missing observer-photo slot cannot be manufactured or relabeled",
            "a photo anchor only selects a geographic neighbourhood, not identical individual plants",
            "sample is conditioned on historical photographic colour classifiability and observational effort",
            "this is a conservative anchored geometry sensitivity and is NOT exact unconstrained 10km diameter feasibility",
            "old and new photo classes cannot serve as a new independent species-disjoint confirmation",
            "no year-specific climate or flower-colour mechanism analysed",
        ]
    }
    records=[]
    for cohort in COHORTS:
        d,info=load_photo_opportunity(paths[cohort],cohort)
        expected_anchor=0
        for tid,g in d.groupby("inat_taxon_id",sort=True):
            candidate=min_observer_gap(g)
            original=best_site_month(g,DIAMETER_KM)
            if (candidate is not None and candidate["gap_min_additional_distinct_observer_photos"]==0)!=(original is not None):
                raise RuntimeError(f"Unexpected original eligibility drift: {cohort}/{tid}")
            expected_anchor+=original is not None
            records.append({
                "cohort":cohort,"inat_taxon_id":str(tid),"species":str(g.species.iloc[0]),
                "n_dated_classifiable_photos":int(len(g)),
                "two_year_month_already_observed":candidate is not None,
                "minimum_additional_observer_photo_slots":(candidate["gap_min_additional_distinct_observer_photos"] if candidate else None),
                "candidate_years":(",".join(map(str,candidate["years"])) if candidate else None),
                "candidate_month":(candidate["month"] if candidate else None),
                "anchor_photo_id":(candidate["anchor_photo_id"] if candidate else None),
                "existing_observer_counts":(",".join(map(str,candidate["existing_distinct_observers_per_year"])) if candidate else None),
                "original_anchor_eligible":original is not None,
            })
        subset=[z for z in records if z["cohort"]==cohort]
        gaps={str(n):sum(z["minimum_additional_observer_photo_slots"]==n for z in subset) for n in range(5)}
        coverage=sum(z["two_year_month_already_observed"] for z in subset)
        n=SOURCE_ELIGIBLE_SPECIES[cohort]
        report["cohorts"][cohort]={
            **info,"source_species":n,
            "n_geographically_year_month_repeatable":coverage,
            "n_no_repeated_year_month_anchored_candidate":n-coverage,
            "minimal_additional_observer_photo_slots_histogram":gaps,
            "n_original_photo_anchor_eligible_gap0":gaps["0"],
            "n_with_nonzero_gap_at_most_1":sum(z["minimum_additional_observer_photo_slots"]==1 for z in subset),
            "n_with_nonzero_gap_at_most_2":sum(z["minimum_additional_observer_photo_slots"] in (1,2) for z in subset),
            "n_with_gap_at_most_2_including_existing":sum(z["minimum_additional_observer_photo_slots"] in (0,1,2) for z in subset),
            "original_10km_30species_gate_pass":bool(gaps["0"]>=30),
        }
        if sum(gaps.values())!=coverage or gaps["0"]!=expected_anchor:
            raise RuntimeError(f"Missing-slot tally inconsistent: {cohort}")
    report["original_all_cohorts_gate_hold"] = not all(
        report["cohorts"][c]["original_10km_30species_gate_pass"] for c in COHORTS
    )
    if not report["original_all_cohorts_gate_hold"]:
        raise RuntimeError("Unexpected recovery changed frozen original 10km HOLD")
    return report,pd.DataFrame(records)


def main():
    ap=argparse.ArgumentParser()
    for c in COHORTS:
        ap.add_argument("--"+c,type=Path,required=True)
    ap.add_argument("--outdir",type=Path,required=True)
    a=ap.parse_args()
    result,rows=audit({c:getattr(a,c) for c in COHORTS})
    a.outdir.mkdir(parents=True,exist_ok=True)
    rows.to_csv(a.outdir/"per_species_label_blind_siteyear_gap.csv",index=False)
    (a.outdir/"result.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__=="__main__":
    main()

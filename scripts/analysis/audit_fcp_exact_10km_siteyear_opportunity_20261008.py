#!/usr/bin/env python3
"""Exact four-photo minimum witness for FCP 10-km site-year-month feasibility.

Exploratory methodology check AFTER anchor-based coverage HOLD.
Does not change the frozen original gate; tests whether restricting site
centres to source photos and half-diameter circles missed valid four-point
sets with max pairwise distance <=10 km. Photo colours are never consulted.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from audit_fcp_siteyear_same_month_opportunity_20261008 import (
    SOURCE_ELIGIBLE_SPECIES,
    best_site_month,
    gc_matrix_km,
    load_photo_opportunity,
)

DIAMETER_KM = 10.0
MIN_YEAR_MONTH_OBSERVERS = 2


def exact_four_photo_witness(g: pd.DataFrame, diameter_km: float = DIAMETER_KM) -> dict | None:
    """Exact YES/NO for >=2 distinct observers in each of two years, same month.

    For >=2 distinct years, >=2 photos/observers per year and pairwise max
    <=diameter, some four-photo subset meets the same requirements; enumerating
    unordered same-year observer-distinct pairs is therefore exhaustive for
    the existence question. Does NOT exhaustively enumerate larger sites or
    maximize the number of repeat years.
    """
    if diameter_km <= 0:
        raise ValueError("Diameter must be positive")
    required = {"latitude", "longitude", "photo_id", "observer", "year", "month"}
    if not required.issubset(g.columns):
        raise ValueError(f"Missing columns: {sorted(required-set(g.columns))}")
    if len(g) < 4:
        return None
    g = g.sort_values("photo_id", kind="stable").reset_index(drop=True)
    lat = g.latitude.to_numpy(float)
    lon = g.longitude.to_numpy(float)
    years = g.year.to_numpy(int)
    months = g.month.to_numpy(int)
    observers = g.observer.fillna("").astype(str).to_numpy()
    dist = gc_matrix_km(lat, lon)
    adjacent = dist <= diameter_km + 1e-8
    for mo in sorted(set(months)):
        mi = np.flatnonzero((months == mo) & (observers != ""))
        if len(mi) < 4:
            continue
        pairs_by_year = {}
        for year in sorted(set(years[mi])):
            ix = mi[years[mi] == year]
            if len(ix) < 2 or len(set(observers[ix])) < MIN_YEAR_MONTH_OBSERVERS:
                continue
            a,b = np.triu_indices(len(ix), k=1)
            keep = adjacent[ix[a], ix[b]] & (observers[ix[a]] != observers[ix[b]])
            if np.any(keep):
                pairs_by_year[year] = np.column_stack([ix[a[keep]], ix[b[keep]]])
        yy = sorted(pairs_by_year)
        for ya, y in enumerate(yy):
            for z in yy[ya+1:]:
                p = pairs_by_year[y]
                q = pairs_by_year[z]
                for i, j in p:
                    possible = (
                        adjacent[i, q[:,0]] & adjacent[i, q[:,1]] &
                        adjacent[j, q[:,0]] & adjacent[j, q[:,1]]
                    )
                    if np.any(possible):
                        k,l = q[int(np.flatnonzero(possible)[0])]
                        witness = np.asarray([i,j,k,l], dtype=int)
                        maxdist = float(dist[np.ix_(witness,witness)].max())
                        if maxdist > diameter_km + 1e-6:
                            raise RuntimeError("Exact site diameter witness incorrect")
                        return {
                            "month": int(mo),
                            "years": [int(y), int(z)],
                            "photo_ids": [str(g.photo_id.iloc[int(a)]) for a in witness],
                            "n_observers_each_year": [
                                int(len(set(observers[[i,j]]))),
                                int(len(set(observers[[k,l]]))),
                            ],
                            "max_diameter_km": maxdist,
                            "n_photo_witness": 4,
                            "is_minimal_witness_only": True,
                        }
    return None


def audit(paths: dict[str, Path]) -> tuple[dict,pd.DataFrame]:
    report = {
        "schema": "fcp_exact_10km_siteyear_opportunity_sensitivity_v1",
        "date_jst": "2026-10-08",
        "status": "post_outcome_geometry_sensitivity_not_new_confirmatory_gate",
        "original_anchor_hold_unchanged": True,
        "confirmatory_decisions_changed": False,
        "diameter_km": DIAMETER_KM,
        "exact_four_photo_existence_check": True,
        "cohorts": {},
        "hard_nonclaims": [
            "exact witness tests existence of at least one valid 4-photo combination, not overall site abundance",
            "the original precommitted anchor-based 10km >=30 per cohort gate remains HOLD",
            "reanalysis after first HOLD is retrospective and cannot be advertised as independent confirmation",
            "two source observers may have photographed the same plant; individual identity unknown",
            "source classifiability and photo selection are not neutral independent natural-population sampling",
            "no same-year specific climate or phenotype association was evaluated",
        ],
    }
    rows=[]
    for cohort in ("discovery","validation","third"):
        d, info = load_photo_opportunity(paths[cohort], cohort)
        n_anchored=n_exact=0
        for tid, g in d.groupby("inat_taxon_id",sort=True):
            anchor=best_site_month(g, DIAMETER_KM)
            exact=exact_four_photo_witness(g, DIAMETER_KM)
            if anchor is not None and exact is None:
                raise RuntimeError(f"{cohort} taxon {tid}: anchored true but exact false")
            n_anchored += anchor is not None
            n_exact += exact is not None
            rows.append({
                "cohort":cohort, "inat_taxon_id":str(tid),"species":str(g.species.iloc[0]),
                "n_dated_classifiable_source_photos":len(g),
                "anchor_10km":anchor is not None,
                "exact_witness_10km":exact is not None,
                "exact_only_10km":anchor is None and exact is not None,
                "witness_month":exact["month"] if exact else None,
                "witness_years":",".join(map(str,exact["years"])) if exact else None,
                "max_witness_diameter_km":exact["max_diameter_km"] if exact else None,
                "witness_photo_ids":",".join(exact["photo_ids"]) if exact else None,
            })
        report["cohorts"][cohort] = {
            **info,
            "n_original_highdepth_species":SOURCE_ELIGIBLE_SPECIES[cohort],
            "n_anchor10":int(n_anchored),
            "n_exact10":int(n_exact),
            "n_extra_exact_geometry":int(n_exact-n_anchored),
            "original_30_species_anchor_gate_pass":bool(n_anchored>=30),
            "exploratory_30_species_exact_gate_pass":bool(n_exact>=30),
        }
    report["original_all_cohort_gate_hold"] = not all(
        report["cohorts"][c]["original_30_species_anchor_gate_pass"]
        for c in ("discovery","validation","third")
    )
    report["any_exact_geometry_witness"] = any(
        report["cohorts"][c]["n_extra_exact_geometry"]>0
        for c in ("discovery","validation","third")
    )
    return report,pd.DataFrame(rows)


def main():
    ap=argparse.ArgumentParser()
    for c in ("discovery","validation","third"):
        ap.add_argument("--"+c,type=Path,required=True)
    ap.add_argument("--outdir",type=Path,required=True)
    a=ap.parse_args()
    result, rows = audit({c:getattr(a,c) for c in ("discovery","validation","third")})
    a.outdir.mkdir(parents=True,exist_ok=True)
    (a.outdir/"result.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    rows.to_csv(a.outdir/"per_species_exact_10km_witness.csv",index=False)
    print(json.dumps(result,indent=2,sort_keys=True))


if __name__=="__main__":
    main()

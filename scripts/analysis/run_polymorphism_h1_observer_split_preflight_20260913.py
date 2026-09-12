#!/usr/bin/env python3
"""Outcome-blind H1 feasibility audit for observer-disjoint D replication.

This script MUST NOT read morph labels, palette counts, or any D outcome.  It uses
only species identity, observer identity, and the frozen global_classifiable flag
to construct deterministic observer-disjoint halves and report opportunity.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "polymorphism_h1_observer_split_preflight_20260913"
OUT.mkdir(parents=True, exist_ok=True)

SOURCES = {
    "discovery": ROOT / "data" / "derived" / "global_monte_carlo_measured_photos_v1.csv",
    "reserve": ROOT / "data" / "derived" / "rgfca_reserve_replication_measured_photos_v1.csv",
}
USECOLS = ["species", "observer_id", "global_classifiable"]
EXPECTED_N40 = {"discovery": 369, "reserve": 363}
THRESHOLDS = [10, 15, 20]


def bool_series(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False).astype(bool)
    return s.fillna("").astype(str).str.strip().str.lower().isin(["true", "1", "yes"])


def stable_key(species: str, observer: str) -> str:
    return hashlib.sha256(f"{species}\0{observer}".encode("utf-8")).hexdigest()


def split_one_species(species: str, obs_counts: pd.Series) -> list[dict[str, Any]]:
    items = [
        (str(obs), int(n), stable_key(species, str(obs)))
        for obs, n in obs_counts.items()
        if str(obs).strip() and int(n) > 0
    ]
    items.sort(key=lambda x: (-x[1], x[2], x[0]))
    totals = {"A": 0, "B": 0}
    nobs = {"A": 0, "B": 0}
    rows: list[dict[str, Any]] = []
    for observer, n, key in items:
        if totals["A"] < totals["B"]:
            half = "A"
        elif totals["B"] < totals["A"]:
            half = "B"
        elif nobs["A"] < nobs["B"]:
            half = "A"
        elif nobs["B"] < nobs["A"]:
            half = "B"
        else:
            # Stable hash tie-break avoids using morph/color/outcome information.
            half = "A" if int(key[-1], 16) % 2 == 0 else "B"
            if not rows:
                half = "A"
        totals[half] += n
        nobs[half] += 1
        rows.append({
            "species": species,
            "observer_id": observer,
            "half": half,
            "n_classifiable": n,
            "stable_hash": key,
        })
    return rows


def audit(name: str, path: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    # Firewall: only these three columns are ever opened.
    df = pd.read_csv(path, usecols=USECOLS)
    df["species"] = df["species"].astype(str)
    df["observer_id"] = df["observer_id"].fillna("").astype(str)
    keep = bool_series(df["global_classifiable"])
    d = df.loc[keep, ["species", "observer_id"]].copy()

    full_counts = d.groupby("species", sort=True).size().rename("n_full_classifiable")
    eligible_species = full_counts[full_counts >= 40].index.tolist()
    if len(eligible_species) != EXPECTED_N40[name]:
        raise RuntimeError(
            f"{name} n>=40 fingerprint mismatch: {len(eligible_species)} != {EXPECTED_N40[name]}"
        )

    d = d[d["species"].isin(eligible_species)].copy()
    grouped = d.groupby(["species", "observer_id"], sort=True).size().rename("n")

    assignment_rows: list[dict[str, Any]] = []
    for species in eligible_species:
        try:
            obs_counts = grouped.loc[species]
        except KeyError:
            continue
        assignment_rows.extend(split_one_species(str(species), obs_counts))

    assign = pd.DataFrame(assignment_rows)
    if assign.empty:
        raise RuntimeError(f"{name}: no observer assignments produced")
    if assign.duplicated(["species", "observer_id"]).any():
        raise RuntimeError(f"{name}: observer assigned more than once within species")

    pivot = (
        assign.groupby(["species", "half"])["n_classifiable"]
        .sum()
        .unstack(fill_value=0)
        .reindex(eligible_species, fill_value=0)
    )
    for half in ["A", "B"]:
        if half not in pivot.columns:
            pivot[half] = 0
    obs_pivot = (
        assign.groupby(["species", "half"])["observer_id"]
        .nunique()
        .unstack(fill_value=0)
        .reindex(eligible_species, fill_value=0)
    )
    for half in ["A", "B"]:
        if half not in obs_pivot.columns:
            obs_pivot[half] = 0

    species_panel = pd.DataFrame({
        "species": eligible_species,
        "n_full_classifiable": [int(full_counts.loc[s]) for s in eligible_species],
        "n_A": [int(pivot.loc[s, "A"]) for s in eligible_species],
        "n_B": [int(pivot.loc[s, "B"]) for s in eligible_species],
        "observers_A": [int(obs_pivot.loc[s, "A"]) for s in eligible_species],
        "observers_B": [int(obs_pivot.loc[s, "B"]) for s in eligible_species],
    })
    species_panel["n_min_half"] = species_panel[["n_A", "n_B"]].min(axis=1)
    species_panel["observer_disjoint"] = True

    threshold_counts = {
        str(t): int((species_panel["n_min_half"] >= t).sum()) for t in THRESHOLDS
    }
    meta = {
        "cohort": name,
        "source": str(path.relative_to(ROOT)),
        "n_rows_source": int(len(df)),
        "n_classifiable_rows": int(keep.sum()),
        "n_species_full_n_ge_40": int(len(eligible_species)),
        "species_with_at_least_two_classifiable_observers": int(
            ((species_panel["observers_A"] > 0) & (species_panel["observers_B"] > 0)).sum()
        ),
        "eligible_by_min_half_threshold": threshold_counts,
        "min_half_quantiles": {
            "q00": float(species_panel["n_min_half"].quantile(0.00)),
            "q25": float(species_panel["n_min_half"].quantile(0.25)),
            "q50": float(species_panel["n_min_half"].quantile(0.50)),
            "q75": float(species_panel["n_min_half"].quantile(0.75)),
            "q100": float(species_panel["n_min_half"].quantile(1.00)),
        },
    }
    assign.insert(0, "cohort", name)
    species_panel.insert(0, "cohort", name)
    return assign, species_panel, meta


def main() -> None:
    all_assign = []
    all_species = []
    cohorts: dict[str, Any] = {}
    for name, path in SOURCES.items():
        assign, species_panel, meta = audit(name, path)
        all_assign.append(assign)
        all_species.append(species_panel)
        cohorts[name] = meta

    assignment = pd.concat(all_assign, ignore_index=True)
    species_panel = pd.concat(all_species, ignore_index=True)
    assignment.to_csv(OUT / "observer_assignment_preoutcome.csv", index=False)
    species_panel.to_csv(OUT / "species_split_opportunity_preoutcome.csv", index=False)

    result = {
        "analysis": "polymorphism_h1_observer_split_preflight",
        "date_jst": "2026-09-13",
        "purpose": "outcome-blind feasibility only",
        "read_columns": USECOLS,
        "split_algorithm": (
            "within species, observer blocks sorted by descending classifiable count then stable SHA256; "
            "greedily assign entire observers to the lower-count half with deterministic tie-breaking"
        ),
        "full_cohort_gate": "n_classifiable >= 40",
        "candidate_min_half_thresholds": THRESHOLDS,
        "cohorts": cohorts,
        "firewall": {
            "morph_column_opened": False,
            "palette_columns_opened": False,
            "D_computed": False,
            "D_unbiased_computed": False,
            "H1_association_or_reproducibility_computed": False,
        },
    }
    (OUT / "result.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

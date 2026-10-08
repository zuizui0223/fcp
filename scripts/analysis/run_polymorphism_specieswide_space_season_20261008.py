#!/usr/bin/env python3
"""FCP cross-species test: geography after conditioning on observed flower season.

Three species-disjoint high-depth photographic cohorts; retrospective analysis.
Outcome: local 4-state visible-colour depletion at <=50 km, conditional on fixed
species colour composition and calendar-season strata. Not fitness or genetics.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

SOURCE_SHA256 = {
    "discovery": "ee854126eed2cfe23e333abe2c28d14df24895a5c52cbc389c060a5a4d6f91f4",
    "validation": "0e2ed349122739eecfc725fb2d5e313d284cf30752da91f9a0429cff0eeaa5e6",
    "third": "57630fc9f281bce94a0c40a70aaf7bce879dde93d6154175adcd021e8f5c1186",
}
COLOURS = ("white", "yellow_orange", "red_pink", "blue_purple")
RADIUS_KM = 50.0
MIN_CLASSIFIABLE = 40
MIN_LOCAL_PAIRS = 30
MIN_EXCHANGEABLE_ROWS = 10
PERMUTATIONS = 199
N_BOOTSTRAP = 1999
SEED = 2026100819
EARTH_KM = 6371.0088
STRATA = ("unconditional", "quarter", "month", "year_month")
PAIR_POLICIES = ("all", "different_observer")
MIN_COHORT_SPECIES = 80


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fp:
        for block in iter(lambda: fp.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def stable_seed(*items: object) -> int:
    msg = "|".join(map(str, (SEED, *items))).encode("utf-8")
    return int.from_bytes(hashlib.sha256(msg).digest()[:8], "little")


def as_bool(s: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False).astype(bool)
    return s.fillna("").astype(str).str.strip().str.casefold().isin(
        ("true", "1", "yes", "y")
    )


def load_source(path: Path, cohort: str) -> tuple[pd.DataFrame, dict]:
    if file_sha256(path) != SOURCE_SHA256[cohort]:
        raise RuntimeError(f"{cohort}: original measured-photo SHA256 failed")
    d = pd.read_csv(path, low_memory=False)
    names = {
        "inat_taxon_id", "species", "photo_id", "observer_id",
        "observed_on", "latitude", "longitude", "morph",
        "global_classifiable",
    }
    if not names.issubset(d):
        raise ValueError(f"Missing source columns: {sorted(names-set(d))}")
    d = d.loc[
        as_bool(d.global_classifiable) & d.morph.astype(str).isin(COLOURS),
        sorted(names),
    ].copy()
    d["latitude"] = pd.to_numeric(d.latitude, errors="coerce")
    d["longitude"] = pd.to_numeric(d.longitude, errors="coerce")
    d = d.loc[
        d.latitude.between(-90, 90) & d.longitude.between(-180, 180)
    ].copy()
    pre_species = int(d.inat_taxon_id.nunique())
    n = d.groupby("inat_taxon_id").size()
    d = d.loc[d.inat_taxon_id.isin(n[n >= MIN_CLASSIFIABLE].index)].copy()
    if d.photo_id.duplicated().any():
        raise RuntimeError("Duplicated photo_id in one frozen cohort")
    date = pd.to_datetime(d.observed_on, errors="coerce")
    good = date.notna() & date.dt.year.between(1990, 2026)
    d["month"] = date.dt.month.where(good).astype("Int64")
    d["year"] = date.dt.year.where(good).astype("Int64")
    d["quarter"] = ((d["month"] - 1)//3 + 1).astype("Int64")
    d["observer"] = d.observer_id.fillna("").astype(str).str.strip()
    d.loc[d.observer.isin(("nan", "None", "<NA>")), "observer"] = ""
    d = d.sort_values(["inat_taxon_id", "photo_id"], kind="stable").reset_index(drop=True)
    species40 = d.inat_taxon_id.nunique()
    rows40 = len(d)
    return d, {
        "n_pre_40_species": pre_species,
        "n_40_classifiable_species": int(species40),
        "n_40_classifiable_rows": int(rows40),
        "n_rows_with_valid_1990_2026_date": int(good.sum()),
        "n_species_with_at_least_40_valid_dates": int(
            (d.groupby("inat_taxon_id").month.count() >= 40).sum()
        ),
    }


def greatcircle_matrix(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    lat = np.deg2rad(np.asarray(lat, float))
    lon = np.deg2rad(np.asarray(lon, float))
    c = np.cos(lat)
    xyz = np.column_stack([c*np.cos(lon), c*np.sin(lon), np.sin(lat)])
    return np.arccos(np.clip(xyz @ xyz.T, -1, 1))*EARTH_KM


def groups_for(d: pd.DataFrame, label: str) -> list[np.ndarray]:
    if label == "unconditional":
        return [np.arange(len(d))]
    if label == "quarter":
        series = d.quarter.astype(int)
    elif label == "month":
        series = d.month.astype(int)
    elif label == "year_month":
        series = d.year.astype(int)*100+d.month.astype(int)
    else:
        raise ValueError(f"Unknown season partition {label}")
    return [
        np.asarray(indices, dtype=int)
        for indices in series.groupby(series, sort=True).indices.values()
    ]


def effective_swap_rows(labels: np.ndarray, groups: list[np.ndarray]) -> int:
    """Count label positions in strata that actually permit different swaps."""
    return int(sum(
        len(index) for index in groups
        if len(index)>=2 and len(np.unique(labels[index]))>=2
    ))


def permute_strata(
    labels: np.ndarray, groups: list[np.ndarray], rng: np.random.Generator
) -> np.ndarray:
    shuffled = labels.copy()
    for index in groups:
        if len(index)>=2:
            shuffled[index] = rng.permutation(shuffled[index])
    return shuffled


def species_test(
    g: pd.DataFrame, cohort: str, stratum: str, pair_policy: str,
    *, permutations: int = PERMUTATIONS
) -> tuple[dict, np.ndarray] | None:
    g = g.loc[g.month.notna()].copy().reset_index(drop=True)
    if len(g) < MIN_CLASSIFIABLE:
        return None
    labels = pd.Categorical(g.morph, categories=COLOURS).codes.astype(np.int8)
    n = len(g)
    ii, jj = np.triu_indices(n, k=1)
    geo = greatcircle_matrix(g.latitude.to_numpy(float), g.longitude.to_numpy(float))
    mask = geo[ii, jj] <= RADIUS_KM
    if pair_policy == "different_observer":
        obs = g.observer.to_numpy(str)
        mask &= (obs[ii] != "") & (obs[jj] != "") & (obs[ii] != obs[jj])
    elif pair_policy != "all":
        raise ValueError("Bad pair policy")
    if int(mask.sum()) < MIN_LOCAL_PAIRS:
        return None
    ua, va = ii[mask], jj[mask]
    groups = groups_for(g, stratum)
    swappable = effective_swap_rows(labels, groups)
    if stratum != "unconditional" and swappable < MIN_EXCHANGEABLE_ROWS:
        return None
    D_all = float(np.mean(labels[ii] != labels[jj]))
    D_local = float(np.mean(labels[ua] != labels[va]))
    observed_depletion = D_all - D_local
    rng = np.random.default_rng(
        stable_seed("species", cohort, int(g.inat_taxon_id.iloc[0]), stratum, pair_policy)
    )
    conditional_local = np.empty(permutations, float)
    for b in range(permutations):
        perm = permute_strata(labels, groups, rng)
        conditional_local[b] = float(np.mean(perm[ua] != perm[va]))
    expected_depletion = D_all-conditional_local
    if not np.isfinite(expected_depletion).all():
        raise ValueError("Non-finite seasonal null")
    if stratum != "unconditional" and np.ptp(expected_depletion) < 1e-12:
        return None
    row = {
        "cohort": cohort,
        "species": str(g.species.iloc[0]),
        "genus": str(g.species.iloc[0]).split()[0],
        "inat_taxon_id": int(g.inat_taxon_id.iloc[0]),
        "stratum": stratum,
        "pair_policy": pair_policy,
        "n_dated_photos": n,
        "n_geographic_local_pairs": int(mask.sum()),
        "n_strata": len(groups),
        "n_effectively_swappable_rows": swappable,
        "D_pair_all": D_all,
        "D_local_pairs": D_local,
        "observed_local_depletion": observed_depletion,
        "null_local_depletion_mean": float(expected_depletion.mean()),
        "excess_over_stratified_null": float(
            observed_depletion - expected_depletion.mean()
        ),
    }
    return row, expected_depletion


def summarize(rows: list[dict], nulls: list[np.ndarray],
              cohort: str, stratum: str, policy: str) -> dict:
    if not rows:
        return {"estimable": False, "n_species": 0, "meets_frozen_80_species_coverage_gate": False, "reason": "no_eligible_photo_coverage", "status": "not_estimable_no_exchangeable_photo_coverage"}
    df = pd.DataFrame(rows)
    mat = np.vstack(nulls)
    excess = df.excess_over_stratified_null.to_numpy(float)
    observed = float(df.observed_local_depletion.mean())
    mean_null = mat.mean(axis=0)
    center = float(mean_null.mean())
    rng = np.random.default_rng(stable_seed("bootstrap", cohort, stratum, policy))
    boot = excess[rng.integers(0, len(excess),
                               size=(N_BOOTSTRAP,len(excess)))].mean(axis=1)
    genus = df.groupby("genus",sort=True).excess_over_stratified_null.mean()
    loo_genus = []
    if len(genus) >= 2:
        for grp in df.genus.unique():
            retained = df.loc[df.genus != grp, "excess_over_stratified_null"]
            if len(retained):
                loo_genus.append(float(retained.mean()))
    return {
        "estimable": True,
        "n_species": int(len(df)),
        "n_genera": int(len(genus)),
        "mean_observed_local_depletion": observed,
        "mean_season_stratified_null_depletion": center,
        "mean_excess_over_season_stratified_null": float(excess.mean()),
        "species_bootstrap_95CI_excess": [
            float(x) for x in np.quantile(boot, [0.025, 0.975])
        ],
        "fraction_species_excess_positive": float(np.mean(excess > 0)),
        "genus_balanced_mean_excess": float(genus.mean()),
        "leave_one_genus_out_min_excess": (
            min(loo_genus) if loo_genus else None
        ),
        "leave_one_genus_out_max_excess": (
            max(loo_genus) if loo_genus else None
        ),
        "permutation_p_upper": float(
            (1+np.count_nonzero(mean_null >= observed)) /
            (len(mean_null)+1)
        ),
        "n_permutations": len(mean_null),
        "meets_frozen_80_species_coverage_gate": bool(len(df) >= MIN_COHORT_SPECIES),
        "null_mode": "shuffle visible photo morphs within fixed species-specific calendar time strata",
        "status": (
            "exploratory_coverage_qualified"
            if len(df) >= MIN_COHORT_SPECIES
            else "coverage_limited_diagnostic_only"
        ),
    }


def cohort_run(d: pd.DataFrame, cohort: str) -> tuple[dict, pd.DataFrame]:
    out = {"scenarios": {}, "n_source_species_40_classifiable": int(d.inat_taxon_id.nunique())}
    metrics = []
    for mode in STRATA:
        for policy in PAIR_POLICIES:
            label = f"{mode}__{policy}"
            rows, nulls = [], []
            for _, g in d.groupby("inat_taxon_id", sort=True):
                t = species_test(g, cohort, mode, policy)
                if t is None:
                    continue
                row, null = t
                rows.append(row)
                nulls.append(null)
            out["scenarios"][label] = summarize(rows,nulls,cohort,mode,policy)
            metrics.extend(rows)
    out["monthly_spatial_excess_qualified_both_pair_policies"] = bool(
        all(
            out["scenarios"][f"month__{policy}"]["meets_frozen_80_species_coverage_gate"]
            and out["scenarios"][f"month__{policy}"]["mean_excess_over_season_stratified_null"]>0
            and out["scenarios"][f"month__{policy}"]["permutation_p_upper"]<0.05
            for policy in PAIR_POLICIES
        )
    )
    return out, pd.DataFrame(metrics)


def main() -> None:
    parser = argparse.ArgumentParser()
    for cohort in SOURCE_SHA256:
        parser.add_argument(f"--{cohort}", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    result = {
        "schema":"fcp_generalizable_space_vs_observed_season_posthoc_v1",
        "date_jst":"2026-10-08",
        "status":"complete",
        "role":"posthoc_comparative_species_balanced_diagnostic",
        "source_sha256":SOURCE_SHA256,
        "confirmatory_decisions_changed":False,
        "cohort_species_are_disjoint":True,
        "same_photo_provider_and_colour_classifier":True,
        "number_of_permutations":PERMUTATIONS,
        "minimum_n_species_per_cohort":MIN_COHORT_SPECIES,
        "primary_spatial_radius_km":RADIUS_KM,
        "declared_primary_strata":"calendar month, repeated year observations pooled",
        "month_year_is_conservative_coverage_diagnostic":True,
        "hard_nonclaims":[
            "not a census of world flower-colour polymorphism prevalence",
            "only species with >=40 classifiable dated photos and >=30 local photo pairs are tested",
            "within-calendar-month label permutation does not eliminate all phenology confounding",
            "calendar-month and month-year temporal groups do not measure flowering fitness or selection",
            "one species may have local groups at very different locations and different habitats",
            "white photographic label has known illumination/clipping coupling",
            "genus and leave-genus diagnostics do not replace complete phylogenetic comparative analysis",
            "repeated observations within species and photo provider do not establish an independent biological replicate of global flora",
            "the residual may reflect neutral gene flow, historical migration, habitat, observer or camera effects",
            "no genotype, pigment chemistry, pollinator agent, plant fitness, adaptive maintenance or causal isolation-by-environment inferred",
        ],
        "cohorts":{},
        "source_coverage":{},
    }
    tables=[]
    for cohort in SOURCE_SHA256:
        frame, coverage=load_source(getattr(args,cohort),cohort)
        output, rows=cohort_run(frame,cohort)
        result["cohorts"][cohort]=output
        result["source_coverage"][cohort]=coverage
        tables.append(rows)
    result["cross_cohort_monthly_qualified_both_policies"] = bool(
        all(
            result["cohorts"][c]["monthly_spatial_excess_qualified_both_pair_policies"]
            for c in SOURCE_SHA256
        )
    )
    pd.concat(tables,ignore_index=True).to_csv(
        args.outdir/"species_space_season_metrics.csv",index=False
    )
    (args.outdir/"result.json").write_text(
        json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8"
    )
    print(json.dumps(result,indent=2,ensure_ascii=False))


if __name__=="__main__":
    main()

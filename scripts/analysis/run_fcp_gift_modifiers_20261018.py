#!/usr/bin/env python3
"""Exploratory GIFT trait contrasts in fitted within-species climate slopes.

Equal-species means of unshrunk Firth slopes are contrasted. Genus-cluster
bootstrap intervals resample both trait groups jointly, without refitting photo
models; they are not causal effects, phylogenetic correction, or molecular tests.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

CONTRASTS = {
    "self_fertilization": ("present", "absent"),
    "lifecycle": ("annual", "perennial"),
    "woodiness": ("woody", "non-woody"),
}
VARIABLES = ("BIO5", "BIO14", "SRADmean")
FAMILY_SIZE = 9


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def admitted_levels(series: pd.Series, first: str, second: str) -> pd.Series:
    """Unknown, variable, biennial or combined labels never become a reference."""
    return series.where(series.isin([first, second]))


def differences_from_multiplicities(counts, n_first, n_second, sums_first, sums_second):
    n1, n0 = counts @ n_first, counts @ n_second
    valid = (n1 > 0) & (n0 > 0)
    values = (counts[valid] @ sums_first) / n1[valid, None] - (counts[valid] @ sums_second) / n0[valid, None]
    return values, valid


def summarize_contrast(wide: pd.DataFrame, trait: str, repetitions: int, seed: int) -> list[dict]:
    first, second = CONTRASTS[trait]
    group = admitted_levels(wide[trait], first, second)
    q = wide.loc[group.notna()].copy()
    q["_first"] = q[trait].eq(first)
    q["genus"] = q.species.str.split().str[0]
    one, zero = q.loc[q._first], q.loc[~q._first]
    n1, n0 = len(one), len(zero)
    g1, g0 = one.genus.nunique(), zero.genus.nunique()
    base = {"trait": trait, "first_level": first, "reference_level": second, "n_first": n1, "n_reference": n0, "genera_first": g1, "genera_reference": g0, "n_excluded_missing_or_other_level": len(wide) - len(q), "bootstrap_repetitions": repetitions, "seed": seed}
    if min(n1, n0) < 20 or min(g1, g0) < 5:
        return [{**base, "variable": v, "status": "insufficient_group_coverage"} for v in VARIABLES]
    genera = sorted(q.genus.unique())
    ng = len(genera)
    def aggregates(sub):
        n = sub.groupby("genus").size().reindex(genera, fill_value=0).to_numpy(float)
        sums = sub.groupby("genus")[list(VARIABLES)].sum().reindex(genera, fill_value=0).to_numpy(float)
        return n, sums
    counts1, sums1 = aggregates(one)
    counts0, sums0 = aggregates(zero)
    rng = np.random.default_rng(seed)
    multiplicities = rng.multinomial(ng, np.full(ng, 1 / ng), size=repetitions)
    bootstrap, valid = differences_from_multiplicities(multiplicities, counts1, counts0, sums1, sums0)
    if valid.mean() < 0.95:
        return [{**base, "variable": v, "status": "too_many_empty_group_bootstraps", "invalid_bootstrap_draws": int((~valid).sum())} for v in VARIABLES]
    result = []
    alpha = 0.05 / FAMILY_SIZE
    for j, variable in enumerate(VARIABLES):
        low, high = np.quantile(bootstrap[:, j], [0.025, 0.975])
        family_low, family_high = np.quantile(bootstrap[:, j], [alpha / 2, 1 - alpha / 2])
        mean1, mean0 = float(one[variable].mean()), float(zero[variable].mean())
        result.append({**base, "variable": variable, "status": "estimated", "genus_clusters_resampled": ng, "invalid_bootstrap_draws": int((~valid).sum()), "first_equal_mean": mean1, "reference_equal_mean": mean0, "difference_first_minus_reference": mean1 - mean0, "percentile_ci_low": float(low), "percentile_ci_high": float(high), "family9_percentile_ci_low": float(family_low), "family9_percentile_ci_high": float(family_high), "first_median_within_species_climate_sd_native": float(one[variable + "_native_sd"].median()), "reference_median_within_species_climate_sd_native": float(zero[variable + "_native_sd"].median())})
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, default=Path("results/fcp_gift_modifiers_20261018"))
    parser.add_argument("--repetitions", type=int, default=5000)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    folder = root / "results/fcp_real_447500_environment_model_20261018"
    receipt = json.loads((folder / "execution_receipt.json").read_text())
    for name, expected in receipt["output_sha256"].items():
        if sha(folder / name) != expected:
            raise ValueError(f"Environmental output hash mismatch: {name}")
    join_path = args.outdir / "species_trait_join.csv"
    extraction = json.loads((args.outdir / "profile.json").read_text())
    if sha(join_path) != extraction["output_sha256"][join_path.name]:
        raise ValueError("Trait join differs from frozen extraction")
    traits = pd.read_csv(join_path)
    slopes = pd.read_csv(folder / "species_slopes.csv")
    primary = slopes.loc[slopes.model.eq("primary")].copy()
    assert len(primary) == 1574 * 3 and not primary.duplicated(["inat_taxon_id", "variable"]).any()
    ids = primary[["inat_taxon_id", "species", "source_cohort"]].drop_duplicates()
    wide = ids.merge(primary.pivot(index="inat_taxon_id", columns="variable", values="slope"), on="inat_taxon_id", validate="one_to_one")
    native_sd = primary.pivot(index="inat_taxon_id", columns="variable", values="primary_within_species_sd_native").add_suffix("_native_sd")
    wide = wide.merge(native_sd, on="inat_taxon_id", validate="one_to_one")
    wide = wide.merge(traits[["inat_taxon_id", *CONTRASTS]], on="inat_taxon_id", how="left", validate="one_to_one")
    assert len(wide) == 1574 and np.isfinite(wide[list(VARIABLES)]).all().all()
    rows = []
    for index, trait in enumerate(CONTRASTS):
        rows.extend(summarize_contrast(wide, trait, args.repetitions, 20261018 + index))
    results = pd.DataFrame(rows)
    results.to_csv(args.outdir / "slope_modifier_contrasts.csv", index=False)
    wide.to_csv(args.outdir / "species_slopes_with_traits.csv", index=False)
    profile = {"completed_utc": datetime.now(timezone.utc).isoformat(), "method": "equal-species unshrunk fitted-slope differences with genus-cluster percentile bootstrap", "scope": "exploratory trait-associated standardized-slope contrasts; no causal or molecular mechanism claim", "contrasts": CONTRASTS, "climate_variables": VARIABLES, "family_size": FAMILY_SIZE, "minimum_species_each": 20, "minimum_genera_each": 5, "bootstrap_repetitions": args.repetitions, "family_interval": "Bonferroni nominal 1-0.05/9 percentile interval per contrast, with bootstrap approximation", "uncertainty_scope": "empirical resampling of fitted slopes by genus; no photo refits, no explicit phylogenetic covariance or within-fit variance resampling", "scaling_limit": "slopes are per each species primary climate SD; trait contrasts may reflect sampled climate-range differences as well as biological sensitivity; native SD medians provided", "source_sha256": {str(join_path): sha(join_path), str(folder / "species_slopes.csv"): sha(folder / "species_slopes.csv")}, "code_sha256": sha(Path(__file__)), "output_sha256": {name: sha(args.outdir / name) for name in ["slope_modifier_contrasts.csv", "species_slopes_with_traits.csv"]}}
    (args.outdir / "modifier_profile.json").write_text(json.dumps(profile, indent=2) + "\n")
    print(results.to_csv(index=False), flush=True)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Post-outcome null-conditioning sensitivity audit for FCP geographic ITV.

Read the previously source-verified specieswide space/season result.json.
Compare fixed-species-cohort equal-species mean geographic effects under
unconditional, month-preserving and year-month-preserving photo-label nulls.
This is NOT a decomposition of causal variance or climatic mediation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

COHORTS = ("discovery", "validation", "third")
POLICIES = ("all", "different_observer")
CONDITIONS = ("unconditional", "month", "year_month")


def summarize(source: dict) -> dict:
    if source.get("schema") != "fcp_generalizable_space_vs_observed_season_posthoc_v1":
        raise ValueError("Unexpected source schema; do not silently reinterpret a new analysis")
    if source.get("status") != "complete" or not source.get("cohort_species_are_disjoint"):
        raise ValueError("Source is incomplete or cohorts are not species-disjoint")
    if source.get("number_of_permutations") != 199:
        raise ValueError("Permutation null has changed")
    out = {
        "schema": "fcp_space_time_null_sensitivity_posthoc_v1",
        "role": "descriptive_retrospective_null_sensitivity_not_causal_partition",
        "source_schema": source["schema"],
        "source_sha256": source.get("source_sha256"),
        "confirmatory_decisions_changed": False,
        "comparisons": {},
        "hard_nonclaims": [
            "attenuation is a comparison between conditional randomization estimands, not the fraction caused by time",
            "month and year-month restrictions alter label exchangeability and n conditionally identifiable species",
            "calendar observation date is not year-specific site climate or measured fitness",
            "no genotype, genetic mating population, selection, or local adaptation is inferred",
            "three species-disjoint cohorts use the same source and flower-colour classifier",
        ],
    }
    for cohort in COHORTS:
        s = source["cohorts"][cohort]["scenarios"]
        out["comparisons"][cohort] = {}
        for policy in POLICIES:
            a, b, c = [s[f"{x}__{policy}"] for x in CONDITIONS]
            n = [x["n_species"] for x in (a, b, c)]
            if min(n) <= 0 or len(set(n)) != 1:
                raise ValueError(f"{cohort}/{policy}: conditioning changes analysis denominator")
            ids = [x["n_identifiable_species"] for x in (a, b, c)]
            if not (ids[0] >= ids[1] >= ids[2] >= 0):
                raise ValueError(f"{cohort}/{policy}: inconsistent identifiable species counts")
            obs = [x["mean_observed_local_depletion"] for x in (a, b, c)]
            if max(obs) - min(obs) > 1e-12:
                raise ValueError(f"{cohort}/{policy}: underlying geographic statistic changed")
            effect = [x["mean_excess_over_season_stratified_null"] for x in (a, b, c)]
            if abs(effect[0]) < 1e-12:
                atten_m = atten_ym = None
            else:
                atten_m = (effect[0] - effect[1]) / effect[0]
                atten_ym = (effect[0] - effect[2]) / effect[0]
            out["comparisons"][cohort][policy] = {
                "n_species_same_across_conditions": n[0],
                "n_identifiable_unconditional_month_yearmonth": ids,
                "mean_observed_local_depletion_fixed": obs[0],
                "mean_excess_unconditional_month_yearmonth": effect,
                "relative_attenuation_unconditional_to_month": atten_m,
                "relative_attenuation_unconditional_to_yearmonth": atten_ym,
                "month_permutation_upper_p": b["permutation_p_upper"],
                "yearmonth_permutation_upper_p": c["permutation_p_upper"],
                "month_bootstrap_95CI": b["species_bootstrap_95CI_excess"],
                "yearmonth_bootstrap_95CI": c["species_bootstrap_95CI_excess"],
                "yearmonth_positive_lower_bound": c["species_bootstrap_95CI_excess"][0] > 0,
            }
    return out


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    original_bytes = args.source.read_bytes()
    out = summarize(json.loads(original_bytes))
    out["source_file_sha256"] = hashlib.sha256(original_bytes).hexdigest()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(out["comparisons"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "polymorphism_42111_h1_h2_eligibility_20260912"

CAPACITY = ROOT / "results" / "rgfca_42111_breadth_depth_step8b_20260911" / "result.json"
BREADTH = ROOT / "results" / "rgfca_42111_breadth_measurement_step8f_20260911" / "result.json"
MISSINGNESS = ROOT / "results" / "rgfca_42111_breadth_missingness_audit_step8h_20260911" / "result.json"
DIRECTION = ROOT / "results" / "polymorphism_directionality_step1_20260909" / "result.json"
MIXED = ROOT / "results" / "polymorphism_mixed_step2_20260909" / "result.json"

N_UNIVERSE = 42111
LEGACY_CLASSIFIABLE_MIN = 40
H1_SECOND_PRIMARY = 0.10
H1_SECOND_STRICT = 0.20
U100_RAW_DEPTH = 100


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def binom_ge(n: int, k: int, p: float) -> float:
    return float(sum(math.comb(n, j) * (p ** j) * ((1.0 - p) ** (n - j)) for j in range(k, n + 1)))


def integer_from_fraction(n: int, frac: float, label: str) -> int:
    x = n * frac
    nearest = int(round(x))
    if abs(x - nearest) > 0.02:
        raise RuntimeError(f"{label} does not map cleanly to an integer count: {x}")
    return nearest


def main() -> None:
    cap = read_json(CAPACITY)
    breadth = read_json(BREADTH)
    missing = read_json(MISSINGNESS)
    direction = read_json(DIRECTION)
    mixed = read_json(MIXED)

    if int(cap["species_universe"]) != N_UNIVERSE:
        raise RuntimeError("capacity universe drift")
    if int(breadth["species_denominator"]) != N_UNIVERSE:
        raise RuntimeError("breadth denominator drift")
    if int(missing["species_denominator"]) != N_UNIVERSE:
        raise RuntimeError("missingness denominator drift")

    depths = {int(k): int(v) for k, v in cap["threshold_counts_after_observer_cap"].items()}
    if U100_RAW_DEPTH not in depths:
        raise RuntimeError("100-photo capacity tier missing")
    u100 = depths[U100_RAW_DEPTH]
    if u100 != 4730:
        raise RuntimeError(f"U100 fingerprint drift: {u100}")

    validation_n = int(direction["discovery_species"])
    if validation_n != 369:
        raise RuntimeError(f"validation fingerprint drift: {validation_n}")
    f10 = float(direction["fingerprint"]["second_ge_0_10"])
    f20 = float(direction["fingerprint"]["second_ge_0_20"])
    v_h1_10 = integer_from_fraction(validation_n, f10, "second_ge_0_10")
    v_h1_20 = integer_from_fraction(validation_n, f20, "second_ge_0_20")
    if (v_h1_10, v_h1_20) != (172, 98):
        raise RuntimeError(f"validation H1 fingerprint drift: {(v_h1_10, v_h1_20)}")

    if mixed["decision"]["verdict"] != "MIXED_REMAINS_STRUCTURAL_MISSINGNESS":
        raise RuntimeError("mixed-image gate status changed; refreeze protocol before proceeding")
    if bool(mixed["decision"]["intermediate_colour_interpretation_supported"]):
        raise RuntimeError("mixed rows unexpectedly qualified as biological intermediates")

    if direction["decision"]["verdict"] != "CLAIM2_REJECT_PRIMARY_NOT_SUPPORTED":
        raise RuntimeError("legacy directionality status changed; refreeze protocol before proceeding")
    if bool(direction["decision"]["claim2_survives"]):
        raise RuntimeError("legacy recurrent-direction claim unexpectedly survives")

    p_center = float(breadth["classifiable_fraction"])
    capdiag = missing["family_selection_diagnostics"]["capacity_bin"]
    p_low = float(capdiag["minimum_level_classifiable_fraction"])
    p_high = float(capdiag["maximum_level_classifiable_fraction"])
    if not (0 < p_low <= p_center <= p_high < 1):
        raise RuntimeError("classifiability planning envelope is inconsistent")

    prob_low = binom_ge(U100_RAW_DEPTH, LEGACY_CLASSIFIABLE_MIN, p_low)
    prob_center = binom_ge(U100_RAW_DEPTH, LEGACY_CLASSIFIABLE_MIN, p_center)
    prob_high = binom_ge(U100_RAW_DEPTH, LEGACY_CLASSIFIABLE_MIN, p_high)

    result = {
        "analysis": "polymorphism_42111_h1_h2_eligibility_audit",
        "date_jst": "2026-09-12",
        "status": "frozen_preexecution_eligibility_and_feasibility_ledger",
        "protocol": "docs/POLYMORPHISM_42111_H1_H2_ELIGIBILITY_PROTOCOL_20260912.md",
        "sampling_frame": {
            "U0_species": N_UNIVERSE,
            "U100_definition": "after_observer_cap >= 100 raw eligible photos",
            "U100_species": u100,
            "legacy_validation_species_nclass_ge40": validation_n,
        },
        "H1": {
            "role": "measurement-valid non-trivial secondary-mode gate; not a prevalence estimator",
            "n_classifiable_min": LEGACY_CLASSIFIABLE_MIN,
            "second_fraction_primary_min": H1_SECOND_PRIMARY,
            "second_fraction_strict_min": H1_SECOND_STRICT,
            "mixed_uncertain_promoted": False,
            "legacy_validation_primary_count": v_h1_10,
            "legacy_validation_strict_count": v_h1_20,
            "U100_observed_eligible_n": None,
            "U100_observed_eligible_n_status": "unknown_until_high_depth_measurement",
        },
        "H2": {
            "role": "sign-invariant geometry of species-specific primary-to-secondary displacement vectors",
            "delta_definition": "Delta_i = mu_i,2 - mu_i,1 using normalized nine-colour biological flower palette centroids",
            "primary_statistic": "largest eigenvalue of mean outer-product matrix of unit Delta_i vectors",
            "orientation": "axial_sign_invariant",
            "null": "10000 independent Haar-random orientations in the 8D zero-sum subspace",
            "legacy_M1_M2_test_reused": False,
            "U100_observed_eligible_n": None,
            "U100_observed_eligible_n_status": "unknown_until_H1_and_delta_construction",
        },
        "frozen_negative_constraints": {
            "mixed_bridge_geometry_p": float(mixed["bridge_geometry"]["p_bridge_upper"]),
            "mixed_verdict": mixed["decision"]["verdict"],
            "legacy_directionality_rho_D_R12": float(direction["primary"]["rho_D_R12"]),
            "legacy_directionality_label_p": float(direction["primary"]["label_permutation_upper_p"]),
            "legacy_directionality_rotation_p": float(direction["direction_specificity"]["p_rotation_upper"]),
            "legacy_directionality_verdict": direction["decision"]["verdict"],
        },
        "planning_only_binomial_depth40_projection": {
            "raw_photos": U100_RAW_DEPTH,
            "classifiable_required": LEGACY_CLASSIFIABLE_MIN,
            "classifiability_probability_low": p_low,
            "classifiability_probability_central": p_center,
            "classifiability_probability_high": p_high,
            "probability_reach_nclass_ge40_low": prob_low,
            "probability_reach_nclass_ge40_central": prob_center,
            "probability_reach_nclass_ge40_high": prob_high,
            "projected_species_low": u100 * prob_low,
            "projected_species_central": u100 * prob_center,
            "projected_species_high": u100 * prob_high,
            "biological_eligibility_claim_allowed": False,
            "assumption": "independent-binomial planning approximation only; species heterogeneity ignored",
        },
        "hard_stops": [
            "one-photo breadth classifiability is not species polymorphism",
            "mixed_uncertain is not a fifth biological morph",
            "planning projection is not an observed eligible N",
            "legacy rejected M1/M2 directionality test is not H2",
            "H3 ecology/phylogeny cannot tune H1/H2 eligibility thresholds",
        ],
    }

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    md = [
        "# 42,111-species H1/H2 eligibility audit",
        "",
        f"- U0 sampling frame: **{N_UNIVERSE:,} species**",
        f"- U100 high-depth opportunity ceiling: **{u100:,} species**",
        f"- legacy validation frame: **{validation_n} species** with >=40 classifiable observations",
        f"- validation H1 primary (second state >=10%): **{v_h1_10} species**",
        f"- validation H1 strict (second state >=20%): **{v_h1_20} species**",
        "- U100 biological H1/H2 eligible N: **unknown until measurement**",
        "",
        "## Planning only",
        "",
        f"At 100 raw photos, the independent-binomial probability of reaching >=40 classifiable rows is {prob_center:.3f} at the observed one-anchor classifiability rate, with a capacity-stratum sensitivity envelope {prob_low:.3f}–{prob_high:.3f}.",
        f"Applied to U100 this is about {u100*prob_center:,.0f} species centrally ({u100*prob_low:,.0f}–{u100*prob_high:,.0f}), but these are resource-planning projections, not biological eligibility estimates.",
        "",
        "## Frozen consequence",
        "",
        "The next empirical operation is high-depth measurement of U100 followed by the H1 gate. Only then may species-specific Delta_i vectors be constructed for H2.",
    ]
    (OUT / "RESULT.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

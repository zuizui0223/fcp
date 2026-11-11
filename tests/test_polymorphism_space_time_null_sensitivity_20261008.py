"""Synthetic schema/inference-boundary checks for the FCP time-null audit."""
from __future__ import annotations
import copy
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "analysis"))
from audit_polymorphism_space_time_null_sensitivity_20261008 import summarize


def example_source():
    scenarios = {}
    for policy in ("all", "different_observer"):
        for mode, n_id, effect in (
            ("unconditional", 9, 0.10),
            ("month", 7, 0.08),
            ("year_month", 4, 0.05),
        ):
            scenarios[f"{mode}__{policy}"] = {
                "n_species": 10,
                "n_identifiable_species": n_id,
                "mean_observed_local_depletion": 0.12,
                "mean_excess_over_season_stratified_null": effect,
                "permutation_p_upper": 0.005,
                "species_bootstrap_95CI_excess": [0.01, 0.09],
            }
    return {
        "schema": "fcp_generalizable_space_vs_observed_season_posthoc_v1",
        "status": "complete",
        "cohort_species_are_disjoint": True,
        "number_of_permutations": 199,
        "cohorts": {
            name: {"scenarios": copy.deepcopy(scenarios)}
            for name in ("discovery", "validation", "third")
        },
    }


def test_attenuation_is_descriptive_and_preserves_denominator():
    result = summarize(example_source())
    v = result["comparisons"]["discovery"]["all"]
    assert v["n_species_same_across_conditions"] == 10
    assert v["n_identifiable_unconditional_month_yearmonth"] == [9, 7, 4]
    assert v["relative_attenuation_unconditional_to_month"] == pytest.approx(0.20)
    assert v["relative_attenuation_unconditional_to_yearmonth"] == pytest.approx(0.50)
    assert result["confirmatory_decisions_changed"] is False
    assert "not the fraction caused by time" in result["hard_nonclaims"][0]


def test_rejects_changed_species_denominator():
    source = example_source()
    source["cohorts"]["third"]["scenarios"]["year_month__all"]["n_species"] = 9
    with pytest.raises(ValueError, match="denominator"):
        summarize(source)


def test_rejects_changed_observed_geographic_statistic():
    source = example_source()
    source["cohorts"]["validation"]["scenarios"]["month__all"]["mean_observed_local_depletion"] = 0.11
    with pytest.raises(ValueError, match="statistic changed"):
        summarize(source)


def test_rejects_unexpected_permutation_protocol():
    source = example_source()
    source["number_of_permutations"] = 99
    with pytest.raises(ValueError, match="Permutation null"):
        summarize(source)


def test_rejects_nonmonotone_conditional_identifiability():
    source = example_source()
    source["cohorts"]["discovery"]["scenarios"]["year_month__all"]["n_identifiable_species"] = 8
    with pytest.raises(ValueError, match="identifiable"):
        summarize(source)

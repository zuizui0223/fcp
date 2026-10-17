"""Post hoc FCP 199700 taxon-leverage report must stay on real frozen estimates."""
from __future__ import annotations

import copy
import numpy as np
import pandas as pd
import pytest

from scripts.analysis import audit_fcp_199700_spatial_species_leverage_20261017 as m


def fixture():
    values = {
        "discovery": [0.2, -0.1],
        "validation": [0.1, 0.05],
        "prospective": [0.0, 0.02],
        "wave4": [0.08, -0.02],
    }
    rows = []
    expected = {}
    for i, label in enumerate(m.COHORTS):
        effect = sum(values[label]) / len(values[label])
        expected[label] = (2, effect, 0.005)
        for j, v in enumerate(values[label]):
            rows.append({
                "report_cohort": label,
                "inat_taxon_id": i*1000+j+1,
                "n_classifiable": 45,
                "local_pairs_50km": 40,
                "D_pair": 0.5,
                "D_local_50km": 0.5 - v,
                "depletion_50km": v,
            })
    expected["pooled"] = (8, sum(sum(v) for v in values.values()) / 8, 0.005)
    original = {
        "results": {
            "pooled_n_spatial_eligible_species": 8,
            "cohorts": {
                name: {
                    "radii": {
                        "50": {"n_species": 2, "mean_depletion": expected[name][1],
                               "mean_depletion_p_upper": 0.005}
                    }
                } for name in m.COHORTS
            },
            "pooled": {
                "50": {"n_species": 8, "mean_depletion": expected["pooled"][1],
                       "mean_depletion_p_upper": 0.005}
            },
        }
    }
    return pd.DataFrame(rows), original, expected


def test_species_heterogeneity_explicit_and_no_inference_upgrade():
    species, old, expected = fixture()
    d = m.analyse(species, old, expected=expected)
    pooled = d["pooled"]
    assert d["status"] == "descriptive_exploratory_post_hoc_only"
    assert d["not_new_statistical_confirmation"] is True
    assert pooled["n_species"] == 8
    assert pooled["n_positive"] == 5
    assert pooled["n_negative"] == 2
    assert pooled["n_zero"] == 1
    assert pooled["fraction_positive"] == pytest.approx(5 / 8)
    assert pooled["mean_depletion"] == pytest.approx(0.04125)
    assert pooled["leave_largest_positive_effects_out"]["10pct"]["n_removed"] == 1
    assert pooled["leave_largest_positive_effects_out"]["10pct"]["remaining_mean_depletion"] == pytest.approx(
        (pooled["mean_depletion"]*8 - 0.2) / 7
    )


def test_true_50km_local_pair_eligibility_is_required():
    species, old, expected = fixture()
    species.loc[0, "local_pairs_50km"] = 29
    with pytest.raises(ValueError, match="local-pair gate drift"):
        m.analyse(species, old, expected=expected)


def test_cohort_means_and_original_null_p_are_immutable():
    species, old, expected = fixture()
    species.loc[0, "depletion_50km"] = 0.22
    species.loc[0, "D_local_50km"] = 0.28
    with pytest.raises(ValueError, match="statistic inconsistent"):
        m.analyse(species, old, expected=expected)
    species, old, expected = fixture()
    old["results"]["pooled"]["50"]["mean_depletion_p_upper"] = 0.010
    with pytest.raises(ValueError, match="pooled old four-cohort effect"):
        m.analyse(species, old, expected=expected)


def test_missing_or_reused_species_and_fake_depletion_are_rejected():
    species, old, expected = fixture()
    species.loc[2, "inat_taxon_id"] = species.loc[0, "inat_taxon_id"]
    with pytest.raises(ValueError, match="same taxon repeated"):
        m.analyse(species, old, expected=expected)
    species, old, expected = fixture()
    species.loc[0, "depletion_50km"] = np.inf
    with pytest.raises(ValueError, match="Dpair-Dlocal"):
        m.analyse(species, old, expected=expected)
    species, old, expected = fixture()
    species = species.drop(columns="D_pair")
    with pytest.raises(ValueError, match="source schema changed"):
        m.analyse(species, old, expected=expected)


def test_source_check_fails_on_nonfinite_or_out_of_range_taxon_values():
    with pytest.raises(ValueError, match="invalid per-species depletion"):
        m.summarise(pd.Series([float("nan")]), pd.Series([1]))
    with pytest.raises(ValueError, match="invalid per-species depletion"):
        m.summarise(pd.Series([1.01]), pd.Series([1]))
    with pytest.raises(ValueError, match="repeated taxon identity"):
        m.summarise(pd.Series([0.1, 0.3]), pd.Series([1, 1]))


def test_ties_have_deterministic_drop_top_decile():
    d=m.summarise(pd.Series([0.2,0.2,0,0,-0.1]),pd.Series([1,5,3,2,4]))
    assert d["n_species"]==5
    assert d["leave_largest_positive_effects_out"]["10pct"]["n_removed"]==1
    assert d["leave_largest_positive_effects_out"]["10pct"]["remaining_mean_depletion"]==pytest.approx(0.025)

"""Reusing original 199 vertex nulls for outcome-selected robust summaries."""
import copy

import numpy as np
import pandas as pd
import pytest

from scripts.analysis import run_fcp_199700_robust_matched_null_20261018 as m
from scripts.analysis.run_polymorphism_distributed_polymorphism_posthoc_20261007 import (
    analyse_species,
)


def _synthetic():
    labels = ("white", "red_pink", "blue_purple", "yellow_orange")
    records, old_rows, expected = [], [], {}
    for i, cohort in enumerate(m.COHORTS):
        obs, nulls = [], []
        for j in range(2):
            taxon = i * 10 + j + 1
            for p in range(50):
                records.append({
                    "fcp_source_cohort": cohort,
                    "inat_taxon_id": taxon,
                    "species": f"synthetic-{cohort}-{j}",
                    "photo_id": 1000000 + taxon * 100 + p,
                    "latitude": float(i + p * .001),
                    "longitude": float(j + p * .001),
                    "morph": labels[(p + j) % len(labels)],
                    "global_classifiable": True,
                })
        subset = pd.DataFrame(records).loc[
            lambda z: z.fcp_source_cohort.eq(cohort)
        ]
        for taxon, g in subset.groupby("inat_taxon_id", sort=True):
            row, per = analyse_species(g, m.SEED_LABELS[cohort])
            row["report_cohort"] = cohort
            old_rows.append(row)
            obs.append(float(row["depletion_50km"]))
            nulls.append(per[50.0])
        # All synthetic photos are within 50km, so global-local depletion=0.
        expected[cohort] = (2, 0.0, 1.0)
    expected["pooled"] = (8, 0.0, 1.0)
    old = pd.DataFrame(old_rows)
    frozen_stat = {
        "schema": "fcp_149900_plus_wave4_original_50km_replayed_all_species_v1",
        "total_photo_rows": 199700,
        "all_original_frozen_cohort_50km_values_reproduced": True,
        "posthoc_primary_50km_remains_posthoc": True,
        "original_method": {
            "permutations": 199,
            "equal_species_weights": True,
            "min_local_pairs": 30,
        },
        "results": {
            "pooled_n_spatial_eligible_species": 8,
            "cohorts": {
                c: {"radii": {"50": {
                    "n_species": 2,
                    "mean_depletion": 0.0,
                    "mean_depletion_p_upper": 1.0,
                }}} for c in m.COHORTS
            },
            "pooled": {"50": {
                "n_species": 8,
                "mean_depletion": 0.0,
                "mean_depletion_p_upper": 1.0,
            }},
        },
    }
    return pd.DataFrame(records), old, frozen_stat, expected


def test_identical_null_selection_on_every_world():
    d = np.array([0.5, 0.2, 0.1, -0.1, 0.0])
    z = np.tile(d.reshape(-1, 1), (1, 199))
    # Opposing originals become the highest null-effect species in half the worlds.
    z[:, :100] = z[:, :100][::-1, :]
    a = m.statistics(d, z)
    ten = a["drop_top_10pct"]
    twenty = a["drop_top_20pct"]
    assert ten["n_removed_highest_effect"] == 1
    assert twenty["n_removed_highest_effect"] == 1  # ceil(5 * 0.2)
    assert ten["observed"] == pytest.approx(np.mean(np.sort(d)[:-1]))
    assert ten["null_mean"] == pytest.approx(np.mean(np.sort(z, axis=0)[:-1, :]))
    assert a["fraction_positive"]["observed"] == pytest.approx(3/5)
    assert all(item["null_selection_recomputed_in_every_permutation"] for item in a.values())
    assert all(item["upper_p_exploratory"] == 1.0 for item in a.values())


def test_frozen_original_mean_null_reproduces_before_extra_stats():
    data, frozen_species, stat, expected = _synthetic()
    result = m.matched_from_photos(data, frozen_species, stat, expected=expected)
    assert result["all_cohorts_old_mean_nulls_replayed"] is True
    assert result["status"] == "post_hoc_exploratory_not_multiple_test_corrected"
    assert result["tests"]["pooled"]["eligible_species"] == 8
    for k in (*m.COHORTS, "pooled"):
        test = result["tests"][k]["exploratory_matched_null"]
        assert test["mean"]["observed"] == pytest.approx(0.0)
        assert test["median"]["upper_p_exploratory"] == 1.0
        assert test["fraction_positive"]["upper_p_exploratory"] == 1.0


def test_unknown_species_or_changed_classification_refused():
    data, frozen, stat, expected = _synthetic()
    with pytest.raises(ValueError, match="taxon absent"):
        bad = frozen.copy()
        bad.loc[0, "inat_taxon_id"] = 9999
        m.matched_from_photos(data, bad, stat, expected=expected)
    with pytest.raises(ValueError, match="colour labels changed"):
        bad = data.copy()
        bad.loc[0, "morph"] = "red_pink"
        m.matched_from_photos(bad, frozen, stat, expected=expected)


def test_source_family_and_original_posthoc_radius_cannot_change():
    data, frozen, stat, expected = _synthetic()
    changed = copy.deepcopy(stat)
    changed["original_method"]["permutations"] = 999
    with pytest.raises(ValueError, match="inference protocol"):
        m.matched_from_photos(data, frozen, changed, expected=expected)
    changed = copy.deepcopy(stat)
    changed["posthoc_primary_50km_remains_posthoc"] = False
    with pytest.raises(ValueError, match="inference protocol"):
        m.matched_from_photos(data, frozen, changed, expected=expected)
    changed = copy.deepcopy(stat)
    changed["results"]["cohorts"]["wave4"]["radii"]["50"]["mean_depletion"] = 0.09
    with pytest.raises(ValueError, match="historical 50km mean null"):
        m.matched_from_photos(data, frozen, changed, expected=expected)


def test_all_summaries_are_calibrated_against_full_shape_nulls():
    values = np.array([0.1, -0.1, 0.1, -0.1], float)
    z = np.zeros((4, 198))
    with pytest.raises(ValueError, match="species/null size mismatch"):
        m.statistics(values, z)
    z = np.zeros((4, 199))
    z[0, 0] = np.nan
    with pytest.raises(ValueError, match="missing, nonfinite"):
        m.statistics(values, z)
    z = np.zeros((4, 199))
    a = m.statistics(values, z)
    assert a["median"]["upper_p_exploratory"] == 1.0
    assert a["fraction_positive"]["upper_p_exploratory"] == pytest.approx(.005)
    assert a["drop_top_10pct"]["observed"] == pytest.approx(-1/30)


def test_top_tail_p_uses_new_selection_for_each_permutation():
    x = np.array([0.9, 0.2, 0.1, -0.1, -0.2, 0.3, 0.0, 0.1, 0.2, 0.1])
    z = np.zeros((len(x), 199))
    for col in range(199):
        z[:, col] = .1 * ((np.arange(len(x)) + col) % 7 - 3)
    row = m.statistics(x, z)["drop_top_20pct"]
    k = 2
    expected = np.sort(z, axis=0)[:-k].mean(axis=0)
    observed = np.sort(x)[:-k].mean()
    assert row["observed"] == pytest.approx(observed)
    assert row["null_mean"] == pytest.approx(expected.mean())
    assert row["upper_p_exploratory"] == pytest.approx(
        (1 + (expected >= observed).sum()) / 200
    )


def test_holm_adjustment_within_four_exploratory_robust_tests():
    values = np.array([.1, -.1, .1, -.1, .2], dtype=float)
    null = np.zeros((5, 199))
    result = m.statistics(values, null)
    tested = ("median", "fraction_positive", "drop_top_10pct", "drop_top_20pct")
    for metric in tested:
        q = result[metric]
        assert q["multiplicity_family_size"] == 4
        assert q["holm_p_within_four_robust_tests_exploratory"] >= q["upper_p_exploratory"]
        assert 0 <= q["holm_p_within_four_robust_tests_exploratory"] <= 1
    assert result["mean"]["holm_p_within_four_robust_tests_exploratory"] is None
    assert result["mean"]["multiplicity_family_size"] == 0

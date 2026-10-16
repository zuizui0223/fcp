"""Synthetic controls for retrospective FCP source-overlap exclusion analysis.

These are not original flowers and must not be reported as an ecological test.
"""
import numpy as np
import pandas as pd

from scripts.analysis.fcp_mois_overlap_exclusion_spatial_20261016 import run_on_df, summarize_valid_50km_species


def example(n=50):
    return pd.DataFrame({
        "inat_taxon_id": [123456] * n,
        "species": ["Testus example"] * n,
        "photo_id": [str(k) for k in range(n)],
        "morph": ["white"] * (n // 2) + ["blue_purple"] * (n - n // 2),
        "latitude": [0.0] * n,
        "longitude": [0.0] * (n // 2) + [10.0] * (n - n // 2),
    })


def test_recomputes_local_depletion_after_identifer_only_exclusion():
    d = example()
    original, first = run_on_df(d, "discovery")
    keep = ~d.photo_id.isin({"0", "1"})
    after, second = run_on_df(d.loc[keep], "discovery")
    assert first["n_species"] == second["n_species"] == 1
    assert first["mean_depletion"] > .4
    assert second["mean_depletion"] > .4
    assert first["mean_depletion_p_upper"] <= .05
    assert second["mean_depletion_p_upper"] <= .05
    assert int(after.iloc[0].n_classifiable) == 48
    assert int(original.iloc[0].n_classifiable) == 50


def test_sources_with_less_than_40_classifiable_rows_become_unevaluable():
    d = example(40)
    _, result = run_on_df(d, "validation")
    assert result["n_species"] == 1
    with np.testing.assert_raises_regex(ValueError, "No taxa"):
        run_on_df(d.loc[d.photo_id != "0"], "validation")


def test_recomputations_are_deterministic_at_fixed_source_contract():
    d = example()
    a, s1 = run_on_df(d, "third")
    b, s2 = run_on_df(d.sample(frac=1, random_state=77), "third")
    # Photo ordering is fixed by the upstream cohort's photo-id sort. A
    # different row order must not be interpreted as a new observer/photo sample:
    # the primary depletion itself is order invariant.
    assert np.isclose(s1["mean_depletion"], s2["mean_depletion"])
    assert np.isclose(a.iloc[0].D_pair, b.iloc[0].D_pair)


def test_paired_bootstrap_restricts_to_species_with_finite_50km_pairs():
    # Species 2 has valid global D but lacks enough local 50-km pairs.
    # Species 3 loses local 50-km support only after photo-ID exclusion.
    original = pd.DataFrame({
        "inat_taxon_id": [1, 2, 3],
        "D_pair": [0.5, 0.4, 0.3],
        "depletion_50km": [0.1, np.nan, 0.3],
    })
    reduced = pd.DataFrame({
        "inat_taxon_id": [1, 2, 3],
        "D_pair": [0.45, 0.4, 0.25],
        "depletion_50km": [0.2, np.nan, np.nan],
    })
    r = summarize_valid_50km_species(original, reduced, "discovery")
    assert r["n_species_evaluable_both"] == 1
    assert np.isclose(r["mean_depletion_change"], 0.1)
    assert np.allclose(r["species_bootstrap_95_change_interval"], [0.1, 0.1])
    assert r["spearman_D_before_vs_after"] is None

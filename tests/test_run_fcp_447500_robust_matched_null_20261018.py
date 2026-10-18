"""Fail-closed source/null and historical-transport checks, using unittest."""
import copy
import unittest

import numpy as np
import pandas as pd

from scripts.analysis.run_fcp_447500_robust_matched_null_20261018 import (
    canonical_photos, verify_historical, verify_mean, verify_species, statistics,
)


class SourceAndNullChecks(unittest.TestCase):
    def photos(self):
        return pd.DataFrame({"inat_taxon_id": [2, 1, 1], "photo_id": [30, 20, 10],
                             "fcp_source_cohort": ["new_U100", "discovery", "discovery"]})

    def test_stable_original_photo_order_and_unique_identity(self):
        p = canonical_photos(self.photos())
        self.assertEqual(p.photo_id.tolist(), [10, 20, 30])
        duplicate = self.photos(); duplicate.loc[0, "photo_id"] = 10
        with self.assertRaises(ValueError):
            canonical_photos(duplicate)
        split = self.photos(); split.loc[1, "fcp_source_cohort"] = "validation"
        with self.assertRaises(ValueError):
            canonical_photos(split)

    def test_original_mean_null_rejects_null_or_count_drift(self):
        d = np.array([.1, .2, -.1])
        z = np.tile(np.array([-.1, .0, .1])[:, None], (1, 199))
        s = statistics(d, z)["mean"]
        frozen = {"n_species": 3, "mean_depletion": s["observed"],
                  "mean_depletion_null_mean": s["null_mean"],
                  "mean_depletion_null_q025": s["null_q025"],
                  "mean_depletion_null_q975": s["null_q975"],
                  "mean_depletion_p_upper": s["upper_p_exploratory"]}
        verify_mean(d, z, frozen, "fixture")
        for key in frozen:
            changed = dict(frozen); changed[key] += .01
            with self.assertRaises(ValueError):
                verify_mean(d, z, changed, "fixture")

    def test_each_species_original_values_and_nan_eligibility(self):
        row = {"n_classifiable": 42, "D_pair": .2, "local_pairs_50km": 40,
               "D_local_50km": .15, "depletion_50km": .05,
               "null_mean_depletion_50km": .001}
        verify_species(row, pd.Series(row), "fixture")
        for key in row:
            changed = dict(row); changed[key] += .01
            with self.assertRaises(ValueError):
                verify_species(row, pd.Series(changed), "fixture")
        missing = dict(row); missing["depletion_50km"] = np.nan
        with self.assertRaises(ValueError):
            verify_species(row, pd.Series(missing), "fixture")

    def test_trim_selection_recomputed_in_every_world(self):
        obs = np.arange(10, dtype=float) / 10 - .5
        base = np.arange(10, dtype=float) / 20 - .3
        z = np.column_stack([np.roll(base, k % 10) for k in range(199)])
        result = statistics(obs, z)
        self.assertAlmostEqual(result["drop_top_20pct"]["null_mean"], np.sort(base)[:8].mean())
        self.assertNotAlmostEqual(result["drop_top_20pct"]["null_mean"], z[:8].mean())
        self.assertEqual(result["drop_top_20pct"]["n_removed_highest_effect"], 2)

    def test_full_historical_summary_is_checked_not_only_p(self):
        obs = np.linspace(-.2, .3, 10)
        z = np.tile(np.linspace(-.1, .1, 10)[:, None], (1, 199))
        computed = {name: {"eligible_species": 10,
                           "frozen_50km_mean_upper_p": statistics(obs, z)["mean"]["upper_p_exploratory"],
                           "exploratory_matched_null": statistics(obs, z)}
                    for name in ("discovery", "validation", "prospective", "wave4", "pooled")}
        frozen = {"tests": copy.deepcopy(computed)}
        verify_historical(computed, frozen)
        frozen["tests"]["wave4"]["exploratory_matched_null"]["median"]["null_q975"] += .01
        with self.assertRaises(ValueError):
            verify_historical(computed, frozen)


if __name__ == "__main__":
    unittest.main()

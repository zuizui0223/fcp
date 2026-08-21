import unittest

import numpy as np
import pandas as pd

from fcp_pipeline.models import analyse_metrics, prepare_model_data, fit_model


class ModelHelperTests(unittest.TestCase):
    def make_data(self):
        rows = []
        for i in range(34):
            rows.append({
                "canonical_name": f"Genus{i:02d} species{i:02d}",
                "family": f"F{i % 25:02d}",
                "spatial_scale": "within_population" if i < 20 else "among_population",
                "n_climate_cells": 20 + i,
                "moisture_breadth": 35 - i + (i % 3) * 0.1,
            })
        return pd.DataFrame(rows)

    def test_prepare_model_data_standardizes_predictors(self):
        d = prepare_model_data(self.make_data(), "moisture_breadth")
        self.assertEqual(len(d), 34)
        self.assertAlmostEqual(float(d.metric_z.mean()), 0.0, places=10)
        self.assertAlmostEqual(float(d.effort_z.mean()), 0.0, places=10)
        self.assertEqual(int(d.among.sum()), 14)

    def test_clustered_glm_is_estimable(self):
        fit, d = fit_model(self.make_data(), "moisture_breadth", clustered=True)
        self.assertIsNotNone(fit)
        self.assertEqual(len(d), 34)
        self.assertTrue(np.isfinite(float(fit.params["metric_z"])))

    def test_permutation_result_is_row_order_invariant(self):
        d = self.make_data()
        shuffled = d.sample(frac=1, random_state=7).reset_index(drop=True)
        r1, _ = analyse_metrics(d, "test", ["moisture_breadth"], 99, np.random.default_rng(123))
        r2, _ = analyse_metrics(shuffled, "test", ["moisture_breadth"], 99, np.random.default_rng(123))
        self.assertEqual(r1[0]["permutation_p_two_sided"], r2[0]["permutation_p_two_sided"])
        self.assertAlmostEqual(r1[0]["odds_ratio"], r2[0]["odds_ratio"], places=12)


if __name__ == "__main__":
    unittest.main()

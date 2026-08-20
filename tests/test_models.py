import unittest

import numpy as np
import pandas as pd

from fcp_pipeline.models import prepare_model_data, fit_model


class ModelHelperTests(unittest.TestCase):
    def make_data(self):
        rows = []
        for i in range(34):
            rows.append({
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


if __name__ == "__main__":
    unittest.main()

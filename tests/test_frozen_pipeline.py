import unittest

import pandas as pd

from fcp_pipeline.constants import METRICS
from fcp_pipeline.validation import validate_frozen_dataset, validate_model_results


class FrozenPipelineTests(unittest.TestCase):
    def make_dataset(self):
        families = [f"F{i:02d}" for i in range(25)]
        rows = []
        for i in range(34):
            row = {
                "canonical_name": f"Genus{i} species{i}",
                "family": families[i % 25],
                "spatial_scale": "within_population" if i < 20 else "among_population",
                "classification_source": "baseline_unambiguous",
                "n_climate_cells": 20 + i,
            }
            for j, metric in enumerate(METRICS):
                row[metric] = i + j / 10
            rows.append(row)
        return pd.DataFrame(rows)

    def test_valid_frozen_dataset(self):
        out = validate_frozen_dataset(self.make_dataset())
        self.assertEqual(len(out), 34)
        self.assertEqual(out.family.nunique(), 25)

    def test_rejects_wrong_class_balance(self):
        d = self.make_dataset()
        d.loc[0, "spatial_scale"] = "among_population"
        with self.assertRaises(ValueError):
            validate_frozen_dataset(d)

    def test_rejects_expanded_classification_source(self):
        d = self.make_dataset()
        d.loc[0, "classification_source"] = "high_confidence_enrichment"
        with self.assertRaises(ValueError):
            validate_frozen_dataset(d)

    def test_model_results_require_five_metrics_and_9999_permutations(self):
        r = pd.DataFrame({
            "metric": METRICS,
            "analysis_status": ["complete"] * 5,
            "permutations_valid": [9999] * 5,
        })
        validate_model_results(r)
        r.loc[0, "permutations_valid"] = 9998
        with self.assertRaises(ValueError):
            validate_model_results(r)


if __name__ == "__main__":
    unittest.main()

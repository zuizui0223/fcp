import importlib.util
import json
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts" / "analysis"
for name in (
    "run_hypervolume_transfer_benchmark",
    "run_hypervolume_actual_support_blocked",
    "run_hypervolume_real_climate_synthetic_qualification",
    "run_hypervolume_real_climate_identifiability_diagnostic",
):
    if name not in sys.modules:
        spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)

import run_hypervolume_real_climate_identifiability_diagnostic as diag


class TestRealClimateIdentifiabilityDiagnostic(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.geometry = Path("inputs/geometry_only.csv")
        cls.climate = Path("climate/real_climate_only.csv")
        cls.parent_scores_path = Path("parent/scores.csv")
        cls.geo = diag.parent.ClimateGeometry(cls.geometry, cls.climate)
        cls.parent_scores = pd.read_csv(cls.parent_scores_path)

    def test_contract_is_post_failure_and_colour_closed(self):
        c = json.loads((ROOT / "docs/supporting/hypervolume_real_climate_identifiability_diagnostic_contract_v1.json").read_text())
        self.assertEqual(c["status"], "post_failure_diagnostic_frozen_before_any_diagnostic_outcome")
        self.assertEqual(c["parent_qualification_decision"], "FAIL")
        self.assertFalse(c["inputs"]["biological_colour_values_allowed"])

    def test_parent_full_score_census(self):
        self.assertEqual(len(self.parent_scores), 45000)
        self.assertEqual(set(self.parent_scores.block), set(diag.parent.BLOCKS))

    def test_exact_shared_label_reconstruction(self):
        s = diag.SCENARIOS[0]
        labels, axis, threshold = diag.shared_params_and_labels(self.geo, "thermal_regime", s, 0)
        ref = diag.parent.labels_for(self.geo, "thermal_regime", "evaluation", s, 0)
        np.testing.assert_array_equal(labels, ref)
        self.assertAlmostEqual(float(np.linalg.norm(axis)), 1.0, places=12)
        self.assertEqual(threshold.shape, (369,))

    def test_parent_score_reconstruction(self):
        block = "water_balance"
        scenario = diag.SCENARIOS[0]
        replicate = 0
        idx = self.geo.schedule("evaluation", replicate)
        x = self.geo.climate[block][idx]
        labels, _axis, _threshold = diag.shared_params_and_labels(self.geo, block, scenario, replicate)
        _bll, cll = diag.parent.axis_log_likelihood_pair(labels[idx][None, ...], diag.parent.likelihood_features(x))
        score, _f = diag.parent.score_from_ll(cll)
        row = diag.lookup_parent_row(self.parent_scores, block, scenario, replicate)
        self.assertLess(abs(float(score.mean()) - float(row.score)), 2e-12)
        for fold in range(4):
            self.assertLess(abs(float(score[0, fold]) - float(row[f"score_fold{fold}"])), 2e-12)

    def test_true_axis_augmentation_is_finite(self):
        block = "atmospheric_energy_dryness"
        scenario = diag.SCENARIOS[2]
        replicate = 1
        idx = self.geo.schedule("evaluation", replicate)
        x = self.geo.climate[block][idx]
        labels, axis, _threshold = diag.shared_params_and_labels(self.geo, block, scenario, replicate)
        y = labels[idx][None, ...]
        _b, fixed = diag.parent.axis_log_likelihood_pair(y, diag.parent.likelihood_features(x))
        true = diag.candidate_ll(y, diag.likelihood_features_axes(x, axis[None, :]))
        augmented = np.concatenate([fixed, true], axis=-1)
        score, mean_f, mass = diag.score_generic(augmented, augmented.shape[-1] - 1)
        self.assertEqual(score.shape, (1, 4))
        self.assertTrue(np.isfinite(score).all())
        self.assertTrue(np.isfinite(mean_f).all())
        self.assertTrue(np.isfinite(mass).all())
        self.assertTrue(((mass >= 0) & (mass <= 1)).all())

    def test_axis_coverage_metrics_bounded(self):
        _labels, axis, _threshold = diag.shared_params_and_labels(self.geo, "thermal_regime", diag.SCENARIOS[0], 2)
        vector = float(np.max(np.abs(diag.parent.AXES @ axis)))
        self.assertGreaterEqual(vector, 0.0)
        self.assertLessEqual(vector, 1.0 + 1e-12)
        idx = self.geo.schedule("evaluation", 2)
        corr = diag.max_projection_corr(self.geo.climate["thermal_regime"][idx][0], axis)
        self.assertGreaterEqual(corr, 0.0)
        self.assertLessEqual(corr, 1.0 + 1e-12)

    def test_species_metrics_are_bounded_and_nonnegative(self):
        block = "water_balance"
        scenario = diag.SCENARIOS[0]
        idx = self.geo.schedule("evaluation", 3)
        x = self.geo.climate[block][idx]
        labels, axis, threshold = diag.shared_params_and_labels(self.geo, block, scenario, 3)
        m = diag.species_metrics(self.geo, idx, x, labels[idx], axis, threshold)
        for key in ("train_transition_exposure", "evaluation_transition_exposure", "train_colour_variability", "evaluation_colour_variability"):
            self.assertGreaterEqual(m[key], 0.0)
            self.assertLessEqual(m[key], 1.0)
        self.assertGreaterEqual(m["train_true_projection_sd"], 0.0)
        self.assertGreaterEqual(m["evaluation_true_projection_sd"], 0.0)

    def test_only_four_prespecified_full_sharing_scenarios(self):
        self.assertEqual(len(diag.SCENARIOS), 4)
        self.assertTrue(all(float(s["shared_fraction"]) == 1.0 for s in diag.SCENARIOS))
        self.assertTrue(all(s["world"].startswith("climate_shared") for s in diag.SCENARIOS))


if __name__ == "__main__":
    unittest.main()

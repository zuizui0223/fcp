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
):
    if name not in sys.modules:
        spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)

import run_hypervolume_real_climate_synthetic_qualification as run


class TestRealClimateSyntheticQualification(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.geometry = Path("inputs/geometry_only.csv")
        cls.climate = Path("climate/real_climate_only.csv")
        cls.geo = run.ClimateGeometry(cls.geometry, cls.climate)

    def test_contract_frozen_before_outcomes(self):
        c = json.loads((ROOT / "docs/supporting/hypervolume_real_climate_synthetic_qualification_contract_v1.json").read_text())
        self.assertEqual(c["status"], "frozen_before_any_real_climate_synthetic_outcome")
        self.assertEqual(c["seed"], run.SEED)
        self.assertEqual(list(c["climate_blocks"]), list(run.BLOCKS))
        self.assertFalse(c["inputs"]["biological_colour_values_allowed"])

    def test_exact_input_hashes(self):
        self.assertEqual(run.digest(self.geometry), run.GEOMETRY_SHA)
        self.assertEqual(run.digest(self.climate), run.CLIMATE_SHA)

    def test_joint_complete_frame_and_capacity(self):
        self.assertEqual(int(self.geo.mask.sum()), 21418)
        self.assertEqual([x["training_species"] for x in self.geo.capacity], [133, 135, 141, 165])
        self.assertEqual([x["evaluation_species"] for x in self.geo.capacity], [60, 62, 59, 26])
        self.assertTrue(all(x["minimum_training_to_sector_photo_km"] >= 500.0 for x in self.geo.capacity))

    def test_all_three_blocks_use_same_photo_frame(self):
        for block in run.BLOCKS:
            self.assertTrue(np.isfinite(self.geo.climate[block][self.geo.mask]).all())

    def test_axis_grid_is_fixed_unit_r4(self):
        self.assertEqual(run.AXES.shape, (32, 4))
        np.testing.assert_allclose(np.linalg.norm(run.AXES, axis=1), 1.0, atol=1e-12)

    def test_schedule_is_deterministic_and_complete(self):
        a = self.geo.schedule("calibration", 0)
        b = self.geo.schedule("calibration", 0)
        np.testing.assert_array_equal(a, b)
        self.assertEqual(a.shape, (4, 40, 20))
        self.assertTrue(self.geo.mask[a].all())

    def test_train_and_evaluation_species_are_disjoint_in_schedule(self):
        idx = self.geo.schedule("evaluation", 3)
        for fold in range(4):
            sids = self.geo.sid[idx[fold, :, 0]]
            self.assertTrue(set(sids[:20]).issubset(self.geo.train_sid))
            self.assertTrue(set(sids[20:]).issubset(self.geo.test_sid))
            self.assertFalse(set(sids[:20]) & set(sids[20:]))

    def test_labels_are_deterministic_synthetic_only(self):
        s = run.POSITIVE[2]
        y1 = run.labels_for(self.geo, "thermal_regime", "evaluation", s, 7)
        y2 = run.labels_for(self.geo, "thermal_regime", "evaluation", s, 7)
        np.testing.assert_array_equal(y1, y2)
        self.assertEqual(y1.dtype, np.dtype(bool))
        self.assertEqual(len(y1), 36900)

    def test_scenario_family_census(self):
        self.assertEqual(len(run.NUISANCE), 12)
        self.assertEqual(len(run.POSITIVE), 6)
        self.assertEqual(len(run.MODELS), 2)
        self.assertEqual(set(run.BLOCKS), {"thermal_regime", "water_balance", "atmospheric_energy_dryness"})

    def test_likelihood_pair_is_finite(self):
        idx = self.geo.schedule("calibration", 1)
        x = self.geo.climate["water_balance"][idx]
        features = run.likelihood_features(x)
        scenarios = run.NUISANCE[:2]
        y = np.stack([run.labels_for(self.geo, "water_balance", "calibration", s, 1)[idx] for s in scenarios])
        bll, cll = run.axis_log_likelihood_pair(y, features)
        self.assertEqual(bll.shape, (2, 4, 40, 32))
        self.assertEqual(cll.shape, (2, 4, 40, 32))
        self.assertTrue(np.isfinite(bll).all())
        self.assertTrue(np.isfinite(cll).all())

    def test_scores_are_finite_and_four_fold(self):
        idx = self.geo.schedule("calibration", 2)
        x = self.geo.climate["atmospheric_energy_dryness"][idx]
        features = run.likelihood_features(x)
        y = np.stack([run.labels_for(self.geo, "atmospheric_energy_dryness", "calibration", run.NUISANCE[0], 2)[idx]])
        bll, cll = run.axis_log_likelihood_pair(y, features)
        for ll in (bll, cll):
            score, mean_f = run.score_from_ll(ll)
            self.assertEqual(score.shape, (1, 4))
            self.assertEqual(mean_f.shape, (1, 4))
            self.assertTrue(np.isfinite(score).all())
            self.assertTrue(np.isfinite(mean_f).all())
            self.assertTrue(((mean_f >= 0) & (mean_f <= 1)).all())

    def test_climate_artifact_has_no_biological_colour_columns(self):
        cols = pd.read_csv(self.climate, nrows=0).columns.tolist()
        self.assertEqual(cols[:5], ["photo_id", "species", "latitude", "longitude", "global_classifiable"])
        self.assertEqual(tuple(cols[5:]), run.ALL_ENV)


if __name__ == "__main__":
    unittest.main()

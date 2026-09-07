"""Independent algebra and geometry checks for threshold robustness."""
from pathlib import Path
import itertools
import sys
import unittest

import numpy as np
from scipy.special import expit, ndtr

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts/analysis"))
import run_hypervolume_threshold_robustness as m


class ThresholdLikelihoodTests(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(92017)
        self.x = self.rng.normal(size=(2, 3, 5, 2))
        self.y = self.rng.integers(0, 2, size=(2, 3, 5)).astype(bool)

    def slow(self, x, y, link, axis):
        proj = x @ m.AXES[axis]
        terms = []
        for offset in m.OFFSETS:
            for slope in m.SLOPES[link]:
                raw = slope * (proj - offset)
                p0 = ndtr(slope * np.tanh((proj - offset) / .5)) if link == "saturating_probit" else expit(raw)
                for sign in (1, -1):
                    p = p0 if sign == 1 else 1 - p0
                    terms.append(np.prod(np.where(y, p, 1 - p)))
        return np.log(np.mean(terms))

    def test_linear_intercept_likelihood_matches_explicit_mixture(self):
        ll = m.axis_log_likelihood_intercept(self.y, self.x, "linear_logit")
        for axis in (0, 13, 31, 47):
            self.assertAlmostEqual(ll[1, 2, axis], self.slow(self.x[1, 2], self.y[1, 2], "linear_logit", axis), places=11)

    def test_probit_intercept_likelihood_matches_explicit_mixture(self):
        ll = m.axis_log_likelihood_intercept(self.y, self.x, "saturating_probit")
        for axis in (0, 13, 31, 47):
            self.assertAlmostEqual(ll[0, 1, axis], self.slow(self.x[0, 1], self.y[0, 1], "saturating_probit", axis), places=11)

    def test_whole_species_colour_swap_invariant(self):
        for link in m.SLOPES:
            a = m.axis_log_likelihood_intercept(self.y, self.x, link)
            b = m.axis_log_likelihood_intercept(~self.y, self.x, link)
            np.testing.assert_allclose(a, b, atol=1e-12, rtol=0)

    def test_zero_coordinates_have_no_axis_information(self):
        for link in m.SLOPES:
            ll = m.axis_log_likelihood_intercept(self.y, np.zeros_like(self.x), link)
            lr = m.pooled.normalized_log_evidence(ll[None, ...])
            np.testing.assert_allclose(lr, 0, atol=1e-12, rtol=0)

    def test_all_binary_labelings_normalize(self):
        x = self.x[:1, :1, :3]
        labels = np.array(list(itertools.product([False, True], repeat=3)))
        for link in m.SLOPES:
            vals = []
            for y in labels:
                vals.append(np.exp(m.axis_log_likelihood_intercept(y.reshape(1, 1, 3), x, link)[0, 0, 7]))
            self.assertAlmostEqual(sum(vals), 1.0, places=11)

    def test_bad_shapes_and_links_rejected(self):
        with self.assertRaises(ValueError):
            m.axis_log_likelihood_intercept(self.y, self.x, "unknown")
        with self.assertRaises(ValueError):
            m.axis_log_likelihood_intercept(self.y[..., :4], self.x, "linear_logit")


class FrozenDesignTests(unittest.TestCase):
    def test_seed_scenarios_and_models(self):
        self.assertEqual(m.SEED, 2026090705)
        self.assertEqual(m.REPS, 250)
        self.assertEqual(len(m.NUISANCE), 12)
        self.assertEqual(len(m.POSITIVE), 6)
        self.assertEqual(len(m.MODELS), 3)
        ids = [m.scenario_id(x) for x in m.NUISANCE + m.POSITIVE]
        self.assertEqual(len(ids), len(set(ids)))
        self.assertEqual(m.SPEC, "0e9b50009b4706c5c6ad5360988e1c9953487ee1")

    def test_parent_source_hashes(self):
        for name, expected in m.PARENT_HASHES.items():
            self.assertEqual(m.digest(Path(m.__file__).with_name(name)), expected)

    def test_stage_seed_separation(self):
        self.assertNotEqual(m.seed_for("schedule", "calibration", 0), m.seed_for("schedule", "evaluation", 0))
        self.assertNotEqual(m.seed_for("labels", "calibration", "a", 0), m.seed_for("labels", "calibration", "b", 0))


class ActualGeometryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = Path(__file__).resolve().parents[1] / "inputs" / "geometry_only.csv"
        if not path.exists():
            raise unittest.SkipTest("coordinate-only input not present")
        cls.geo = m.support.Geometry(path)

    def test_capacity_and_support_axes(self):
        self.assertEqual(len(self.geo.species), 369)
        self.assertEqual(int(self.geo.mask.sum()), 21424)
        axes = m.support_aligned_axes(self.geo)
        self.assertEqual(axes.shape, (369, 2))
        np.testing.assert_allclose(np.linalg.norm(axes, axis=1), 1, atol=1e-12)

    def test_schedule_preserves_roles_region_buffer_and_unique_photos(self):
        idx = m.schedule(self.geo, "unit", 0)
        for fold in range(4):
            train = set(self.geo.sid[idx[fold, :20, 0]])
            test = set(self.geo.sid[idx[fold, 20:, 0]])
            self.assertTrue(train <= self.geo.train_sid)
            self.assertTrue(test <= self.geo.test_sid)
            self.assertFalse(train & test)
            self.assertTrue(self.geo.mask[idx[fold]].all())
            self.assertTrue((self.geo.sector[idx[fold, 20:]] == fold).all())
            for j in range(40):
                self.assertEqual(len(set(idx[fold, j])), 20)

    def test_generators_are_reproducible_and_scenario_distinct(self):
        aligned = m.support_aligned_axes(self.geo)
        a = m.labels_for(self.geo, "unit", m.NUISANCE[7], 0, aligned)
        b = m.labels_for(self.geo, "unit", m.NUISANCE[7], 0, aligned)
        c = m.labels_for(self.geo, "unit", m.NUISANCE[8], 0, aligned)
        np.testing.assert_array_equal(a, b)
        self.assertEqual(a.dtype, np.bool_)
        self.assertEqual(a.shape, (36900,))
        self.assertFalse(np.array_equal(a, c))


if __name__ == "__main__":
    unittest.main()

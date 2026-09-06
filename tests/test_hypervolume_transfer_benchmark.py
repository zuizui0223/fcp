import importlib.util
from pathlib import Path
import unittest
import numpy as np
from sklearn.metrics import adjusted_rand_score

PATH = Path(__file__).resolve().parents[1] / 'scripts/analysis/run_hypervolume_transfer_benchmark.py'
spec = importlib.util.spec_from_file_location('hypervolume_transfer', PATH)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

class TestHypervolumeTransfer(unittest.TestCase):
    def test_ari_matches_sklearn(self):
        rng = np.random.default_rng(8)
        for _ in range(100):
            a, b = rng.integers(0, 2, (2, 20)).astype(bool)
            self.assertAlmostEqual(float(m.adjusted_rand_binary(a, b)), adjusted_rand_score(a, b), places=12)

    def test_perfect_and_reverse_partition(self):
        a = np.array([0]*10+[1]*10, dtype=bool)
        self.assertEqual(float(m.adjusted_rand_binary(a, a)), 1)
        self.assertEqual(float(m.adjusted_rand_binary(a, ~a)), 1)
        self.assertEqual(float(m.adjusted_rand_binary(a, np.ones_like(a))), 0)

    def test_constant_truth_is_no_information(self):
        a = np.zeros(20, dtype=bool)
        self.assertEqual(float(m.adjusted_rand_binary(a, a)), 0)

    def test_affinity_identity(self):
        mu = np.zeros((1, 2, 3))
        c = np.broadcast_to(np.eye(3), (1, 2, 3, 3)).copy()
        self.assertAlmostEqual(float(m.affinity(mu, c)[0]), 1.)
        mu[0, 1, 0] = 2
        self.assertAlmostEqual(float(m.affinity(mu, c)[0]), np.exp(-0.5))

    def test_environment_even_and_finite(self):
        r = np.random.default_rng(2)
        x = m.unit_vectors(r, (100, 3))
        np.testing.assert_allclose(m.environment(x), m.environment(-x), rtol=1e-14, atol=1e-14)
        self.assertTrue(np.isfinite(m.environment(x)).all())

    def test_generated_reproducible_and_disjoint_streams(self):
        a = m.generate_world('unit_test', 'harmonic', 'geographic_shared', 1., 0)
        b = m.generate_world('unit_test', 'harmonic', 'geographic_shared', 1., 0)
        for x, y in zip(a, b):
            np.testing.assert_array_equal(x, y)
        self.assertNotEqual(m.seed_for('calibration', 0), m.seed_for('evaluation', 0))

    def test_geometry_independent_of_outcome(self):
        a = m.generate_world('unit_test', 'harmonic', 'geographic_shared', 1., 0)
        b = m.generate_world('unit_test', 'harmonic', 'no_structure', 2., 0)
        for i in (0, 1):
            np.testing.assert_array_equal(a[i], b[i])
            np.testing.assert_array_equal(m.log_total_volume(a[i]), m.log_total_volume(b[i]))

    def test_score_label_swap_invariant(self):
        g, e, y = m.generate_world('unit_test', 'harmonic', 'geographic_shared', 1., 7)
        a = m.score_domain(g[None], y[None])
        flips = np.arange(60) % 3 == 0
        y2 = y.copy(); y2[flips] = ~y2[flips]
        b = m.score_domain(g[None], y2[None])
        for key in a:
            np.testing.assert_allclose(a[key], b[key], rtol=1e-12, atol=1e-12)

    def test_target_labels_do_not_enter_training_fit(self):
        g, e, y = m.generate_world('unit_test', 'harmonic', 'geographic_shared', 1., 2)
        mu, c, _ = m.fit_classes(g[:30], y[:30])
        y2 = y.copy(); y2[30:] = ~y2[30:]
        mu2, c2, _ = m.fit_classes(g[:30], y2[:30])
        np.testing.assert_array_equal(mu, mu2)
        np.testing.assert_array_equal(c, c2)

    def test_low_information_remains_in_denominator(self):
        g, e, y = m.generate_world('unit_test', 'harmonic', 'geographic_shared', 1., 9)
        score = m.score_domain(g[None], np.zeros_like(y)[None])
        self.assertEqual(float(score['transfer'][0]), 0)
        self.assertEqual(float(score['valid_species_fraction'][0]), 0)

    def test_covariance_positive_definite(self):
        x = np.ones((3, 20, 3))
        y = np.tile(np.arange(20)%2 == 0, (3, 1))
        mu, c, valid = m.fit_classes(x, y)
        self.assertTrue((np.linalg.eigvalsh(c) > 0).all())
        self.assertTrue(valid.all())

    def test_batch_matches_single(self):
        worlds = [m.generate_world('unit_test', 'harmonic', 'mixed_specific', .5, i) for i in range(2)]
        batch = m.score_domain(np.stack([z[0] for z in worlds]), np.stack([z[2] for z in worlds]))
        for i, (g, e, y) in enumerate(worlds):
            single = m.score_domain(g[None], y[None])
            for key in single:
                np.testing.assert_allclose(single[key][0], batch[key][i], atol=1e-12)

    def test_invalid_mapping_rejected(self):
        with self.assertRaises(ValueError):
            m.environment(np.zeros((20,3)), 'other')
        with self.assertRaises(ValueError):
            m.generate_world('evaluation','harmonic','confounded_shared',1.,0)

    def test_transfer_against_independent_scipy_logpdf(self):
        from scipy.stats import multivariate_normal
        g, e, y = m.generate_world('unit_test', 'harmonic', 'mixed_specific', 1., 23)
        mu, cov, valid = m.fit_classes(g, y)
        values = []
        for i in range(30):
            for j in range(30, 60):
                if not (valid[i] and valid[j]):
                    values.append(0.); continue
                l0 = multivariate_normal.logpdf(g[j], mean=mu[i, 0], cov=cov[i, 0])
                l1 = multivariate_normal.logpdf(g[j], mean=mu[i, 1], cov=cov[i, 1])
                pred = l1 > l0
                values.append(0. if pred.all() or not pred.any() else adjusted_rand_score(y[j], pred))
        self.assertAlmostEqual(m.score_domain(g[None], y[None])['transfer'][0], np.mean(values), places=12)

    def test_wilson_non_degenerate_at_zero_and_one(self):
        self.assertGreater(m.wilson(0, 500)[1], 0)
        self.assertLess(m.wilson(500, 500)[0], 1)

if __name__ == '__main__':
    unittest.main()

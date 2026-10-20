"""Numerical and leakage checks for the prespecified transport analyses."""
import unittest

import numpy as np
import pandas as pd

try:
    from scripts.analysis import run_fcp_environment_transport_20261019 as m
except ImportError:
    m = None


class TransportTests(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(m, "transport implementation is not present")

    def test_units_preserve_log_odds_and_covariance_quadratic_form(self):
        beta = np.array([.3, -.7, 1.2])
        cov = np.array([[.4, .1, -.02], [.1, .8, .03], [-.02, .03, .2]])
        sd = np.array([3., 27., 1700.])
        inc = np.array([1., 10., 1000.])
        native_difference = np.array([4., 9., -350.])
        b, v = m.convert_units(beta, cov, sd, inc)
        old_x = native_difference / sd
        new_x = native_difference / inc
        self.assertAlmostEqual(float(old_x @ beta), float(new_x @ b), places=13)
        self.assertAlmostEqual(float(old_x @ cov @ old_x), float(new_x @ v @ new_x), places=13)

    def test_heldout_genus_slopes_cannot_change_training_coefficient(self):
        records = []
        for i in range(32):
            for j, variable in enumerate(('BIO5', 'BIO14', 'SRADmean')):
                records.append({'inat_taxon_id':i, 'species':f'Genus{i // 2} species{i}',
                                'variable':variable, 'slope':.02 * i + .1 * j,
                                'variance_approx':.4 + .01 * i,
                                'primary_within_species_sd_native':[2., 20., 2000.][j]})
        rows = pd.DataFrame(records)
        fold = m.stable_fold('Genus0', m.GENUS_SEED, 10)
        before = m.learn_fold_coefficients(rows, fold)
        changed = rows.copy()
        held = changed.species.str.split().str[0].map(lambda g:m.stable_fold(g, m.GENUS_SEED, 10)).eq(fold)
        changed.loc[held, 'slope'] = 1e6
        changed.loc[held, 'variance_approx'] = 1e-9
        after = m.learn_fold_coefficients(changed, fold)
        np.testing.assert_array_equal(before['coefficient'], after['coefficient'])
        self.assertTrue(all(before.n_training_species == rows.loc[~held].inat_taxon_id.nunique()))

    def test_cell_grouping_and_antimeridian_canonicalization(self):
        lat = np.array([-5.2, -5.4, 4.1, 4.1])
        lon = np.array([12.2, 12.7, -180., 180.])
        cells, folds = m.cell_folds(lat, lon)
        self.assertEqual(cells[0], cells[1])
        self.assertEqual(folds[0], folds[1])
        self.assertEqual(cells[2], cells[3])
        self.assertEqual(folds[2], folds[3])

    def test_great_circle_distance_crosses_antimeridian(self):
        a = np.array([[0., 179.9], [89.9, 0.]])
        b = np.array([[0., -179.9], [89.9, 180.]])
        distance = m.great_circle_km(a, b)
        self.assertAlmostEqual(distance[0, 0], 22.2390160467, places=7)
        self.assertAlmostEqual(distance[1, 1], 22.2390160467, places=7)
        np.testing.assert_allclose(distance, m.great_circle_km(b, a).T, atol=1e-9)
        np.testing.assert_allclose(np.diag(m.great_circle_km(a, a)), 0.)

    def test_buffer_excludes_every_nearby_training_photo(self):
        coords = np.array([[0., 179.9], [0., -179.9], [0., 179.], [2., 179.9]])
        folds = np.array([0, 1, 1, 2])
        train0, test0, _ = m.spatial_split_indices(coords, folds, 0, 0.)
        train50, test50, distances = m.spatial_split_indices(coords, folds, 0, 50.)
        np.testing.assert_array_equal(train0, [1, 2, 3])
        np.testing.assert_array_equal(train50, [2, 3])
        np.testing.assert_array_equal(test0, test50)
        self.assertTrue(np.all(distances > 50.))
        self.assertEqual(set(train50).intersection(test50), set())

    def test_scaling_and_longitude_origin_are_fitted_only_to_training(self):
        rng = np.random.default_rng(5)
        train_climate = rng.normal(size=(30, 3))
        train_coords = np.column_stack([rng.normal(40., 2., 30), rng.normal(179., .4, 30)])
        basis = m.fit_feature_basis(train_climate, train_coords, 'geoclimate')
        self.assertEqual(basis['status'], 'ok')
        x = m.apply_feature_basis(basis, train_climate, train_coords)
        np.testing.assert_allclose(x[:, 1:].mean(axis=0), 0., atol=1e-12)
        np.testing.assert_allclose(x[:, 1:].std(axis=0), 1., atol=1e-12)
        old_mean = basis['mean'].copy()
        old_scale = basis['scale'].copy()
        old_origin = basis['longitude_origin_radians']
        extreme = m.apply_feature_basis(basis, np.ones((2, 3)) * 1000., np.array([[80., -170.], [-80., 0.]]))
        self.assertGreater(np.abs(extreme[:, 1:]).max(), 100.)
        np.testing.assert_array_equal(basis['mean'], old_mean)
        np.testing.assert_array_equal(basis['scale'], old_scale)
        self.assertEqual(basis['longitude_origin_radians'], old_origin)

    def test_longitude_unidentifiable_does_not_drop_predictors(self):
        rng = np.random.default_rng(8)
        climate = rng.normal(size=(24, 3))
        coords = np.column_stack([rng.normal(size=24), np.tile([0., 90., 180., 270.], 6)])
        self.assertEqual(m.fit_feature_basis(climate, coords, 'geography')['status'], 'geographic_basis_not_identifiable')
        self.assertEqual(m.fit_feature_basis(climate, coords, 'geoclimate')['status'], 'geographic_basis_not_identifiable')
        self.assertEqual(m.fit_feature_basis(climate, coords, 'climate')['status'], 'ok')

    def test_constant_response_training_remains_eligible(self):
        rng = np.random.default_rng(27)
        climate = rng.normal(size=(30, 3))
        coords = rng.normal(size=(30, 2))
        result = m.fit_one_model(climate, coords, np.ones(30), climate[:4], coords[:4], 'climate')
        self.assertEqual(result['status'], 'ok')
        self.assertTrue(result['constant_training_response'])
        self.assertTrue(np.isfinite(result['test_eta']).all())
        intercept = m.fit_one_model(climate, coords, np.zeros(30), climate[:4], coords[:4], 'intercept')
        self.assertEqual(intercept['status'], 'ok')
        np.testing.assert_allclose(m.expit(intercept['test_eta']), .5 / 31., atol=1e-7)

    def test_count_and_feature_gates_do_not_use_test_labels(self):
        rng = np.random.default_rng(18)
        climate = rng.normal(size=(30, 3))
        coords = rng.normal(size=(30, 2))
        result = m.fit_one_model(climate[:19], coords[:19], np.ones(19), climate[19:], coords[19:], 'intercept')
        self.assertEqual(result['status'], 'training_n_below_20')
        constant = climate.copy()
        constant[:, 1] = 0.
        result = m.fit_one_model(constant, coords, np.zeros(30), climate[:4], coords[:4], 'climate')
        self.assertEqual(result['status'], 'zero_or_near_zero_training_predictor_sd')

    def test_pairwise_loss_and_auc_use_all_pairs_and_half_ties(self):
        white = np.array([0., 2.])
        chromatic = np.array([0., 1.])
        result = m.pairwise_metrics(white, chromatic)
        exact_loss = np.mean([np.logaddexp(0., -float(a - b)) for a in white for b in chromatic])
        self.assertAlmostEqual(result['pair_logloss'], exact_loss, places=14)
        self.assertAlmostEqual(result['pair_logloss_gain'], np.log(2.) - exact_loss, places=14)
        self.assertEqual(result['pair_auc'], .625)
        zero = m.pairwise_metrics(np.zeros(3), np.zeros(5))
        self.assertEqual(zero['pair_auc'], .5)
        self.assertAlmostEqual(zero['pair_logloss_gain'], 0., places=14)
        self.assertEqual(zero['n_pairs'], 15)

    def test_paired_metrics_use_only_matched_finite_predictions(self):
        y = np.array([0., 1., 0., 1.])
        baseline = np.array([0., 0., np.nan, 1.])
        augmented = np.array([-1., 1., -2., np.nan])
        keep, loss_gain, brier_gain = m.paired_gains(y, baseline, augmented)
        np.testing.assert_array_equal(keep, [True, True, False, False])
        np.testing.assert_allclose(loss_gain, np.log(2.) - np.logaddexp(0., -1.))
        np.testing.assert_allclose(brier_gain, .25 - (1. / (1. + np.exp(1.))) ** 2)

    def test_genus_bootstrap_weights_species_not_genus_means(self):
        values = np.array([[10., np.nan], [20., 2.], [30., 8.]])
        codes = np.array([0, 0, 1])
        weights = np.array([[1, 1], [2, 0], [0, 2]])
        means = m.cluster_weighted_means(values, codes, weights)
        np.testing.assert_allclose(means[:, 0], [20., 15., 30.])
        np.testing.assert_allclose(means[:, 1], [5., 2., 8.])

    def test_smoothed_intercept_recalibration_reaches_target_at_extremes(self):
        self.assertTrue(hasattr(m, 'recalibrate_intercept'), 'smoothed recalibration is not implemented')
        for y in (np.zeros(30), np.ones(30), np.tile([0., 1.], 15)):
            eta = np.linspace(-1000., 1000., 30)
            result = m.recalibrate_intercept(eta, y)
            self.assertTrue(np.isfinite(result['shift']))
            self.assertAlmostEqual(float(m.expit(eta + result['shift']).mean()),
                                   float((y.sum() + .5) / (len(y) + 1.)), places=12)

    def test_recalibration_does_not_change_intercept_only_firth(self):
        self.assertTrue(hasattr(m, 'recalibrate_intercept'), 'smoothed recalibration is not implemented')
        y = np.r_[np.ones(3), np.zeros(27)]
        exact = np.log((3. + .5) / (27. + .5))
        result = m.recalibrate_intercept(np.repeat(exact, 30), y)
        self.assertAlmostEqual(result['shift'], 0., places=12)

    def test_recalibration_and_slopes_do_not_depend_on_test_features(self):
        rng = np.random.default_rng(101)
        climate = rng.normal(size=(40, 3))
        coords = rng.normal(size=(40, 2))
        y = (climate[:, 0] > 0.).astype(float)
        first = m.fit_one_model(climate, coords, y, climate[:3], coords[:3], 'climate')
        second = m.fit_one_model(climate, coords, y, climate[:3] * 1000., coords[:3] * 10., 'climate')
        self.assertIn('intercept_recalibration_shift', first)
        self.assertEqual(first['intercept_recalibration_shift'], second['intercept_recalibration_shift'])
        np.testing.assert_array_equal(first['beta'], second['beta'])
        self.assertAlmostEqual(first['recalibrated_training_mean_probability'], (y.sum() + .5) / 41., places=12)

    def test_exact_training_feature_boundaries_are_not_extrapolation(self):
        rng = np.random.default_rng(7)
        climate = rng.normal(size=(30, 3))
        coords = rng.normal(size=(30, 2))
        y = np.tile([0., 1.], 15)
        result = m.fit_one_model(climate, coords, y, climate, coords, 'climate')
        self.assertEqual(result['status'], 'ok')
        self.assertFalse(result['test_outside_train_feature_range'].any())


if __name__ == '__main__':
    unittest.main()

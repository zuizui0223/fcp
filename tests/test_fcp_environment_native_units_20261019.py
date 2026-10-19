"""Analytic checks for native-unit re-expression and cluster arithmetic."""
import tempfile
from pathlib import Path
import unittest

import numpy as np
import pandas as pd

from scripts.analysis.run_fcp_environment_native_units_20261019 import (
    PHYSICAL_INCREMENTS,
    convert_primary_slopes,
    means_from_multiplicities,
    reml_in_physical_units,
    verify_file,
)


def example_slopes():
    covariance = np.array([[4., 1., -.6], [1., 9., .3], [-.6, .3, 1.]])
    rows = []
    for tid, species, scales in [(1, "Alpha one", [2., 5., 200.]),
                                  (2, "Beta two", [.5, 20., 800.])]:
        for j, variable in enumerate(PHYSICAL_INCREMENTS):
            b = [1., -2., .3][j]
            rows.append({"inat_taxon_id": tid, "species": species,
                         "source_cohort": "test", "model": "primary",
                         "status": "ok", "slope": b,
                         "se_approx": np.sqrt(covariance[j, j]),
                         "variance_approx": covariance[j, j],
                         "ci_low_wald": b - 1.96 * np.sqrt(covariance[j, j]),
                         "ci_high_wald": b + 1.96 * np.sqrt(covariance[j, j]),
                         "variable": variable, "intercept": -.7,
                         "primary_within_species_sd_native": scales[j],
                         "slope_per_native_unit": b / scales[j],
                         "primary_predictor_min": 0.,
                         "primary_predictor_max": 4. * scales[j],
                         **{"cov_with_" + v: covariance[j, k]
                            for k, v in enumerate(PHYSICAL_INCREMENTS)}})
    return pd.DataFrame(rows)


class NativeUnitChecks(unittest.TestCase):
    def test_joint_covariance_and_predictions_match_analytic_conversion(self):
        old = example_slopes()
        new = convert_primary_slopes(old.sample(frac=1, random_state=2))
        variables = list(PHYSICAL_INCREMENTS)
        increment = np.array(list(PHYSICAL_INCREMENTS.values()))
        for tid in [1, 2]:
            a = old.loc[old.inat_taxon_id.eq(tid)].set_index("variable").loc[variables]
            b = new.loc[new.inat_taxon_id.eq(tid)].set_index("variable").loc[variables]
            factor = increment / a.primary_within_species_sd_native.to_numpy()
            covariance = a[["cov_with_" + v for v in variables]].to_numpy()
            expected_covariance = np.diag(factor) @ covariance @ np.diag(factor)
            np.testing.assert_allclose(b.slope, a.slope * factor)
            np.testing.assert_allclose(b[["cov_with_" + v for v in variables]], expected_covariance)
            np.testing.assert_allclose(b.variance_approx, np.diag(expected_covariance))
            centered_native = np.array([[1., -3., 100.], [-2., 9., -400.]])
            old_design = centered_native / a.primary_within_species_sd_native.to_numpy()
            new_design = centered_native / increment
            np.testing.assert_allclose(old_design @ a.slope, new_design @ b.slope)
            np.testing.assert_allclose(
                np.einsum("ij,jk,ik->i", old_design, covariance, old_design),
                np.einsum("ij,jk,ik->i", new_design, expected_covariance, new_design))
            np.testing.assert_allclose(b.slope / b.se_approx, a.slope / a.se_approx)
            np.testing.assert_allclose(b.slope / increment, a.slope_per_native_unit)

    def test_zero_sd_missing_axis_and_invalid_covariance_fail_explicitly(self):
        for bad in [example_slopes().drop(index=0),
                    example_slopes().assign(primary_within_species_sd_native=0),
                    example_slopes().assign(cov_with_BIO5=-4.)]:
            with self.assertRaises(ValueError):
                convert_primary_slopes(bad)

    def test_common_reml_rescaling_preserves_analytic_tau_and_profile(self):
        y = np.array([-2., -1., 0., 1., 3.])
        variance = np.repeat(.25, len(y))
        base, *_ = reml_in_physical_units(y, variance)
        self.assertAlmostEqual(base["tau2"], y.var(ddof=1) - .25, places=6)
        self.assertAlmostEqual(base["reml_mean"], y.mean(), places=8)
        for multiplier in [1e-7, 1e7]:
            result, *_ = reml_in_physical_units(y * multiplier, variance * multiplier**2)
            for name in ["reml_mean", "reml_mean_ci_low", "reml_mean_ci_high", "tau"]:
                self.assertAlmostEqual(result[name] / multiplier, base[name], places=6)
            for name in ["tau2", "tau2_profile_ci_low", "tau2_profile_ci_high"]:
                self.assertAlmostEqual(result[name] / multiplier**2, base[name], places=5)

    def test_zero_heterogeneity_retains_modified_hksj_and_unknown_sign_fraction(self):
        result, *_ = reml_in_physical_units(np.zeros(6), np.full(6, .04))
        self.assertEqual(result["tau2"], 0.)
        self.assertEqual(result["hksj_q_used"], 1.)
        self.assertAlmostEqual(result["reml_mean_se_mhksj"], .2 / np.sqrt(6))
        self.assertIsNone(result["working_normal_latent_positive_fraction"])

    def test_genus_cluster_means_keep_species_weight_and_shared_axes(self):
        # Alpha has values 1,3; Beta has 10. Resampling a genus retains its size.
        multiplicities = np.array([[1, 1], [2, 0], [0, 2]])
        sums = np.array([[4., 8.], [10., 20.]])
        means = means_from_multiplicities(multiplicities, np.array([2, 1]), sums)
        np.testing.assert_allclose(means[:, 0], [14 / 3, 2., 10.])
        np.testing.assert_allclose(means[:, 1], means[:, 0] * 2)

    def test_frozen_hash_check_rejects_changed_file(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "frozen.txt"
            path.write_text("original")
            import hashlib
            expected = hashlib.sha256(b"original").hexdigest()
            self.assertTrue(verify_file(path, expected, "test")["verified"])
            path.write_text("changed")
            with self.assertRaisesRegex(ValueError, "SHA256 mismatch"):
                verify_file(path, expected, "test")


if __name__ == "__main__":
    unittest.main()

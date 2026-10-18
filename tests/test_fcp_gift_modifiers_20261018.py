import unittest
import numpy as np
import pandas as pd

from scripts.analysis.run_fcp_gift_modifiers_20261018 import admitted_levels, differences_from_multiplicities


class TraitContrastChecks(unittest.TestCase):
    def test_unknown_and_combined_traits_are_not_absence(self):
        x = pd.Series(["present", "absent", None, "variable", "annual/perennial"])
        actual = admitted_levels(x, "present", "absent")
        self.assertEqual(actual.notna().tolist(), [True, True, False, False, False])

    def test_joint_genus_draws_preserve_both_group_denominators(self):
        # Genus A has two first-group species (1,3) and one reference (10);
        # genus B one first (5) and two reference (0,2).
        counts = np.array([[1, 1], [2, 0], [0, 2]])
        actual, valid = differences_from_multiplicities(counts, np.array([2., 1.]), np.array([1., 2.]), np.array([[4.], [5.]]), np.array([[10.], [2.]]))
        np.testing.assert_allclose(actual[:, 0], [-1., -8., 4.])
        self.assertTrue(valid.all())

    def test_empty_trait_group_is_invalid_not_zero_effect(self):
        actual, valid = differences_from_multiplicities(np.array([[2, 0], [1, 1], [0, 2]]), np.array([1., 0.]), np.array([0., 1.]), np.array([[4.], [0.]]), np.array([[0.], [2.]]))
        self.assertEqual(valid.tolist(), [False, True, False])
        np.testing.assert_allclose(actual, [[2.]])


if __name__ == "__main__":
    unittest.main()

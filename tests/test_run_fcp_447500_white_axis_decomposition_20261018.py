"""Exact pair decomposition, source-null agreement and bounded descriptive gates."""
import unittest

import numpy as np
import pandas as pd

from scripts.analysis.run_fcp_447500_white_axis_decomposition_20261018 import (
    component_worlds, describe_composition, check_reconstruction, component_summary,
)


class DecompositionChecks(unittest.TestCase):
    def test_white_and_chromatic_pairs_partition_total_in_each_world(self):
        labels = np.array([0, 0, 1, 2])
        perm = np.array([[0, 1, 0, 2], [2, 1, 0, 0]])
        u, v = np.triu_indices(4, 1)
        local = np.array([True, True, False, False, True, False])
        g, loc, delta, null_loc, null_delta = component_worlds(labels, perm, u, v, local)
        np.testing.assert_allclose(g, [5 / 6, 4 / 6, 1 / 6])
        self.assertAlmostEqual(loc[0], loc[1] + loc[2])
        np.testing.assert_allclose(null_loc[:, 0], null_loc[:, 1] + null_loc[:, 2])
        np.testing.assert_allclose(delta[0], delta[1] + delta[2])
        np.testing.assert_allclose(null_delta[:, 0], null_delta[:, 1] + null_delta[:, 2])
        # Label permutations preserve global component proportions exactly.
        for p in perm:
            gg, *_ = component_worlds(p, perm, u, v, local)
            np.testing.assert_allclose(gg, g)

    def test_original_null_gate_rejects_a_changed_world(self):
        labels = np.array([0, 0, 1, 2])
        perm = np.tile(labels, (199, 1))
        u, v = np.triu_indices(4, 1)
        g, loc, d, nl, z = component_worlds(labels, perm, u, v, np.ones(6, bool))
        check_reconstruction(d, z, d[0], z[:, 0], 'fixture')
        bad = z[:, 0].copy(); bad[30] += .001
        with self.assertRaises(ValueError):
            check_reconstruction(d, z, d[0], bad, 'fixture')

    def test_composition_denominators_include_zero_classifiable_species(self):
        table = pd.DataFrame({"n_classifiable": [0, 10, 40, 50],
                              "n_states": [0, 2, 2, 2], "n_white": [0, 5, 20, 0],
                              "n_nonwhite": [0, 5, 20, 50],
                              "white_nonwhite_pairs": [0, 25, 400, 0],
                              "all_unequal_pairs": [0, 25, 400, 625]})
        result = describe_composition(table)
        self.assertEqual(result['all_measured_taxa']['n_species'], 4)
        self.assertEqual(result['all_measured_taxa']['white_in_observed_polymorphism_species'], 2)
        self.assertEqual(result['ge40_classifiable']['n_species'], 2)
        self.assertEqual(result['ge40_classifiable']['white_in_observed_polymorphism_species'], 1)
        self.assertAlmostEqual(result['ge40_classifiable']['white_in_observed_polymorphism_fraction'], .5)

    def test_two_directional_tests_have_holm_two_and_no_test_on_total(self):
        global_values = np.tile([.4, .3, .1], (3, 1))
        observed_local = global_values - [.06, .04, .02]
        null_local = np.tile(global_values[:, None, :], (1, 199, 1))
        result = component_summary(global_values, observed_local, null_local)
        self.assertIsNone(result['total']['p_upper_posthoc'])
        self.assertEqual(result['white_nonwhite']['p_upper_posthoc'], .005)
        self.assertEqual(result['white_nonwhite']['holm_two_component_p_posthoc'], .01)
        self.assertEqual(result['nonwhite_nonwhite']['holm_two_component_p_posthoc'], .01)


if __name__ == '__main__':
    unittest.main()

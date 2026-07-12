import math
import unittest

from multi_patch import (
    MultiPatchParameters,
    classify_phase,
    complete_graph,
    invasion_rates,
    path_graph,
    spectral_components,
)


class MultiPatchTests(unittest.TestCase):
    def test_two_patch_reduction(self):
        q = 2.0
        m = 0.5
        params = MultiPatchParameters(
            b=0.3,
            d=0.2,
            h=(q, -q),
            adjacency=((0.0, 1.0), (1.0, 0.0)),
            m=m,
        )
        h_a, h_b, h_bal, d_asym = spectral_components(params)
        expected = math.sqrt(q * q + m * m) - m
        self.assertAlmostEqual(h_a, expected, places=10)
        self.assertAlmostEqual(h_b, expected, places=10)
        self.assertAlmostEqual(h_bal, expected, places=10)
        self.assertAlmostEqual(d_asym, 0.0, places=10)

    def test_homogeneous_landscape_reduces_to_minimal_model(self):
        params = MultiPatchParameters(
            b=0.8,
            d=0.4,
            h=(0.0, 0.0, 0.0),
            adjacency=path_graph(3),
            m=3.0,
        )
        h_a, h_b, h_bal, d_asym = spectral_components(params)
        self.assertAlmostEqual(h_a, 0.0, places=10)
        self.assertAlmostEqual(h_b, 0.0, places=10)
        self.assertAlmostEqual(h_bal, 0.0, places=10)
        self.assertAlmostEqual(d_asym, 0.0, places=10)
        self.assertEqual(classify_phase(params), "protected_polymorphism")

    def test_asymmetric_landscape_generates_directional_shift(self):
        params = MultiPatchParameters(
            b=0.0,
            d=0.0,
            h=(3.0, -1.0, -1.0),
            adjacency=complete_graph(3),
            m=0.2,
        )
        h_a, h_b, h_bal, d_asym = spectral_components(params)
        self.assertGreater(h_bal, 0.0)
        self.assertNotAlmostEqual(d_asym, 0.0, places=8)
        a_into_b, b_into_a = invasion_rates(params)
        self.assertAlmostEqual(a_into_b, h_bal + d_asym, places=10)
        self.assertAlmostEqual(b_into_a, h_bal - d_asym, places=10)

    def test_stronger_migration_erodes_balancing_bonus(self):
        common = dict(
            b=0.0,
            d=0.0,
            h=(2.0, -2.0, 2.0, -2.0),
            adjacency=path_graph(4),
        )
        low = MultiPatchParameters(m=0.1, **common)
        high = MultiPatchParameters(m=10.0, **common)
        h_bal_low = spectral_components(low)[2]
        h_bal_high = spectral_components(high)[2]
        self.assertGreater(h_bal_low, h_bal_high)
        self.assertGreaterEqual(h_bal_high, -1e-10)


if __name__ == "__main__":
    unittest.main()

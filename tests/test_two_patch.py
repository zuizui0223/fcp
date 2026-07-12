import math
import unittest

from two_patch import (
    TwoPatchParameters,
    classify_phase,
    effective_heterogeneity,
    invasion_rates,
    simulate,
)


class TestTwoPatchModel(unittest.TestCase):
    def test_effective_heterogeneity_limits(self):
        self.assertAlmostEqual(effective_heterogeneity(2.0, 0.0), 2.0)
        self.assertGreater(effective_heterogeneity(2.0, 1.0), 0.0)
        self.assertLess(effective_heterogeneity(2.0, 1000.0), 0.01)

    def test_closed_form_invasion_rates(self):
        params = TwoPatchParameters(b=0.4, d=0.2, h=0.7, m=0.3)
        bonus = math.sqrt(0.7**2 + 0.3**2) - 0.3
        a_into_b, b_into_a = invasion_rates(params)
        self.assertAlmostEqual(a_into_b, 0.2 + 0.4 + bonus)
        self.assertAlmostEqual(b_into_a, 0.4 - 0.2 + bonus)

    def test_spatial_heterogeneity_can_rescue_polymorphism(self):
        no_space = TwoPatchParameters(b=0.2, d=0.6, h=0.0, m=0.1)
        with_space = TwoPatchParameters(b=0.2, d=0.6, h=1.0, m=0.1)
        self.assertNotEqual(classify_phase(no_space), "protected_polymorphism")
        self.assertEqual(classify_phase(with_space), "protected_polymorphism")

    def test_gene_flow_erodes_spatial_contribution(self):
        low_m = TwoPatchParameters(b=0.1, d=0.7, h=1.0, m=0.05)
        high_m = TwoPatchParameters(b=0.1, d=0.7, h=1.0, m=10.0)
        self.assertEqual(classify_phase(low_m), "protected_polymorphism")
        self.assertNotEqual(classify_phase(high_m), "protected_polymorphism")

    def test_numerical_interior_state_under_protected_polymorphism(self):
        params = TwoPatchParameters(b=0.3, d=0.0, h=0.8, m=0.2)
        p1, p2 = simulate(params, p1_0=0.01, p2_0=0.01, dt=0.01, steps=50_000)
        self.assertGreater(p1, 1e-4)
        self.assertGreater(p2, 1e-4)
        self.assertLess(p1, 1.0 - 1e-4)
        self.assertLess(p2, 1.0 - 1e-4)


if __name__ == "__main__":
    unittest.main()

import math
import unittest

from finite_population import (
    branching_establishment_probability,
    detection_probability,
    polymorphism_detection_probability,
    probability_at_least_one_established_origin,
)
from temporal import classify_temporal_phase, invasion_rates_from_series
from spatiotemporal import compress_invasion_exponents, protected_polymorphism


class TemporalTheoryTests(unittest.TestCase):
    def test_variance_does_not_change_additive_temporal_invasion(self):
        low_variance = [0.2, 0.2, 0.2, 0.2]
        high_variance = [-1.8, 2.2, -1.8, 2.2]
        self.assertEqual(invasion_rates_from_series(0.5, low_variance),
                         invasion_rates_from_series(0.5, high_variance))

    def test_temporal_mean_recovers_static_boundary(self):
        self.assertEqual(
            classify_temporal_phase(0.6, [0.2, 0.2, 0.2]),
            "protected_polymorphism",
        )
        self.assertEqual(
            classify_temporal_phase(0.1, [0.4, 0.4, 0.4]),
            "A_fixation",
        )

    def test_general_compression(self):
        b_eff, d_eff = compress_invasion_exponents(0.7, 0.3)
        self.assertAlmostEqual(b_eff, 0.5)
        self.assertAlmostEqual(d_eff, 0.2)
        self.assertTrue(protected_polymorphism(0.7, 0.3))
        self.assertGreater(b_eff, abs(d_eff))


class ObservationBridgeTests(unittest.TestCase):
    def test_detection_probability(self):
        self.assertAlmostEqual(detection_probability(0.1, 10), 1.0 - 0.9 ** 10)

    def test_polymorphism_detection_is_zero_at_fixation(self):
        self.assertEqual(polymorphism_detection_probability(0.0, 20), 0.0)
        self.assertEqual(polymorphism_detection_probability(1.0, 20), 0.0)

    def test_polymorphism_detection_increases_with_effort(self):
        self.assertLess(
            polymorphism_detection_probability(0.05, 10),
            polymorphism_detection_probability(0.05, 100),
        )

    def test_establishment_pipeline(self):
        self.assertEqual(branching_establishment_probability(-0.1), 0.0)
        p = probability_at_least_one_established_origin(1000, 1e-5, 1000, 0.1)
        self.assertGreater(p, 0.0)
        self.assertLess(p, 1.0)


if __name__ == "__main__":
    unittest.main()

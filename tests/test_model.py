import math
import unittest

from model import classify_phase, simulate


class TestPhaseModel(unittest.TestCase):
    def test_protected_polymorphism_condition(self):
        r = classify_phase(1.0, 0.2)
        self.assertEqual(r.phase, "protected_polymorphism")
        self.assertTrue(math.isclose(r.equilibrium_p_A, 0.6))

    def test_fixation_regimes(self):
        self.assertEqual(classify_phase(0.5, 1.0).phase, "A_fixation")
        self.assertEqual(classify_phase(0.5, -1.0).phase, "B_fixation")

    def test_bistability(self):
        self.assertEqual(classify_phase(-1.0, 0.0).phase, "bistability")

    def test_simulation_converges_to_internal_equilibrium(self):
        traj = simulate(
            0.1,
            balancing_strength=1.0,
            directional_bias=0.2,
            dt=0.01,
            steps=5000,
        )
        self.assertAlmostEqual(traj[-1], 0.6, places=3)


if __name__ == "__main__":
    unittest.main()

import unittest

from outcome_structure import classify_spatial_outcome, has_metapopulation_polymorphism


class TestSpatialOutcome(unittest.TestCase):
    def test_monomorphic_states(self):
        self.assertEqual(classify_spatial_outcome([1.0, 1.0]).label, "A_monomorphic")
        self.assertEqual(classify_spatial_outcome([0.0, 0.0]).label, "B_monomorphic")

    def test_geographic_mosaic_is_not_local_coexistence(self):
        out = classify_spatial_outcome([0.99, 0.01])
        self.assertEqual(out.label, "geographic_mosaic")
        self.assertTrue(has_metapopulation_polymorphism([0.99, 0.01]))
        self.assertEqual(out.locally_mixed_patches, 0)

    def test_local_coexistence(self):
        out = classify_spatial_outcome([0.5, 0.4, 0.6])
        self.assertEqual(out.label, "local_coexistence")
        self.assertEqual(out.locally_mixed_patches, 3)

    def test_mixed_spatial_polymorphism(self):
        out = classify_spatial_outcome([0.99, 0.5, 0.01])
        self.assertEqual(out.label, "mixed_spatial_polymorphism")
        self.assertEqual(out.locally_mixed_patches, 1)

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):
            classify_spatial_outcome([])
        with self.assertRaises(ValueError):
            classify_spatial_outcome([1.2])


if __name__ == "__main__":
    unittest.main()

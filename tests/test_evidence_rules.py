import unittest

from fcp_pipeline.evidence import rule_label


class EvidenceRuleTests(unittest.TestCase):
    def test_within_population(self):
        label, within, geographic = rule_label(
            "White and purple morphs coexist within populations."
        )
        self.assertEqual(label, "within_population")
        self.assertEqual((within, geographic), (1, 0))

    def test_among_population(self):
        label, within, geographic = rule_label(
            "Flower colour shows geographic variation among populations."
        )
        self.assertEqual(label, "among_population")
        self.assertEqual((within, geographic), (0, 1))

    def test_mixed(self):
        label, _, _ = rule_label(
            "Colour morphs coexist within populations and also show geographic variation."
        )
        self.assertEqual(label, "mixed")

    def test_unclear(self):
        label, _, _ = rule_label("This species has variable flowers.")
        self.assertEqual(label, "unclear")


if __name__ == "__main__":
    unittest.main()

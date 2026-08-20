import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "docs" / "jbi_manuscript.md"
SI_INDEX = ROOT / "docs" / "jbi_supporting_information_index.md"
DATA_SHA = "bdc06dd671f41ce062ebf4ba687437909d9617b268657504c1c6c5e991d417ed"


class ManuscriptConsistencyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.text = MANUSCRIPT.read_text(encoding="utf-8")
        cls.si = SI_INDEX.read_text(encoding="utf-8")

    def test_canonical_frozen_scope_is_explicit(self):
        for token in (
            "34 species",
            "25 families",
            "20 within-population",
            "14 geographically structured",
            "9,999 label permutations",
            DATA_SHA,
        ):
            self.assertIn(token, self.text)

    def test_current_primary_permutation_results_are_reported(self):
        for token in (
            "permutation p = 0.6131",
            "permutation p = 0.0423",
            "permutation p = 0.3567",
            "permutation p = 0.3859",
            "permutation p = 0.2372",
            "permutation p-value was 0.212",
        ):
            self.assertIn(token, self.text)

    def test_stale_transient_results_are_absent(self):
        stale = (
            "permutation p = 0.0475",
            "permutation p = 0.6029",
            "Holm permutation p = 0.2375",
            "workflow run `31142541223`",
            "artifact ID `8980386463`",
            "### Candidate-versus-control climatic niches",
            "### Candidate-versus-control and coarse spatial sensitivity analyses",
            "### Coarse occurrence-cloud alternatives",
            "### Broader evidence and occurrence-sampling sensitivity",
        )
        for token in stale:
            self.assertNotIn(token, self.text)

    def test_uncertainty_boundary_is_retained(self):
        for token in (
            "source-traceable, rule-derived classifications",
            "not evidence for a uniquely established moisture mechanism",
            "do not estimate fundamental physiological tolerance",
            "CR2/Satterthwaite",
        ):
            self.assertIn(token, self.text)

    def test_si_index_points_to_durable_freeze(self):
        self.assertIn("data/frozen/frozen_34species_five_metric_dataset.csv", self.si)
        self.assertIn(DATA_SHA, self.si)
        self.assertIn("Historical S1–S19 numbering", self.si)
        self.assertIn("not part of the active submission analysis", self.si)


if __name__ == "__main__":
    unittest.main()

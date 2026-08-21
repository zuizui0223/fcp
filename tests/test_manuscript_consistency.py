import csv
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
README = ROOT / "README.md"
MANUSCRIPT = ROOT / "docs" / "jbi_manuscript.md"
OLD_MANUSCRIPT = ROOT / "docs" / "jbi_manuscript_editorial_revision_v2.md"
SI_INDEX = ROOT / "docs" / "jbi_supporting_information_index.md"
CR2_SUMMARY = ROOT / "docs" / "supporting" / "cr2_satterthwaite_summary.csv"
DATA_SHA = "bdc06dd671f41ce062ebf4ba687437909d9617b268657504c1c6c5e991d417ed"


class ManuscriptConsistencyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.readme = README.read_text(encoding="utf-8")
        cls.text = MANUSCRIPT.read_text(encoding="utf-8")
        cls.si = SI_INDEX.read_text(encoding="utf-8")

    def test_readme_points_to_canonical_paper_files(self):
        for token in (
            "docs/jbi_manuscript.md",
            "docs/PIPELINE_34SPECIES.md",
            "docs/jbi_supporting_information_index.md",
            "docs/jbi_submission_completion_checklist.md",
            "data/frozen/frozen_34species_five_metric_dataset.csv",
            ".github/workflows/34species-paper.yml",
        ):
            self.assertIn(token, self.readme)

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
        self.assertFalse(OLD_MANUSCRIPT.exists(), "Superseded manuscript must not remain in the active tree")

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
        self.assertIn("cr2_satterthwaite_summary.csv", self.si)
        self.assertIn("not part of the active submission analysis", self.si)

    def test_verified_cr2_summary_is_frozen(self):
        expected = {
            "temperature_breadth": (0.814440941821187, 7.2857484381671, 0.624516849429256),
            "moisture_breadth": (0.407005114097468, 11.4800993443417, 0.0632256137089587),
            "climatic_heterogeneity": (0.677093756471007, 10.4615034395777, 0.396201860489011),
            "pca_dispersion": (0.708498853726923, 11.0853401094357, 0.455026933950585),
            "pca_hull_area": (0.572197950065939, 8.77626044207433, 0.117422659202637),
        }
        with CR2_SUMMARY.open(newline="", encoding="utf-8") as handle:
            rows = {row["metric"]: row for row in csv.DictReader(handle)}
        self.assertEqual(set(rows), set(expected))
        for metric, (odds_ratio, df, p_value) in expected.items():
            row = rows[metric]
            self.assertEqual(int(row["n_species"]), 34)
            self.assertEqual(int(row["n_families"]), 25)
            self.assertAlmostEqual(float(row["odds_ratio"]), odds_ratio, places=12)
            self.assertAlmostEqual(float(row["satterthwaite_df"]), df, places=10)
            self.assertAlmostEqual(float(row["cr2_satterthwaite_p"]), p_value, places=12)


if __name__ == "__main__":
    unittest.main()

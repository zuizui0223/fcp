from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "docs" / "jbi_manuscript.md"
SI_INDEX = ROOT / "docs" / "jbi_supporting_information_index.md"
FIGURE_PLAN = ROOT / "docs" / "FIGURE_PLAN.md"
FIGDIR = ROOT / "docs" / "figures"

EXPECTED_STEMS = (
    "figure1_geographic_context",
    "figure2_five_metric_forest",
    "figure3_raw_species_metrics",
    "figure4_leave_one_family_out",
    "figure5_inference_method_sensitivity",
    "figureS1_34_species_distribution_context",
    "figureS2_power_precision_design",
)


class FigureIntegrationTests(unittest.TestCase):
    def test_canonical_figure_files_exist(self):
        for stem in EXPECTED_STEMS:
            for ext in ("png", "pdf"):
                path = FIGDIR / f"{stem}.{ext}"
                self.assertTrue(path.exists(), f"Missing figure: {path}")
                self.assertGreater(path.stat().st_size, 1000, f"Figure appears empty: {path}")

    def test_manuscript_calls_and_captions_exist(self):
        text = MANUSCRIPT.read_text(encoding="utf-8")
        for token in (
            "Figure 1",
            "Figure 2",
            "Figure 3",
            "Figure 4",
            "Figure 5",
            "Supporting Figure S1",
            "Supporting Figure S2",
            "## Figure captions",
            "broader exact GBIF occurrence subset",
        ):
            self.assertIn(token, text)

    def test_supporting_index_tracks_canonical_figures(self):
        text = SI_INDEX.read_text(encoding="utf-8")
        self.assertIn("## 6. Canonical figures", text)
        for token in (
            "figure4_leave_one_family_out",
            "figure5_inference_method_sensitivity",
            "figureS1_34_species_distribution_context",
            "figureS2_power_precision_design",
            "broader exact GBIF",
        ):
            self.assertIn(token, text)

    def test_figure_selection_is_question_led(self):
        text = FIGURE_PLAN.read_text(encoding="utf-8")
        for token in (
            "selected from the paper question backward",
            "Figure 2 — Five-metric forest plot",
            "Figure 4 — Leave-one-family-out stability",
            "Figure 5 — Inference-method sensitivity",
            "Supporting Figure S2 — Finite-sample power/precision design diagnostic",
            "This is the central result figure",
            "What is deliberately not promoted to a main figure",
        ):
            self.assertIn(token, text)


if __name__ == "__main__":
    unittest.main()

from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "docs" / "JBI_UPSTREAM_REAUDIT_PROTOCOL.md"
SCRIPT = ROOT / "scripts" / "literature" / "build_systematic_spatial_evidence_axes.py"
SOURCE_REVIEW = ROOT / "scripts" / "literature" / "prepare_systematic_source_review.py"
HISTORICAL_RESCUE = ROOT / "scripts" / "literature" / "prepare_historical_34_source_rescue.py"
WORKFLOW = ROOT / ".github" / "workflows" / "jbi-upstream-reaudit.yml"


class UpstreamSpatialReauditTests(unittest.TestCase):
    def test_protocol_retains_mixed_and_two_axes(self):
        text = PROTOCOL.read_text(encoding="utf-8")
        self.assertIn("local_coexistence_documented", text)
        self.assertIn("geographic_structure_documented", text)
        self.assertIn("mixed_evidence", text)
        self.assertIn("two orthogonal outcomes", text)
        self.assertIn("Legacy binary analysis", text)

    def test_builder_does_not_drop_mixed(self):
        text = SCRIPT.read_text(encoding="utf-8")
        self.assertIn('return "mixed_evidence"', text)
        self.assertNotRegex(text, re.compile(r"mixed.*continue", re.I))
        self.assertIn("adjudicated_local_coexistence", text)
        self.assertIn("adjudicated_geographic_structure", text)

    def test_taxon_validation_precedes_species_state_aggregation(self):
        text = SCRIPT.read_text(encoding="utf-8")
        self.assertIn("resolve_gbif_name", text)
        self.assertIn("accepted_by_input", text)
        self.assertIn("species_records", text)
        self.assertLess(text.index("accepted_by_input"), text.index("species_records"))

    def test_source_review_separates_blind_sheet_and_coordinator_key(self):
        text = SOURCE_REVIEW.read_text(encoding="utf-8")
        self.assertIn("--review-out", text)
        self.assertIn("--key-out", text)
        self.assertIn("FORBIDDEN_REVIEW_COLUMNS", text)
        self.assertIn("automated_within_signal", text)
        self.assertIn("reviewer-facing sheet", text)

    def test_historical_rescue_hides_old_label_from_reviewer_sheet(self):
        text = HISTORICAL_RESCUE.read_text(encoding="utf-8")
        self.assertIn("--review-out", text)
        self.assertIn("--key-out", text)
        self.assertIn("historical_spatial_scale", text)
        self.assertIn("reviewer_facing_historical_label_columns", text)
        self.assertIn("coordinator-only", text.lower())

    def test_workflow_enforces_blinding_and_stops_before_climate_modeling(self):
        text = WORKFLOW.read_text(encoding="utf-8")
        self.assertIn("build_systematic_spatial_evidence_axes.py", text)
        self.assertIn("systematic_source_review_blind.csv", text)
        self.assertIn("systematic_source_review_coordinator_key.csv", text)
        self.assertIn("historical_34_source_review_blind.csv", text)
        self.assertIn("historical_34_source_review_coordinator_key.csv", text)
        self.assertIn("reviewer_facing_automated_signal_columns", text)
        self.assertNotIn("run_34species_models.py", text)
        self.assertNotIn("compute_climatic_niche_metrics.py", text)
        self.assertNotIn("WorldClim", text)


if __name__ == "__main__":
    unittest.main()

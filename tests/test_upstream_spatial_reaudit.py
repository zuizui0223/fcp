from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "docs" / "JBI_UPSTREAM_REAUDIT_PROTOCOL.md"
SCRIPT = ROOT / "scripts" / "literature" / "build_systematic_spatial_evidence_axes.py"
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

    def test_workflow_stops_before_climate_modeling(self):
        text = WORKFLOW.read_text(encoding="utf-8")
        self.assertIn("build_systematic_spatial_evidence_axes.py", text)
        self.assertNotIn("run_34species_models.py", text)
        self.assertNotIn("compute_climatic_niche_metrics.py", text)
        self.assertNotIn("WorldClim", text)


if __name__ == "__main__":
    unittest.main()

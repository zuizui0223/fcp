from pathlib import Path
import csv
import hashlib
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "docs" / "JBI_UPSTREAM_REAUDIT_PROTOCOL.md"
LEGACY_BUILDER = ROOT / "scripts" / "literature" / "build_systematic_spatial_evidence_axes.py"
SOURCE_REVIEW = ROOT / "scripts" / "literature" / "prepare_systematic_source_review.py"
HISTORICAL_RESCUE = ROOT / "scripts" / "literature" / "prepare_historical_34_source_rescue.py"
LEGACY_WORKFLOW = ROOT / ".github" / "workflows" / "jbi-upstream-reaudit.yml"
CS_WORKFLOW = ROOT / ".github" / "workflows" / "jbi-v22-coexistence-segregation-refined.yml"
CLIMATE_WORKFLOW = ROOT / ".github" / "workflows" / "jbi-cs-climate-rebuild.yml"
MODEL_WORKFLOW = ROOT / ".github" / "workflows" / "jbi-cs-models.yml"
EVIDENCE_FREEZE = ROOT / "data" / "frozen" / "jbi_cs_evidence_freeze_v22.csv"
CLIMATE_FREEZE = ROOT / "data" / "frozen" / "jbi_cs_climate_analysis_v22.csv"


class UpstreamSpatialReauditTests(unittest.TestCase):
    def test_protocol_defines_independent_positive_C_and_S_axes(self):
        text = PROTOCOL.read_text(encoding="utf-8")
        self.assertIn("C_local_coexistence_documented", text)
        self.assertIn("S_spatial_segregation_documented", text)
        self.assertIn("coexistence_and_segregation", text)
        self.assertIn("two independent positive-evidence axes", text)
        self.assertIn("C=0", text)
        self.assertIn("not documented", text)
        self.assertIn("historical binary analysis", text.lower())
        self.assertNotIn("two orthogonal outcomes", text)

    def test_legacy_builder_does_not_drop_mixed_navigation_cases(self):
        text = LEGACY_BUILDER.read_text(encoding="utf-8")
        self.assertIn('return "mixed_evidence"', text)
        self.assertNotRegex(text, re.compile(r"mixed.*continue", re.I))
        self.assertIn("adjudicated_local_coexistence", text)
        self.assertIn("adjudicated_geographic_structure", text)

    def test_taxon_validation_precedes_legacy_species_state_aggregation(self):
        text = LEGACY_BUILDER.read_text(encoding="utf-8")
        self.assertIn("resolve_gbif_name", text)
        self.assertIn("accepted_by_input", text)
        self.assertIn("species_records", text)
        self.assertLess(text.index("accepted_by_input"), text.index("species_records"))

    def test_archived_blind_review_scaffold_remains_separated(self):
        source = SOURCE_REVIEW.read_text(encoding="utf-8")
        rescue = HISTORICAL_RESCUE.read_text(encoding="utf-8")
        self.assertIn("--review-out", source)
        self.assertIn("--key-out", source)
        self.assertIn("FORBIDDEN_REVIEW_COLUMNS", source)
        self.assertIn("reviewer-facing sheet", source)
        self.assertIn("--review-out", rescue)
        self.assertIn("--key-out", rescue)
        self.assertIn("historical_spatial_scale", rescue)
        self.assertIn("coordinator-only", rescue.lower())

    def test_legacy_upstream_workflow_remains_non_modeling_provenance(self):
        text = LEGACY_WORKFLOW.read_text(encoding="utf-8")
        self.assertIn("build_systematic_spatial_evidence_axes.py", text)
        self.assertIn("systematic_source_review_blind.csv", text)
        self.assertIn("systematic_source_review_coordinator_key.csv", text)
        self.assertNotIn("run_34species_models.py", text)
        self.assertNotIn("compute_climatic_niche_metrics.py", text)

    def test_canonical_CS_workflow_uses_final_refined_chain(self):
        text = CS_WORKFLOW.read_text(encoding="utf-8")
        self.assertIn("build_v22_coexistence_segregation_refined_v5.py", text)
        self.assertIn("C_local_coexistence_documented", text)
        self.assertIn("S_spatial_segregation_documented", text)
        self.assertIn("coexistence_and_segregation", text)
        self.assertIn("Epimedium pubescens", text)

    def test_evidence_freeze_has_expected_CS_counts(self):
        with EVIDENCE_FREEZE.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 34)
        states = [row["organization_state"] for row in rows]
        self.assertEqual(states.count("local_coexistence_only"), 11)
        self.assertEqual(states.count("spatial_segregation_only"), 8)
        self.assertEqual(states.count("coexistence_and_segregation"), 15)
        self.assertEqual(sum(int(row["C_local_coexistence_documented"]) for row in rows), 26)
        self.assertEqual(sum(int(row["S_spatial_segregation_documented"]) for row in rows), 23)
        self.assertNotIn("spatial_scale", rows[0])

    def test_climate_rebuild_and_primary_models_use_frozen_CS_boundary(self):
        climate_workflow = CLIMATE_WORKFLOW.read_text(encoding="utf-8")
        model_workflow = MODEL_WORKFLOW.read_text(encoding="utf-8")
        self.assertIn("wc2.1_10m_bio.zip", climate_workflow)
        self.assertIn("--min-cells 20", climate_workflow)
        self.assertIn("jbi_cs_climate_analysis_dataset.csv", climate_workflow)
        self.assertIn("run_cs_models_v2.py", model_workflow)
        self.assertIn("9999", model_workflow)
        self.assertIn("161adbe80ee3b38a60b17cd0ad1e048eb9d454ae1cfea5825192995ef39a9a42", model_workflow)

    def test_durable_CS_climate_freeze_matches_artifact_hash_and_counts(self):
        digest = hashlib.sha256(CLIMATE_FREEZE.read_bytes()).hexdigest()
        self.assertEqual(digest, "161adbe80ee3b38a60b17cd0ad1e048eb9d454ae1cfea5825192995ef39a9a42")
        with CLIMATE_FREEZE.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 34)
        self.assertTrue(all(int(float(row["n_climate_cells"])) >= 20 for row in rows))
        self.assertEqual(sum(int(row["C_local_coexistence_documented"]) for row in rows), 26)
        self.assertEqual(sum(int(row["S_spatial_segregation_documented"]) for row in rows), 23)


if __name__ == "__main__":
    unittest.main()

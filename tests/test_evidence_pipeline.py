import unittest

from build_review_queue import build_queue, candidate_names
from ingest_evidence import merge, validate_evidence


class EvidencePipelineTests(unittest.TestCase):
    def test_candidate_is_prioritized(self):
        census = [
            {"gbif_key": "1", "canonical_name": "Background species", "family": "Poaceae", "evidence_status": "unreviewed"},
            {"gbif_key": "2", "canonical_name": "Dactylorhiza sambucina", "family": "Orchidaceae", "evidence_status": "unreviewed"},
        ]
        candidates = candidate_names([{"species": "Dactylorhiza sambucina"}])
        queue = build_queue(census, candidates, seed="x", batch_size=10, max_batches=None)
        self.assertEqual(queue[0]["canonical_name"], "Dactylorhiza sambucina")
        self.assertEqual(queue[0]["review_priority"], "0")

    def test_unknown_is_not_converted_to_monomorphic(self):
        census = [{
            "gbif_key": "1", "canonical_name": "Example species",
            "outcome_class": "unknown", "evidence_status": "unreviewed",
            "pollination_mode": "unknown", "display_opportunity": "unknown",
            "sampling_effort": "unknown", "notes": "",
        }]
        merged = merge(census, [])
        self.assertEqual(merged[0]["outcome_class"], "unknown")
        self.assertEqual(merged[0]["evidence_status"], "unreviewed")

    def test_unreviewed_monomorphism_is_rejected(self):
        with self.assertRaises(ValueError):
            validate_evidence({
                "outcome_class": "monomorphic_confirmed",
                "evidence_status": "unreviewed",
            })

    def test_verified_polymorphism_merges(self):
        census = [{
            "gbif_key": "1", "canonical_name": "Example species",
            "outcome_class": "unknown", "evidence_status": "unreviewed",
            "pollination_mode": "unknown", "display_opportunity": "unknown",
            "sampling_effort": "unknown", "notes": "",
        }]
        evidence = [{
            "gbif_key": "1", "canonical_name": "Example species",
            "outcome_class": "local_coexistence", "evidence_status": "verified",
            "pollination_mode": "insect", "display_opportunity": "showy_perianth",
            "sampling_effort": "primary_population_study", "notes": "ok",
        }]
        merged = merge(census, evidence)
        self.assertEqual(merged[0]["outcome_class"], "local_coexistence")
        self.assertEqual(merged[0]["evidence_status"], "verified")


if __name__ == "__main__":
    unittest.main()

import csv
import tempfile
import unittest
from pathlib import Path

from build_background_manifest import build_manifest


class BackgroundManifestTests(unittest.TestCase):
    def test_expands_target_counts_without_species_preassignment(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            strata = tmp_path / "strata.csv"
            out = tmp_path / "manifest.csv"
            strata.write_text(
                "stratum_id,broad_clade,family,expected_pollination_axis,target_n\n"
                "A,eudicot,FamilyA,biotic_insect,2\n"
                "B,monocot,FamilyB,abiotic_wind,3\n",
                encoding="utf-8",
            )
            n = build_manifest(strata, out)
            self.assertEqual(n, 5)
            with out.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual([r["sample_id"] for r in rows], ["A-001", "A-002", "B-001", "B-002", "B-003"])
            self.assertTrue(all(r["species"] == "" for r in rows))
            self.assertTrue(all(r["evidence_status"] == "pending_random_draw" for r in rows))

    def test_rejects_nonpositive_target(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            strata = tmp_path / "strata.csv"
            out = tmp_path / "manifest.csv"
            strata.write_text(
                "stratum_id,broad_clade,family,expected_pollination_axis,target_n\n"
                "A,eudicot,FamilyA,biotic_insect,0\n",
                encoding="utf-8",
            )
            with self.assertRaises(ValueError):
                build_manifest(strata, out)


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

import copy
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "data" / "audit_jbi_ch1_scaleup_inat_feasibility.py"
SPEC = importlib.util.spec_from_file_location("jbi_ch1_scaleup_inat_audit", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


class FakeCore:
    def __init__(self, *, selection_failure: str | None = None, qc_failure: str | None = None):
        self.selection_failure = selection_failure
        self.qc_failure = qc_failure

    def fetch_candidates(self, config, species):
        rows = [
            {
                "species": species,
                "observation_id": index + 1,
                "photo_id": f"{species}-{index + 1}",
                "photo_url": f"https://example.invalid/{species}/{index + 1}.jpg",
                "latitude": float(index % 30),
                "longitude": float(index % 60),
                "observed_month": index % 6 + 1,
                "observer": f"observer-{index % 20}",
                "spatial_cell": f"cell-{index % 10}",
                "selection_hash": f"{index:064x}",
            }
            for index in range(260)
        ]
        return {"id": 1000 + len(species), "name": species, "rank": "species"}, rows

    def select_rows(self, config, species, candidates):
        if species == self.selection_failure:
            raise RuntimeError("only 173 photographs satisfy frozen balance caps")
        return candidates[:200]

    def species_qc(self, config, species, taxon, candidates, selected):
        failed = species == self.qc_failure
        return {
            "gate_pass": not failed,
            "unique_observers": 20,
            "unique_spatial_cells": 10,
            "unique_calendar_months": 6,
            "maximum_observer_fraction": 0.05,
            "maximum_spatial_cell_fraction": 0.10,
            "maximum_month_fraction": 0.17,
            "gate_failures": ["synthetic_failure"] if failed else [],
            "observer_counts_top10": [],
            "month_counts": [],
        }


def cohort_rows():
    return [
        {
            "cohort_order": str(index + 1),
            "rank": str(index + 10),
            "canonical_name": f"Species {index + 1}",
            "family": f"Family {index % 8}",
        }
        for index in range(12)
    ]


def contract():
    return json.loads(
        (ROOT / "docs" / "supporting" / "jbi_ch1_scaleup_inat_feasibility_contract_v1.json").read_text(
            encoding="utf-8"
        )
    )


class ScaleupInatFeasibilityTests(unittest.TestCase):
    def test_all_species_pass_without_image_download(self):
        reports, rows = AUDIT.audit_cohort(FakeCore(), cohort_rows(), contract())
        self.assertEqual(len(reports), 12)
        self.assertEqual(len(rows), 2400)
        self.assertTrue(all(report["gate_pass"] for report in reports))
        self.assertTrue(all(report["candidate_images_downloaded"] is False for report in reports))
        self.assertTrue(all(report["flower_colour_pixels_inspected"] is False for report in reports))
        self.assertEqual(len({row["photo_id"] for row in rows}), 2400)

    def test_selection_failure_is_reported_not_replaced(self):
        failing = "Species 4"
        reports, rows = AUDIT.audit_cohort(
            FakeCore(selection_failure=failing), cohort_rows(), contract()
        )
        failure = next(report for report in reports if report["species"] == failing)
        self.assertEqual(failure["status"], "insufficient_balanced_metadata")
        self.assertFalse(failure["gate_pass"])
        self.assertEqual(len(rows), 2200)
        self.assertNotIn(failing, {row["species"] for row in rows})

    def test_qc_failure_is_withheld(self):
        failing = "Species 9"
        reports, rows = AUDIT.audit_cohort(
            FakeCore(qc_failure=failing), cohort_rows(), contract()
        )
        failure = next(report for report in reports if report["species"] == failing)
        self.assertEqual(failure["status"], "metadata_qc_failed")
        self.assertEqual(failure["gate_failures"], ["synthetic_failure"])
        self.assertNotIn(failing, {row["species"] for row in rows})

    def test_forbidden_colour_input_fails_closed(self):
        bad = copy.deepcopy(contract())
        bad["audit_rules"]["flower_colour_pixels_inspected"] = True
        with self.assertRaisesRegex(ValueError, "forbidden feasibility inputs"):
            AUDIT.audit_cohort(FakeCore(), cohort_rows(), bad)

    def test_candidate_manifest_contains_metadata_only_rows(self):
        _, rows = AUDIT.audit_cohort(FakeCore(), cohort_rows(), contract())
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "candidate.csv"
            AUDIT.write_candidate_manifest(path, rows)
            text = path.read_text(encoding="utf-8")
            self.assertIn("photo_url", text.splitlines()[0])
            self.assertNotIn("flower_colour", text.splitlines()[0])
            self.assertEqual(len(text.splitlines()), 2401)


if __name__ == "__main__":
    unittest.main()

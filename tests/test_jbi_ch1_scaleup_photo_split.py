from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "data" / "freeze_jbi_ch1_scaleup_photo_split.py"
SPEC = importlib.util.spec_from_file_location("jbi_ch1_scaleup_photo_freeze", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
FREEZE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(FREEZE)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_frame() -> pd.DataFrame:
    rows = []
    for species_index in range(12):
        species = f"Species {species_index + 1}"
        for photo_index in range(200):
            rows.append(
                {
                    "species": species,
                    "photo_id": f"p{species_index + 1:02d}_{photo_index + 1:03d}",
                    "observation_id": f"o{species_index + 1:02d}_{photo_index + 1:03d}",
                    "photo_url": f"https://example.invalid/{species_index + 1}/{photo_index + 1}.jpg",
                    "observer": f"observer_{photo_index % 20}",
                    "observed_on": f"2025-{photo_index % 12 + 1:02d}-01",
                    "latitude": species_index + photo_index / 1000,
                    "longitude": -100 + species_index + photo_index / 1000,
                }
            )
    return pd.DataFrame(rows)


def report(candidate: Path) -> dict:
    return {
        "protocol": "jbi-ch1-scaleup-inat-feasibility-v1",
        "status": "pass_all_12",
        "candidate_manifest_valid_for_final_freeze": True,
        "candidate_manifest_rows": 2400,
        "species_passed": 12,
        "failed_species": [],
        "candidate_manifest_sha256": digest(candidate),
        "candidate_images_downloaded": False,
        "flower_colour_pixels_inspected": False,
        "stage_a_effects_used": False,
        "stage_b_surfaces_used": False,
        "environmental_layers_used": False,
    }


class ScaleupPhotoFreezeTests(unittest.TestCase):
    def setUp(self):
        self.contract_source = (
            ROOT / "docs" / "supporting" / "jbi_ch1_scaleup_photo_split_contract_v1.json"
        )

    def run_freeze(self, tmp: Path, frame: pd.DataFrame, report_patch: dict | None = None):
        candidate = tmp / "candidate.csv"
        frame.to_csv(candidate, index=False, lineterminator="\n")
        contract = tmp / "contract.json"
        contract.write_bytes(self.contract_source.read_bytes())
        report_data = report(candidate)
        if report_patch:
            report_data.update(report_patch)
        report_path = tmp / "report.json"
        report_path.write_text(json.dumps(report_data), encoding="utf-8")
        source = tmp / "frozen_source.csv"
        split = tmp / "split.csv"
        manifest = tmp / "manifest.json"
        result = FREEZE.freeze_scaleup(
            contract_path=contract,
            report_path=report_path,
            candidate_path=candidate,
            source_output=source,
            split_output=split,
            manifest_output=manifest,
        )
        return result, source, split, manifest

    def test_freezes_exact_960_1440_split(self):
        with tempfile.TemporaryDirectory() as tmp_name:
            result, source, split_path, manifest = self.run_freeze(
                Path(tmp_name), source_frame()
            )
            split = pd.read_csv(split_path, dtype={"photo_id": str})
            self.assertEqual(result["status"], "scaleup_photo_source_and_split_frozen")
            self.assertEqual(len(split), 2400)
            self.assertEqual(int((split["split"] == "calibration").sum()), 960)
            self.assertEqual(int((split["split"] == "evaluation").sum()), 1440)
            counts = split.groupby(["species", "split"]).size().unstack()
            self.assertTrue((counts["calibration"] == 80).all())
            self.assertTrue((counts["evaluation"] == 120).all())
            self.assertEqual(source.read_bytes(), (Path(tmp_name) / "candidate.csv").read_bytes())
            self.assertEqual(json.loads(manifest.read_text())["assignment_sha256"], result["assignment_sha256"])
            self.assertFalse(result["image_pixels_read"])
            self.assertFalse(result["evaluation_opened_for_rule_tuning"])

    def test_assignment_is_independent_of_row_order_and_metadata(self):
        base = source_frame()
        with tempfile.TemporaryDirectory() as a_name, tempfile.TemporaryDirectory() as b_name:
            first, _, _, _ = self.run_freeze(Path(a_name), base)
            changed = base.sample(frac=1, random_state=42).reset_index(drop=True)
            changed["observer"] = "changed"
            changed["observed_on"] = "1900-01-01"
            changed["latitude"] = -changed["latitude"]
            second, _, _, _ = self.run_freeze(Path(b_name), changed)
            self.assertEqual(first["assignment_sha256"], second["assignment_sha256"])
            self.assertEqual(
                first["source_species_photo_id_sha256"],
                second["source_species_photo_id_sha256"],
            )

    def test_failed_species_report_cannot_freeze(self):
        with tempfile.TemporaryDirectory() as tmp_name:
            with self.assertRaisesRegex(ValueError, "feasibility precondition failed"):
                self.run_freeze(
                    Path(tmp_name),
                    source_frame(),
                    {
                        "status": "audit_complete_replacement_required",
                        "candidate_manifest_valid_for_final_freeze": False,
                        "candidate_manifest_rows": 2200,
                        "species_passed": 11,
                        "failed_species": ["Species 12"],
                    },
                )

    def test_candidate_hash_mismatch_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp_name:
            with self.assertRaisesRegex(ValueError, "SHA256"):
                self.run_freeze(
                    Path(tmp_name),
                    source_frame(),
                    {"candidate_manifest_sha256": "0" * 64},
                )

    def test_post_outcome_candidate_manifest_is_rejected(self):
        frame = source_frame()
        frame["flower_colour_state"] = "purple"
        with tempfile.TemporaryDirectory() as tmp_name:
            with self.assertRaisesRegex(ValueError, "downstream measurement outcome"):
                self.run_freeze(Path(tmp_name), frame)

    def test_contract_cannot_enable_evaluation_tuning(self):
        with tempfile.TemporaryDirectory() as tmp_name:
            tmp = Path(tmp_name)
            candidate = tmp / "candidate.csv"
            source_frame().to_csv(candidate, index=False, lineterminator="\n")
            contract_data = json.loads(self.contract_source.read_text(encoding="utf-8"))
            contract_data["status"] = "not_frozen"
            contract_path = tmp / "contract.json"
            contract_path.write_text(json.dumps(contract_data), encoding="utf-8")
            report_path = tmp / "report.json"
            report_path.write_text(json.dumps(report(candidate)), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "not frozen"):
                FREEZE.freeze_scaleup(
                    contract_path=contract_path,
                    report_path=report_path,
                    candidate_path=candidate,
                    source_output=tmp / "source.csv",
                    split_output=tmp / "split.csv",
                    manifest_output=tmp / "manifest.json",
                )


if __name__ == "__main__":
    unittest.main()

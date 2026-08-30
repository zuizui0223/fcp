from __future__ import annotations

import copy
import importlib.util
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "data" / "select_jbi_ch1_scaleup_cohort.py"
SPEC = importlib.util.spec_from_file_location("jbi_ch1_scaleup_selector", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
SELECTOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SELECTOR)


EXPECTED_SPECIES = [
    "Anemone coronaria",
    "Hesperis matronalis",
    "Digitalis purpurea",
    "Phlox drummondii",
    "Iris lutescens",
    "Silene littorea",
    "Lotus corniculatus",
    "Lobelia siphilitica",
    "Ipomopsis aggregata",
    "Anemone pavonina",
    "Platystemon californicus",
    "Castilleja coccinea",
]


class ScaleupSelectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ranking_path = ROOT / "data" / "global_flower_colour_species_ranked.csv"
        cls.ledger_path = ROOT / "docs" / "supporting" / "jbi_ch1_scaleup_species_ledger_v1.csv"
        cls.contract_path = ROOT / "docs" / "supporting" / "jbi_ch1_scaleup_contract_v1.json"
        cls.ranking = SELECTOR.read_csv(
            cls.ranking_path, SELECTOR.RANKING_REQUIRED, "ranking"
        )
        cls.ledger = SELECTOR.read_csv(
            cls.ledger_path, SELECTOR.LEDGER_REQUIRED, "ledger"
        )
        cls.contract = json.loads(cls.contract_path.read_text(encoding="utf-8"))

    def test_repository_ledger_selects_frozen_twelve(self):
        selected, skipped = SELECTOR.select_cohort(
            self.ranking, self.ledger, self.contract
        )
        self.assertEqual(
            [row["canonical_name"] for row in selected], EXPECTED_SPECIES
        )
        self.assertEqual([row["cohort_order"] for row in selected], list(range(1, 13)))
        self.assertEqual(max(row["rank"] for row in selected), 34)
        self.assertEqual(len({row["canonical_name"] for row in selected}), 12)
        self.assertLessEqual(
            max(
                sum(row["family"] == family for row in selected)
                for family in {row["family"] for row in selected}
            ),
            2,
        )
        self.assertEqual(
            [row["canonical_name"] for row in skipped], ["Linanthus parryae"]
        )

    def test_completed_development_species_cannot_reenter(self):
        ledger = copy.deepcopy(self.ledger)
        ledger.append(
            {
                "rank": "1",
                "canonical_name": "Ipomoea purpurea",
                "family": "Convolvulaceae",
                "best_doi": "10.1111/j.1558-5646.1987.tb02468.x",
                "decision": "eligible",
                "evidence_class": "natural_discrete_polymorphism",
                "photo_state_risk": "low",
                "decision_basis": "synthetic invalid test row",
            }
        )
        with self.assertRaisesRegex(ValueError, "cannot re-enter"):
            SELECTOR.select_cohort(self.ranking, ledger, self.contract)

    def test_eligible_record_requires_low_photo_state_risk(self):
        ledger = copy.deepcopy(self.ledger)
        target = next(row for row in ledger if row["canonical_name"] == "Anemone coronaria")
        target["photo_state_risk"] = "unclear"
        with self.assertRaisesRegex(ValueError, "photo_state_risk"):
            SELECTOR.select_cohort(self.ranking, ledger, self.contract)

    def test_rank_mismatch_fails_closed(self):
        ledger = copy.deepcopy(self.ledger)
        target = next(row for row in ledger if row["canonical_name"] == "Hesperis matronalis")
        target["rank"] = "999"
        with self.assertRaisesRegex(ValueError, "rank mismatch"):
            SELECTOR.select_cohort(self.ranking, ledger, self.contract)

    def test_forbidden_colour_dependent_contract_fails_closed(self):
        contract = copy.deepcopy(self.contract)
        contract["selection_used_stage_b_surfaces"] = True
        with self.assertRaisesRegex(ValueError, "forbidden post-colour input"):
            SELECTOR.select_cohort(self.ranking, self.ledger, contract)

    def test_family_cap_is_enforced_before_target(self):
        ranking = [
            {
                "rank": str(i + 1),
                "canonical_name": f"Species {i + 1}",
                "family": "Family A" if i < 4 else f"Family {i}",
                "best_doi": f"10.test/{i + 1}",
            }
            for i in range(6)
        ]
        ledger = [
            {
                **row,
                "decision": "eligible",
                "evidence_class": "natural_discrete_polymorphism",
                "photo_state_risk": "low",
                "decision_basis": "synthetic supported evidence",
            }
            for row in ranking
        ]
        contract = copy.deepcopy(self.contract)
        contract["target_new_species"] = 4
        contract["maximum_species_per_family"] = 2
        selected, skipped = SELECTOR.select_cohort(ranking, ledger, contract)
        self.assertEqual(
            [row["canonical_name"] for row in selected],
            ["Species 1", "Species 2", "Species 5", "Species 6"],
        )
        self.assertEqual(
            [row["canonical_name"] for row in skipped],
            ["Species 3", "Species 4"],
        )


if __name__ == "__main__":
    unittest.main()

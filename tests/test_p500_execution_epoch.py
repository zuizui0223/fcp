import copy
import json
import unittest
from pathlib import Path

from fcp_pipeline.p500_execution_epoch import (
    EXPECTED_FREEZE,
    advance_epoch,
    canonical_json_sha256,
    validate_preopen_freeze,
)


ROOT = Path(__file__).resolve().parents[1]
RECEIPT_PATH = ROOT / "docs/POLYMORPHISM_H2_P500_EXECUTION_EPOCH_20260915.json"


def load_receipt():
    return json.loads(RECEIPT_PATH.read_text(encoding="utf-8"))


def dispatch_event(receipt, event_id="dispatch-1"):
    return {
        "kind": "first_dispatch_record",
        "event_id": event_id,
        "timestamp_utc": "2026-09-15T06:00:00Z",
        "prior_receipt_sha256": canonical_json_sha256(receipt),
        "image_requests_started_before_record": False,
        "image_bytes_opened_before_record": False,
        "measurement_started_before_record": False,
        "H2_opened_before_record": False,
        "frozen_rows": 49900,
        "frozen_species": 499,
    }


class TestP500ExecutionEpoch(unittest.TestCase):
    def test_committed_preopen_receipt_is_valid(self):
        receipt = load_receipt()
        validate_preopen_freeze(receipt)
        self.assertEqual(receipt["state"], "PREOPEN_FROZEN")
        self.assertFalse(receipt["historical_nonaccess_claimed"])

    def test_all_scientific_constants_are_frozen(self):
        receipt = load_receipt()
        for key, expected in EXPECTED_FREEZE.items():
            self.assertEqual(receipt[key], expected)

    def test_preopen_rejects_any_opened_surface(self):
        for flag in load_receipt()["opening_flags"]:
            with self.subTest(flag=flag):
                mutated = load_receipt()
                mutated["opening_flags"][flag] = True
                with self.assertRaisesRegex(RuntimeError, "already opened"):
                    validate_preopen_freeze(mutated)

    def test_preopen_rejects_scientific_mutation(self):
        mutations = {
            "axis_refit_allowed": True,
            "primary_threshold": 0.20,
            "strict_sensitivity_threshold": 0.10,
            "structured_null_replicates": 1999,
            "primary_seed": 1,
            "frozen_rows": 49899,
            "frozen_species": 498,
            "species_replacement_allowed": True,
            "row_replacement_allowed": True,
            "target_relaxation_allowed": True,
            "rerun_after_outcome_allowed": True,
            "historical_nonaccess_claimed": True,
        }
        for field, value in mutations.items():
            with self.subTest(field=field):
                mutated = load_receipt()
                mutated[field] = value
                with self.assertRaisesRegex(RuntimeError, "frozen field mismatch"):
                    validate_preopen_freeze(mutated)

    def test_only_first_transition_is_dispatch(self):
        receipt = load_receipt()
        event = dispatch_event(receipt)
        with self.assertRaisesRegex(RuntimeError, "illegal state transition"):
            advance_epoch(receipt, "MEASUREMENT_COMPLETE", event)
        with self.assertRaisesRegex(RuntimeError, "illegal state transition"):
            advance_epoch(receipt, "INFERENCE_OPENED", event)

    def test_dispatch_must_be_recorded_before_any_access(self):
        receipt = load_receipt()
        for field in (
            "image_requests_started_before_record",
            "image_bytes_opened_before_record",
            "measurement_started_before_record",
            "H2_opened_before_record",
        ):
            with self.subTest(field=field):
                event = dispatch_event(receipt)
                event[field] = True
                with self.assertRaisesRegex(RuntimeError, "chronology violation"):
                    advance_epoch(receipt, "DISPATCH_RECORDED", event)

    def test_dispatch_is_hash_bound_to_exact_prior_receipt(self):
        receipt = load_receipt()
        event = dispatch_event(receipt)
        event["prior_receipt_sha256"] = "0" * 64
        with self.assertRaisesRegex(RuntimeError, "digest mismatch"):
            advance_epoch(receipt, "DISPATCH_RECORDED", event)

    def test_duplicate_dispatch_cannot_be_replayed(self):
        receipt = load_receipt()
        event = dispatch_event(receipt)
        dispatched = advance_epoch(receipt, "DISPATCH_RECORDED", event)
        event2 = copy.deepcopy(event)
        event2["prior_receipt_sha256"] = canonical_json_sha256(dispatched)
        with self.assertRaisesRegex(RuntimeError, "illegal state transition"):
            advance_epoch(dispatched, "DISPATCH_RECORDED", event2)

    def test_measurement_completion_requires_full_fixed_census(self):
        receipt = load_receipt()
        dispatched = advance_epoch(receipt, "DISPATCH_RECORDED", dispatch_event(receipt))
        completion = {
            "kind": "measurement_completion_record",
            "event_id": "complete-1",
            "timestamp_utc": "2026-09-15T07:00:00Z",
            "prior_receipt_sha256": canonical_json_sha256(dispatched),
            "terminal_rows": 49900,
            "unique_measurement_ids": 49900,
            "terminal_partitions": 256,
            "duplicate_measurement_ids": 0,
            "missing_measurement_ids": 0,
            "early_stopping_used": False,
            "replacement_used": False,
        }
        completed = advance_epoch(dispatched, "MEASUREMENT_COMPLETE", completion)
        self.assertEqual(completed["state"], "MEASUREMENT_COMPLETE")

        for field, bad in (
            ("terminal_rows", 49899),
            ("unique_measurement_ids", 49899),
            ("terminal_partitions", 255),
            ("duplicate_measurement_ids", 1),
            ("missing_measurement_ids", 1),
            ("early_stopping_used", True),
            ("replacement_used", True),
        ):
            with self.subTest(field=field):
                bad_event = copy.deepcopy(completion)
                bad_event[field] = bad
                with self.assertRaises(RuntimeError):
                    advance_epoch(dispatched, "MEASUREMENT_COMPLETE", bad_event)

    def test_inference_cannot_open_before_completion_receipt(self):
        receipt = load_receipt()
        dispatched = advance_epoch(receipt, "DISPATCH_RECORDED", dispatch_event(receipt))
        event = {
            "kind": "inference_opening_record",
            "event_id": "infer-1",
            "timestamp_utc": "2026-09-15T08:00:00Z",
            "prior_receipt_sha256": canonical_json_sha256(dispatched),
            "complete_terminal_census_verified": True,
            "metadata_colour_join_opened_before_record": False,
            "H2_opened_before_record": False,
        }
        with self.assertRaisesRegex(RuntimeError, "illegal state transition"):
            advance_epoch(dispatched, "INFERENCE_OPENED", event)


if __name__ == "__main__":
    unittest.main()

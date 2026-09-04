from __future__ import annotations

import pandas as pd
import pytest

from fcp_pipeline.global_capacity_recovery import (
    RETRYABLE_429,
    frozen_recovery_rows,
    merge_transport_recovery,
    nonempty_error_mask,
    retryable_429_mask,
)


def _audit() -> pd.DataFrame:
    rows = []
    for i in range(5):
        error = RETRYABLE_429 if i in {1, 3} else ""
        rows.append({
            "global_row_index": i,
            "shard_index": i % 2,
            "species": f"sp_{i}",
            "inat_taxon_id": 100 + i,
            "raw_results": 0 if error else 200,
            "locally_eligible": 0 if error else 120,
            "prior_excluded": 0,
            "wrong_taxon": 0,
            "after_observer_cap": 0 if error else 110,
            "maximum_span_km": 0.0 if error else 1000.0,
            "capped_row_id_sha256": "" if error else f"hash{i}",
            "request_error": error,
            "eligible_raw_100": False if error else True,
            "eligible_raw_80": False if error else True,
            "eligible_raw_60": False if error else True,
        })
    return pd.DataFrame(rows)


def test_retry_set_is_every_and_only_exact_429_in_global_order():
    audit = _audit().sample(frac=1.0, random_state=7).reset_index(drop=True)
    retry = frozen_recovery_rows(audit)
    assert retry["global_row_index"].tolist() == [1, 3]
    assert retryable_429_mask(audit).sum() == 2
    assert nonempty_error_mask(audit).sum() == 2


def test_transport_merge_changes_only_frozen_429_rows():
    original = _audit()
    recovered = frozen_recovery_rows(original).copy()
    recovered["request_error"] = ""
    recovered["raw_results"] = 200
    recovered["locally_eligible"] = 100
    recovered["after_observer_cap"] = 90
    recovered["maximum_span_km"] = 500.0
    recovered["capped_row_id_sha256"] = ["r1", "r3"]
    recovered["eligible_raw_100"] = False
    recovered["eligible_raw_80"] = True
    recovered["eligible_raw_60"] = True
    merged = merge_transport_recovery(original, recovered)
    assert len(merged) == len(original)
    assert nonempty_error_mask(merged).sum() == 0
    assert merged.loc[merged.global_row_index == 0, "capped_row_id_sha256"].iloc[0] == "hash0"
    assert merged.loc[merged.global_row_index == 1, "capped_row_id_sha256"].iloc[0] == "r1"
    assert bool(merged.loc[merged.global_row_index == 1, "eligible_raw_80"].iloc[0]) is True


def test_transport_merge_refuses_partial_retry_set():
    original = _audit()
    recovered = frozen_recovery_rows(original).iloc[:1].copy()
    with pytest.raises(ValueError, match="exactly the frozen 429 set"):
        merge_transport_recovery(original, recovered)


def test_transport_merge_refuses_non429_replacement():
    original = _audit()
    recovered = frozen_recovery_rows(original).copy()
    extra = original.loc[original.global_row_index == 0].copy()
    recovered = pd.concat([recovered, extra], ignore_index=True)
    with pytest.raises(ValueError, match="exactly the frozen 429 set"):
        merge_transport_recovery(original, recovered)

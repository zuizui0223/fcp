import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from fcp_pipeline.photo_first_measurement_execution import (
    ACQUISITION_FIELDS,
    WORKER_FIELDS,
    build_measurement_firewall,
    canonical_candidate_sha256,
    reassemble_complete_measurement,
    select_measurement_partition,
    validate_candidate_pixel_gate,
    validate_terminal_partition_results,
)


ROOT = Path(__file__).resolve().parents[1]
MEASUREMENT_CONTRACT = json.loads(
    (ROOT / "docs/supporting/random_photo_first_measurement_contract_v1.json").read_text()
)
EXECUTION_CONTRACT = json.loads(
    (ROOT / "docs/supporting/random_photo_first_measurement_execution_v1.json").read_text()
)


def candidate_table(n=40):
    rows = []
    for index in range(n):
        rows.append(
            {
                "cell_id": index % 8,
                "observation_id": 10000 + index,
                "photo_id": 20000 + index,
                "photo_url_large": f"https://example.org/{20000 + index}/large.jpg",
                "photo_license": "cc-by",
                "attribution": "synthetic",
                "species": f"Genus species{index % 7}",
                "inat_taxon_id": 30000 + index % 7,
                "latitude": -40.0 + index,
                "longitude": -160.0 + 5.0 * index,
                "positional_accuracy_m": 100.0,
                "observed_on": "2026-06-01",
                "observer_id": str(40000 + index),
                "observer": f"observer{index}",
            }
        )
    return pd.DataFrame(rows)


def candidate_manifest(candidate, *, pixels_may_open=True):
    return {
        "protocol": "random-photo-first-candidate-pool-v1",
        "status": "metadata_pool_frozen_before_candidate_image_pixels",
        "outcome_firewall": {
            "candidate_image_pixels_opened": False,
            "morph_used_for_selection": False,
            "species_list_fixed_or_targeted": False,
            "legacy_pr21_terminal_records_used": False,
        },
        "counts": {"observations": len(candidate)},
        "candidate_table_sha256": canonical_candidate_sha256(candidate),
        "premeasurement_h1_gate": {
            "pixels_may_open": pixels_may_open,
            "candidate_pool_can_reach_fixed_replicate_size": pixels_may_open,
        },
    }


def terminal_frame(ids):
    return pd.DataFrame(
        {
            "measurement_id": list(ids),
            "morph": ["red_pink"] * len(ids),
            "measurement_status": ["classified_four_state_morph"] * len(ids),
            "roi_status": ["automated_colour_state_admitted"] * len(ids),
            "flower_effective_pixels": [150] * len(ids),
        }
    )


def test_capacity_gate_fails_closed_before_firewall_or_pixels():
    candidate = candidate_table(10)
    manifest = candidate_manifest(candidate, pixels_may_open=False)
    with pytest.raises(
        RuntimeError, match="not_evaluable_candidate_sampling_capacity_before_pixels"
    ):
        validate_candidate_pixel_gate(
            candidate, manifest, MEASUREMENT_CONTRACT, EXECUTION_CONTRACT
        )


def test_firewall_worker_has_only_blind_fields_and_source_url_is_sealed():
    candidate = candidate_table(20)
    firewall = build_measurement_firewall(
        candidate,
        candidate_manifest(candidate),
        MEASUREMENT_CONTRACT,
        EXECUTION_CONTRACT,
    )
    assert tuple(firewall.worker_manifest.columns) == WORKER_FIELDS
    assert tuple(firewall.acquisition_key.columns) == ACQUISITION_FIELDS
    assert "species" not in firewall.worker_manifest.columns
    assert "latitude" not in firewall.worker_manifest.columns
    assert "longitude" not in firewall.worker_manifest.columns
    assert "photo_url_large" not in firewall.worker_manifest.columns
    assert "photo_url_large" in firewall.acquisition_key.columns
    assert "photo_url_large" not in firewall.metadata_join_key.columns
    assert firewall.worker_manifest.measurement_id.str.startswith("FCPR-").all()
    assert firewall.worker_manifest.measurement_id.nunique() == len(candidate)
    assert firewall.manifest["candidate_pixels_opened"] is False
    assert firewall.manifest["coordinate_colour_join_opened"] is False


def test_128_partitions_cover_each_measurement_exactly_once():
    candidate = candidate_table(500)
    firewall = build_measurement_firewall(
        candidate,
        candidate_manifest(candidate),
        MEASUREMENT_CONTRACT,
        EXECUTION_CONTRACT,
    )
    seen = []
    for shard in range(32):
        for partition in range(4):
            selected = select_measurement_partition(
                firewall.worker_manifest,
                firewall.acquisition_key,
                semantic_shard_index=shard,
                compute_partition_index=partition,
            )
            seen.extend(selected.measurement_id.astype(str).tolist())
    assert len(seen) == len(candidate)
    assert len(set(seen)) == len(candidate)
    assert set(seen) == set(firewall.worker_manifest.measurement_id.astype(str))


def test_terminal_partition_rejects_metadata_leak_and_nonterminal_state():
    ids = ["FCPR-AAA", "FCPR-BBB"]
    good = terminal_frame(ids)
    validate_terminal_partition_results(good, ids)

    leaked = good.copy()
    leaked["latitude"] = [1.0, 2.0]
    with pytest.raises(ValueError, match="leaked metadata"):
        validate_terminal_partition_results(leaked, ids)

    bad = good.copy()
    bad.loc[0, "measurement_status"] = "retry_later"
    with pytest.raises(ValueError, match="unknown morph/status"):
        validate_terminal_partition_results(bad, ids)


def test_reassembly_refuses_127_partition_receipts_before_coordinate_join():
    candidate = candidate_table(20)
    firewall = build_measurement_firewall(
        candidate,
        candidate_manifest(candidate),
        MEASUREMENT_CONTRACT,
        EXECUTION_CONTRACT,
    )
    empty = terminal_frame([])
    with pytest.raises(ValueError, match="incomplete_measurement_partitions"):
        reassemble_complete_measurement(
            [empty.copy() for _ in range(127)],
            firewall.worker_manifest,
            firewall.metadata_join_key,
        )


def test_complete_128_partition_reassembly_opens_join_only_after_exact_coverage():
    candidate = candidate_table(250)
    firewall = build_measurement_firewall(
        candidate,
        candidate_manifest(candidate),
        MEASUREMENT_CONTRACT,
        EXECUTION_CONTRACT,
    )
    partition_results = []
    for shard in range(32):
        for partition in range(4):
            selected = select_measurement_partition(
                firewall.worker_manifest,
                firewall.acquisition_key,
                semantic_shard_index=shard,
                compute_partition_index=partition,
            )
            frame = terminal_frame(selected.measurement_id.astype(str).tolist())
            validate_terminal_partition_results(
                frame, selected.measurement_id.astype(str).tolist()
            )
            partition_results.append(frame)

    result = reassemble_complete_measurement(
        partition_results,
        firewall.worker_manifest,
        firewall.metadata_join_key,
    )
    assert len(result.joined_photos) == len(candidate)
    assert result.joined_photos.measurement_id.nunique() == len(candidate)
    assert {"species", "latitude", "longitude", "morph"}.issubset(
        result.joined_photos.columns
    )
    assert result.result_manifest[
        "coordinate_colour_join_opened_after_complete_measurement"
    ] is True
    assert result.result_manifest["terminal_result_rows"] == len(candidate)
    assert result.result_manifest["legacy_pr21_terminal_records_used"] is False

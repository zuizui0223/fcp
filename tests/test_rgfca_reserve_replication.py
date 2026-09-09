"""Fail-closed reserve selection and measurement tests; no network or pixels."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from fcp_pipeline.global_measurement_budget import select_measurement_rows
from fcp_pipeline.rgfca_reserve_replication import (
    AUDIT, CONTRACT, PROTOCOL, audit_reserve, build_firewall, geometry_audit,
    overlap_counts, palette_vectors, reassemble, reserve_rows, sha,
)


def small_candidate():
    return pd.DataFrame({"inat_taxon_id": np.repeat([1, 2, 3, 4], 3),
                         "photo_id": range(12), "observation_id": range(100, 112)})


def test_entire_complement_is_order_invariant():
    frame = small_candidate()
    discovery = select_measurement_rows(frame, target_photos_per_species=3, maximum_species=2)
    a = reserve_rows(frame, discovery, target=3, budget=2)
    b = reserve_rows(frame.sample(frac=1, random_state=8), discovery, target=3, budget=2)
    pd.testing.assert_frame_equal(a, b)
    assert len(a) == 6 and not set(a.inat_taxon_id) & set(discovery.inat_taxon_id)
    assert overlap_counts(a, discovery) == {"photo_id": 0, "observation_id": 0}


def test_wrong_discovery_cannot_define_a_favorable_reserve():
    c = small_candidate()
    d = select_measurement_rows(c, target_photos_per_species=3, maximum_species=2)
    d.loc[d.index[0], "observation_id"] = 999999
    with pytest.raises(ValueError, match="original outcome-blind"):
        reserve_rows(c, d, target=3, budget=2)


def test_duplicate_observation_is_rejected():
    c = small_candidate()
    d = select_measurement_rows(c, target_photos_per_species=3, maximum_species=2)
    c.loc[1, "observation_id"] = c.loc[0, "observation_id"]
    with pytest.raises(ValueError, match="duplicated"):
        reserve_rows(c, d, target=3, budget=2)


def test_palette_mapping_and_invalid_classified_row():
    frame = pd.DataFrame({f"flower_fraction_{p}": [v] for p, v in
                          zip(("white", "yellow", "orange", "bronze", "red", "pink", "magenta", "blue", "purple"),
                              (.6, .05, .02, .03, .04, .02, .04, .1, .1))})
    frame["morph"] = "white"
    frame["measurement_status"] = "classified_four_state_morph"
    from fcp_pipeline.photo_first_measurement import REFERENCE_RGB
    for colour in REFERENCE_RGB:
        frame[f"background_palette_count_{colour}"] = 100 if colour == "green" else 0
    frame["background_effective_pixels"] = 100
    out = palette_vectors(frame)
    assert np.allclose(out[["colour_white", "colour_yellow_orange", "colour_red_pink", "colour_blue_purple"]], [[.6, .1, .1, .2]])
    frame.loc[0, "flower_fraction_white"] = np.nan
    with pytest.raises(ValueError, match="invalid"):
        palette_vectors(frame)


def test_background_keeps_green_and_excludes_flower():
    from fcp_pipeline.photo_first_measurement import REFERENCE_RGB
    from fcp_pipeline.rgfca_reserve_image_control import background_palette_counts
    rgb = np.zeros((20, 20, 3), dtype=np.uint8)
    flower = np.zeros((20, 20), dtype=bool)
    flower[:10] = True
    background = ~flower
    rgb[flower] = REFERENCE_RGB["pink"]
    rgb[background] = REFERENCE_RGB["green"]
    m = {"flower_mask": flower, "background_mask": background, "background_effective_pixels": 200}
    counts = background_palette_counts(rgb, m)
    assert counts["green"] == 200 and sum(counts.values()) == 200 and counts["pink"] == 0
    m["background_mask"] = np.ones_like(background)
    with pytest.raises(ValueError, match="includes flower"):
        background_palette_counts(rgb, m)


def test_background_receipt_mismatch_fails():
    from fcp_pipeline.rgfca_reserve_image_control import background_palette_counts
    rgb = np.zeros((10, 20, 3), dtype=np.uint8)
    m = {"flower_mask": np.zeros((10, 20), dtype=bool), "background_mask": np.ones((10, 20), dtype=bool), "background_effective_pixels": 199}
    with pytest.raises(ValueError, match="pixel count"):
        background_palette_counts(rgb, m)


def test_geometry_requires_precision_observers_dates():
    c = pd.DataFrame({"inat_taxon_id": [1, 1, 1], "species": ["A"] * 3,
                      "latitude": [0., 1., 2.], "longitude": [0., 1., 2.],
                      "positional_accuracy_m": [5., 5., 5.], "observer_id": [1, 2, 3],
                      "observed_on": ["2020-01-01", "2021-04-01", "2022-07-01"]})
    assert geometry_audit(c).iloc[0].observers == 3
    c.loc[0, "positional_accuracy_m"] = 5001
    with pytest.raises(ValueError, match="5-km"):
        geometry_audit(c)
    c.loc[0, "positional_accuracy_m"] = 5
    c["observer_id"] = 1
    with pytest.raises(ValueError, match="observer cap"):
        geometry_audit(c)


@pytest.fixture(scope="module")
def audited():
    return audit_reserve()


def test_exact_full_metadata_audit(audited):
    reserve, _, report = audited
    frozen = json.loads(AUDIT.read_text())
    assert report == {k: v for k, v in frozen.items() if k != "geometry_table_sha256"}
    assert reserve.groupby("inat_taxon_id").size().eq(100).all()
    assert report["discovery_taxon_overlap"] == 0
    assert report["colour_fields_parsed"] is False
    assert report["measurement_authorized"] is False


def test_actual_firewall_has_no_context_and_full_256_coverage(audited, tmp_path):
    reserve, _, report = audited
    contract = json.loads(CONTRACT.read_text())
    result = build_firewall(reserve, report, tmp_path, contract)
    assignments = pd.read_csv(tmp_path / "partition_assignments.csv")
    assert assignments.measurement_id.nunique() == len(assignments) == 50000
    assert len(assignments.groupby(["measurement_batch", "semantic_shard", "compute_partition"])) == 256
    for b in range(2):
        worker = pd.read_csv(tmp_path / f"batch_{b}/worker_packet/measurement_manifest.csv")
        assert list(worker.columns) == ["measurement_id", "image_filename", "photo_license"]
        assert worker.image_filename.str.fullmatch(r"FCPG-[A-F0-9]{24}\.jpg").all()
    assert result["candidate_pixels_opened"] is False
    assert result["coordinate_colour_join_opened"] is False
    empty = tmp_path / "empty-results"
    empty.mkdir()
    with pytest.raises(ValueError, match="incomplete_or_unexpected"):
        reassemble(tmp_path, empty)
    # Tampering a sealed context file must fail even before inspecting receipts.
    key = tmp_path / "sealed_keys/metadata_join_key.csv"
    key.write_text("tampered")
    with pytest.raises(ValueError, match="firewall file changed"):
        reassemble(tmp_path, empty)


def test_primary_and_robustness_are_fixed_before_pixels():
    c = json.loads(CONTRACT.read_text())
    assert c["protocol"] == PROTOCOL and c["status"] == "frozen_before_reserve_pixels"
    assert c["cohort"]["rows"] == 50000 and c["cohort"]["taxa"] == 500
    assert c["primary"]["permutations"] == 999
    assert "intersection rule" in c["required_robustness"]["gate"]
    assert c["cohort"]["prior_experiment_photo_and_observation_overlap_allowed"] is False

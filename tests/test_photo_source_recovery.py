from pathlib import Path

import pandas as pd
import pytest

from scripts.data.recover_jbi_ch1_photo_source_manifest import (
    canonicalize_source,
    choose_candidate,
    inspect_candidate,
)


def write_manifest(path: Path, *, offset: int = 0, colour: bool = False) -> None:
    rows = []
    for s in range(6):
        for i in range(200):
            row = {
                "taxon_name": f"Species {s+1}",
                "photo_id": f"p{s+1:02d}_{i+1+offset:03d}",
                "observer": f"u{i % 7}",
                "latitude": s + i / 1000,
            }
            if colour:
                row["flower_colour_state"] = "purple"
            rows.append(row)
    pd.DataFrame(rows).to_csv(path, index=False)


def test_inspect_accepts_exact_6x200_pre_measurement_csv(tmp_path):
    path = tmp_path / "candidate.csv"
    write_manifest(path)
    result = inspect_candidate(path)
    assert result["status"] == "eligible"
    assert result["species_count"] == 6
    assert set(result["per_species_counts"].values()) == {200}


def test_inspect_rejects_post_outcome_manifest(tmp_path):
    path = tmp_path / "post.csv"
    write_manifest(path, colour=True)
    result = inspect_candidate(path)
    assert result["status"] == "rejected"
    assert result["reason"] == "contains_downstream_outcome_columns"


def test_identical_duplicate_manifests_are_not_ambiguous(tmp_path):
    a = tmp_path / "a.csv"
    b = tmp_path / "nested"
    b.mkdir()
    b = b / "copy.csv"
    write_manifest(a)
    write_manifest(b)
    selected, duplicates = choose_candidate([inspect_candidate(a), inspect_candidate(b)])
    assert len(duplicates) == 2
    assert Path(selected["path"]).name == "a.csv"


def test_conflicting_eligible_manifests_fail_closed(tmp_path):
    a = tmp_path / "a.csv"
    b = tmp_path / "b.csv"
    write_manifest(a, offset=0)
    write_manifest(b, offset=1000)
    with pytest.raises(RuntimeError, match="multiple conflicting"):
        choose_candidate([inspect_candidate(a), inspect_candidate(b)])


def test_canonicalize_renames_species_and_preserves_acquisition_metadata(tmp_path):
    path = tmp_path / "candidate.csv"
    write_manifest(path)
    frame = canonicalize_source(path, "taxon_name", "photo_id")
    assert list(frame.columns[:2]) == ["species", "photo_id"]
    assert "observer" in frame.columns
    assert "latitude" in frame.columns
    assert len(frame) == 1200

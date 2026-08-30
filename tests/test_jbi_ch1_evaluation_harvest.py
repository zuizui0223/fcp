import json
from pathlib import Path

import pandas as pd
import pytest

from scripts.data.harvest_jbi_ch1_evaluation_features import (
    EXPECTED_SHARDS,
    join_coordinates,
    read_jsonl_files,
    validate_against_split,
)


SPECIES = [
    "Antirrhinum majus",
    "Dactylorhiza sambucina",
    "Gentiana lutea",
    "Ipomoea purpurea",
    "Lysimachia arvensis",
    "Raphanus sativus",
]


def make_split() -> pd.DataFrame:
    rows = []
    for species_index, species in enumerate(SPECIES):
        for i in range(200):
            split = "calibration" if i < 80 else "evaluation"
            rows.append(
                {
                    "species": species,
                    "photo_id": f"p{species_index:02d}_{i:03d}",
                    "split": split,
                    "latitude": -40.0 + species_index * 10.0 + i / 1000.0,
                    "longitude": -150.0 + species_index * 20.0 + i / 1000.0,
                }
            )
    return pd.DataFrame(rows)


def make_artifact_rows() -> list[dict]:
    rows = []
    for species_index, species in enumerate(SPECIES):
        evaluation_ids = [f"p{species_index:02d}_{i:03d}" for i in range(80, 200)]
        for local_index, photo_id in enumerate(evaluation_ids):
            rows.append(
                {
                    "species": species,
                    "photo_id": photo_id,
                    "blind_id": f"blind-{species_index:02d}-{local_index:03d}",
                    "feature_status": "ok",
                    "feature_method": "florence_open_vocab_box",
                    "evaluation_row": True,
                    "evaluation_feature_measurement": True,
                    "calibration_only": False,
                    "final_label": False,
                    "compute_shard_index": local_index // 20,
                    "compute_shard_count": 6,
                    "feature_vector": [float(species_index), float(local_index)],
                }
            )
    return rows


def write_shards(root: Path, rows: list[dict]) -> None:
    by_cell: dict[tuple[str, int], list[dict]] = {}
    for row in rows:
        key = (row["species"], int(row["compute_shard_index"]))
        by_cell.setdefault(key, []).append(row)
    assert len(by_cell) == EXPECTED_SHARDS
    for (species, shard), group in by_cell.items():
        slug = species.lower().replace(" ", "_")
        path = root / f"{slug}_{shard}.jsonl"
        path.write_text(
            "\n".join(json.dumps(row, sort_keys=True) for row in group) + "\n",
            encoding="utf-8",
        )


def test_exact_720_row_harvest_passes_and_joins_frozen_coordinates(tmp_path):
    rows = make_artifact_rows()
    write_shards(tmp_path, rows)
    parsed, files = read_jsonl_files(tmp_path)
    split = make_split()

    result = validate_against_split(parsed, split)
    joined = join_coordinates(parsed, split)

    assert len(files) == 36
    assert len(parsed) == 720
    assert result["frozen_evaluation_id_match"] is True
    assert result["calibration_overlap_count"] == 0
    assert set(result["species_counts"].values()) == {120}
    assert len(joined) == 720
    assert all("latitude" in row and "longitude" in row for row in joined)


def test_duplicate_evaluation_photo_id_fails_closed():
    rows = make_artifact_rows()
    rows[1]["photo_id"] = rows[0]["photo_id"]
    with pytest.raises(ValueError, match="duplicate species/photo IDs"):
        validate_against_split(rows, make_split())


def test_calibration_id_substitution_fails_exact_set_check():
    rows = make_artifact_rows()
    rows[0]["photo_id"] = "p00_000"
    with pytest.raises(ValueError, match="differ from frozen evaluation set"):
        validate_against_split(rows, make_split())


def test_incorrect_evaluation_flag_fails_closed():
    rows = make_artifact_rows()
    rows[0]["evaluation_feature_measurement"] = False
    with pytest.raises(ValueError, match="evaluation_feature_measurement must be true"):
        validate_against_split(rows, make_split())


def test_wrong_shard_count_fails_before_concatenation(tmp_path):
    rows = make_artifact_rows()
    write_shards(tmp_path, rows)
    next(tmp_path.glob("*.jsonl")).unlink()
    with pytest.raises(ValueError, match="expected 36 JSONL shards"):
        read_jsonl_files(tmp_path)


def test_wrong_rows_per_shard_fails_before_concatenation(tmp_path):
    rows = make_artifact_rows()
    write_shards(tmp_path, rows)
    path = next(tmp_path.glob("*.jsonl"))
    lines = path.read_text(encoding="utf-8").splitlines()
    path.write_text("\n".join(lines[:-1]) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="expected 20 rows"):
        read_jsonl_files(tmp_path)


def test_cross_species_photo_id_swap_fails_exact_species_conditioning():
    rows = make_artifact_rows()
    rows[0]["species"], rows[120]["species"] = rows[120]["species"], rows[0]["species"]
    with pytest.raises(ValueError, match="differ from frozen evaluation set"):
        validate_against_split(rows, make_split())

import pandas as pd
import pytest

from fcp_pipeline.photo_split import (
    SplitSpec,
    assignment_hash,
    canonical_id_hash,
    freeze_photo_split,
    validate_source_manifest,
)


def make_manifest() -> pd.DataFrame:
    rows = []
    for species_index in range(6):
        species = f"Species {species_index + 1}"
        for photo_index in range(200):
            rows.append(
                {
                    "species": species,
                    "photo_id": f"p{species_index + 1:02d}_{photo_index + 1:03d}",
                    "observation_id": f"o{species_index + 1:02d}_{photo_index + 1:03d}",
                    "observer": f"user_{photo_index % 13}",
                    "observed_on": f"2025-{(photo_index % 12) + 1:02d}-{(photo_index % 27) + 1:02d}",
                    "latitude": species_index * 10 + photo_index / 1000,
                    "longitude": -120 + species_index * 5 + photo_index / 1000,
                }
            )
    return pd.DataFrame(rows)


def assignment_map(frame: pd.DataFrame) -> dict[str, str]:
    return dict(zip(frame["photo_id"], frame["split"], strict=True))


def test_freeze_produces_exact_480_720_gate():
    frozen = freeze_photo_split(
        make_manifest(), species_col="species", photo_id_col="photo_id"
    )
    assert len(frozen) == 1200
    assert (frozen["split"] == "calibration").sum() == 480
    assert (frozen["split"] == "evaluation").sum() == 720
    counts = frozen.groupby(["species", "split"]).size().unstack()
    assert (counts["calibration"] == 80).all()
    assert (counts["evaluation"] == 120).all()


def test_assignment_is_row_order_invariant():
    source = make_manifest()
    a = freeze_photo_split(source, species_col="species", photo_id_col="photo_id")
    b = freeze_photo_split(
        source.sample(frac=1, random_state=991).reset_index(drop=True),
        species_col="species",
        photo_id_col="photo_id",
    )
    assert assignment_map(a) == assignment_map(b)
    assert assignment_hash(a, species_col="species", photo_id_col="photo_id") == assignment_hash(
        b, species_col="species", photo_id_col="photo_id"
    )


def test_assignment_ignores_geography_observer_and_date():
    source = make_manifest()
    a = freeze_photo_split(source, species_col="species", photo_id_col="photo_id")
    changed = source.copy()
    changed["observer"] = "completely_changed"
    changed["observed_on"] = "1900-01-01"
    changed["latitude"] = -changed["latitude"]
    changed["longitude"] = 0.0
    b = freeze_photo_split(changed, species_col="species", photo_id_col="photo_id")
    assert assignment_map(a) == assignment_map(b)


def test_different_salt_changes_assignment_but_not_counts():
    source = make_manifest()
    a = freeze_photo_split(source, species_col="species", photo_id_col="photo_id")
    b = freeze_photo_split(
        source,
        species_col="species",
        photo_id_col="photo_id",
        spec=SplitSpec(salt="different-predeclared-salt"),
    )
    assert assignment_map(a) != assignment_map(b)
    assert (b["split"] == "calibration").sum() == 480
    assert (b["split"] == "evaluation").sum() == 720


def test_rejects_duplicate_photo_ids():
    source = make_manifest()
    source.loc[1, "photo_id"] = source.loc[0, "photo_id"]
    with pytest.raises(ValueError, match="globally unique"):
        validate_source_manifest(source, species_col="species", photo_id_col="photo_id")


def test_rejects_wrong_species_sample_size():
    source = make_manifest().iloc[:-1].copy()
    with pytest.raises(ValueError, match="exactly 200"):
        validate_source_manifest(source, species_col="species", photo_id_col="photo_id")


def test_rejects_manifest_that_already_contains_colour_outcomes():
    source = make_manifest()
    source["flower_colour_state"] = "purple"
    with pytest.raises(ValueError, match="downstream measurement outcome"):
        validate_source_manifest(source, species_col="species", photo_id_col="photo_id")


def test_canonical_id_hash_ignores_row_order_and_metadata():
    source = make_manifest()
    changed = source.sample(frac=1, random_state=44).reset_index(drop=True)
    changed["observer"] = "x"
    assert canonical_id_hash(source, species_col="species", photo_id_col="photo_id") == canonical_id_hash(
        changed, species_col="species", photo_id_col="photo_id"
    )

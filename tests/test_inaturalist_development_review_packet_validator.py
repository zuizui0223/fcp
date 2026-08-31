import csv
import json
from pathlib import Path

from PIL import Image

from scripts.data.build_inaturalist_photo_development_sample import sha256
from scripts.data.validate_inaturalist_development_review_packet import validate_packet


def write_csv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def test_validator_rejects_public_local_path_before_dataset_checks(tmp_path: Path):
    artifact = tmp_path / "artifact"
    artifact.mkdir()
    private_manifest = artifact / "artifact_manifest.json"
    private_manifest.write_text(json.dumps({"private_files_sha256": {}}), encoding="utf-8")
    public = tmp_path / "public.json"
    public.write_text(
        json.dumps(
            {
                "artifact_manifest_sha256": sha256(private_manifest),
                "source_artifact": "C:/Users/example/private",
            }
        ),
        encoding="utf-8",
    )
    for name in [
        "private_photo_provenance.csv",
        "technical_image_profile.csv",
        "reviewer_A_annotation_sheet.csv",
        "reviewer_B_annotation_sheet.csv",
        "species_codebook_template.csv",
    ]:
        write_csv(artifact / name, ["placeholder"], [{"placeholder": "x"}])
    try:
        validate_packet(artifact, public)
    except RuntimeError as error:
        assert "public manifest exposes a local absolute user path" in str(error)
    else:
        raise AssertionError("validator accepted an exposed local path")

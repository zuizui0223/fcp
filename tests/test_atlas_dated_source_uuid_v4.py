from __future__ import annotations

import io
import json
from pathlib import Path
import tarfile

import pytest

from fcp_pipeline.atlas_dated_source_uuid_v4 import (
    resolve_uuid_bucket_rows,
    scan_snapshot_uuid_one_pass,
    validate_frozen_open_data_url,
    validate_live_uuid_results,
    validate_uuid_bucket_amendment,
)


AMENDMENT = Path("docs/supporting/jbi_atlas_dated_source_uuid_bucket_amendment_v4.json")


def _contract():
    return json.loads(AMENDMENT.read_text(encoding="utf-8"))


def _selected(photo_id: str = "20"):
    return {
        "cohort_id": "C01",
        "species": "Plantus alpha",
        "inat_taxon_id": "11",
        "inat_genus_id": "22",
        "observation_id": "10",
        "photo_id": photo_id,
        "photo_url_large": f"https://inaturalist-open-data.s3.amazonaws.com/photos/{photo_id}/large.jpeg",
        "photo_license": "cc-by",
        "attribution": "Frozen Attribution",
        "latitude": "1.0",
        "longitude": "2.0",
        "positional_accuracy_m": "10.0",
        "observed_on": "2024-01-02",
        "observed_month": "1",
        "local_solar_quarter": "1",
        "observer_id": "30",
        "observer": "u30",
        "primary_thinning_cell": "1,2",
        "sensitivity_thinning_cell": "1,2",
        "selection_hash": "abc",
        "candidate_image_pixels_opened": "false",
    }


def _api(photo_id: str = "20", *, license_code: str = "cc-by"):
    return {
        "id": 10,
        "uuid": "obs-uuid-10",
        "user": {"id": 30},
        # Taxon drift is intentionally allowed here: live selection state is already frozen.
        "taxon": {"id": 999},
        "photos": [
            {
                "id": int(photo_id),
                "license_code": license_code,
                "attribution": "Current Attribution",
            }
        ],
    }


def test_v4_amendment_preserves_v3_stop_and_stable_uuid_role() -> None:
    contract = _contract()
    validate_uuid_bucket_amendment(contract)
    assert contract["trigger"]["candidate_image_pixels_opened"] is False
    assert contract["trigger"]["replacement_or_resampling_permitted"] is False
    assert contract["schema_evidence"]["documented_observation_key"] == "observation_uuid"
    assert contract["snapshot_resolution"]["primary_exact_association_key"] == [
        "photo_id",
        "observation_uuid",
    ]


def test_live_uuid_enrichment_uses_identity_not_mutable_taxon_or_coordinates() -> None:
    contract = _contract()
    enriched = validate_live_uuid_results([_selected()], [_api()], contract)
    assert enriched[0]["observation_uuid"] == "obs-uuid-10"
    assert enriched[0]["live_current_taxon_id"] == "999"
    assert enriched[0]["uuid_enrichment_used_pixels"] is False

    changed = _api(license_code="cc-by-nc")
    with pytest.raises(ValueError, match="licence changed"):
        validate_live_uuid_results([_selected()], [changed], contract)

    detached = _api(photo_id="21")
    with pytest.raises(ValueError, match="no longer attached"):
        validate_live_uuid_results([_selected()], [detached], contract)


def _add(bundle: tarfile.TarFile, name: str, text: str) -> None:
    payload = text.encode("utf-8")
    info = tarfile.TarInfo(name)
    info.size = len(payload)
    bundle.addfile(info, io.BytesIO(payload))


def test_one_pass_snapshot_scan_uses_known_uuid_directly(tmp_path: Path) -> None:
    archive = tmp_path / "snapshot.tar.gz"
    with tarfile.open(archive, "w:gz") as bundle:
        _add(
            bundle,
            "observations.csv",
            "observation_uuid\tobserver_id\tlatitude\tlongitude\tpositional_accuracy\ttaxon_id\tquality_grade\tobserved_on\n"
            # Deliberately different taxon/coordinates from the frozen live row.
            "obs-uuid-10\t30\t9.0\t8.0\t99\t999\tresearch\t2024-01-02\n",
        )
        _add(bundle, "observers.csv", "observer_id\tlogin\tname\n30\tu30\tUser 30\n")
        _add(
            bundle,
            "photos.csv",
            "photo_uuid\tphoto_id\tobservation_uuid\tobserver_id\textension\tlicense\twidth\theight\tposition\n"
            "photo-uuid-20\t20\tobs-uuid-10\t30\tjpeg\tCC-BY\t1000\t800\t0\n",
        )
        _add(
            bundle,
            "taxa.csv",
            "taxon_id\tancestry\trank_level\trank\tname\tactive\n"
            "11\t1\\47125\\22\t10\tspecies\tPlantus alpha\tt\n"
            "22\t1\\47125\t20\tgenus\tPlantus\tt\n",
        )
    with archive.open("rb") as handle:
        scanned = scan_snapshot_uuid_one_pass(
            handle,
            observation_uuids={"obs-uuid-10"},
            photo_ids={"20"},
            observer_ids={"30"},
            taxon_ids={"11"},
            genus_ids={"22"},
        )
    assert scanned["observations"]["obs-uuid-10"]["taxon_id"] == "999"
    assert scanned["photos"]["20"][0]["observation_uuid"] == "obs-uuid-10"
    assert set(scanned["taxa"]) == {"11", "22"}


def test_exact_snapshot_link_wins_and_metadata_omission_uses_head_only() -> None:
    contract = _contract()
    first = validate_live_uuid_results([_selected("20")], [_api("20")], contract)[0]
    second_raw = dict(_selected("21"))
    second_raw["observation_id"] = "11"
    second_raw["observer_id"] = "31"
    second_api = {
        "id": 11,
        "uuid": "obs-uuid-11",
        "user": {"id": 31},
        "taxon": {"id": 11},
        "photos": [{"id": 21, "license_code": "cc-by", "attribution": "A"}],
    }
    second = validate_live_uuid_results([second_raw], [second_api], contract)[0]
    scanned = {
        "photos": {
            "20": [
                {
                    "photo_uuid": "p20",
                    "photo_id": "20",
                    "observation_uuid": "obs-uuid-10",
                    "observer_id": "30",
                    "extension": "jpeg",
                    "license": "CC-BY",
                    "width": "1000",
                    "height": "800",
                    "position": "0",
                }
            ]
        },
        "observations": {"obs-uuid-10": {}, "obs-uuid-11": {}},
        "observers": {"30": {}},
        "taxa": {"11": {}, "22": {}},
    }
    headed: list[str] = []

    def head(row):
        headed.append(str(row["photo_id"]))
        return {
            "status": 200,
            "content_length_bytes": 12345,
            "etag": "etag-21",
            "content_type": "image/jpeg",
        }

    audit, rows = resolve_uuid_bucket_rows([first, second], scanned, contract, bucket_head=head)
    # A two-row unit fixture cannot satisfy the immutable production denominator of 60,000,
    # but resolution classes are still fully audited before the final denominator check.
    assert audit["resolution_class_counts"] == {
        "snapshot_metadata_associated": 1,
        "snapshot_metadata_unrepresented_bucket_verified": 1,
    }
    assert audit["bucket_fallback_failure_count"] == 0
    assert headed == ["21"]
    assert rows == []


def test_bucket_fallback_rejects_non_open_data_url() -> None:
    contract = _contract()
    with pytest.raises(ValueError, match="not frozen to the Open Data bucket"):
        validate_frozen_open_data_url("https://example.org/photos/20/large.jpeg", "20", contract)

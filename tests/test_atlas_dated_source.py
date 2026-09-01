from __future__ import annotations

import io
import json
from pathlib import Path
import tarfile

import pytest

from fcp_pipeline.atlas_dated_source import (
    scan_snapshot,
    selected_tsv_rows,
    table_name,
    validate_dated_source_amendment,
)
from fcp_pipeline.atlas_dated_source_m2m import (
    scan_snapshot_m2m,
    selected_tsv_multimap,
    validate_m2m_amendment,
)


AMENDMENT_PATH = Path(
    "docs/supporting/jbi_atlas_dated_source_amendment_v1.json"
)
M2M_AMENDMENT_PATH = Path(
    "docs/supporting/jbi_atlas_dated_source_m2m_amendment_v2.json"
)


def test_dated_source_amendment_is_frozen_before_outcomes_and_pixels() -> None:
    amendment = json.loads(AMENDMENT_PATH.read_text(encoding="utf-8"))
    validate_dated_source_amendment(amendment)
    broken = json.loads(json.dumps(amendment))
    broken["authorization"]["live_feasibility_can_authorize_images"] = True
    with pytest.raises(ValueError, match="cannot independently authorize"):
        validate_dated_source_amendment(broken)

    m2m = json.loads(M2M_AMENDMENT_PATH.read_text(encoding="utf-8"))
    validate_m2m_amendment(m2m)
    m2m["trigger"]["selected_association_rows_inspected"] = True
    with pytest.raises(ValueError, match="before association rows"):
        validate_m2m_amendment(m2m)


def test_table_names_and_exact_tsv_selection_fail_closed() -> None:
    assert table_name("nested/photos.csv.gz") == "photos"
    assert table_name("observations.csv") == "observations"
    assert table_name("readme.txt") is None
    rows = selected_tsv_rows(
        io.StringIO("id\tvalue\n1\ta\n2\tb\n"),
        key="id",
        wanted={"2"},
        required_fields=("id", "value"),
    )
    assert rows == {"2": {"id": "2", "value": "b"}}
    with pytest.raises(ValueError, match="duplicates"):
        selected_tsv_rows(
            io.StringIO("id\tvalue\n2\ta\n2\tb\n"),
            key="id",
            wanted={"2"},
            required_fields=("id", "value"),
        )


def test_m2m_selector_preserves_links_but_rejects_duplicate_associations() -> None:
    text = (
        "photo_uuid\tphoto_id\tobservation_uuid\tobserver_id\n"
        "p-1\t99\tobs-1\t7\n"
        "p-1\t99\tobs-2\t7\n"
    )
    rows = selected_tsv_multimap(
        io.StringIO(text),
        key="photo_id",
        wanted={"99"},
        required_fields=("photo_uuid", "photo_id", "observation_uuid", "observer_id"),
        composite_key=("photo_uuid", "observation_uuid"),
    )
    assert [row["observation_uuid"] for row in rows["99"]] == ["obs-1", "obs-2"]
    with pytest.raises(ValueError, match="duplicates selected association key"):
        selected_tsv_multimap(
            io.StringIO(text + "p-1\t99\tobs-2\t7\n"),
            key="photo_id",
            wanted={"99"},
            required_fields=(
                "photo_uuid",
                "photo_id",
                "observation_uuid",
                "observer_id",
            ),
            composite_key=("photo_uuid", "observation_uuid"),
        )


def _add_tsv(bundle: tarfile.TarFile, name: str, text: str) -> None:
    payload = text.encode("utf-8")
    info = tarfile.TarInfo(name=name)
    info.size = len(payload)
    bundle.addfile(info, io.BytesIO(payload))


def test_snapshot_scanner_handles_observations_before_photos(tmp_path: Path) -> None:
    archive = tmp_path / "snapshot.tar.gz"
    with tarfile.open(archive, "w:gz") as bundle:
        _add_tsv(
            bundle,
            "observations.csv",
            "observation_uuid\tobserver_id\tlatitude\tlongitude\tpositional_accuracy\t"
            "taxon_id\tquality_grade\tobserved_on\tanomaly_score\n"
            "obs-1\t7\t1.0\t2.0\t10\t11\tresearch\t2020-01-02\t\n",
        )
        _add_tsv(bundle, "observers.csv", "observer_id\tlogin\tname\n7\tu7\tUser 7\n")
        _add_tsv(
            bundle,
            "photos.csv",
            "photo_uuid\tphoto_id\tobservation_uuid\tobserver_id\textension\tlicense\t"
            "width\theight\tposition\n"
            "p-1\t99\tobs-1\t7\tjpg\tCC-BY\t1000\t800\t0\n",
        )
        _add_tsv(
            bundle,
            "taxa.csv",
            "taxon_id\tancestry\trank_level\trank\tname\tactive\n"
            "11\t1/47125/22\t10\tspecies\tPlant one\tt\n"
            "22\t1/47125\t20\tgenus\tPlantus\tt\n",
        )
    result = scan_snapshot(
        archive,
        taxon_ids={"11"},
        photo_ids={"99"},
        observer_ids={"7"},
        genus_ids={"22"},
    )
    assert result["second_pass_for_observations"] is True
    assert result["photos"]["99"]["observation_uuid"] == "obs-1"
    assert result["observations"]["obs-1"]["taxon_id"] == "11"
    assert set(result["taxa"]) == {"11", "22"}


def test_m2m_snapshot_scanner_retains_all_selected_photo_links(tmp_path: Path) -> None:
    archive = tmp_path / "snapshot-m2m.tar.gz"
    with tarfile.open(archive, "w:gz") as bundle:
        _add_tsv(
            bundle,
            "observations.csv",
            "observation_uuid\tobserver_id\tlatitude\tlongitude\tpositional_accuracy\t"
            "taxon_id\tquality_grade\tobserved_on\n"
            "obs-1\t7\t1.0\t2.0\t10\t11\tresearch\t2020-01-02\n"
            "obs-2\t7\t3.0\t4.0\t10\t11\tresearch\t2020-01-03\n",
        )
        _add_tsv(bundle, "observers.csv", "observer_id\tlogin\tname\n7\tu7\tUser 7\n")
        _add_tsv(
            bundle,
            "photos.csv",
            "photo_uuid\tphoto_id\tobservation_uuid\tobserver_id\textension\tlicense\t"
            "width\theight\tposition\n"
            "p-1\t99\tobs-1\t7\tjpg\tCC-BY\t1000\t800\t0\n"
            "p-1\t99\tobs-2\t7\tjpg\tCC-BY\t1000\t800\t1\n",
        )
        _add_tsv(
            bundle,
            "taxa.csv",
            "taxon_id\tancestry\trank_level\trank\tname\tactive\n"
            "11\t1/47125/22\t10\tspecies\tPlant one\tt\n"
            "22\t1/47125\t20\tgenus\tPlantus\tt\n",
        )
    result = scan_snapshot_m2m(
        archive,
        taxon_ids={"11"},
        photo_ids={"99"},
        observer_ids={"7"},
        genus_ids={"22"},
    )
    assert len(result["photos"]["99"]) == 2
    assert set(result["observations"]) == {"obs-1", "obs-2"}

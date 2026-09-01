from __future__ import annotations

import hashlib
import io
import json
from pathlib import Path
import tarfile

import pytest

from fcp_pipeline.atlas_dated_source_m2m import (
    _association_matches,
    scan_snapshot_m2m,
)
from fcp_pipeline.atlas_dated_source_streaming import (
    HashingCountingReader,
    NONMATCH_SENTINEL,
    drain_to_eof,
    scan_snapshot_m2m_one_pass,
    validate_streaming_amendment,
)


M2M_PATH = Path("docs/supporting/jbi_atlas_dated_source_m2m_amendment_v2.json")
STREAM_PATH = Path(
    "docs/supporting/jbi_atlas_dated_source_streaming_amendment_v3.json"
)


def _contracts() -> tuple[dict, dict]:
    return (
        json.loads(M2M_PATH.read_text(encoding="utf-8")),
        json.loads(STREAM_PATH.read_text(encoding="utf-8")),
    )


def _add_tsv(bundle: tarfile.TarFile, name: str, text: str) -> None:
    payload = text.encode("utf-8")
    info = tarfile.TarInfo(name=name)
    info.size = len(payload)
    bundle.addfile(info, io.BytesIO(payload))


def _archive(observation_rows: str, photo_rows: str) -> bytes:
    sink = io.BytesIO()
    with tarfile.open(fileobj=sink, mode="w:gz") as bundle:
        _add_tsv(
            bundle,
            "observations.csv",
            "observation_uuid\tobserver_id\tlatitude\tlongitude\tpositional_accuracy\t"
            "taxon_id\tquality_grade\tobserved_on\n"
            + observation_rows,
        )
        _add_tsv(
            bundle,
            "observers.csv",
            "observer_id\tlogin\tname\n7\tu7\tUser 7\n",
        )
        _add_tsv(
            bundle,
            "photos.csv",
            "photo_uuid\tphoto_id\tobservation_uuid\tobserver_id\textension\tlicense\t"
            "width\theight\tposition\n"
            + photo_rows,
        )
        _add_tsv(
            bundle,
            "taxa.csv",
            "taxon_id\tancestry\trank_level\trank\tname\tactive\n"
            "11\t1/47125/22\t10\tspecies\tPlant one\tt\n"
            "22\t1/47125\t20\tgenus\tPlantus\tt\n",
        )
    return sink.getvalue()


def _selected() -> dict[str, str]:
    return {
        "inat_taxon_id": "11",
        "inat_genus_id": "22",
        "photo_id": "99",
        "photo_license": "cc-by",
        "observer_id": "7",
        "observed_on": "2020-01-02",
        "latitude": "1.0",
        "longitude": "2.0",
        "positional_accuracy_m": "10",
    }


def _stream_scan(payload: bytes, selected_rows: list[dict[str, str]]) -> tuple[dict, HashingCountingReader]:
    m2m, _stream = _contracts()
    reader = HashingCountingReader(io.BytesIO(payload))
    scanned = scan_snapshot_m2m_one_pass(
        reader,
        selected_rows=selected_rows,
        taxon_ids={"11"},
        photo_ids={"99"},
        observer_ids={"7"},
        genus_ids={"22"},
        m2m=m2m,
    )
    drain_to_eof(reader)
    return scanned, reader


def test_streaming_amendment_is_frozen_before_rows_and_pixels() -> None:
    m2m, streaming = _contracts()
    validate_streaming_amendment(streaming, m2m)
    broken = json.loads(json.dumps(streaming))
    broken["trigger"]["selected_snapshot_association_rows_inspected"] = True
    with pytest.raises(ValueError, match="before source rows"):
        validate_streaming_amendment(broken, m2m)


def test_hashing_reader_covers_exact_compressed_stream() -> None:
    payload = _archive(
        "obs-1\t7\t1.0\t2.0\t10\t11\tresearch\t2020-01-02\n",
        "p-1\t99\tobs-1\t7\tjpg\tCC-BY\t1000\t800\t0\n",
    )
    _scanned, reader = _stream_scan(payload, [_selected()])
    assert reader.bytes_read == len(payload)
    assert reader.hexdigest == hashlib.sha256(payload).hexdigest()


def test_one_pass_prefilter_preserves_exact_match_count_with_extra_nonmatching_link(
    tmp_path: Path,
) -> None:
    m2m, _stream = _contracts()
    payload = _archive(
        "obs-1\t7\t1.0\t2.0\t10\t11\tresearch\t2020-01-02\n"
        "obs-2\t7\t8.0\t9.0\t10\t11\tresearch\t2020-01-03\n",
        "p-1\t99\tobs-1\t7\tjpg\tCC-BY\t1000\t800\t0\n"
        "p-1\t99\tobs-2\t7\tjpg\tCC-BY\t1000\t800\t1\n",
    )
    archive = tmp_path / "snapshot.tar.gz"
    archive.write_bytes(payload)
    full = scan_snapshot_m2m(
        archive,
        taxon_ids={"11"},
        photo_ids={"99"},
        observer_ids={"7"},
        genus_ids={"22"},
    )
    one_pass, _reader = _stream_scan(payload, [_selected()])
    selected = _selected()

    full_matches = sum(
        _association_matches(selected, photo, full["observations"][photo["observation_uuid"]], m2m)
        for photo in full["photos"]["99"]
    )
    stream_matches = sum(
        _association_matches(
            selected,
            photo,
            one_pass["observations"][photo["observation_uuid"]],
            m2m,
        )
        for photo in one_pass["photos"]["99"]
    )
    assert full_matches == stream_matches == 1
    assert one_pass["prefiltered_nonmatching_linked_observations"] == 1
    assert one_pass["observations"]["obs-2"]["observer_id"] == NONMATCH_SENTINEL
    assert one_pass["second_pass_for_observations"] is False


def test_one_pass_prefilter_preserves_ambiguity_when_two_links_can_match() -> None:
    m2m, _stream = _contracts()
    payload = _archive(
        "obs-1\t7\t1.0\t2.0\t10\t11\tresearch\t2020-01-02\n"
        "obs-2\t7\t1.0\t2.0\t10\t11\tresearch\t2020-01-02\n",
        "p-1\t99\tobs-1\t7\tjpg\tCC-BY\t1000\t800\t0\n"
        "p-1\t99\tobs-2\t7\tjpg\tCC-BY\t1000\t800\t1\n",
    )
    one_pass, _reader = _stream_scan(payload, [_selected()])
    selected = _selected()
    matches = sum(
        _association_matches(
            selected,
            photo,
            one_pass["observations"][photo["observation_uuid"]],
            m2m,
        )
        for photo in one_pass["photos"]["99"]
    )
    assert matches == 2
    assert one_pass["linked_candidate_observations"] == 2
    assert one_pass["prefiltered_nonmatching_linked_observations"] == 0

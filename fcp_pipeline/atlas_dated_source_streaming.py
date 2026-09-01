"""Bounded-storage one-pass execution of the frozen dated-source M:M resolver.

The biological/source rules live in ``atlas_dated_source_m2m`` and its v2
amendment.  This module changes only I/O: observations occur before photos in
the fixed tarball, so rows that could possibly satisfy the already-frozen
observation-side exact-match conditions are retained prospectively.  Later
photo links whose observation UUID was necessarily incapable of matching are
represented by deterministic nonmatch sentinels before delegating final
resolution to the unchanged v2 resolver.
"""

from __future__ import annotations

import hashlib
import tarfile
from typing import Any, BinaryIO, Mapping, Sequence, TextIO

from .atlas_dated_source import _close, selected_tsv_rows, table_name
from .atlas_dated_source_m2m import (
    TABLES,
    _text_member,
    selected_tsv_multimap,
    validate_m2m_amendment,
)


PROTOCOL = "jbi-atlas-dated-source-streaming-execution-amendment-v3"
NONMATCH_SENTINEL = "__fcp_prefilter_nonmatch__"


class HashingCountingReader:
    """Transparent sequential reader that hashes the exact compressed bytes."""

    def __init__(self, raw: BinaryIO) -> None:
        self.raw = raw
        self._digest = hashlib.sha256()
        self.bytes_read = 0

    def read(self, size: int = -1) -> bytes:
        data = self.raw.read(size)
        if data:
            self._digest.update(data)
            self.bytes_read += len(data)
        return data

    def readinto(self, buffer: bytearray | memoryview) -> int:
        data = self.read(len(buffer))
        n = len(data)
        buffer[:n] = data
        return n

    def readable(self) -> bool:
        return True

    @property
    def hexdigest(self) -> str:
        return self._digest.hexdigest()


def drain_to_eof(reader: HashingCountingReader, block_size: int = 16 * 1024 * 1024) -> None:
    """Consume any compressed bytes buffered beyond tar end markers."""

    while reader.read(block_size):
        pass


def validate_streaming_amendment(
    streaming: Mapping[str, Any], m2m: Mapping[str, Any]
) -> None:
    """Fail closed if the pre-row bounded-storage execution freeze drifts."""

    validate_m2m_amendment(m2m)
    if streaming.get("protocol") != PROTOCOL:
        raise ValueError("unexpected dated-source streaming amendment protocol")
    trigger = streaming.get("trigger", {})
    if (
        trigger.get("candidate_image_pixels_opened") is not False
        or trigger.get("continuous_colour_used") is not False
        or trigger.get("selected_snapshot_association_rows_inspected") is not False
    ):
        raise ValueError("streaming execution was not frozen before source rows")
    source = streaming.get("exact_snapshot_stream", {})
    parent_snapshot = m2m.get("immutable_parents", {}).get("snapshot", {})
    if (
        source.get("official_url")
        != "https://inaturalist-open-data.s3.amazonaws.com/metadata/inaturalist-open-data-20260827.tar.gz"
        or source.get("content_length_bytes") != parent_snapshot.get("content_length_bytes")
        or source.get("sha256") != parent_snapshot.get("sha256")
        or source.get("stream_must_be_drained_to_eof_for_hash") is not True
        or source.get("archive_must_not_be_persisted") is not True
    ):
        raise ValueError("streaming snapshot identity or storage firewall changed")
    prefilter = streaming.get("one_pass_observation_prefilter", {})
    if (
        prefilter.get("photo_fields_are_not_used_for_prefiltering") is not True
        or "deterministic nonmatch sentinel"
        not in str(prefilter.get("linked_observation_that_fails_prefilter", ""))
    ):
        raise ValueError("one-pass observation prefilter semantics changed")
    final = streaming.get("unchanged_final_resolution", {})
    if (
        final.get("association_row_unique_key") != ["photo_uuid", "observation_uuid"]
        or final.get("all_photo_links_for_each_selected_photo_id_are_retained") is not True
        or final.get("asset_field_conflict_rule_unchanged") is not True
        or final.get("required_exact_matches_per_selected_photo") != 1
        or final.get("replacement_or_resampling_permitted") is not False
        or final.get("selected_species") != 200
        or final.get("selected_photos") != 60000
    ):
        raise ValueError("streaming execution changed final M:M resolution")
    authorization = streaming.get("authorization", {})
    if authorization.get("streaming_v3_can_authorize_images_alone") is not False:
        raise ValueError("streaming execution cannot independently authorize images")


def _observation_side_match_is_possible(
    selected: Mapping[str, Any],
    observation: Mapping[str, str],
    m2m: Mapping[str, Any],
) -> bool:
    """Necessary (not sufficient) condition for the frozen v2 exact match."""

    source = m2m["unchanged_source_rules"]
    resolution = m2m["association_resolution"]
    exact = (
        str(selected["observer_id"]) == observation["observer_id"]
        and str(selected["inat_taxon_id"]) == observation["taxon_id"]
        and str(selected["observed_on"]) == observation["observed_on"]
        and str(source["quality_grade"]).casefold()
        == observation["quality_grade"].casefold()
    )
    return exact and _close(
        selected["latitude"],
        observation["latitude"],
        float(resolution["coordinate_tolerance_degrees"]),
    ) and _close(
        selected["longitude"],
        observation["longitude"],
        float(resolution["coordinate_tolerance_degrees"]),
    ) and _close(
        selected["positional_accuracy_m"],
        observation["positional_accuracy"],
        float(resolution["positional_accuracy_tolerance_m"]),
    )


def selected_observation_candidates(
    handle: TextIO,
    *,
    selected_rows: Sequence[Mapping[str, Any]],
    required_fields: Sequence[str],
    m2m: Mapping[str, Any],
) -> dict[str, dict[str, str]]:
    """Retain every observation that could possibly be an exact v2 match.

    The index intentionally excludes photo information because photos occur
    later in the stream.  The predicate is a logical subset of the unchanged
    v2 final match: a row rejected here cannot become an exact match later.
    """

    header_line = handle.readline()
    if not header_line:
        raise ValueError("dated-source table is empty")
    header = header_line.rstrip("\r\n").split("\t")
    if len(header) != len(set(header)) or "observation_uuid" not in header:
        raise ValueError("dated-source observations header is invalid")
    missing = set(required_fields) - set(header)
    if missing:
        raise ValueError(f"dated-source table lacks fields: {sorted(missing)}")

    index: dict[tuple[str, str, str], list[Mapping[str, Any]]] = {}
    for selected in selected_rows:
        key = (
            str(selected["observer_id"]),
            str(selected["inat_taxon_id"]),
            str(selected["observed_on"]),
        )
        index.setdefault(key, []).append(selected)

    key_indices = (
        header.index("observer_id"),
        header.index("taxon_id"),
        header.index("observed_on"),
    )
    retained: dict[str, dict[str, str]] = {}
    for raw_line in handle:
        values = raw_line.rstrip("\r\n").split("\t")
        if len(values) != len(header):
            raise ValueError("dated-source TSV row width changed")
        key = tuple(values[i] for i in key_indices)
        candidates = index.get(key)
        if not candidates:
            continue
        row = dict(zip(header, values, strict=True))
        if not any(
            _observation_side_match_is_possible(selected, row, m2m)
            for selected in candidates
        ):
            continue
        uuid = row["observation_uuid"]
        if uuid in retained:
            raise ValueError(
                f"dated-source observations duplicates selected candidate UUID: {uuid}"
            )
        retained[uuid] = row
    return retained


def _nonmatch_observation(observation_uuid: str) -> dict[str, str]:
    """Represent a linked row proven incapable of satisfying the v2 match."""

    return {
        "observation_uuid": observation_uuid,
        "observer_id": NONMATCH_SENTINEL,
        "latitude": "0",
        "longitude": "0",
        "positional_accuracy": "0",
        "taxon_id": NONMATCH_SENTINEL,
        "quality_grade": NONMATCH_SENTINEL,
        "observed_on": NONMATCH_SENTINEL,
    }


def scan_snapshot_m2m_one_pass(
    fileobj: BinaryIO,
    *,
    selected_rows: Sequence[Mapping[str, Any]],
    taxon_ids: set[str],
    photo_ids: set[str],
    observer_ids: set[str],
    genus_ids: set[str],
    m2m: Mapping[str, Any],
) -> dict[str, Any]:
    """Scan the exact tar.gz once while retaining all selected photo links."""

    validate_m2m_amendment(m2m)
    photos: dict[str, list[dict[str, str]]] = {}
    observations: dict[str, dict[str, str]] = {}
    observers: dict[str, dict[str, str]] = {}
    taxa: dict[str, dict[str, str]] = {}
    members: list[str] = []
    required = {
        "photos": (
            "photo_uuid",
            "photo_id",
            "observation_uuid",
            "observer_id",
            "extension",
            "license",
            "width",
            "height",
            "position",
        ),
        "observations": (
            "observation_uuid",
            "observer_id",
            "latitude",
            "longitude",
            "positional_accuracy",
            "taxon_id",
            "quality_grade",
            "observed_on",
        ),
        "observers": ("observer_id", "login", "name"),
        "taxa": ("taxon_id", "ancestry", "rank_level", "rank", "name", "active"),
    }
    with tarfile.open(fileobj=fileobj, mode="r|gz") as bundle:
        for member in bundle:
            table = table_name(member.name)
            if table is None or not member.isfile():
                continue
            members.append(member.name)
            extracted = bundle.extractfile(member)
            if extracted is None:
                raise ValueError(f"cannot read dated-source member: {member.name}")
            with _text_member(extracted, member.name) as handle:
                if table == "observations":
                    observations = selected_observation_candidates(
                        handle,
                        selected_rows=selected_rows,
                        required_fields=required[table],
                        m2m=m2m,
                    )
                elif table == "observers":
                    observers = selected_tsv_rows(
                        handle,
                        key="observer_id",
                        wanted=observer_ids,
                        required_fields=required[table],
                    )
                elif table == "photos":
                    photos = selected_tsv_multimap(
                        handle,
                        key="photo_id",
                        wanted=photo_ids,
                        required_fields=required[table],
                        composite_key=("photo_uuid", "observation_uuid"),
                    )
                elif table == "taxa":
                    taxa = selected_tsv_rows(
                        handle,
                        key="taxon_id",
                        wanted=taxon_ids | genus_ids,
                        required_fields=required[table],
                    )

    if {table_name(name) for name in members} != set(TABLES):
        raise ValueError("dated-source archive table inventory changed")

    expected_observations = {
        row["observation_uuid"] for links in photos.values() for row in links
    }
    linked_candidates = {
        uuid: row for uuid, row in observations.items() if uuid in expected_observations
    }
    unlinked_candidate_count = len(observations) - len(linked_candidates)
    sentinel_uuids = expected_observations - set(linked_candidates)
    final_observations = dict(linked_candidates)
    final_observations.update(
        {uuid: _nonmatch_observation(uuid) for uuid in sentinel_uuids}
    )
    return {
        "members": members,
        "photos": photos,
        "observations": final_observations,
        "observers": observers,
        "taxa": taxa,
        "second_pass_for_observations": False,
        "one_pass_observation_prefilter": True,
        "linked_candidate_observations": len(linked_candidates),
        "prefiltered_nonmatching_linked_observations": len(sentinel_uuids),
        "unlinked_prefilter_candidates_discarded": unlinked_candidate_count,
    }

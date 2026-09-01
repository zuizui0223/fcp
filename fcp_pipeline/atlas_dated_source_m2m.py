"""Versioned resolver for iNaturalist's many-to-many photo-observation table."""

from __future__ import annotations

import codecs
import gzip
import math
from pathlib import Path
import tarfile
from typing import Any, BinaryIO, Mapping, Sequence, TextIO

from .atlas_dated_source import (
    _bool,
    _close,
    _selected_ids,
    selected_tsv_rows,
    table_name,
)


PROTOCOL = "jbi-atlas-dated-source-m2m-amendment-v2"
TABLES = ("observations", "observers", "photos", "taxa")
ASSET_FIELDS = (
    "photo_uuid",
    "photo_id",
    "observer_id",
    "extension",
    "license",
    "width",
    "height",
)


def validate_m2m_amendment(amendment: Mapping[str, Any]) -> None:
    """Fail closed if the post-v1, pre-pixel schema correction drifts."""

    if amendment.get("protocol") != PROTOCOL:
        raise ValueError("unexpected dated-source M:M amendment protocol")
    trigger = amendment.get("trigger", {})
    if (
        trigger.get("candidate_image_pixels_opened") is not False
        or trigger.get("continuous_colour_used") is not False
        or trigger.get("selected_association_rows_inspected") is not False
    ):
        raise ValueError("M:M resolver was not frozen before association rows")
    parents = amendment.get("immutable_parents", {})
    snapshot = parents.get("snapshot", {})
    live = parents.get("live_selection", {})
    if (
        parents.get("v1_amendment", {}).get("sha256_lf_canonical_v1")
        != "e9b81da27890c5cb8aa8afb7e712cb2da725626b1aee7adac6676ddf3844c523"
        or snapshot.get("snapshot_date") != "2026-08-27"
        or snapshot.get("content_length_bytes") != 35093052336
        or snapshot.get("sha256")
        != "c98202c07796b275fe41fc1518fc394ac09caf2dede370a4ee64ce6d68b0c50d"
        or live.get("selected_species") != 200
        or live.get("selected_photos") != 60000
        or live.get("replacement_or_resampling_permitted") is not False
    ):
        raise ValueError("M:M resolver immutable parents changed")
    resolution = amendment.get("association_resolution", {})
    if (
        resolution.get("photo_asset_lookup_key") != "photo_id"
        or resolution.get("association_row_unique_key")
        != ["photo_uuid", "observation_uuid"]
        or resolution.get("multiple_observation_links_per_photo_asset_permitted")
        is not True
        or resolution.get("asset_fields_required_constant_across_links")
        != list(ASSET_FIELDS)
        or resolution.get("coordinate_tolerance_degrees") != 1e-7
        or resolution.get("positional_accuracy_tolerance_m") != 0.0
        or resolution.get("required_exact_matches_per_selected_photo") != 1
    ):
        raise ValueError("M:M association-resolution rule changed")
    source = amendment.get("unchanged_source_rules", {})
    if (
        source.get("quality_grade") != "research"
        or source.get("taxon_rank") != "species"
        or source.get("taxon_active") is not True
        or source.get("required_ancestor_taxon_id") != 47125
    ):
        raise ValueError("M:M source rules changed")
    authorization = amendment.get("authorization", {})
    if (
        authorization.get("v1_stop_can_authorize_images") is not False
        or authorization.get("v2_reconciliation_can_authorize_images_alone")
        is not False
    ):
        raise ValueError("M:M reconciliation cannot independently authorize images")


def _text_member(extracted: BinaryIO, member_name: str) -> TextIO:
    raw: BinaryIO = extracted
    if member_name.casefold().endswith(".gz"):
        raw = gzip.GzipFile(fileobj=extracted)
    return codecs.getreader("utf-8")(raw, errors="strict")


def selected_tsv_multimap(
    handle: TextIO,
    *,
    key: str,
    wanted: set[str],
    required_fields: Sequence[str],
    composite_key: Sequence[str],
) -> dict[str, list[dict[str, str]]]:
    """Retain all association rows while enforcing the documented row key."""

    header_line = handle.readline()
    if not header_line:
        raise ValueError("dated-source table is empty")
    header = header_line.rstrip("\r\n").split("\t")
    if len(header) != len(set(header)) or key not in header:
        raise ValueError("dated-source table header is invalid")
    missing = (set(required_fields) | set(composite_key)) - set(header)
    if missing:
        raise ValueError(f"dated-source table lacks fields: {sorted(missing)}")
    key_index = header.index(key)
    selected: dict[str, list[dict[str, str]]] = {}
    seen_composites: set[tuple[str, ...]] = set()
    for raw_line in handle:
        values = raw_line.rstrip("\r\n").split("\t")
        if len(values) != len(header):
            raise ValueError("dated-source TSV row width changed")
        value = values[key_index]
        if value not in wanted:
            continue
        row = dict(zip(header, values, strict=True))
        composite = tuple(row[field] for field in composite_key)
        if composite in seen_composites:
            raise ValueError(
                "dated-source photos table duplicates selected association key: "
                + "|".join(composite)
            )
        seen_composites.add(composite)
        selected.setdefault(value, []).append(row)
    return selected


def scan_snapshot_m2m(
    archive: Path,
    *,
    taxon_ids: set[str],
    photo_ids: set[str],
    observer_ids: set[str],
    genus_ids: set[str],
) -> dict[str, Any]:
    """Stream the fixed archive while preserving all selected photo links."""

    photos: dict[str, list[dict[str, str]]] = {}
    observations: dict[str, dict[str, str]] = {}
    observers: dict[str, dict[str, str]] = {}
    taxa: dict[str, dict[str, str]] = {}
    members: list[str] = []
    observations_preceded_photos = False
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
    with tarfile.open(archive, mode="r|gz") as bundle:
        for member in bundle:
            table = table_name(member.name)
            if table is None or not member.isfile():
                continue
            members.append(member.name)
            extracted = bundle.extractfile(member)
            if extracted is None:
                raise ValueError(f"cannot read dated-source member: {member.name}")
            with _text_member(extracted, member.name) as handle:
                if table == "photos":
                    photos = selected_tsv_multimap(
                        handle,
                        key="photo_id",
                        wanted=photo_ids,
                        required_fields=required[table],
                        composite_key=("photo_uuid", "observation_uuid"),
                    )
                elif table == "observations":
                    if not photos:
                        observations_preceded_photos = True
                        continue
                    wanted = {
                        row["observation_uuid"]
                        for links in photos.values()
                        for row in links
                    }
                    observations = selected_tsv_rows(
                        handle,
                        key="observation_uuid",
                        wanted=wanted,
                        required_fields=required[table],
                    )
                elif table == "observers":
                    observers = selected_tsv_rows(
                        handle,
                        key="observer_id",
                        wanted=observer_ids,
                        required_fields=required[table],
                    )
                elif table == "taxa":
                    taxa = selected_tsv_rows(
                        handle,
                        key="taxon_id",
                        wanted=taxon_ids | genus_ids,
                        required_fields=required[table],
                    )
    if observations_preceded_photos:
        wanted = {
            row["observation_uuid"]
            for links in photos.values()
            for row in links
        }
        with tarfile.open(archive, mode="r|gz") as bundle:
            for member in bundle:
                if table_name(member.name) != "observations" or not member.isfile():
                    continue
                extracted = bundle.extractfile(member)
                if extracted is None:
                    raise ValueError("cannot read dated-source observations")
                with _text_member(extracted, member.name) as handle:
                    observations = selected_tsv_rows(
                        handle,
                        key="observation_uuid",
                        wanted=wanted,
                        required_fields=required["observations"],
                    )
                break
    if {table_name(name) for name in members} != set(TABLES):
        raise ValueError("dated-source archive table inventory changed")
    return {
        "members": members,
        "photos": photos,
        "observations": observations,
        "observers": observers,
        "taxa": taxa,
        "second_pass_for_observations": observations_preceded_photos,
    }


def _asset_conflict(links: Sequence[Mapping[str, str]]) -> list[str]:
    return [field for field in ASSET_FIELDS if len({row[field] for row in links}) != 1]


def _association_matches(
    selected: Mapping[str, Any],
    photo: Mapping[str, str],
    observation: Mapping[str, str],
    amendment: Mapping[str, Any],
) -> bool:
    source = amendment["unchanged_source_rules"]
    resolution = amendment["association_resolution"]
    exact = (
        str(selected["observer_id"]) == photo["observer_id"]
        and str(selected["observer_id"]) == observation["observer_id"]
        and str(selected["inat_taxon_id"]) == observation["taxon_id"]
        and str(selected["observed_on"]) == observation["observed_on"]
        and str(source["quality_grade"]).casefold()
        == observation["quality_grade"].casefold()
        and str(selected["photo_license"]).casefold() == photo["license"].casefold()
        and photo["license"].casefold()
        in {str(value).casefold() for value in source["allowed_photo_licenses"]}
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


def reconcile_rows_m2m(
    panels: Sequence[Mapping[str, Any]],
    selected_rows: Sequence[Mapping[str, Any]],
    scanned: Mapping[str, Any],
    amendment: Mapping[str, Any],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Resolve exactly one matching association for every frozen photo asset."""

    validate_m2m_amendment(amendment)
    taxon_ids, photo_ids, observer_ids, genus_ids = _selected_ids(panels, selected_rows)
    photos = scanned["photos"]
    observations = scanned["observations"]
    observers = scanned["observers"]
    taxa = scanned["taxa"]
    missing = {
        "photo_assets": sorted(photo_ids - set(photos)),
        "observers": sorted(observer_ids - set(observers)),
        "taxa": sorted((taxon_ids | genus_ids) - set(taxa)),
    }
    expected_observations = {
        row["observation_uuid"] for links in photos.values() for row in links
    }
    missing["linked_observations"] = sorted(expected_observations - set(observations))
    asset_conflicts = {
        photo_id: fields
        for photo_id, links in photos.items()
        if (fields := _asset_conflict(links))
    }
    association_failures: list[dict[str, Any]] = []
    association_failure_total = 0
    output: list[dict[str, Any]] = []
    source = amendment["unchanged_source_rules"]
    for raw in selected_rows:
        selected = dict(raw)
        photo_id = str(selected["photo_id"])
        links = photos.get(photo_id, [])
        if not links or photo_id in asset_conflicts:
            continue
        matches = []
        for photo in links:
            observation = observations.get(photo["observation_uuid"])
            if observation is not None and _association_matches(
                selected, photo, observation, amendment
            ):
                matches.append((photo, observation))
        if len(matches) != 1:
            association_failure_total += 1
            if len(association_failures) < 100:
                association_failures.append(
                    {
                        "photo_id": photo_id,
                        "association_rows": len(links),
                        "exact_matches": len(matches),
                    }
                )
            continue
        photo, observation = matches[0]
        observer = observers.get(str(selected["observer_id"]))
        taxon_id = str(selected["inat_taxon_id"])
        genus_id = str(selected["inat_genus_id"])
        taxon = taxa.get(taxon_id)
        genus = taxa.get(genus_id)
        if observer is None or taxon is None or genus is None:
            continue
        ancestry = set(str(taxon["ancestry"]).split("/"))
        taxon_valid = (
            taxon["rank"].casefold() == str(source["taxon_rank"]).casefold()
            and _bool(taxon["active"]) is bool(source["taxon_active"])
            and str(source["required_ancestor_taxon_id"]) in ancestry
            and genus["rank"].casefold() == "genus"
            and genus_id in ancestry
        )
        if not taxon_valid:
            association_failure_total += 1
            if len(association_failures) < 100:
                association_failures.append(
                    {
                        "photo_id": photo_id,
                        "association_rows": len(links),
                        "exact_matches": 1,
                        "failure": "taxon_or_genus_status",
                    }
                )
            continue
        display_name = observer["name"].strip() or observer["login"].strip()
        license_code = photo["license"].upper()
        attribution = (
            f"{display_name}, no rights reserved (CC0)"
            if photo["license"].casefold() == "cc0"
            else f"© {display_name}, some rights reserved ({license_code})"
        )
        extension = photo["extension"].strip().casefold()
        output.append(
            {
                **selected,
                "snapshot_date": amendment["immutable_parents"]["snapshot"][
                    "snapshot_date"
                ],
                "snapshot_photo_uuid": photo["photo_uuid"],
                "snapshot_observation_uuid": photo["observation_uuid"],
                "photo_extension": extension,
                "photo_width": photo["width"],
                "photo_height": photo["height"],
                "photo_position": photo["position"],
                "photo_url_large": (
                    "https://inaturalist-open-data.s3.amazonaws.com/"
                    f"photos/{photo_id}/large.{extension}"
                ),
                "photo_license": photo["license"].casefold(),
                "observer": observer["login"],
                "observer_name": observer["name"],
                "attribution": attribution,
                "candidate_image_pixels_opened": False,
            }
        )
    missing_counts = {name: len(values) for name, values in missing.items()}
    passed = (
        not any(missing_counts.values())
        and not asset_conflicts
        and association_failure_total == 0
        and len(output) == 60000
    )
    multi_link_counts = [len(links) for links in photos.values() if len(links) > 1]
    audit = {
        "protocol": PROTOCOL,
        "status": (
            "pass_dated_source_m2m_scaleout_freeze"
            if passed
            else "not_evaluable_dated_source_m2m_reconciliation"
        ),
        "candidate_image_pixels_opened": False,
        "continuous_colour_used": False,
        "selected_species": len(taxon_ids),
        "selected_photo_assets": len(selected_rows),
        "resolved_photo_assets": len(photos),
        "selected_photo_assets_with_multiple_observation_links": len(multi_link_counts),
        "maximum_observation_links_per_selected_photo_asset": max(
            multi_link_counts, default=1
        ),
        "selected_association_rows": sum(len(links) for links in photos.values()),
        "resolved_linked_observations": len(observations),
        "resolved_observers": len(observers),
        "resolved_taxa_and_genera": len(taxa),
        "missing_counts": missing_counts,
        "missing_examples": {name: values[:20] for name, values in missing.items()},
        "asset_field_conflict_count": len(asset_conflicts),
        "asset_field_conflict_examples": [
            {"photo_id": key, "fields": value}
            for key, value in list(asset_conflicts.items())[:20]
        ],
        "association_resolution_failure_count": association_failure_total,
        "association_resolution_failure_examples": association_failures,
        "frozen_observations": len(output) if passed else 0,
        "replacement_permitted": False,
        "image_acquisition_authorized": False,
        "claim_ceiling": amendment["claim_ceiling"],
    }
    return audit, output if passed else []

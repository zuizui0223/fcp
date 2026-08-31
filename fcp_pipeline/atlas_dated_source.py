"""Fail-closed reconciliation of a pre-colour atlas with dated iNat metadata."""

from __future__ import annotations

import codecs
import gzip
import io
import math
from pathlib import Path
import tarfile
from typing import Any, BinaryIO, Mapping, Sequence, TextIO


PROTOCOL = "jbi-atlas-dated-source-amendment-v1"
TABLES = ("observations", "observers", "photos", "taxa")


def validate_dated_source_amendment(amendment: Mapping[str, Any]) -> None:
    if amendment.get("protocol") != PROTOCOL:
        raise ValueError("unexpected dated-source amendment protocol")
    trigger = amendment.get("trigger", {})
    if (
        trigger.get("candidate_image_pixels_opened") is not False
        or trigger.get("continuous_colour_used") is not False
        or trigger.get("result_available_when_frozen") is not False
    ):
        raise ValueError("dated-source amendment was not frozen before outcomes")
    snapshot = amendment.get("snapshot", {})
    if (
        snapshot.get("snapshot_date") != "2026-08-27"
        or snapshot.get("content_length_bytes") != 35093052336
        or snapshot.get("moving_latest_pointer_permitted") is not False
        or not str(snapshot.get("object_key", "")).endswith("20260827.tar.gz")
    ):
        raise ValueError("dated snapshot identity changed")
    rules = amendment.get("exact_reconciliation", {})
    if (
        rules.get("selected_species") != 200
        or rules.get("selected_photos") != 60000
        or rules.get("coordinate_tolerance_degrees") != 1e-7
        or rules.get("positional_accuracy_tolerance_m") != 0.0
        or rules.get("quality_grade") != "research"
        or rules.get("taxon_rank") != "species"
        or rules.get("taxon_active") is not True
        or rules.get("required_ancestor_taxon_id") != 47125
    ):
        raise ValueError("dated-source reconciliation thresholds changed")
    authorization = amendment.get("authorization", {})
    if (
        authorization.get("live_feasibility_can_authorize_images") is not False
        or authorization.get("dated_source_reconciliation_can_authorize_images_alone")
        is not False
    ):
        raise ValueError("dated-source gate cannot independently authorize images")
    failure = str(amendment.get("failure_rule", ""))
    if "No replacement" not in failure or "STOP" not in failure:
        raise ValueError("dated-source failure rule changed")


def table_name(member_name: str) -> str | None:
    name = Path(member_name).name.casefold()
    if name.endswith(".gz"):
        name = name[:-3]
    if name.endswith(".csv"):
        name = name[:-4]
    return name if name in TABLES else None


def text_member(extracted: BinaryIO, member_name: str) -> TextIO:
    raw: BinaryIO = extracted
    if member_name.casefold().endswith(".gz"):
        raw = gzip.GzipFile(fileobj=extracted)
    # tarfile's streaming _FileInFile wraps a non-seekable _Stream that does not
    # implement seekable().  The incremental codec reader only needs read(), so
    # it preserves true one-pass operation for the 35 GB compressed snapshot.
    return codecs.getreader("utf-8")(raw, errors="strict")


def selected_tsv_rows(
    handle: TextIO,
    *,
    key: str,
    wanted: set[str],
    required_fields: Sequence[str],
) -> dict[str, dict[str, str]]:
    """Scan an unquoted TSV once and retain only exact requested keys."""

    header_line = handle.readline()
    if not header_line:
        raise ValueError("dated-source table is empty")
    header = header_line.rstrip("\r\n").split("\t")
    if len(header) != len(set(header)) or key not in header:
        raise ValueError("dated-source table header is invalid")
    missing = set(required_fields) - set(header)
    if missing:
        raise ValueError(f"dated-source table lacks fields: {sorted(missing)}")
    key_index = header.index(key)
    selected: dict[str, dict[str, str]] = {}
    for raw_line in handle:
        values = raw_line.rstrip("\r\n").split("\t")
        if len(values) != len(header):
            raise ValueError("dated-source TSV row width changed")
        value = values[key_index]
        if value not in wanted:
            continue
        if value in selected:
            raise ValueError(f"dated-source table duplicates selected {key}: {value}")
        selected[value] = dict(zip(header, values, strict=True))
    return selected


def _selected_ids(
    panels: Sequence[Mapping[str, Any]],
    observations: Sequence[Mapping[str, Any]],
) -> tuple[set[str], set[str], set[str], set[str]]:
    if len(panels) != 200 or len(observations) != 60000:
        raise ValueError("dated-source input denominator changed")
    taxon_ids = {str(row.get("taxon_id", "")).strip() for row in panels}
    if "" in taxon_ids or len(taxon_ids) != 200:
        raise ValueError("dated-source panels do not contain 200 unique taxa")
    photo_ids = {str(row.get("photo_id", "")).strip() for row in observations}
    observer_ids = {str(row.get("observer_id", "")).strip() for row in observations}
    genus_ids = {
        str(row.get("inat_genus_id", "")).strip() for row in observations
    }
    if "" in photo_ids or len(photo_ids) != 60000:
        raise ValueError("dated-source observations do not contain 60000 unique photos")
    if "" in observer_ids or "" in genus_ids:
        raise ValueError("dated-source observations lack observer or genus identity")
    observed_taxa = {
        str(row.get("inat_taxon_id", "")).strip() for row in observations
    }
    if observed_taxa != taxon_ids:
        raise ValueError("panel and observation taxon identities disagree")
    return taxon_ids, photo_ids, observer_ids, genus_ids


def scan_snapshot(
    archive: Path,
    *,
    taxon_ids: set[str],
    photo_ids: set[str],
    observer_ids: set[str],
    genus_ids: set[str],
) -> dict[str, Any]:
    """Stream the archive, using a second pass only when observations precede photos."""

    photos: dict[str, dict[str, str]] = {}
    observations: dict[str, dict[str, str]] = {}
    observers: dict[str, dict[str, str]] = {}
    taxa: dict[str, dict[str, str]] = {}
    members: list[str] = []
    observation_member_seen_before_photos = False

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
            with text_member(extracted, member.name) as handle:
                if table == "photos":
                    photos = selected_tsv_rows(
                        handle,
                        key="photo_id",
                        wanted=photo_ids,
                        required_fields=required[table],
                    )
                elif table == "observations":
                    if not photos:
                        observation_member_seen_before_photos = True
                        continue
                    wanted = {row["observation_uuid"] for row in photos.values()}
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

    if observation_member_seen_before_photos:
        wanted = {row["observation_uuid"] for row in photos.values()}
        with tarfile.open(archive, mode="r|gz") as bundle:
            for member in bundle:
                if table_name(member.name) != "observations" or not member.isfile():
                    continue
                extracted = bundle.extractfile(member)
                if extracted is None:
                    raise ValueError("cannot read dated-source observations")
                with text_member(extracted, member.name) as handle:
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
        "second_pass_for_observations": observation_member_seen_before_photos,
    }


def _bool(value: str) -> bool:
    folded = value.strip().casefold()
    if folded in {"t", "true", "1"}:
        return True
    if folded in {"f", "false", "0"}:
        return False
    raise ValueError(f"invalid dated-source boolean: {value!r}")


def _close(first: object, second: object, tolerance: float) -> bool:
    try:
        a, b = float(first), float(second)
    except (TypeError, ValueError):
        return False
    return math.isfinite(a) and math.isfinite(b) and abs(a - b) <= tolerance


def reconcile_rows(
    panels: Sequence[Mapping[str, Any]],
    selected_rows: Sequence[Mapping[str, Any]],
    scanned: Mapping[str, Any],
    amendment: Mapping[str, Any],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Reconcile exact selected rows; mismatches are evidence, never replacements."""

    validate_dated_source_amendment(amendment)
    taxon_ids, photo_ids, observer_ids, genus_ids = _selected_ids(panels, selected_rows)
    photos = scanned["photos"]
    observations = scanned["observations"]
    observers = scanned["observers"]
    taxa = scanned["taxa"]
    missing = {
        "photos": sorted(photo_ids - set(photos)),
        "observers": sorted(observer_ids - set(observers)),
        "taxa": sorted((taxon_ids | genus_ids) - set(taxa)),
    }
    expected_observation_uuids = {
        row["observation_uuid"] for row in photos.values()
    }
    missing["observations"] = sorted(expected_observation_uuids - set(observations))
    duplicate_observation_links = len(photos) - len(expected_observation_uuids)

    rules = amendment["exact_reconciliation"]
    tolerance = float(rules["coordinate_tolerance_degrees"])
    accuracy_tolerance = float(rules["positional_accuracy_tolerance_m"])
    allowed_licenses = {
        str(value).casefold() for value in rules["allowed_photo_licenses"]
    }
    mismatches: list[dict[str, str]] = []
    mismatch_total = 0

    def mismatch(photo_id: str, field: str, expected: object, observed: object) -> None:
        nonlocal mismatch_total
        mismatch_total += 1
        if len(mismatches) < 100:
            mismatches.append(
                {
                    "photo_id": photo_id,
                    "field": field,
                    "expected": str(expected),
                    "observed": str(observed),
                }
            )

    output: list[dict[str, Any]] = []
    for raw in selected_rows:
        row = dict(raw)
        photo_id = str(row["photo_id"])
        photo = photos.get(photo_id)
        if photo is None:
            continue
        observation = observations.get(photo["observation_uuid"])
        observer = observers.get(str(row["observer_id"]))
        taxon_id = str(row["inat_taxon_id"])
        genus_id = str(row["inat_genus_id"])
        taxon = taxa.get(taxon_id)
        genus = taxa.get(genus_id)
        if observation is None or observer is None or taxon is None or genus is None:
            continue
        comparisons = {
            "photo_observer_id": (row["observer_id"], photo["observer_id"]),
            "observation_observer_id": (row["observer_id"], observation["observer_id"]),
            "taxon_id": (taxon_id, observation["taxon_id"]),
            "observed_on": (row["observed_on"], observation["observed_on"]),
            "quality_grade": (rules["quality_grade"], observation["quality_grade"].casefold()),
            "photo_license": (str(row["photo_license"]).casefold(), photo["license"].casefold()),
        }
        for field, (expected, observed) in comparisons.items():
            if str(expected) != str(observed):
                mismatch(photo_id, field, expected, observed)
        if photo["license"].casefold() not in allowed_licenses:
            mismatch(photo_id, "allowed_photo_license", sorted(allowed_licenses), photo["license"])
        if not _close(row["latitude"], observation["latitude"], tolerance):
            mismatch(photo_id, "latitude", row["latitude"], observation["latitude"])
        if not _close(row["longitude"], observation["longitude"], tolerance):
            mismatch(photo_id, "longitude", row["longitude"], observation["longitude"])
        if not _close(
            row["positional_accuracy_m"],
            observation["positional_accuracy"],
            accuracy_tolerance,
        ):
            mismatch(
                photo_id,
                "positional_accuracy",
                row["positional_accuracy_m"],
                observation["positional_accuracy"],
            )
        ancestry = set(str(taxon["ancestry"]).split("/")) | set(
            str(taxon["ancestry"]).split("\\")
        )
        if (
            taxon["rank"].casefold() != str(rules["taxon_rank"]).casefold()
            or _bool(taxon["active"]) is not bool(rules["taxon_active"])
            or str(rules["required_ancestor_taxon_id"]) not in ancestry
        ):
            mismatch(photo_id, "taxon_status_or_ancestry", "active angiosperm species", taxon)
        if genus["rank"].casefold() != "genus" or genus_id not in ancestry:
            mismatch(photo_id, "genus_identity", genus_id, genus)

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
                **row,
                "snapshot_date": amendment["snapshot"]["snapshot_date"],
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

    missing_counts = {key: len(value) for key, value in missing.items()}
    passed = (
        not any(missing_counts.values())
        and duplicate_observation_links == 0
        and mismatch_total == 0
        and len(output) == 60000
    )
    audit = {
        "protocol": PROTOCOL,
        "status": (
            "pass_dated_source_scaleout_freeze"
            if passed
            else "not_evaluable_dated_source_reconciliation"
        ),
        "candidate_image_pixels_opened": False,
        "continuous_colour_used": False,
        "selected_species": len(taxon_ids),
        "selected_photos": len(selected_rows),
        "resolved_photos": len(photos),
        "resolved_observations": len(observations),
        "resolved_observers": len(observers),
        "resolved_taxa_and_genera": len(taxa),
        "missing_counts": missing_counts,
        "missing_examples": {key: value[:20] for key, value in missing.items()},
        "duplicate_selected_photo_links_to_same_observation": duplicate_observation_links,
        "mismatch_count": mismatch_total,
        "mismatch_examples": mismatches,
        "frozen_observations": len(output) if passed else 0,
        "replacement_permitted": False,
        "image_acquisition_authorized": False,
        "claim_ceiling": amendment["claim_ceiling"],
    }
    return audit, output if passed else []

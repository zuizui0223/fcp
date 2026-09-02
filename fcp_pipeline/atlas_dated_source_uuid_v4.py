"""Stable-UUID and bucket-aware dated-source provenance for terminal atlas v4.

This module deliberately separates two roles that the stopped v3 resolver mixed:

* live API metadata defines the already-frozen scientific selection;
* the dated Open Data snapshot and bucket establish source/association provenance.

No function in this module accepts image pixels or changes the frozen 60k selection.
"""

from __future__ import annotations

from collections import Counter
import tarfile
from typing import Any, BinaryIO, Callable, Mapping, Sequence
from urllib.parse import urlparse

from .atlas_dated_source import selected_tsv_rows, table_name
from .atlas_dated_source_m2m import TABLES, _text_member, selected_tsv_multimap


PROTOCOL = "jbi-atlas-dated-source-uuid-bucket-amendment-v4"
PASS_LABEL = "pass_dated_source_uuid_bucket_scaleout_freeze"
ALLOWED_LICENSES = {"cc0", "cc-by", "cc-by-sa", "cc-by-nc", "cc-by-nc-sa"}


def validate_uuid_bucket_amendment(amendment: Mapping[str, Any]) -> None:
    if amendment.get("protocol") != PROTOCOL:
        raise ValueError("unexpected UUID/bucket dated-source protocol")
    if amendment.get("status") != (
        "frozen_after_v3_stop_and_before_any_terminal_scaleout_candidate_pixel"
    ):
        raise ValueError("UUID/bucket amendment was not frozen before pixels")
    trigger = amendment.get("trigger", {})
    if (
        trigger.get("candidate_image_pixels_opened") is not False
        or trigger.get("continuous_colour_used") is not False
        or trigger.get("replacement_or_resampling_permitted") is not False
        or trigger.get("v3_stop_result_git_blob_sha")
        != "0742186c425cb506d16b1cda4c274534cd8d5a8d"
    ):
        raise ValueError("v4 trigger or immutable v3 STOP changed")
    evidence = amendment.get("schema_evidence", {})
    if (
        evidence.get("documented_observation_key") != "observation_uuid"
        or evidence.get("open_data_readme_git_blob_sha")
        != "b90afdfff6dd0e150774f32eba8f71bdd6cddc1d"
        or evidence.get("metadata_readme_git_blob_sha")
        != "89362939543dc4b8bec75997586038c61b7737e7"
        or evidence.get("documented_taxon_ancestry_separator") != "backslash"
    ):
        raise ValueError("official Open Data schema evidence changed")
    geometry = amendment.get("immutable_parents", {}).get(
        "terminal_geometry_artifact", {}
    )
    if (
        geometry.get("selected_species") != 200
        or geometry.get("selected_observations") != 60000
        or geometry.get("selected_photos") != 60000
        or geometry.get("artifact_id") != 9769165047
    ):
        raise ValueError("terminal 200-species denominator changed")
    snapshot = amendment.get("immutable_parents", {}).get("snapshot", {})
    if (
        snapshot.get("snapshot_date") != "2026-08-27"
        or snapshot.get("content_length_bytes") != 35093052336
        or snapshot.get("sha256")
        != "c98202c07796b275fe41fc1518fc394ac09caf2dede370a4ee64ce6d68b0c50d"
    ):
        raise ValueError("exact dated snapshot changed")
    enrich = amendment.get("live_uuid_enrichment", {})
    if (
        enrich.get("lookup_key") != "frozen integer observation_id"
        or enrich.get("required_exact_results_per_selected_observation") != 1
        or enrich.get("required_unique_observation_uuid") is not True
        or enrich.get("selected_photo_must_still_be_attached_to_returned_observation")
        is not True
        or enrich.get("selected_photo_current_license_must_equal_frozen_license")
        is not True
        or {str(x).casefold() for x in enrich.get("allowed_photo_licenses", [])}
        != ALLOWED_LICENSES
    ):
        raise ValueError("live UUID enrichment rule changed")
    resolution = amendment.get("snapshot_resolution", {})
    if (
        resolution.get("primary_exact_association_key")
        != ["photo_id", "observation_uuid"]
        or resolution.get("selected_observation_uuid_set_known_before_snapshot_scan")
        is not True
        or resolution.get("observer_table_role", "").startswith("diagnostic") is False
        or resolution.get("taxon_table_role", "").startswith("diagnostic") is False
    ):
        raise ValueError("snapshot association roles changed")
    fallback = amendment.get("metadata_omission_fallback", {})
    if (
        fallback.get("permitted_only_when_exact_snapshot_association_is_absent")
        is not True
        or fallback.get("verification_method")
        != "HTTP HEAD only; do not read image response body"
        or fallback.get("required_http_status") != 200
        or fallback.get("required_positive_content_length") is not True
    ):
        raise ValueError("bucket fallback rule changed")
    passed = amendment.get("pass_rule", {})
    if (
        passed.get("pass_label") != PASS_LABEL
        or passed.get("all_60000_rows_retained") is not True
        or passed.get("replacement_or_resampling_permitted") is not False
        or passed.get("image_acquisition_authorized_by_this_gate_alone") is not False
    ):
        raise ValueError("v4 pass rule changed")


def _result_id(row: Mapping[str, Any]) -> str:
    value = row.get("id")
    if value in (None, ""):
        raise ValueError("live API result lacks id")
    return str(value)


def _photo_by_id(observation: Mapping[str, Any], photo_id: str) -> Mapping[str, Any] | None:
    matches = [
        photo
        for photo in observation.get("photos", ()) or ()
        if str(photo.get("id", "")) == photo_id
    ]
    if len(matches) > 1:
        raise ValueError(f"live API repeats selected photo id {photo_id}")
    return matches[0] if matches else None


def validate_live_uuid_results(
    selected_rows: Sequence[Mapping[str, Any]],
    api_results: Sequence[Mapping[str, Any]],
    amendment: Mapping[str, Any],
) -> list[dict[str, Any]]:
    """Attach stable UUIDs without reselecting or re-ranking observations."""

    validate_uuid_bucket_amendment(amendment)
    selected = {str(row["observation_id"]): dict(row) for row in selected_rows}
    if len(selected) != len(selected_rows):
        raise ValueError("frozen observation IDs are not unique")
    results: dict[str, Mapping[str, Any]] = {}
    for result in api_results:
        oid = _result_id(result)
        if oid not in selected:
            raise ValueError(f"live UUID response contains unrequested observation {oid}")
        if oid in results:
            raise ValueError(f"live UUID response duplicates observation {oid}")
        results[oid] = result
    missing = set(selected) - set(results)
    if missing:
        raise ValueError(f"live UUID enrichment missing {len(missing)} frozen observations")

    output: list[dict[str, Any]] = []
    seen_uuid: set[str] = set()
    for raw in selected_rows:
        row = dict(raw)
        oid = str(row["observation_id"])
        result = results[oid]
        uuid = str(result.get("uuid") or "").strip()
        if not uuid or uuid in seen_uuid:
            raise ValueError("live observation UUID is missing or duplicated")
        seen_uuid.add(uuid)
        user = result.get("user") or {}
        if str(user.get("id") or "") != str(row["observer_id"]):
            raise ValueError(f"live observer identity changed for observation {oid}")
        photo_id = str(row["photo_id"])
        photo = _photo_by_id(result, photo_id)
        if photo is None:
            raise ValueError(f"selected photo {photo_id} is no longer attached to observation {oid}")
        current_license = str(photo.get("license_code") or "").casefold()
        frozen_license = str(row["photo_license"]).casefold()
        if current_license != frozen_license or current_license not in ALLOWED_LICENSES:
            raise ValueError(f"selected photo licence changed for photo {photo_id}")
        output.append(
            {
                **row,
                "observation_uuid": uuid,
                "live_current_taxon_id": str((result.get("taxon") or {}).get("id") or ""),
                "live_current_photo_attribution": str(photo.get("attribution") or ""),
                "uuid_enrichment_used_pixels": False,
            }
        )
    return output


def scan_snapshot_uuid_one_pass(
    fileobj: BinaryIO,
    *,
    observation_uuids: set[str],
    photo_ids: set[str],
    observer_ids: set[str],
    taxon_ids: set[str],
    genus_ids: set[str],
) -> dict[str, Any]:
    """Scan the exact tar stream once using already-known stable UUIDs."""

    observations: dict[str, dict[str, str]] = {}
    photos: dict[str, list[dict[str, str]]] = {}
    observers: dict[str, dict[str, str]] = {}
    taxa: dict[str, dict[str, str]] = {}
    members: list[str] = []
    required = {
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
                    observations = selected_tsv_rows(
                        handle,
                        key="observation_uuid",
                        wanted=observation_uuids,
                        required_fields=required[table],
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
    return {
        "members": members,
        "observations": observations,
        "photos": photos,
        "observers": observers,
        "taxa": taxa,
    }


def validate_frozen_open_data_url(url: str, photo_id: str, amendment: Mapping[str, Any]) -> None:
    fallback = amendment["metadata_omission_fallback"]
    parsed = urlparse(str(url))
    if parsed.scheme != "https" or parsed.netloc != fallback["required_frozen_url_host"]:
        raise ValueError(f"selected photo {photo_id} is not frozen to the Open Data bucket")
    prefix = f"/photos/{photo_id}/large."
    if not parsed.path.startswith(prefix) or len(parsed.path) <= len(prefix):
        raise ValueError(f"selected photo {photo_id} has an unexpected Open Data path")


def validate_bucket_head(head: Mapping[str, Any], photo_id: str, amendment: Mapping[str, Any]) -> None:
    fallback = amendment["metadata_omission_fallback"]
    if int(head.get("status", 0)) != int(fallback["required_http_status"]):
        raise ValueError(f"bucket HEAD failed for photo {photo_id}")
    if int(head.get("content_length_bytes") or 0) <= 0:
        raise ValueError(f"bucket HEAD lacks positive content length for photo {photo_id}")


def resolve_uuid_bucket_rows(
    enriched_rows: Sequence[Mapping[str, Any]],
    scanned: Mapping[str, Any],
    amendment: Mapping[str, Any],
    *,
    bucket_head: Callable[[Mapping[str, Any]], Mapping[str, Any]],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Resolve every frozen row by exact snapshot link or documented bucket fallback."""

    validate_uuid_bucket_amendment(amendment)
    photos = scanned["photos"]
    observations = scanned["observations"]
    observers = scanned["observers"]
    taxa = scanned["taxa"]
    output: list[dict[str, Any]] = []
    classes: Counter[str] = Counter()
    fallback_failures: list[dict[str, Any]] = []
    ambiguous_exact: list[dict[str, Any]] = []

    for raw in enriched_rows:
        row = dict(raw)
        photo_id = str(row["photo_id"])
        uuid = str(row["observation_uuid"])
        links = photos.get(photo_id, [])
        exact = [link for link in links if link["observation_uuid"] == uuid]
        if len(exact) > 1:
            ambiguous_exact.append(
                {"photo_id": photo_id, "observation_uuid": uuid, "exact_links": len(exact)}
            )
            continue
        if len(exact) == 1:
            photo = exact[0]
            resolution_class = "snapshot_metadata_associated"
            source_fields = {
                "snapshot_photo_uuid": photo["photo_uuid"],
                "snapshot_observation_uuid": photo["observation_uuid"],
                "snapshot_photo_extension": photo["extension"],
                "snapshot_photo_width": photo["width"],
                "snapshot_photo_height": photo["height"],
                "snapshot_photo_position": photo["position"],
                "snapshot_photo_license": photo["license"].casefold(),
                "bucket_head_status": "not_required_for_snapshot_association",
                "bucket_head_etag": "",
                "bucket_head_content_length_bytes": "",
                "bucket_head_content_type": "",
            }
        else:
            validate_frozen_open_data_url(str(row["photo_url_large"]), photo_id, amendment)
            try:
                head = dict(bucket_head(row))
                validate_bucket_head(head, photo_id, amendment)
            except (OSError, ValueError, TypeError) as exc:
                fallback_failures.append(
                    {"photo_id": photo_id, "observation_uuid": uuid, "failure": f"{type(exc).__name__}: {exc}"}
                )
                continue
            resolution_class = (
                "snapshot_metadata_unrepresented_bucket_verified"
                if not links
                else "snapshot_metadata_association_changed_bucket_verified"
            )
            source_fields = {
                "snapshot_photo_uuid": "",
                "snapshot_observation_uuid": "",
                "snapshot_photo_extension": "",
                "snapshot_photo_width": "",
                "snapshot_photo_height": "",
                "snapshot_photo_position": "",
                "snapshot_photo_license": "",
                "bucket_head_status": str(head.get("status")),
                "bucket_head_etag": str(head.get("etag") or ""),
                "bucket_head_content_length_bytes": str(head.get("content_length_bytes") or ""),
                "bucket_head_content_type": str(head.get("content_type") or ""),
            }
        classes[resolution_class] += 1
        output.append(
            {
                **row,
                **source_fields,
                "source_resolution_class": resolution_class,
                "snapshot_observation_metadata_present": uuid in observations,
                "snapshot_observer_metadata_present": str(row["observer_id"]) in observers,
                "snapshot_selected_taxon_metadata_present": str(row["inat_taxon_id"]) in taxa,
                "snapshot_selected_genus_metadata_present": str(row["inat_genus_id"]) in taxa,
                "candidate_image_pixels_opened": False,
            }
        )

    passed = (
        not ambiguous_exact
        and not fallback_failures
        and len(output) == len(enriched_rows)
        and len(output) == 60000
    )
    audit = {
        "protocol": PROTOCOL,
        "status": PASS_LABEL if passed else "not_evaluable_dated_source_uuid_bucket_reconciliation",
        "candidate_image_pixels_opened": False,
        "continuous_colour_used": False,
        "selected_species": 200,
        "selected_photo_assets": len(enriched_rows),
        "frozen_observations": len(output) if passed else 0,
        "replacement_permitted": False,
        "image_acquisition_authorized": False,
        "resolution_class_counts": dict(sorted(classes.items())),
        "snapshot_observation_metadata_present": sum(
            str(row["observation_uuid"]) in observations for row in enriched_rows
        ),
        "snapshot_observer_metadata_present": len(observers),
        "snapshot_taxa_and_genera_metadata_present": len(taxa),
        "ambiguous_exact_snapshot_association_count": len(ambiguous_exact),
        "ambiguous_exact_snapshot_association_examples": ambiguous_exact[:20],
        "bucket_fallback_failure_count": len(fallback_failures),
        "bucket_fallback_failure_examples": fallback_failures[:20],
        "claim_ceiling": amendment["claim_ceiling"],
    }
    return audit, output if passed else []


__all__ = [
    "PASS_LABEL",
    "PROTOCOL",
    "resolve_uuid_bucket_rows",
    "scan_snapshot_uuid_one_pass",
    "validate_bucket_head",
    "validate_frozen_open_data_url",
    "validate_live_uuid_results",
    "validate_uuid_bucket_amendment",
]

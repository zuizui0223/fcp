#!/usr/bin/env python3
"""Stream and reconcile the exact 2026-08-27 iNaturalist snapshot once.

No image URL is fetched by this program. The 35.09 GB metadata tarball is read
sequentially from its frozen official S3 URL, never persisted, and SHA-256 is
computed over the exact compressed bytes while the frozen 60k association rows
are resolved under the v2 many-to-many rules.
"""

from __future__ import annotations

import argparse
import csv
from datetime import timezone
from email.utils import parsedate_to_datetime
import hashlib
import json
from pathlib import Path
import sys
import tarfile
from typing import Any, Mapping, Sequence
from urllib.request import Request, urlopen


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.atlas_dated_source import _selected_ids
from fcp_pipeline.atlas_dated_source_m2m import reconcile_rows_m2m, validate_m2m_amendment
from fcp_pipeline.atlas_dated_source_streaming import (
    HashingCountingReader,
    drain_to_eof,
    scan_snapshot_m2m_one_pass,
    validate_streaming_amendment,
)


DEFAULT_M2M = ROOT / "docs/supporting/jbi_atlas_dated_source_m2m_amendment_v2.json"
DEFAULT_STREAM = (
    ROOT / "docs/supporting/jbi_atlas_dated_source_streaming_amendment_v3.json"
)
DEFAULT_RECEIPT = ROOT / "docs/supporting/jbi_atlas_inaturalist_snapshot_receipt_v1.json"


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path} must contain one JSON object")
    return value


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        return
    fields = list(rows[0])
    union = {key for row in rows for key in row}
    fields.extend(sorted(union - set(fields)))
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_blob_sha(path: Path) -> str:
    payload = path.read_bytes()
    header = f"blob {len(payload)}\0".encode("ascii")
    return hashlib.sha1(header + payload).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def normalize_etag(value: str | None) -> str:
    return (value or "").strip().strip('"')


def http_identity(url: str, *, method: str) -> tuple[Any, dict[str, Any]]:
    request = Request(
        url,
        method=method,
        headers={
            "User-Agent": "fcp-jbi-atlas-dated-source-v3/1",
            "Accept-Encoding": "identity",
        },
    )
    response = urlopen(request, timeout=120)
    content_length = response.headers.get("Content-Length")
    last_modified = response.headers.get("Last-Modified")
    identity = {
        "status": int(getattr(response, "status", response.getcode())),
        "content_length_bytes": int(content_length) if content_length else None,
        "etag": normalize_etag(response.headers.get("ETag")),
        "last_modified_utc": (
            parsedate_to_datetime(last_modified)
            .astimezone(timezone.utc)
            .strftime("%Y-%m-%dT%H:%M:%SZ")
            if last_modified
            else None
        ),
        "content_encoding": response.headers.get("Content-Encoding"),
    }
    return response, identity


def validate_http_identity(actual: Mapping[str, Any], expected: Mapping[str, Any]) -> None:
    require(
        actual.get("status") == 200,
        f"snapshot HTTP status changed: {actual.get('status')}",
    )
    require(
        actual.get("content_length_bytes") == expected.get("content_length_bytes"),
        "snapshot Content-Length changed",
    )
    require(actual.get("etag") == expected.get("etag"), "snapshot ETag changed")
    require(
        actual.get("last_modified_utc") == expected.get("last_modified_utc"),
        "snapshot Last-Modified changed",
    )
    # Content-Encoding is recorded, not used as an identity gate. S3 may store
    # gzip as object metadata even when the wire bytes are the exact tar.gz.
    # The authoritative protection is the full streamed byte count + SHA-256.


def verify_inputs(
    metadata_dir: Path,
    streaming: Mapping[str, Any],
    m2m_path: Path,
    receipt_path: Path,
) -> tuple[dict[str, Any], list[dict[str, str]], list[dict[str, str]]]:
    parents = streaming["immutable_parents"]
    require(
        git_blob_sha(m2m_path) == parents["m2m_v2"]["git_blob_sha"],
        "v2 M:M amendment Git blob changed",
    )
    require(
        git_blob_sha(receipt_path) == parents["snapshot_receipt"]["git_blob_sha"],
        "snapshot receipt Git blob changed",
    )

    geometry = parents["terminal_geometry_artifact"]
    for filename, expected_sha in geometry["files"].items():
        path = metadata_dir / filename
        require(path.exists(), f"missing terminal geometry file: {filename}")
        require(
            sha256_file(path) == expected_sha,
            f"terminal geometry SHA changed: {filename}",
        )

    manifest_path = metadata_dir / "scaleout_metadata_manifest.json"
    require(manifest_path.exists(), "missing scaleout metadata manifest")
    manifest = read_json(manifest_path)
    require(
        manifest.get("candidate_image_pixels_opened") is False,
        "terminal metadata manifest says image pixels were opened",
    )
    for filename, expected_sha in geometry["files"].items():
        require(
            manifest.get("files", {}).get(filename) == expected_sha,
            f"terminal metadata manifest SHA changed: {filename}",
        )

    feasibility = read_json(metadata_dir / "scaleout_metadata_feasibility.json")
    require(
        feasibility.get("status") == "pass_live_api_scaleout_feasibility",
        "live scaleout feasibility did not pass",
    )
    require(
        feasibility.get("candidate_image_pixels_opened") is False,
        "feasibility opened candidate pixels",
    )
    require(
        feasibility.get("continuous_colour_used") is False,
        "feasibility used colour",
    )
    require(
        int(feasibility.get("frozen_species", -1)) == 200,
        "frozen species count changed",
    )
    require(
        int(feasibility.get("frozen_observations", -1)) == 60000,
        "frozen observation count changed",
    )

    panels = read_csv(metadata_dir / "scaleout_species_panels.csv")
    observations = read_csv(metadata_dir / "scaleout_observation_manifest.csv")
    require(len(panels) == 200, "terminal panel must contain 200 species")
    require(len(observations) == 60000, "terminal manifest must contain 60,000 rows")
    require(
        len({row["photo_id"] for row in observations}) == 60000,
        "photo IDs are not unique",
    )
    require(
        len({row["observation_id"] for row in observations}) == 60000,
        "live observation IDs are not unique",
    )
    require(
        all(
            str(row.get("candidate_image_pixels_opened", "")).casefold() == "false"
            for row in observations
        ),
        "terminal observation manifest contains opened-image rows",
    )
    taxon_ids, photo_ids, observer_ids, genus_ids = _selected_ids(panels, observations)
    require(len(taxon_ids) == 200, "terminal taxon set changed")
    require(len(photo_ids) == 60000, "terminal photo set changed")
    require(len(genus_ids) == 200, "terminal genus-distinct selection changed")
    require(observer_ids, "terminal observer set is empty")
    return feasibility, panels, observations


def failure_audit(
    *,
    streaming: Mapping[str, Any],
    technical_failure: str,
    stream_state: Mapping[str, Any],
) -> dict[str, Any]:
    return {
        "protocol": streaming["protocol"],
        "status": "not_evaluable_dated_source_m2m_streaming_reconciliation",
        "candidate_image_pixels_opened": False,
        "continuous_colour_used": False,
        "selected_species": 200,
        "selected_photo_assets": 60000,
        "frozen_observations": 0,
        "technical_failure": technical_failure,
        "replacement_permitted": False,
        "image_acquisition_authorized": False,
        "snapshot_stream": dict(stream_state),
        "claim_ceiling": streaming["claim_ceiling"],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata-dir", type=Path, required=True)
    parser.add_argument("--m2m-amendment", type=Path, default=DEFAULT_M2M)
    parser.add_argument("--streaming-amendment", type=Path, default=DEFAULT_STREAM)
    parser.add_argument("--snapshot-receipt", type=Path, default=DEFAULT_RECEIPT)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    m2m = read_json(args.m2m_amendment)
    streaming = read_json(args.streaming_amendment)
    receipt = read_json(args.snapshot_receipt)
    validate_m2m_amendment(m2m)
    validate_streaming_amendment(streaming, m2m)
    expected = streaming["exact_snapshot_stream"]
    require(
        receipt.get("status") == "pass_exact_archive_identity_and_tar_integrity",
        "snapshot receipt did not pass",
    )
    require(
        receipt.get("source", {}).get("official_url") == expected["official_url"],
        "snapshot URL changed",
    )
    require(
        receipt.get("source", {}).get("sha256") == expected["sha256"],
        "snapshot receipt SHA changed",
    )

    feasibility, panels, observations = verify_inputs(
        args.metadata_dir,
        streaming,
        args.m2m_amendment,
        args.snapshot_receipt,
    )
    taxon_ids, photo_ids, observer_ids, genus_ids = _selected_ids(panels, observations)

    stream_state: dict[str, Any] = {
        "official_url": expected["official_url"],
        "archive_persisted": False,
        "candidate_image_pixels_opened": False,
    }
    frozen_rows: list[dict[str, Any]] = []
    audit: dict[str, Any]
    try:
        head_response, head_identity = http_identity(
            expected["official_url"], method="HEAD"
        )
        stream_state["head_identity"] = head_identity
        try:
            validate_http_identity(head_identity, expected)
        finally:
            head_response.close()

        get_response, get_identity = http_identity(
            expected["official_url"], method="GET"
        )
        stream_state["get_identity"] = get_identity
        reader = HashingCountingReader(get_response)
        try:
            validate_http_identity(get_identity, expected)
            scanned = scan_snapshot_m2m_one_pass(
                reader,
                selected_rows=observations,
                taxon_ids=taxon_ids,
                photo_ids=photo_ids,
                observer_ids=observer_ids,
                genus_ids=genus_ids,
                m2m=m2m,
            )
            drain_to_eof(reader)
        finally:
            get_response.close()

        stream_state.update(
            {
                "bytes_read": reader.bytes_read,
                "computed_sha256": reader.hexdigest,
                "archive_members": scanned["members"],
                "one_pass_observation_prefilter": True,
                "linked_candidate_observations": scanned[
                    "linked_candidate_observations"
                ],
                "prefiltered_nonmatching_linked_observations": scanned[
                    "prefiltered_nonmatching_linked_observations"
                ],
                "unlinked_prefilter_candidates_discarded": scanned[
                    "unlinked_prefilter_candidates_discarded"
                ],
            }
        )
        require(
            reader.bytes_read == int(expected["content_length_bytes"]),
            "stream byte count changed",
        )
        require(reader.hexdigest == expected["sha256"], "stream SHA-256 changed")

        audit, frozen_rows = reconcile_rows_m2m(panels, observations, scanned, m2m)
        audit["execution_protocol"] = streaming["protocol"]
        audit["snapshot_stream"] = stream_state
        audit["one_pass_equivalence_rule"] = streaming[
            "one_pass_observation_prefilter"
        ]["role"]
        audit["image_acquisition_authorized"] = False
        audit["claim_ceiling"] = streaming["claim_ceiling"]
    except (OSError, UnicodeError, ValueError, tarfile.TarError) as exc:
        audit = failure_audit(
            streaming=streaming,
            technical_failure=f"{type(exc).__name__}: {exc}",
            stream_state=stream_state,
        )
        frozen_rows = []

    audit["parents"] = {
        "streaming_amendment_git_blob_sha": git_blob_sha(args.streaming_amendment),
        "m2m_amendment_git_blob_sha": git_blob_sha(args.m2m_amendment),
        "snapshot_receipt_git_blob_sha": git_blob_sha(args.snapshot_receipt),
        "metadata_feasibility_sha256": sha256_file(
            args.metadata_dir / "scaleout_metadata_feasibility.json"
        ),
        "species_panels_sha256": sha256_file(
            args.metadata_dir / "scaleout_species_panels.csv"
        ),
        "selected_observations_sha256": sha256_file(
            args.metadata_dir / "scaleout_observation_manifest.csv"
        ),
    }
    audit["live_metadata_feasibility_status"] = feasibility["status"]

    args.output_dir.mkdir(parents=True, exist_ok=True)
    reconciliation_path = (
        args.output_dir / "dated_source_m2m_streaming_reconciliation.json"
    )
    write_json(reconciliation_path, audit)
    frozen_path = args.output_dir / "dated_source_m2m_observation_manifest.csv"
    if frozen_rows:
        write_csv(frozen_path, frozen_rows)

    manifest = {
        "protocol": streaming["protocol"],
        "status": audit["status"],
        "candidate_image_pixels_opened": False,
        "archive_persisted": False,
        "files": {reconciliation_path.name: sha256_file(reconciliation_path)},
    }
    if frozen_path.exists():
        manifest["files"][frozen_path.name] = sha256_file(frozen_path)
    manifest_path = args.output_dir / "dated_source_m2m_streaming_manifest.json"
    write_json(manifest_path, manifest)
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0 if audit["status"] == "pass_dated_source_m2m_scaleout_freeze" else 2


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Resolve the frozen 60k atlas with stable live UUIDs plus the exact Open Data snapshot.

The program is metadata-only.  It never GETs an image body.  Bucket fallback uses
HTTP HEAD only, and only when the exact monthly photos.csv association is absent.
"""

from __future__ import annotations

import argparse
import csv
from email.utils import parsedate_to_datetime
import hashlib
import json
from pathlib import Path
import sys
import time
from typing import Any, Mapping, Sequence
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.atlas_dated_source_streaming import HashingCountingReader, drain_to_eof
from fcp_pipeline.atlas_dated_source_uuid_v4 import (
    PASS_LABEL,
    resolve_uuid_bucket_rows,
    scan_snapshot_uuid_one_pass,
    validate_live_uuid_results,
    validate_uuid_bucket_amendment,
)


DEFAULT_AMENDMENT = ROOT / "docs/supporting/jbi_atlas_dated_source_uuid_bucket_amendment_v4.json"
USER_AGENT = "fcp-jbi-atlas-dated-source-uuid-v4/1 (metadata-only provenance)"


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
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    union = {key for row in rows for key in row}
    fields.extend(sorted(union - set(fields)))
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
    return hashlib.sha1(f"blob {len(payload)}\0".encode("ascii") + payload).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def get_json(url: str, *, attempts: int = 4, pause_seconds: float = 0.1) -> dict[str, Any]:
    error: Exception | None = None
    for attempt in range(attempts):
        try:
            request = Request(
                url,
                headers={"User-Agent": USER_AGENT, "Accept": "application/json"},
            )
            with urlopen(request, timeout=120) as response:
                payload = json.load(response)
            if pause_seconds > 0:
                time.sleep(pause_seconds)
            if not isinstance(payload, dict):
                raise ValueError("live UUID endpoint did not return a JSON object")
            return payload
        except (HTTPError, URLError, TimeoutError, ValueError) as exc:
            error = exc
            if isinstance(exc, HTTPError) and exc.code not in {429, 500, 502, 503, 504}:
                raise
            if attempt + 1 < attempts:
                time.sleep(max(1.0, pause_seconds) * (2**attempt))
    raise RuntimeError(f"live UUID request failed after {attempts} attempts: {error}")


def fetch_live_uuid_results(
    selected_rows: Sequence[Mapping[str, Any]],
    *,
    endpoint: str,
    batch_size: int,
    pause_seconds: float,
) -> list[dict[str, Any]]:
    if batch_size < 1 or batch_size > 200:
        raise ValueError("UUID API batch_size must lie in 1..200")
    observation_ids = [str(row["observation_id"]) for row in selected_rows]
    if len(observation_ids) != len(set(observation_ids)):
        raise ValueError("selected observation IDs are not unique")
    output: list[dict[str, Any]] = []
    for start in range(0, len(observation_ids), batch_size):
        batch = observation_ids[start : start + batch_size]
        params = {
            "id": ",".join(batch),
            "per_page": len(batch),
            "page": 1,
            "order_by": "id",
            "order": "asc",
            "verifiable": "any",
        }
        payload = get_json(
            f"{endpoint}?{urlencode(params)}",
            pause_seconds=pause_seconds,
        )
        results = payload.get("results") or []
        if not isinstance(results, list):
            raise ValueError("live UUID endpoint results changed type")
        output.extend(dict(row) for row in results)
        completed = min(start + len(batch), len(observation_ids))
        if completed % 5000 == 0 or completed == len(observation_ids):
            print(f"uuid_metadata_checked={completed}/{len(observation_ids)}", flush=True)
    return output


def normalize_etag(value: str | None) -> str:
    return (value or "").strip().strip('"')


def http_snapshot(url: str, *, method: str) -> tuple[Any, dict[str, Any]]:
    request = Request(
        url,
        method=method,
        headers={"User-Agent": USER_AGENT, "Accept-Encoding": "identity"},
    )
    response = urlopen(request, timeout=180)
    modified = response.headers.get("Last-Modified")
    identity = {
        "status": int(getattr(response, "status", response.getcode())),
        "content_length_bytes": (
            int(response.headers["Content-Length"])
            if response.headers.get("Content-Length")
            else None
        ),
        "etag": normalize_etag(response.headers.get("ETag")),
        "last_modified_utc": (
            parsedate_to_datetime(modified).strftime("%Y-%m-%dT%H:%M:%SZ")
            if modified
            else None
        ),
        "content_encoding": response.headers.get("Content-Encoding"),
    }
    return response, identity


def validate_snapshot_http(identity: Mapping[str, Any], snapshot: Mapping[str, Any]) -> None:
    require(identity.get("status") == 200, "snapshot HTTP status changed")
    require(
        identity.get("content_length_bytes") == snapshot.get("content_length_bytes"),
        "snapshot Content-Length changed",
    )
    require(identity.get("etag") == snapshot.get("etag"), "snapshot ETag changed")
    require(
        identity.get("last_modified_utc") == snapshot.get("last_modified_utc"),
        "snapshot Last-Modified changed",
    )


def bucket_head(row: Mapping[str, Any]) -> dict[str, Any]:
    request = Request(
        str(row["photo_url_large"]),
        method="HEAD",
        headers={"User-Agent": USER_AGENT, "Accept-Encoding": "identity"},
    )
    try:
        with urlopen(request, timeout=90) as response:
            return {
                "status": int(getattr(response, "status", response.getcode())),
                "content_length_bytes": (
                    int(response.headers["Content-Length"])
                    if response.headers.get("Content-Length")
                    else None
                ),
                "etag": normalize_etag(response.headers.get("ETag")),
                "content_type": str(response.headers.get("Content-Type") or ""),
            }
    except HTTPError as exc:
        return {
            "status": int(exc.code),
            "content_length_bytes": None,
            "etag": "",
            "content_type": "",
        }


def verify_terminal_metadata(metadata_dir: Path, amendment: Mapping[str, Any]) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    parent = amendment["immutable_parents"]["terminal_geometry_artifact"]
    for filename, expected in parent["files"].items():
        path = metadata_dir / filename
        require(path.is_file(), f"terminal metadata file missing: {filename}")
        require(sha256_file(path) == expected, f"terminal metadata SHA changed: {filename}")
    manifest_path = metadata_dir / "scaleout_metadata_manifest.json"
    require(manifest_path.is_file(), "terminal metadata manifest missing")
    manifest = read_json(manifest_path)
    require(
        manifest.get("candidate_image_pixels_opened") is False,
        "terminal metadata says candidate pixels were opened",
    )
    panels = read_csv(metadata_dir / "scaleout_species_panels.csv")
    rows = read_csv(metadata_dir / "scaleout_observation_manifest.csv")
    require(len(panels) == 200, "terminal panel count changed")
    require(len(rows) == 60000, "terminal observation denominator changed")
    require(len({row["photo_id"] for row in rows}) == 60000, "terminal photo IDs are not unique")
    require(len({row["observation_id"] for row in rows}) == 60000, "terminal observation IDs are not unique")
    require(
        all(str(row.get("candidate_image_pixels_opened", "")).casefold() == "false" for row in rows),
        "terminal metadata contains opened candidate rows",
    )
    return panels, rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata-dir", type=Path, required=True)
    parser.add_argument("--amendment", type=Path, default=DEFAULT_AMENDMENT)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--uuid-batch-size", type=int, default=200)
    parser.add_argument("--api-pause-seconds", type=float, default=0.1)
    args = parser.parse_args()

    amendment = read_json(args.amendment)
    validate_uuid_bucket_amendment(amendment)
    panels, selected_rows = verify_terminal_metadata(args.metadata_dir, amendment)
    enrich_rule = amendment["live_uuid_enrichment"]

    args.output_dir.mkdir(parents=True, exist_ok=True)
    enrichment_path = args.output_dir / "dated_source_uuid_enrichment.csv"
    try:
        api_results = fetch_live_uuid_results(
            selected_rows,
            endpoint=str(enrich_rule["endpoint"]),
            batch_size=args.uuid_batch_size,
            pause_seconds=args.api_pause_seconds,
        )
        enriched = validate_live_uuid_results(selected_rows, api_results, amendment)
    except (OSError, ValueError, RuntimeError) as exc:
        failure = {
            "protocol": amendment["protocol"],
            "status": "not_evaluable_live_uuid_enrichment",
            "candidate_image_pixels_opened": False,
            "continuous_colour_used": False,
            "selected_photo_assets": 60000,
            "technical_failure": f"{type(exc).__name__}: {exc}",
            "replacement_permitted": False,
            "image_acquisition_authorized": False,
            "claim_ceiling": amendment["claim_ceiling"],
        }
        write_json(args.output_dir / "dated_source_uuid_bucket_reconciliation.json", failure)
        return 2
    write_csv(enrichment_path, enriched)

    snapshot = amendment["immutable_parents"]["snapshot"]
    stream_state: dict[str, Any] = {
        "official_url": snapshot["official_url"],
        "archive_persisted": False,
        "candidate_image_pixels_opened": False,
    }
    frozen_rows: list[dict[str, Any]] = []
    try:
        head_response, head_identity = http_snapshot(str(snapshot["official_url"]), method="HEAD")
        try:
            validate_snapshot_http(head_identity, snapshot)
        finally:
            head_response.close()
        stream_state["head_identity"] = head_identity

        response, get_identity = http_snapshot(str(snapshot["official_url"]), method="GET")
        reader = HashingCountingReader(response)
        try:
            validate_snapshot_http(get_identity, snapshot)
            scanned = scan_snapshot_uuid_one_pass(
                reader,
                observation_uuids={str(row["observation_uuid"]) for row in enriched},
                photo_ids={str(row["photo_id"]) for row in enriched},
                observer_ids={str(row["observer_id"]) for row in enriched},
                taxon_ids={str(row["inat_taxon_id"]) for row in enriched},
                genus_ids={str(row["inat_genus_id"]) for row in enriched},
            )
            drain_to_eof(reader)
        finally:
            response.close()
        require(reader.bytes_read == snapshot["content_length_bytes"], "snapshot byte count changed")
        require(reader.hexdigest == snapshot["sha256"], "snapshot SHA-256 changed")
        stream_state.update(
            {
                "get_identity": get_identity,
                "bytes_read": reader.bytes_read,
                "computed_sha256": reader.hexdigest,
                "archive_members": scanned["members"],
            }
        )
        audit, frozen_rows = resolve_uuid_bucket_rows(
            enriched,
            scanned,
            amendment,
            bucket_head=bucket_head,
        )
        audit["snapshot_stream"] = stream_state
        audit["uuid_enrichment_rows"] = len(enriched)
        audit["uuid_enrichment_sha256"] = sha256_file(enrichment_path)
        audit["parents"] = {
            "amendment_git_blob_sha": git_blob_sha(args.amendment),
            "metadata_feasibility_sha256": sha256_file(args.metadata_dir / "scaleout_metadata_feasibility.json"),
            "species_panels_sha256": sha256_file(args.metadata_dir / "scaleout_species_panels.csv"),
            "selected_observations_sha256": sha256_file(args.metadata_dir / "scaleout_observation_manifest.csv"),
        }
    except (OSError, ValueError, RuntimeError, tarfile.TarError) as exc:  # type: ignore[name-defined]
        audit = {
            "protocol": amendment["protocol"],
            "status": "not_evaluable_dated_source_uuid_bucket_reconciliation",
            "candidate_image_pixels_opened": False,
            "continuous_colour_used": False,
            "selected_species": 200,
            "selected_photo_assets": 60000,
            "frozen_observations": 0,
            "technical_failure": f"{type(exc).__name__}: {exc}",
            "replacement_permitted": False,
            "image_acquisition_authorized": False,
            "snapshot_stream": stream_state,
            "uuid_enrichment_rows": len(enriched),
            "uuid_enrichment_sha256": sha256_file(enrichment_path),
            "claim_ceiling": amendment["claim_ceiling"],
        }
        frozen_rows = []

    reconciliation_path = args.output_dir / "dated_source_uuid_bucket_reconciliation.json"
    write_json(reconciliation_path, audit)
    frozen_path = args.output_dir / "dated_source_uuid_bucket_observation_manifest.csv"
    if frozen_rows:
        write_csv(frozen_path, frozen_rows)
    files = {
        reconciliation_path.name: sha256_file(reconciliation_path),
        enrichment_path.name: sha256_file(enrichment_path),
    }
    if frozen_path.is_file():
        files[frozen_path.name] = sha256_file(frozen_path)
    manifest = {
        "protocol": amendment["protocol"],
        "status": audit["status"],
        "candidate_image_pixels_opened": False,
        "archive_persisted": False,
        "files": files,
    }
    write_json(args.output_dir / "dated_source_uuid_bucket_manifest.json", manifest)
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0 if audit["status"] == PASS_LABEL else 2


if __name__ == "__main__":
    import tarfile
    raise SystemExit(main())

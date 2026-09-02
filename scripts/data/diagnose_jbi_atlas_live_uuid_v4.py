#!/usr/bin/env python3
"""Diagnose live UUID enrichment for the frozen terminal 60k without opening pixels."""

from __future__ import annotations

import argparse
import csv
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

USER_AGENT = "fcp-jbi-atlas-live-uuid-diagnostic-v4/1 (metadata-only)"
DEFAULT_AMENDMENT = ROOT / "docs/supporting/jbi_atlas_dated_source_uuid_bucket_amendment_v4.json"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path} must contain one JSON object")
    return value


def write_json(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def get_json(url: str, *, attempts: int = 4, pause_seconds: float = 0.1) -> dict[str, Any]:
    last: Exception | None = None
    for attempt in range(attempts):
        try:
            request = Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
            with urlopen(request, timeout=120) as response:
                payload = json.load(response)
            if pause_seconds > 0:
                time.sleep(pause_seconds)
            if not isinstance(payload, dict):
                raise ValueError("API response is not a JSON object")
            return payload
        except (HTTPError, URLError, TimeoutError, ValueError) as exc:
            last = exc
            if isinstance(exc, HTTPError) and exc.code not in {429, 500, 502, 503, 504}:
                raise
            if attempt + 1 < attempts:
                time.sleep(max(1.0, pause_seconds) * (2**attempt))
    raise RuntimeError(f"API request failed after {attempts} attempts: {last}")


def fetch_results(rows: Sequence[Mapping[str, Any]], endpoint: str, batch_size: int, pause: float) -> list[dict[str, Any]]:
    ids = [str(row["observation_id"]) for row in rows]
    if len(ids) != 60000 or len(set(ids)) != 60000:
        raise RuntimeError("diagnostic requires the exact 60,000 unique frozen observation IDs")
    output: list[dict[str, Any]] = []
    for start in range(0, len(ids), batch_size):
        batch = ids[start:start + batch_size]
        params = {
            "id": ",".join(batch),
            "per_page": len(batch),
            "page": 1,
            "order_by": "id",
            "order": "asc",
        }
        payload = get_json(f"{endpoint}?{urlencode(params)}", pause_seconds=pause)
        results = payload.get("results") or []
        if not isinstance(results, list):
            raise ValueError("API results changed type")
        output.extend(dict(row) for row in results)
        completed = min(start + len(batch), len(ids))
        if completed % 5000 == 0 or completed == len(ids):
            print(f"uuid_diagnostic_checked={completed}/{len(ids)}", flush=True)
    return output


def sample(values: Sequence[str], limit: int = 20) -> list[str]:
    return sorted(set(values))[:limit]


def diagnose(selected_rows: Sequence[Mapping[str, Any]], results: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    selected = {str(row["observation_id"]): dict(row) for row in selected_rows}
    result_lists: dict[str, list[Mapping[str, Any]]] = {}
    unexpected: list[str] = []
    for result in results:
        oid = str(result.get("id") or "")
        if oid not in selected:
            unexpected.append(oid)
            continue
        result_lists.setdefault(oid, []).append(result)

    missing_ids = sorted(set(selected) - set(result_lists))
    duplicate_ids = sorted(oid for oid, values in result_lists.items() if len(values) != 1)
    observer_mismatch: list[str] = []
    missing_uuid: list[str] = []
    duplicate_uuid_ids: list[str] = []
    photo_detached: list[str] = []
    license_changed: list[str] = []
    taxon_changed: list[str] = []
    uuid_owner: dict[str, str] = {}

    for oid, values in result_lists.items():
        if len(values) != 1:
            continue
        result = values[0]
        selected_row = selected[oid]
        user_id = str((result.get("user") or {}).get("id") or "")
        if user_id != str(selected_row["observer_id"]):
            observer_mismatch.append(oid)
        uuid = str(result.get("uuid") or "").strip()
        if not uuid:
            missing_uuid.append(oid)
        else:
            previous = uuid_owner.setdefault(uuid, oid)
            if previous != oid:
                duplicate_uuid_ids.extend([previous, oid])
        photo_id = str(selected_row["photo_id"])
        photos = [p for p in (result.get("photos") or []) if str(p.get("id") or "") == photo_id]
        if len(photos) != 1:
            photo_detached.append(oid)
        else:
            current_license = str(photos[0].get("license_code") or "").casefold()
            if current_license != str(selected_row["photo_license"]).casefold():
                license_changed.append(oid)
        current_taxon = str((result.get("taxon") or {}).get("id") or "")
        if current_taxon and current_taxon != str(selected_row["inat_taxon_id"]):
            taxon_changed.append(oid)

    hard_fail_ids = set(missing_ids) | set(duplicate_ids) | set(observer_mismatch) | set(missing_uuid) | set(duplicate_uuid_ids) | set(photo_detached) | set(license_changed)
    return {
        "protocol": "jbi-atlas-live-uuid-diagnostic-v4",
        "status": "diagnostic_complete_no_pixels",
        "candidate_image_pixels_opened": False,
        "selected_observations": len(selected_rows),
        "api_result_rows": len(results),
        "unique_api_observation_ids": len(result_lists),
        "counts": {
            "missing_api_observation_ids": len(missing_ids),
            "duplicate_api_observation_ids": len(duplicate_ids),
            "unexpected_api_observation_ids": len(unexpected),
            "observer_identity_mismatch": len(observer_mismatch),
            "missing_observation_uuid": len(missing_uuid),
            "duplicate_observation_uuid": len(set(duplicate_uuid_ids)),
            "selected_photo_no_longer_uniquely_attached": len(photo_detached),
            "selected_photo_license_changed": len(license_changed),
            "current_taxon_changed_diagnostic_only": len(taxon_changed),
            "v4_hard_failure_observations": len(hard_fail_ids),
        },
        "samples": {
            "missing_api_observation_ids": sample(missing_ids),
            "duplicate_api_observation_ids": sample(duplicate_ids),
            "observer_identity_mismatch": sample(observer_mismatch),
            "missing_observation_uuid": sample(missing_uuid),
            "duplicate_observation_uuid": sample(duplicate_uuid_ids),
            "selected_photo_no_longer_uniquely_attached": sample(photo_detached),
            "selected_photo_license_changed": sample(license_changed),
            "current_taxon_changed_diagnostic_only": sample(taxon_changed),
        },
        "replacement_or_resampling_permitted": False,
        "biological_inference": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--observation-manifest", type=Path, required=True)
    parser.add_argument("--amendment", type=Path, default=DEFAULT_AMENDMENT)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--batch-size", type=int, default=200)
    parser.add_argument("--pause-seconds", type=float, default=0.1)
    args = parser.parse_args()
    if args.batch_size < 1 or args.batch_size > 200:
        raise ValueError("batch-size must lie in 1..200")
    amendment = read_json(args.amendment)
    rows = read_csv(args.observation_manifest)
    results = fetch_results(rows, str(amendment["live_uuid_enrichment"]["endpoint"]), args.batch_size, args.pause_seconds)
    audit = diagnose(rows, results)
    write_json(args.output, audit)
    print(json.dumps(audit, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

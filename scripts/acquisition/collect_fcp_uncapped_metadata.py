#!/usr/bin/env python3
"""Checkpointed metadata-only iNaturalist observation/photo census (NO photo cap).

One invocation handles one frozen species. It has no dependency on colours,
climate, geography endpoints or observed polymorphism. It records every licensed
photo of every eligible observed flowering taxon until the stable API cursor is
exhausted. The open API is NOT suited to unthrottled world-scale scraping: use
bulk exports or approved bounded batches and obey service recommendations.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import time
from pathlib import Path
from typing import Callable
from urllib.parse import urlencode
from urllib.request import Request, urlopen

API = "https://api.inaturalist.org/v1/observations"
PAGE_SIZE = 200
MIN_PAUSE_SECONDS = 1.0
# Preserve the photo license and attribution. Only photos permitting derivative
# scientific processing appear in the candidate ledger.
ACCEPTED_LICENSES = {"cc0", "cc-by", "cc-by-sa", "cc-by-nc", "cc-by-nc-sa"}
FIELDS = (
    "inat_taxon_id", "species", "observation_id", "photo_id", "photo_url_large",
    "photo_license", "attribution", "observer_id", "observed_on", "latitude",
    "longitude", "positional_accuracy_m", "source_created_d2", "source_query"
)


def query_parameters(taxon_id: int, created_before: str, cursor: int) -> dict:
    params = {"taxon_id": taxon_id, "quality_grade": "research", "photos": "true",
              "geo": "true", "captive": "false", "created_d2": created_before,
              "per_page": PAGE_SIZE, "order_by": "id", "order": "asc"}
    if cursor:
        params["id_above"] = cursor
    return params


def network_fetch(params: dict, *, pause_seconds: float) -> dict:
    url = f"{API}?{urlencode(params)}"
    req = Request(url, headers={"Accept": "application/json", "User-Agent": "FCP uncapped flower-colour measurement planning/1.0 (zuizui0223/fcp)"})
    with urlopen(req, timeout=50) as resp:
        payload = json.load(resp)
    time.sleep(max(MIN_PAUSE_SECONDS, pause_seconds))
    return payload


def _coordinates(obs: dict) -> tuple[float | None, float | None]:
    geo = obs.get("geojson") or {}
    coords = geo.get("coordinates") or []
    if len(coords) >= 2:
        try:
            return float(coords[1]), float(coords[0])
        except (ValueError, TypeError):
            pass
    raw = obs.get("location") or ""
    try:
        lat, lon = map(float, raw.split(",", 1))
        return lat, lon
    except (ValueError, TypeError):
        return None, None


def _large_photo_url(s: str) -> str:
    # Mirror the established FCP iNat large-photo substitution, do not synthesize
    # URLs from photo IDs that were not returned by the provider.
    for variant in ("square", "small", "medium", "thumb"):
        s = s.replace(f"/{variant}.", "/large.")
    return s


def observation_rows(obs: dict, *, species: str, taxon_id: int,
                     freeze: str, query_label: str) -> list[dict]:
    taxon = obs.get("taxon") or {}
    this_taxon = taxon.get("id")
    ancestors = taxon.get("ancestor_ids") or []
    if this_taxon != taxon_id and taxon_id not in ancestors:
        return []
    try:
        obs_id = int(obs["id"])
        user = int((obs["user"] or {})["id"])
    except (KeyError, ValueError, TypeError):
        return []
    lat, lon = _coordinates(obs)
    if lat is None or lon is None or not (math.isfinite(lat) and math.isfinite(lon)):
        return []
    if not (-90 <= lat <= 90 and -180 <= lon <= 180):
        return []
    accuracy = obs.get("positional_accuracy")
    try:
        accuracy = float(accuracy) if accuracy is not None else None
        if accuracy is not None and (not math.isfinite(accuracy) or accuracy < 0):
            accuracy = None
    except (TypeError, ValueError):
        accuracy = None
    results = []
    for photo in obs.get("photos") or []:
        try:
            pid = int(photo["id"])
        except (KeyError, ValueError, TypeError):
            continue
        lic = str(photo.get("license_code") or "").lower()
        url = str(photo.get("url") or "")
        if lic not in ACCEPTED_LICENSES or not url:
            continue
        results.append(dict(zip(FIELDS, (
            taxon_id, species, obs_id, pid, _large_photo_url(url), lic,
            str(photo.get("attribution") or ""), user,
            str(obs.get("observed_on") or ""), lat, lon, accuracy, freeze, query_label
        ))))
    return results


def _write_json(path: Path, obj: dict) -> None:
    temp = path.with_suffix(path.suffix + ".partial")
    temp.write_text(json.dumps(obj, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temp, path)


def collect(taxon_id: int, species: str, created_before: str, out_dir: Path,
            request: Callable[[dict], dict], max_batches: int | None = None) -> dict:
    if taxon_id < 1 or not species.strip():
        raise ValueError("Canonical positive taxon ID and nonempty species required")
    if max_batches is not None and max_batches < 1:
        raise ValueError("max_batches must be positive")
    out_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_path = out_dir / "checkpoint.json"
    ledger_path = out_dir / "photo_manifest.ndjson"
    ledger_path.touch(exist_ok=True)  # A complete zero-photo species is valid metadata.
    initial = {"schema": "fcp_uncapped_metadata_census_v1", "taxon_id": taxon_id,
               "species": species, "created_d2": created_before,
               "cursor": 0, "batches": 0, "status": "in_progress", "page_size": PAGE_SIZE,
               "colour_pixels_opened": False, "photo_cap": None}
    checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8")) if checkpoint_path.exists() else initial
    for key, exp in (("taxon_id",taxon_id),("species",species),("created_d2",created_before)):
        if checkpoint.get(key) != exp:
            raise ValueError("Checkpoint does not describe current frozen species/date")
    if checkpoint.get("status") == "complete":
        return checkpoint
    cursor = int(checkpoint["cursor"])
    batches = int(checkpoint["batches"])
    used = 0
    while True:
        params = query_parameters(taxon_id, created_before, cursor)
        response = request(params)
        result = response.get("results")
        if not isinstance(result, list) or len(result) > PAGE_SIZE:
            raise ValueError("Bad/incomplete provider response: not an observations page")
        ids = [int(item["id"]) for item in result]
        if any(i <= cursor for i in ids) or ids != sorted(set(ids)):
            raise ValueError("Provider pagination cursor not strictly monotonic")
        rows = []
        for item in result:
            rows.extend(observation_rows(item, species=species, taxon_id=taxon_id,
                                         freeze=created_before, query_label="inat_v1_observations_id_above"))
        # Durable rows first; failed checkpoint update can replay only the most
        # recent page. Dedupe by photo ID at the later assembly gate.
        if rows:
            with ledger_path.open("a", encoding="utf-8") as handle:
                for row in rows:
                    handle.write(json.dumps(row, ensure_ascii=False, sort_keys=True) + "\n")
        if ids:
            cursor = ids[-1]
        batches += 1;used += 1
        terminal = len(result) < PAGE_SIZE
        checkpoint.update({"cursor": cursor, "batches": batches,
                           "status": "complete" if terminal else "in_progress",
                           "last_page_observations": len(result),
                           "last_page_photo_rows": len(rows),
                           "photo_cap": None, "colour_pixels_opened": False})
        _write_json(checkpoint_path, checkpoint)
        if terminal or (max_batches is not None and used >= max_batches):
            return checkpoint


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--taxon-id", type=int, required=True)
    p.add_argument("--species", required=True)
    p.add_argument("--created-before", required=True, help="Frozen YYYY-MM-DD metadata upper date")
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--pause-seconds", type=float, default=1.2)
    p.add_argument("--max-batches", type=int, default=None,
                   help="Diagnostics only: output remains incomplete, never a completed census")
    a = p.parse_args(argv)
    result = collect(a.taxon_id,a.species,a.created_before,a.out_dir,
                     lambda params: network_fetch(params,pause_seconds=a.pause_seconds),
                     max_batches=a.max_batches)
    print(json.dumps({"taxon_id":a.taxon_id,"batches":result["batches"],
                      "status":result["status"],"cursor":result["cursor"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

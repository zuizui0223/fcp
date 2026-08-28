#!/usr/bin/env python3
"""Acquire the Chapter 1 6 x 200 iNaturalist development-photo manifest.

The workflow records URLs and metadata only; it does not redistribute image binaries.
Selection is deterministic given the API snapshot defined by created_d2 and the frozen
configuration. One photo is retained per observation.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import time
from typing import Any
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import pandas as pd


USER_AGENT = "zuizui0223-fcp-jbi-ch1/1.0 (research reproducibility)"


def get_json(url: str, *, pause: float) -> dict[str, Any]:
    req = Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
    with urlopen(req, timeout=60) as response:
        payload = json.load(response)
    if pause > 0:
        time.sleep(pause)
    return payload


def exact_taxon(base_url: str, name: str, *, pause: float) -> dict[str, Any]:
    url = f"{base_url}/taxa?" + urlencode({"q": name, "rank": "species", "per_page": 30})
    payload = get_json(url, pause=pause)
    exact = [
        row
        for row in payload.get("results", [])
        if str(row.get("name", "")).casefold() == name.casefold()
        and str(row.get("rank", "")).casefold() == "species"
    ]
    if len(exact) != 1:
        raise RuntimeError(f"expected one exact species taxon for {name!r}, found {len(exact)}")
    return exact[0]


def large_photo_url(url: str | None) -> str | None:
    if not url:
        return None
    for token in ("square", "small", "medium"):
        marker = f"/{token}."
        if marker in url:
            return url.replace(marker, "/large.")
    return url


def coordinates(obs: dict[str, Any]) -> tuple[float | None, float | None]:
    geojson = obs.get("geojson") or {}
    coords = geojson.get("coordinates") or []
    if isinstance(coords, list) and len(coords) >= 2:
        try:
            return float(coords[1]), float(coords[0])
        except (TypeError, ValueError):
            pass
    raw = obs.get("location")
    if isinstance(raw, str) and "," in raw:
        try:
            lat, lon = raw.split(",", 1)
            return float(lat), float(lon)
        except ValueError:
            pass
    return None, None


def observed_month(obs: dict[str, Any]) -> int | None:
    value = str(obs.get("observed_on") or "")
    try:
        return int(value[5:7]) if len(value) >= 7 else None
    except ValueError:
        return None


def stable_hash(salt: str, *parts: object) -> str:
    payload = "\x1f".join([salt, *[str(p) for p in parts]]).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def spatial_cell(lat: float, lon: float, degrees: float) -> str:
    lat_bin = math.floor(lat / degrees) * degrees
    lon_bin = math.floor(lon / degrees) * degrees
    return f"{lat_bin:.6f},{lon_bin:.6f}"


def fetch_candidates(config: dict[str, Any], species: str) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    api = config["inat_api"]
    base_url = api["base_url"].rstrip("/")
    pause = float(api.get("request_pause_seconds", 1.0))
    taxon = exact_taxon(base_url, species, pause=pause)

    per_page = int(api.get("per_page", 200))
    max_candidates = int(api.get("maximum_candidates_per_species", 3000))
    max_pages = math.ceil(max_candidates / per_page)
    params = {
        "taxon_id": taxon["id"],
        "quality_grade": api.get("quality_grade", "research"),
        "photos": "true" if api.get("photos", True) else "false",
        "geo": "true" if api.get("geo", True) else "false",
        "captive": "true" if api.get("captive", False) else "false",
        "created_d2": api["created_d2"],
        "per_page": per_page,
        "order_by": "id",
        "order": "asc",
    }

    rows: list[dict[str, Any]] = []
    seen_obs: set[int] = set()
    seen_photo: set[int] = set()
    max_accuracy = api.get("maximum_positional_accuracy_m")

    for page in range(1, max_pages + 1):
        params["page"] = page
        url = f"{base_url}/observations?" + urlencode(params)
        payload = get_json(url, pause=pause)
        results = payload.get("results", [])
        if not results:
            break
        for obs in results:
            try:
                obs_id = int(obs["id"])
            except (KeyError, TypeError, ValueError):
                continue
            if obs_id in seen_obs:
                continue
            seen_obs.add(obs_id)

            accuracy = obs.get("positional_accuracy")
            if max_accuracy is not None and accuracy is not None:
                try:
                    if float(accuracy) > float(max_accuracy):
                        continue
                except (TypeError, ValueError):
                    pass

            lat, lon = coordinates(obs)
            if lat is None or lon is None or not (-90 <= lat <= 90) or not (-180 <= lon <= 180):
                continue

            photos = obs.get("photos") or []
            if not photos:
                continue
            photo = photos[0]
            try:
                photo_id = int(photo["id"])
            except (KeyError, TypeError, ValueError):
                continue
            if photo_id in seen_photo:
                continue
            seen_photo.add(photo_id)

            user = obs.get("user") or {}
            month = observed_month(obs)
            rows.append(
                {
                    "species": species,
                    "inat_taxon_id": int(taxon["id"]),
                    "inat_taxon_name": taxon.get("name"),
                    "observation_id": obs_id,
                    "photo_id": str(photo_id),
                    "photo_url": large_photo_url(photo.get("url")),
                    "photo_url_api": photo.get("url"),
                    "photo_license": photo.get("license_code"),
                    "attribution": photo.get("attribution"),
                    "observation_license": obs.get("license_code"),
                    "quality_grade": obs.get("quality_grade"),
                    "latitude": lat,
                    "longitude": lon,
                    "positional_accuracy_m": accuracy,
                    "geoprivacy": obs.get("geoprivacy"),
                    "obscured": bool(obs.get("obscured")),
                    "observed_on": obs.get("observed_on"),
                    "observed_month": month,
                    "time_observed_at": obs.get("time_observed_at"),
                    "created_at": obs.get("created_at"),
                    "observer_id": user.get("id"),
                    "observer": user.get("login"),
                    "place_guess": obs.get("place_guess"),
                }
            )
            if len(rows) >= max_candidates:
                break
        if len(rows) >= max_candidates:
            break
        if len(results) < per_page:
            break
    return taxon, rows


def select_rows(config: dict[str, Any], species: str, candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    selection = config["selection"]
    target = int(selection["target_photographs_per_species"])
    salt = selection["stable_hash_salt"]
    cell_degrees = float(selection["spatial_cell_degrees"])
    caps = selection["hard_caps"]
    max_observer = int(caps["maximum_per_observer"])
    max_cell = int(caps["maximum_per_spatial_cell"])
    max_month = int(caps["maximum_per_calendar_month"])

    prepared: list[dict[str, Any]] = []
    for row in candidates:
        row = dict(row)
        row["spatial_cell"] = spatial_cell(float(row["latitude"]), float(row["longitude"]), cell_degrees)
        row["selection_hash"] = stable_hash(salt, species, row["observation_id"], row["photo_id"])
        prepared.append(row)

    # Maximize spatial coverage first: iterate deterministically through cells in cycles,
    # applying observer/month caps. A second deterministic fill pass is used only when
    # the round-robin stage cannot reach 200 while obeying the same caps.
    by_cell: dict[str, list[dict[str, Any]]] = {}
    for row in prepared:
        by_cell.setdefault(row["spatial_cell"], []).append(row)
    for rows in by_cell.values():
        rows.sort(key=lambda r: (r["selection_hash"], int(r["observation_id"])))

    cell_order = sorted(by_cell, key=lambda c: stable_hash(salt, species, "cell", c))
    selected: list[dict[str, Any]] = []
    selected_ids: set[str] = set()
    observer_counts: Counter[str] = Counter()
    cell_counts: Counter[str] = Counter()
    month_counts: Counter[int | None] = Counter()

    def can_take(row: dict[str, Any]) -> bool:
        observer = str(row.get("observer") or row.get("observer_id") or "unknown")
        month = row.get("observed_month")
        cell = row["spatial_cell"]
        return (
            observer_counts[observer] < max_observer
            and cell_counts[cell] < max_cell
            and month_counts[month] < max_month
        )

    def take(row: dict[str, Any]) -> None:
        observer = str(row.get("observer") or row.get("observer_id") or "unknown")
        month = row.get("observed_month")
        cell = row["spatial_cell"]
        selected.append(row)
        selected_ids.add(str(row["photo_id"]))
        observer_counts[observer] += 1
        cell_counts[cell] += 1
        month_counts[month] += 1

    cursors = {cell: 0 for cell in cell_order}
    while len(selected) < target:
        progress = False
        for cell in cell_order:
            rows = by_cell[cell]
            while cursors[cell] < len(rows):
                row = rows[cursors[cell]]
                cursors[cell] += 1
                if str(row["photo_id"]) in selected_ids or not can_take(row):
                    continue
                take(row)
                progress = True
                break
            if len(selected) >= target:
                break
        if not progress:
            break

    if len(selected) < target:
        remaining = sorted(prepared, key=lambda r: (r["selection_hash"], int(r["observation_id"])))
        for row in remaining:
            if str(row["photo_id"]) in selected_ids or not can_take(row):
                continue
            take(row)
            if len(selected) >= target:
                break

    if len(selected) != target:
        raise RuntimeError(
            f"{species}: only {len(selected)} photographs satisfy frozen balance caps; target={target}, candidates={len(candidates)}"
        )
    return selected


def species_qc(config: dict[str, Any], species: str, taxon: dict[str, Any], candidates: list[dict[str, Any]], selected: list[dict[str, Any]]) -> dict[str, Any]:
    observers = Counter(str(r.get("observer") or r.get("observer_id") or "unknown") for r in selected)
    cells = Counter(str(r["spatial_cell"]) for r in selected)
    months = Counter(r.get("observed_month") for r in selected if r.get("observed_month") is not None)
    n = len(selected)
    qc = {
        "species": species,
        "inat_taxon_id": int(taxon["id"]),
        "candidate_count": len(candidates),
        "selected_count": n,
        "unique_observers": len(observers),
        "unique_spatial_cells": len(cells),
        "unique_calendar_months": len(months),
        "maximum_observer_fraction": max(observers.values()) / n if n else None,
        "maximum_spatial_cell_fraction": max(cells.values()) / n if n else None,
        "maximum_month_fraction": max(months.values()) / n if months and n else None,
        "observer_counts_top10": observers.most_common(10),
        "month_counts": sorted((str(k), int(v)) for k, v in months.items()),
    }
    gates = config["selection"]["qc_gates"]
    failures: list[str] = []
    if qc["unique_observers"] < int(gates["minimum_unique_observers"]):
        failures.append("minimum_unique_observers")
    if qc["unique_spatial_cells"] < int(gates["minimum_unique_spatial_cells"]):
        failures.append("minimum_unique_spatial_cells")
    if qc["unique_calendar_months"] < int(gates["minimum_unique_calendar_months"]):
        failures.append("minimum_unique_calendar_months")
    if qc["maximum_observer_fraction"] > float(gates["maximum_observer_fraction"]):
        failures.append("maximum_observer_fraction")
    qc["gate_failures"] = failures
    qc["gate_pass"] = not failures
    return qc


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("docs/supporting/jbi_ch1_inat_acquisition_v1.json"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/frozen/jbi_ch1_photo_source_manifest.csv"),
    )
    parser.add_argument(
        "--qc",
        type=Path,
        default=Path("docs/supporting/jbi_ch1_inat_acquisition_qc_v1.json"),
    )
    args = parser.parse_args()

    config = json.loads(args.config.read_text(encoding="utf-8"))
    species_list = list(config["target_selection"]["species"])
    if len(species_list) != 6 or len(set(species_list)) != 6:
        raise RuntimeError("acquisition config must contain exactly six unique target species")

    all_selected: list[dict[str, Any]] = []
    qc_rows: list[dict[str, Any]] = []
    for species in species_list:
        print(f"Acquiring {species}...", flush=True)
        taxon, candidates = fetch_candidates(config, species)
        selected = select_rows(config, species, candidates)
        qc = species_qc(config, species, taxon, candidates, selected)
        qc_rows.append(qc)
        if not qc["gate_pass"]:
            raise RuntimeError(f"{species}: acquisition QC gate failed: {qc['gate_failures']}")
        all_selected.extend(selected)
        print(
            f"{species}: candidates={len(candidates)} selected={len(selected)} observers={qc['unique_observers']} cells={qc['unique_spatial_cells']} months={qc['unique_calendar_months']}",
            flush=True,
        )

    frame = pd.DataFrame(all_selected)
    if len(frame) != 1200:
        raise RuntimeError(f"expected exactly 1200 selected photographs, found {len(frame)}")
    if frame["photo_id"].astype(str).duplicated().any():
        dup = frame.loc[frame["photo_id"].astype(str).duplicated(keep=False), "photo_id"].iloc[0]
        raise RuntimeError(f"photo ID duplicated across final manifest: {dup}")

    counts = frame.groupby("species").size()
    if len(counts) != 6 or not (counts == 200).all():
        raise RuntimeError(f"final per-species counts invalid: {counts.to_dict()}")

    frame = frame.sort_values(["species", "selection_hash", "photo_id"], kind="mergesort").reset_index(drop=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.output, index=False, lineterminator="\n")

    qc_payload = {
        "protocol": config["protocol"],
        "status": "pass",
        "config_sha256": hashlib.sha256(args.config.read_bytes()).hexdigest(),
        "source_manifest_sha256": hashlib.sha256(args.output.read_bytes()).hexdigest(),
        "target_selection": config["target_selection"],
        "inat_api": config["inat_api"],
        "selection": config["selection"],
        "species_qc": qc_rows,
        "total_selected": int(len(frame)),
        "all_species_gate_pass": all(row["gate_pass"] for row in qc_rows),
    }
    args.qc.parent.mkdir(parents=True, exist_ok=True)
    args.qc.write_text(json.dumps(qc_payload, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(qc_payload, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

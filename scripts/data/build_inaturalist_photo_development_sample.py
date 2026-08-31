#!/usr/bin/env python3
"""Build a morph-blind iNaturalist API development sample.

This is deliberately not the frozen inferential pilot.  The scientific pilot
must be reconstructed from a dated iNaturalist Open Dataset snapshot.  This
script uses the moving public API only to exercise the fixed sampling,
provenance, image-availability, and annotation workflow before that large
snapshot is acquired.

Exact coordinates, observer identifiers, image URLs, and attribution remain in
the local artifact directory.  The repository receives only an aggregate,
claim-limited manifest.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import mimetypes
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Any, Iterable

from PIL import Image
from pyproj import Transformer


API_ROOT = "https://api.inaturalist.org/v1"
CSV_ROOT = "https://www.inaturalist.org/observations.csv"
USER_AGENT = "FCP-morph-blind-photo-development/1.0 (research workflow audit)"
REUSABLE_LICENSES = {"cc0", "cc-by"}
GRID_M = 50_000
DEFAULT_SEED = "fcp-inaturalist-photo-development-v1"
SPECIES = {
    "Digitalis purpurea": 53983,
    "Erythranthe lewisii": 777190,
    "Hepatica nobilis": 639660,
    "Hesperis matronalis": 47697,
    "Orchis mascula": 132621,
    "Protea repens": 355849,
}
TRANSFORMER = Transformer.from_crs("EPSG:4326", "EPSG:6933", always_xy=True)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_rank(seed: str, *values: object) -> str:
    token = "\x1f".join([seed, *(str(value) for value in values)])
    return hashlib.sha256(token.encode("utf-8")).hexdigest()


def request_bytes(url: str, retries: int = 4) -> tuple[bytes, dict[str, str]]:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                headers = {key.lower(): value for key, value in response.headers.items()}
                return response.read(), headers
        except Exception:
            if attempt + 1 == retries:
                raise
            time.sleep(2**attempt)
    raise AssertionError("unreachable")


def api_params(taxon_id: int) -> dict[str, Any]:
    return {
        "taxon_id": taxon_id,
        "quality_grade": "research",
        "photos": "true",
        "captive": "false",
        "geo": "true",
        "geoprivacy": "open",
        "taxon_geoprivacy": "open",
        "acc_below": 1000,
        "d1": "2016-01-01",
        "d2": "2025-12-31",
        "photo_license": "cc0,cc-by",
        "order_by": "id",
        "order": "asc",
    }


def parse_csv(payload: bytes) -> list[dict[str, str]]:
    text = payload.decode("utf-8-sig")
    return list(csv.DictReader(text.splitlines()))


def fetch_species_rows(
    species: str,
    taxon_id: int,
    raw_dir: Path,
    request_interval: float,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    page = 1
    while True:
        params = api_params(taxon_id) | {"per_page": 200, "page": page}
        url = f"{CSV_ROOT}?{urllib.parse.urlencode(params)}"
        page_path = raw_dir / f"{taxon_id}_page_{page:03d}.csv"
        headers_path = raw_dir / f"{taxon_id}_page_{page:03d}.headers.json"
        if page_path.exists() and headers_path.exists():
            payload = page_path.read_bytes()
            source = "reused"
        else:
            payload, headers = request_bytes(url)
            page_path.write_bytes(payload)
            headers_path.write_text(
                json.dumps(headers, indent=2, sort_keys=True) + "\n", encoding="utf-8"
            )
            source = "downloaded"
        page_rows = parse_csv(payload)
        rows.extend(page_rows)
        print(
            f"{species}: metadata page {page}, n={len(page_rows)}, source={source}",
            flush=True,
        )
        if len(page_rows) < 200:
            break
        if page >= 50:
            raise RuntimeError(f"{species}: API pagination reached the 10,000-row ceiling")
        page += 1
        if source == "downloaded":
            time.sleep(request_interval)
    return rows


def spatial_fields(row: dict[str, str]) -> dict[str, Any]:
    lon = float(row["longitude"])
    lat = float(row["latitude"])
    x, y = TRANSFORMER.transform(lon, lat)
    observed = date.fromisoformat(row["observed_on"])
    week = observed.isocalendar().week
    return {
        "longitude": lon,
        "latitude": lat,
        "grid_x": math.floor(x / GRID_M),
        "grid_y": math.floor(y / GRID_M),
        "flowering_week": week,
    }


def validate_row(row: dict[str, str], species: str, taxon_id: int) -> dict[str, Any]:
    if row.get("scientific_name") != species:
        raise ValueError(f"{species}: unexpected scientific_name={row.get('scientific_name')}")
    if int(row["taxon_id"]) != taxon_id:
        raise ValueError(f"{species}: unexpected taxon_id={row.get('taxon_id')}")
    if row.get("quality_grade", "").casefold() != "research":
        raise ValueError(f"{species}: non-research-grade row")
    if row.get("captive_cultivated", "").casefold() not in {"false", "f", "0", ""}:
        raise ValueError(f"{species}: captive/cultivated row")
    accuracy = float(row["positional_accuracy"])
    if not (0 <= accuracy < 1000):
        raise ValueError(f"{species}: positional_accuracy outside gate: {accuracy}")
    if row.get("geoprivacy", "").casefold() not in {"", "open"}:
        raise ValueError(f"{species}: non-open user geoprivacy")
    if row.get("taxon_geoprivacy", "").casefold() not in {"", "open"}:
        raise ValueError(f"{species}: non-open taxon geoprivacy")
    observed = date.fromisoformat(row["observed_on"])
    if not (date(2016, 1, 1) <= observed <= date(2025, 12, 31)):
        raise ValueError(f"{species}: date outside gate: {observed}")
    if not row.get("image_url"):
        raise ValueError(f"{species}: missing image_url despite photo gate")
    return spatial_fields(row)


def round_robin_sample(
    rows: Iterable[dict[str, Any]], species: str, seed: str, target: int
) -> list[dict[str, Any]]:
    """Select one encounter per observer-cell-week and spread effort among cells."""

    deduplicated: dict[tuple[str, int, int, int], dict[str, Any]] = {}
    for row in rows:
        key = (
            str(row["user_id"]),
            int(row["grid_x"]),
            int(row["grid_y"]),
            int(row["flowering_week"]),
        )
        previous = deduplicated.get(key)
        if previous is None or stable_rank(seed, species, row["id"]) < stable_rank(
            seed, species, previous["id"]
        ):
            deduplicated[key] = row

    by_cell: dict[tuple[int, int], list[dict[str, Any]]] = defaultdict(list)
    for row in deduplicated.values():
        by_cell[(int(row["grid_x"]), int(row["grid_y"]))].append(row)
    for cell, candidates in by_cell.items():
        candidates.sort(
            key=lambda row: (
                stable_rank(seed, species, cell, row["flowering_week"]),
                stable_rank(seed, species, row["id"]),
            )
        )

    cell_order = sorted(by_cell, key=lambda cell: stable_rank(seed, species, *cell))
    selected: list[dict[str, Any]] = []
    depth = 0
    while len(selected) < target:
        added = 0
        for cell in cell_order:
            candidates = by_cell[cell]
            if depth < len(candidates):
                selected.append(candidates[depth])
                added += 1
                if len(selected) == target:
                    break
        if added == 0:
            break
        depth += 1
    return selected


def fetch_selected_details(
    species: str,
    selected: list[dict[str, Any]],
    raw_dir: Path,
    request_interval: float,
) -> dict[int, dict[str, Any]]:
    details: dict[int, dict[str, Any]] = {}
    for chunk_index in range(0, len(selected), 100):
        ids = [str(row["id"]) for row in selected[chunk_index : chunk_index + 100]]
        params = {"id": ",".join(ids), "per_page": len(ids), "order_by": "id", "order": "asc"}
        url = f"{API_ROOT}/observations?{urllib.parse.urlencode(params)}"
        part = chunk_index // 100 + 1
        json_path = raw_dir / f"selected_{species.replace(' ', '_')}_{part:02d}.json"
        headers_path = raw_dir / f"selected_{species.replace(' ', '_')}_{part:02d}.headers.json"
        downloaded = False
        if json_path.exists() and headers_path.exists():
            payload = json_path.read_bytes()
        else:
            payload, headers = request_bytes(url)
            downloaded = True
            json_path.write_bytes(payload)
            headers_path.write_text(
                json.dumps(headers, indent=2, sort_keys=True) + "\n", encoding="utf-8"
            )
        decoded = json.loads(payload)
        for result in decoded.get("results", []):
            details[int(result["id"])] = result
        if downloaded:
            time.sleep(request_interval)
    if len(details) != len(selected):
        raise RuntimeError(
            f"{species}: enriched {len(details)} of {len(selected)} selected observations"
        )
    return details


def reusable_photo(observation: dict[str, Any]) -> dict[str, Any]:
    photos = [
        photo
        for photo in observation.get("photos", [])
        if str(photo.get("license_code", "")).casefold() in REUSABLE_LICENSES
    ]
    if not photos:
        raise ValueError(f"observation {observation.get('id')} has no reusable photo")
    return min(photos, key=lambda photo: int(photo["id"]))


def medium_url(photo: dict[str, Any]) -> str:
    url = str(photo["url"])
    return url.replace("/square.", "/medium.")


def assign_partition(rows: list[dict[str, Any]], species: str, seed: str) -> None:
    def identity(row: dict[str, Any]) -> int:
        return int(row["id"] if "id" in row else row["observation_id"])

    ranked = sorted(rows, key=lambda row: stable_rank(seed, species, "partition", identity(row)))
    n_development = round(0.4 * len(ranked))
    development_ids = {identity(row) for row in ranked[:n_development]}
    for row in rows:
        row["annotation_partition"] = (
            "development_40" if identity(row) in development_ids else "locked_60"
        )


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"cannot write empty table: {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact-dir", type=Path, required=True)
    parser.add_argument(
        "--public-manifest",
        type=Path,
        default=Path(
            "docs/supporting/jbi_inaturalist_api_development_sample_manifest_v1.json"
        ),
    )
    parser.add_argument("--target-per-species", type=int, default=200)
    parser.add_argument("--seed", default=DEFAULT_SEED)
    parser.add_argument("--request-interval", type=float, default=1.05)
    parser.add_argument("--skip-images", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    artifact_dir = args.artifact_dir.resolve()
    raw_dir = artifact_dir / "raw_api"
    image_dir = artifact_dir / "images"
    raw_dir.mkdir(parents=True, exist_ok=True)
    image_dir.mkdir(parents=True, exist_ok=True)

    acquired_at = datetime.now(timezone.utc).isoformat()
    private_rows: list[dict[str, Any]] = []
    blinded_rows: list[dict[str, Any]] = []
    species_summary: list[dict[str, Any]] = []

    for species, taxon_id in SPECIES.items():
        raw_rows = fetch_species_rows(species, taxon_id, raw_dir, args.request_interval)
        validated: list[dict[str, Any]] = []
        descendant_rows_excluded = 0
        for row in raw_rows:
            if row.get("scientific_name") != species or int(row["taxon_id"]) != taxon_id:
                descendant_rows_excluded += 1
                continue
            validated.append(row | spatial_fields(row) | validate_row(row, species, taxon_id))
        reserve_target = min(
            len(validated), args.target_per_species + max(20, args.target_per_species // 4)
        )
        candidates = round_robin_sample(
            validated, species=species, seed=args.seed, target=reserve_target
        )
        details = fetch_selected_details(species, candidates, raw_dir, args.request_interval)

        species_private: list[dict[str, Any]] = []
        photo_gate_excluded = 0
        image_gate_excluded = 0
        for row in candidates:
            if len(species_private) == args.target_per_species:
                break
            observation = details[int(row["id"])]
            try:
                photo = reusable_photo(observation)
            except ValueError:
                photo_gate_excluded += 1
                continue
            extension = Path(urllib.parse.urlparse(str(photo["url"])).path).suffix or ".jpg"
            blind_id = "FCP-" + stable_rank(args.seed, species, row["id"], photo["id"])[:16].upper()
            image_path = image_dir / f"{blind_id}{extension.lower()}"
            image_sha = ""
            image_bytes = ""
            image_width = ""
            image_height = ""
            if not args.skip_images:
                try:
                    if not image_path.exists():
                        payload, _headers = request_bytes(medium_url(photo))
                        image_path.write_bytes(payload)
                    image_sha = sha256(image_path)
                    image_bytes = image_path.stat().st_size
                    with Image.open(image_path) as image:
                        image_width, image_height = image.size
                except Exception:
                    image_gate_excluded += 1
                    if image_path.exists():
                        image_path.unlink()
                    continue

            private = {
                "blind_id": blind_id,
                "canonical_name": species,
                "taxon_id": taxon_id,
                "observation_id": int(row["id"]),
                "observation_uuid": observation["uuid"],
                "observer_id": int(row["user_id"]),
                "observed_on": row["observed_on"],
                "flowering_week": int(row["flowering_week"]),
                "latitude": row["latitude"],
                "longitude": row["longitude"],
                "positional_accuracy_m": row["positional_accuracy"],
                "grid_x_50km_epsg6933": int(row["grid_x"]),
                "grid_y_50km_epsg6933": int(row["grid_y"]),
                "observation_updated_at": observation.get("updated_at", ""),
                "photo_id": int(photo["id"]),
                "photo_license": photo["license_code"],
                "photo_attribution": photo.get("attribution", ""),
                "photo_medium_url": medium_url(photo),
                "image_file": image_path.name if not args.skip_images else "",
                "image_sha256": image_sha,
                "image_bytes": image_bytes,
                "image_width": image_width,
                "image_height": image_height,
                "annotation_partition": "pending_partition_assignment",
                "acquired_at_utc": acquired_at,
            }
            species_private.append(private)

        assign_partition(species_private, species, args.seed)
        for private in species_private:
            blinded_rows.append(
                {
                    "blind_id": private["blind_id"],
                    "image_file": private["image_file"],
                    "annotation_partition": private["annotation_partition"],
                    "flower_visible": "",
                    "single_encounter_usable": "",
                    "anonymous_morph_code": "",
                    "classification_confidence": "",
                    "exclusion_reason": "",
                    "reviewer_id": "",
                    "reviewed_at_utc": "",
                }
            )

        private_rows.extend(species_private)
        observer_counts = Counter(row["observer_id"] for row in species_private)
        cell_counts = Counter(
            (row["grid_x_50km_epsg6933"], row["grid_y_50km_epsg6933"])
            for row in species_private
        )
        top_observer_n = max(observer_counts.values(), default=0)
        top_cell_n = max(cell_counts.values(), default=0)
        structure_gate_pass = (
            len(species_private) >= 120
            and len(cell_counts) >= 20
            and top_observer_n / max(len(species_private), 1) <= 0.20
        )
        species_summary.append(
            {
                "canonical_name": species,
                "taxon_id": taxon_id,
                "api_taxon_search_rows": len(raw_rows),
                "exact_species_metadata_rows": len(validated),
                "descendant_taxon_rows_excluded": descendant_rows_excluded,
                "reserve_candidates_enriched": len(candidates),
                "detail_photo_license_gate_excluded": photo_gate_excluded,
                "image_download_or_decode_gate_excluded": image_gate_excluded,
                "selected_encounters": len(species_private),
                "unique_50km_cells": len(
                    {(row["grid_x_50km_epsg6933"], row["grid_y_50km_epsg6933"]) for row in species_private}
                ),
                "unique_observers": len({row["observer_id"] for row in species_private}),
                "top_observer_n": top_observer_n,
                "top_observer_share": top_observer_n / max(len(species_private), 1),
                "top_50km_cell_n": top_cell_n,
                "top_50km_cell_share": top_cell_n / max(len(species_private), 1),
                "development_structure_gate_status": (
                    "pass_pending_flower_visibility_and_morph_scoring"
                    if structure_gate_pass
                    else "fail_or_not_evaluable_structure_gate"
                ),
                "development_n": sum(
                    row["annotation_partition"] == "development_40" for row in species_private
                ),
                "locked_n": sum(
                    row["annotation_partition"] == "locked_60" for row in species_private
                ),
                "photo_license_counts": dict(
                    sorted(Counter(row["photo_license"] for row in species_private).items())
                ),
            }
        )
        print(f"{species}: selected {len(species_private)}", flush=True)

    private_path = artifact_dir / "private_provenance.csv"
    blinded_path = artifact_dir / "blinded_annotation_sheet.csv"
    summary_path = artifact_dir / "species_summary.csv"
    write_csv(private_path, private_rows)
    write_csv(blinded_path, blinded_rows)
    write_csv(summary_path, species_summary)

    artifact_manifest = {
        "status": "complete_api_development_sample_not_inferential",
        "protocol": "fcp-inaturalist-photo-development-v1",
        "acquired_at_utc": acquired_at,
        "source": "moving iNaturalist API and observations.csv endpoint",
        "dated_open_dataset_required_for_frozen_pilot": True,
        "dated_open_dataset_selected_at_development_time": "20260827",
        "dated_open_dataset_bundle_bytes": 35093052336,
        "seed": args.seed,
        "target_per_species": args.target_per_species,
        "grid": {"crs": "EPSG:6933", "cell_size_m": GRID_M},
        "observation_window": ["2016-01-01", "2025-12-31"],
        "gate": api_params(0) | {"taxon_id": "species-specific exact active ID"},
        "taxon_query_semantics": (
            "The iNaturalist taxon_id query includes descendant taxa. Rows are admitted only "
            "when scientific_name and row taxon_id both equal the frozen species concept; "
            "descendant exclusions are reported per species."
        ),
        "selection": (
            "one encounter per observer x 50-km cell x ISO flowering week; "
            "deterministic round-robin across cells; up to 200 per species"
        ),
        "annotation_split": "deterministic 40% development and 60% locked within species",
        "images_downloaded": not args.skip_images,
        "n_selected_total": len(private_rows),
        "species_summary": species_summary,
        "private_files_sha256": {
            private_path.name: sha256(private_path),
            blinded_path.name: sha256(blinded_path),
            summary_path.name: sha256(summary_path),
        },
        "claim_ceiling": (
            "Workflow development and image-availability evidence only. The moving API sample "
            "must not estimate morph frequencies, C, G, boundary concordance, organizer effects, "
            "or universality. Frozen inference requires reconstruction from the named dated Open "
            "Dataset snapshot followed by blinded human scoring and all preregistered gates."
        ),
    }
    artifact_manifest_path = artifact_dir / "artifact_manifest.json"
    artifact_manifest_path.write_text(
        json.dumps(artifact_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    public_manifest = {
        key: value
        for key, value in artifact_manifest.items()
        if key not in {"private_files_sha256"}
    }
    public_manifest["artifact_manifest_sha256"] = sha256(artifact_manifest_path)
    public_manifest["public_data_exclusions"] = [
        "exact coordinates",
        "observer identifiers",
        "observation and photo identifiers",
        "image URLs and files",
        "attribution ledger",
    ]
    args.public_manifest.parent.mkdir(parents=True, exist_ok=True)
    args.public_manifest.write_text(
        json.dumps(public_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(public_manifest, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Run the live metadata-only admission and geometry freeze for the FCP atlas."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
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

from fcp_pipeline.image_first_atlas import (
    AtlasMetadataAdapter,
    freeze_atlas_geometry,
    freeze_atlas_metadata,
)


USER_AGENT = "zuizui0223-fcp-image-first-atlas/1.0 (metadata feasibility; research reproducibility)"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def get_json(url: str, *, pause: float, attempts: int = 4) -> dict[str, Any]:
    last_error: Exception | None = None
    for attempt in range(attempts):
        try:
            request = Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
            with urlopen(request, timeout=90) as response:
                payload = json.load(response)
            if pause > 0:
                time.sleep(pause)
            return payload
        except (HTTPError, URLError, TimeoutError) as exc:
            last_error = exc
            if isinstance(exc, HTTPError) and exc.code not in {429, 500, 502, 503, 504}:
                raise
            if attempt + 1 < attempts:
                time.sleep(max(2.0, pause) * (2**attempt))
    raise RuntimeError(f"iNaturalist metadata request failed after {attempts} attempts: {last_error}")


class InaturalistMetadataAdapter(AtlasMetadataAdapter):
    """Read-only iNaturalist v1 adapter with explicit throttling and pagination."""

    def __init__(self, *, base_url: str, pause_seconds: float) -> None:
        self.base_url = base_url.rstrip("/")
        self.pause_seconds = float(pause_seconds)

    @staticmethod
    def _params(query: Mapping[str, Any]) -> dict[str, Any]:
        params: dict[str, Any] = {}
        for key, value in query.items():
            if key == "limit":
                continue
            if isinstance(value, bool):
                params[key] = "true" if value else "false"
            elif isinstance(value, (list, tuple)):
                params[key] = ",".join(str(item) for item in value)
            else:
                params[key] = value
        return params

    def _paged(self, endpoint: str, query: Mapping[str, Any]) -> list[dict[str, Any]]:
        limit = int(query["limit"])
        per_page = min(200, limit)
        pages = math.ceil(limit / per_page)
        params = self._params(query)
        params["per_page"] = per_page
        rows: list[dict[str, Any]] = []
        for page in range(1, pages + 1):
            params["page"] = page
            url = f"{self.base_url}/{endpoint}?{urlencode(params)}"
            payload = get_json(url, pause=self.pause_seconds)
            results = payload.get("results") or []
            rows.extend(dict(row) for row in results)
            if len(rows) >= limit or len(results) < per_page:
                break
        return rows[:limit]

    def species_counts(self, query: Mapping[str, Any]) -> Sequence[Mapping[str, Any]]:
        return self._paged("observations/species_counts", query)

    def observations(
        self,
        taxon_id: int,
        query: Mapping[str, Any],
    ) -> Sequence[Mapping[str, Any]]:
        adjusted = dict(query)
        adjusted["taxon_id"] = int(taxon_id)
        adjusted.pop("rank", None)
        return self._paged("observations", adjusted)


def write_csv(path: Path, rows: Sequence[Mapping[str, Any]], preferred: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    available = {key for row in rows for key in row}
    fields = [key for key in preferred if key in available]
    fields.extend(sorted(available - set(fields)))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--contract",
        type=Path,
        default=Path("docs/supporting/jbi_image_first_atlas_contract_v1.json"),
    )
    parser.add_argument(
        "--cohort",
        type=Path,
        default=Path("data/atlas/jbi_image_first_atlas_cohort_v1.csv"),
    )
    parser.add_argument(
        "--observations",
        type=Path,
        default=Path("data/atlas/jbi_image_first_atlas_observation_manifest_v1.csv"),
    )
    parser.add_argument(
        "--feasibility",
        type=Path,
        default=Path("docs/supporting/jbi_image_first_atlas_metadata_feasibility_v1.json"),
    )
    parser.add_argument(
        "--geometry",
        type=Path,
        default=Path("docs/supporting/jbi_image_first_atlas_geometry_freeze_v1.json"),
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("docs/supporting/jbi_image_first_atlas_freeze_manifest_v1.json"),
    )
    args = parser.parse_args()

    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    source = contract["metadata_source"]
    adapter = InaturalistMetadataAdapter(
        base_url=str(source["base_url"]),
        pause_seconds=float(source["request_pause_seconds"]),
    )
    freeze = freeze_atlas_metadata(contract, adapter)

    write_csv(
        args.cohort,
        freeze.cohort,
        (
            "cohort_order",
            "metadata_rank",
            "species",
            "inat_taxon_id",
            "inat_genus_id",
            "flowering_annotated_observation_count",
            "selected_photographs",
            "candidate_image_pixels_opened",
            "flower_colour_used",
            "literature_classification_used_for_admission",
        ),
    )
    write_csv(
        args.observations,
        freeze.observations,
        (
            "cohort_order",
            "species",
            "inat_taxon_id",
            "inat_genus_id",
            "observation_id",
            "photo_id",
            "photo_url_large",
            "photo_license",
            "attribution",
            "latitude",
            "longitude",
            "positional_accuracy_m",
            "observed_on",
            "observed_month",
            "local_solar_quarter",
            "observer_id",
            "observer",
            "primary_thinning_cell",
            "sensitivity_thinning_cell",
            "selection_hash",
        ),
    )
    write_json(args.feasibility, freeze.audit)

    geometry = freeze_atlas_geometry(freeze.observations, contract) if len(freeze.cohort) == 50 else {
        "protocol": contract["protocol"],
        "status": "not_run_metadata_gate_failed",
        "selection_used_image_pixels": False,
        "selection_used_continuous_colour": False,
        "next_gate": "STOP: metadata cohort incomplete",
    }
    write_json(args.geometry, geometry)

    manifest = {
        "protocol": contract["protocol"],
        "status": (
            "metadata_and_geometry_frozen_before_image_pixels"
            if freeze.audit["status"] == "pass_50_species_metadata_only"
            and geometry["status"] == "geometry_scale_frozen"
            else "not_evaluable"
        ),
        "candidate_image_pixels_opened": False,
        "flower_roi_used": False,
        "continuous_colour_used": False,
        "literature_classification_used_for_admission": False,
        "files": {
            str(args.contract).replace("\\", "/"): sha256(args.contract),
            str(args.cohort).replace("\\", "/"): sha256(args.cohort),
            str(args.observations).replace("\\", "/"): sha256(args.observations),
            str(args.feasibility).replace("\\", "/"): sha256(args.feasibility),
            str(args.geometry).replace("\\", "/"): sha256(args.geometry),
        },
        "species_admitted": freeze.audit["species_admitted"],
        "selected_observations": freeze.audit["selected_observations"],
        "selected_primary_scale_km": geometry.get("selected_primary_scale_km"),
        "next_gate": (
            "flower ROI extraction under the separately frozen measurement contract"
            if geometry["status"] == "geometry_scale_frozen"
            else "STOP before image download or colour measurement"
        ),
    }
    write_json(args.manifest, manifest)
    print(json.dumps(manifest, indent=2, ensure_ascii=False))
    return 0 if manifest["status"] == "metadata_and_geometry_frozen_before_image_pixels" else 2


if __name__ == "__main__":
    raise SystemExit(main())

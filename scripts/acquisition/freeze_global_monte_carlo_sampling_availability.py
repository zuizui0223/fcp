#!/usr/bin/env python3
"""Freeze metadata-only global sampling-availability surfaces for RGFCA bias audits."""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from fcp_pipeline.random_photo_pool import InaturalistObservationClient, equal_area_cell_bounds
from fcp_pipeline.shared_transition_surface import EqualAreaGrid

ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "docs/supporting/global_monte_carlo_observation_bias_contract_v1.json"
OUT_CELLS = ROOT / "data/frozen/global_monte_carlo_sampling_availability_cells_v1.csv"
OUT_MANIFEST = ROOT / "docs/supporting/global_monte_carlo_sampling_availability_manifest_v1.json"
OUTPUTS = (OUT_CELLS, OUT_MANIFEST)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def base_query(census: dict[str, object], grid: EqualAreaGrid, cell_id: int) -> dict[str, object]:
    params: dict[str, object] = {
        "taxon_id": int(census["taxon_id"]),
        "quality_grade": str(census["quality_grade"]),
        "photos": "true",
        "geo": "true",
        "rank": str(census["rank"]),
        "acc_below": int(census["maximum_positional_accuracy_m"]),
        "obscuration": str(census["obscuration"]),
        "order_by": str(census["order_by"]),
        "order": str(census["order"]),
        "per_page": int(census["per_page"]),
        "page": int(census["page"]),
    }
    params.update(equal_area_cell_bounds(grid, cell_id))
    return params


def layer_query(census: dict[str, object], grid: EqualAreaGrid, cell_id: int, layer: str) -> dict[str, object]:
    params = base_query(census, grid, cell_id)
    if layer in {"license_eligible_photo_records", "flowering_annotated_eligible_records"}:
        params["photo_license"] = ",".join(str(x) for x in census["allowed_photo_licenses"])
    if layer == "flowering_annotated_eligible_records":
        params["term_id"] = int(census["flowering_term_id"])
        params["term_value_id"] = int(census["flowering_term_value_id"])
    return params


def safe_fraction(num: float, den: float) -> float:
    return float(num / den) if den > 0 else float("nan")


def main() -> None:
    existing = [str(path.relative_to(ROOT)) for path in OUTPUTS if path.exists()]
    if existing:
        raise RuntimeError("sampling-availability outputs already exist; refusing rerun: " + ", ".join(existing))

    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    if contract.get("status") != "frozen_before_any_rgfca_global_flower_colour_field":
        raise RuntimeError("observation-bias contract is not frozen before RGFCA colour outcomes")
    census = contract["metadata_only_sampling_availability_census"]
    if census.get("candidate_image_pixels_opened") is not False or census.get("flower_colour_used") is not False:
        raise RuntimeError("sampling-availability census unexpectedly permits pixels or colour")

    grid_spec = census["grid"]
    grid = EqualAreaGrid(n_lon=int(grid_spec["n_lon"]), n_sinlat=int(grid_spec["n_sinlat"]))
    if grid.n_cells != int(grid_spec["cells"]):
        raise RuntimeError("sampling-availability grid drift")
    layer_names = [str(item["name"]) for item in census["layers"]]
    expected_layers = [
        "all_research_photo_records",
        "license_eligible_photo_records",
        "flowering_annotated_eligible_records",
    ]
    if layer_names != expected_layers:
        raise RuntimeError(f"sampling-availability layer drift: {layer_names}")
    expected_attempts = int(census["fixed_request_attempts"])
    if expected_attempts != grid.n_cells * len(layer_names):
        raise RuntimeError("sampling-availability request-count contract is inconsistent")

    client = InaturalistObservationClient(
        request_interval_seconds=float(census["request_interval_seconds"]),
        timeout_seconds=45.0,
        max_retries=int(census["request_retries"]),
        user_agent="fcp-rgfca-sampling-availability/1.0 (github.com/zuizui0223/fcp)",
    )

    rows: list[dict[str, object]] = []
    attempts = 0
    errors = 0
    for cell_id in range(grid.n_cells):
        row: dict[str, object] = {"cell_id": int(cell_id), **equal_area_cell_bounds(grid, cell_id)}
        for layer in layer_names:
            attempts += 1
            params = layer_query(census, grid, cell_id, layer)
            try:
                payload = client.observations(params)
                total = int(payload.get("total_results") or 0)
                error = ""
            except Exception as exc:
                total = -1
                error = f"{type(exc).__name__}:{str(exc)[:180]}"
                errors += 1
            row[layer] = total
            row[f"{layer}_request_error"] = error
        all_n = float(row["all_research_photo_records"])
        lic_n = float(row["license_eligible_photo_records"])
        flow_n = float(row["flowering_annotated_eligible_records"])
        row["log1p_all_research_photo_records"] = math.log1p(all_n) if all_n >= 0 else float("nan")
        row["log1p_license_eligible_photo_records"] = math.log1p(lic_n) if lic_n >= 0 else float("nan")
        row["log1p_flowering_annotated_eligible_records"] = math.log1p(flow_n) if flow_n >= 0 else float("nan")
        row["licence_eligible_fraction"] = safe_fraction(lic_n, all_n) if all_n >= 0 and lic_n >= 0 else float("nan")
        row["flowering_annotation_fraction_given_licence"] = safe_fraction(flow_n, lic_n) if lic_n >= 0 and flow_n >= 0 else float("nan")
        rows.append(row)
        if (cell_id + 1) % 20 == 0 or cell_id + 1 == grid.n_cells:
            print(json.dumps({"processed_cells": cell_id + 1, "cells": grid.n_cells, "request_errors": errors}), flush=True)

    if attempts != expected_attempts:
        raise RuntimeError(f"sampling-availability request count drift: {attempts} != {expected_attempts}")

    frame = pd.DataFrame(rows).sort_values("cell_id", kind="mergesort").reset_index(drop=True)
    error_fraction = float(errors / attempts)
    failure = census["failure_policy"]
    status = (
        "complete_metadata_only_sampling_availability_surface"
        if error_fraction <= float(failure["maximum_error_fraction"])
        else str(failure["if_exceeded"])
    )

    valid_effort = frame.loc[frame["all_research_photo_records"] >= 0, "log1p_all_research_photo_records"].to_numpy(dtype=float)
    high_effort_cut = float(np.quantile(valid_effort, 0.90)) if len(valid_effort) else float("nan")
    frame["top_decile_all_photo_effort"] = frame["log1p_all_research_photo_records"] >= high_effort_cut

    OUT_CELLS.parent.mkdir(parents=True, exist_ok=True)
    OUT_MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(OUT_CELLS, index=False, lineterminator="\n")
    manifest = {
        "protocol": contract["protocol"],
        "status": status,
        "candidate_image_pixels_opened": False,
        "flower_colour_used": False,
        "cells": int(len(frame)),
        "layers": layer_names,
        "request_attempts": int(attempts),
        "request_errors": int(errors),
        "request_error_fraction": error_fraction,
        "maximum_error_fraction": float(failure["maximum_error_fraction"]),
        "top_decile_log1p_all_photo_effort_threshold": high_effort_cut,
        "top_decile_cells": int(frame["top_decile_all_photo_effort"].sum()),
        "interpretation": "Counts are target-group sampling-availability proxies that combine biological availability and observer/platform effort. They are frozen robustness covariates, not direct estimates of true sampling probability.",
        "lineage": {"contract_sha256": sha256_file(CONTRACT)},
        "files": {"cell_surface": {"path": str(OUT_CELLS.relative_to(ROOT)), "sha256": sha256_file(OUT_CELLS)}},
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2), flush=True)


if __name__ == "__main__":
    main()

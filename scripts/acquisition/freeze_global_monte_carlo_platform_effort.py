#!/usr/bin/env python3
"""Freeze an all-taxa platform-activity surface for RGFCA observation-bias diagnostics."""
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
CONTRACT = ROOT / "docs/supporting/global_monte_carlo_platform_effort_contract_v1.json"
TARGET = ROOT / "data/frozen/global_monte_carlo_sampling_availability_cells_v1.csv"
OUT = ROOT / "data/frozen/global_monte_carlo_platform_effort_cells_v1.csv"
MANIFEST = ROOT / "docs/supporting/global_monte_carlo_platform_effort_manifest_v1.json"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def gini(values: np.ndarray) -> float:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x) & (x >= 0)]
    if len(x) == 0 or x.sum() <= 0:
        return float("nan")
    x = np.sort(x)
    n = len(x)
    return float((2 * np.sum(np.arange(1, n + 1) * x) / (n * x.sum())) - (n + 1) / n)


def main() -> int:
    if OUT.exists() or MANIFEST.exists():
        raise RuntimeError("platform-effort output already exists; refusing rerun")
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    if contract.get("status") != "frozen_after_target_group_availability_census_and_before_any_rgfca_global_flower_colour_field":
        raise RuntimeError("platform-effort contract status drift")
    if contract["outcome_firewall"]["candidate_image_pixels_opened"] is not False or contract["outcome_firewall"]["flower_colour_used"] is not False:
        raise RuntimeError("platform-effort contract permits forbidden outcomes")

    grid = EqualAreaGrid(int(contract["grid"]["n_lon"]), int(contract["grid"]["n_sinlat"]))
    if grid.n_cells != int(contract["grid"]["cells"]):
        raise RuntimeError("platform-effort grid drift")
    target = pd.read_csv(TARGET).sort_values("cell_id", kind="mergesort").reset_index(drop=True)
    if len(target) != grid.n_cells or target["cell_id"].nunique() != grid.n_cells:
        raise RuntimeError("target-group availability cell frame drift")

    q = contract["query"]
    client = InaturalistObservationClient(
        request_interval_seconds=float(q["request_interval_seconds"]),
        timeout_seconds=45.0,
        max_retries=int(q["request_retries"]),
        user_agent="fcp-rgfca-platform-effort/1.0 (github.com/zuizui0223/fcp)",
    )
    rows: list[dict[str, object]] = []
    errors = 0
    for cell_id in range(grid.n_cells):
        params: dict[str, object] = {
            "quality_grade": str(q["quality_grade"]),
            "photos": "true",
            "geo": "true",
            "acc_below": int(q["maximum_positional_accuracy_m"]),
            "obscuration": str(q["obscuration"]),
            "order_by": str(q["order_by"]),
            "order": str(q["order"]),
            "per_page": int(q["per_page"]),
            "page": int(q["page"]),
        }
        params.update(equal_area_cell_bounds(grid, cell_id))
        try:
            payload = client.observations(params)
            total = int(payload.get("total_results") or 0)
            error = ""
        except Exception as exc:
            total = -1
            error = f"{type(exc).__name__}:{str(exc)[:180]}"
            errors += 1
        target_n = int(target.loc[target["cell_id"] == cell_id, "all_research_photo_records"].iloc[0])
        rows.append({
            "cell_id": cell_id,
            **equal_area_cell_bounds(grid, cell_id),
            "all_taxa_research_photo_records": total,
            "request_error": error,
            "log1p_all_taxa_research_photo_records": (math.log1p(total) if total >= 0 else float("nan")),
            "target_group_all_research_photo_records": target_n,
            "target_group_photo_fraction_of_all_taxa_records": (float(target_n / total) if total > 0 else float("nan")),
        })
        if (cell_id + 1) % 20 == 0 or cell_id + 1 == grid.n_cells:
            print(json.dumps({"processed_cells": cell_id + 1, "cells": grid.n_cells, "errors": errors}), flush=True)

    frame = pd.DataFrame(rows).sort_values("cell_id", kind="mergesort").reset_index(drop=True)
    attempts = len(frame)
    if attempts != int(contract["fixed_request_attempts"]):
        raise RuntimeError("platform-effort request-count drift")
    error_fraction = float(errors / attempts)
    fail = contract["failure_policy"]
    status = "complete_metadata_only_platform_effort_surface" if error_fraction <= float(fail["maximum_error_fraction"]) else str(fail["if_exceeded"])

    valid = frame.loc[frame["all_taxa_research_photo_records"] >= 0].copy()
    counts = valid["all_taxa_research_photo_records"].to_numpy(dtype=float)
    target_log = np.log1p(valid["target_group_all_research_photo_records"].to_numpy(dtype=float))
    platform_log = valid["log1p_all_taxa_research_photo_records"].to_numpy(dtype=float)
    ranks = pd.DataFrame({"platform": platform_log, "target": target_log}).rank(method="average")
    rho = float(ranks["platform"].corr(ranks["target"])) if len(ranks) >= 3 else float("nan")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(OUT, index=False, lineterminator="\n")
    manifest = {
        "protocol": contract["protocol"],
        "status": status,
        "candidate_image_pixels_opened": False,
        "flower_colour_used": False,
        "cells": int(len(frame)),
        "request_attempts": attempts,
        "request_errors": int(errors),
        "request_error_fraction": error_fraction,
        "all_taxa_total_records": int(counts.sum()) if len(counts) else 0,
        "all_taxa_zero_cells": int(np.count_nonzero(counts == 0)),
        "all_taxa_gini_across_equal_area_cells": gini(counts),
        "all_taxa_maximum_cell_share": (float(counts.max() / counts.sum()) if counts.sum() > 0 else None),
        "spearman_log_platform_vs_log_target_group_availability": (rho if np.isfinite(rho) else None),
        "interpretation": "All-taxa record density is a platform-wide activity proxy, not a pure observer-effort probability. It is frozen as a negative-control covariate before any RGFCA global flower-colour field.",
        "lineage": {
            "contract_sha256": sha256_file(CONTRACT),
            "target_group_surface_sha256": sha256_file(TARGET),
            "cell_surface_sha256": sha256_file(OUT),
        },
        "files": {"cell_surface": str(OUT.relative_to(ROOT))},
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

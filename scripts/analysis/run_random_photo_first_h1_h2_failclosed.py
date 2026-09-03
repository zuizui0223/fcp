#!/usr/bin/env python3
"""Fail-closed orchestration for the frozen random photo-first H1 -> H2 run.

This module does not alter any H1/H2 statistic, null, seed, predictor, threshold,
or environmental input. It only prevents an H2 support/sensitivity insufficiency
from erasing an already completed H1 result.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

import run_random_photo_first_h1_h2 as frozen


def _not_evaluable_h2_payload(
    h1_result: dict[str, Any],
    h2_contract: dict[str, Any],
    reason: str,
    *,
    climate_join_opened: bool,
) -> dict[str, Any]:
    h1_p = float(h1_result["primary"]["p_upper"])
    h1_alpha = float(h1_result["alpha"])
    return {
        "protocol": h2_contract["protocol"],
        "status": "not_evaluable_h2_after_evaluable_h1",
        "reason": str(reason),
        "h1_primary_p_upper": h1_p,
        "h1_primary_supported": bool(h1_p < h1_alpha),
        "climate_colour_join_opened": bool(climate_join_opened),
        "confirmatory_climate_concordance_claim": False,
        "claim_ceiling": h2_contract["claim_ceiling"],
    }


def robust_run_h2(
    h1_result: dict[str, Any],
    primary_maps: Any,
    grid: Any,
    h2_contract: dict[str, Any],
) -> dict[str, Any]:
    """Run frozen H2, making primary/secondary non-evaluability explicit."""

    source_spec = h2_contract["environment_source"]["primary_grid_source"]
    climate_join_opened = False
    try:
        climate_source = frozen.verified_climate_source(
            source_spec["path"], source_spec["content_sha256"]
        )
        climate_cells = frozen.aggregate_climate_to_h1_grid(climate_source, grid=grid)
        climate_edges = frozen.build_edge_climate_contrasts(climate_cells, grid=grid)
        frozen.H2_CLIMATE_CELLS.parent.mkdir(parents=True, exist_ok=True)
        climate_cells.to_csv(frozen.H2_CLIMATE_CELLS, index=False, lineterminator="\n")
        climate_edges.to_csv(frozen.H2_CLIMATE_EDGES, index=False, lineterminator="\n")
        climate_join_opened = True

        minimum_edges = int(h2_contract["primary_test"]["minimum_supported_edges"])
        primary = frozen.climate_concordance_test(
            primary_maps.observed.edge_table,
            primary_maps.null_persistence,
            primary_maps.edge_ids,
            climate_edges,
            predictor="multivariate_climate_distance",
            subset="global",
            minimum_supported_edges=minimum_edges,
        )
    except ValueError as exc:
        payload = _not_evaluable_h2_payload(
            h1_result,
            h2_contract,
            str(exc),
            climate_join_opened=climate_join_opened,
        )
        frozen.write_json(frozen.H2_RESULT, payload)
        return payload

    # Secondary decomposition may be unavailable without invalidating primary H2.
    try:
        drivers = frozen.climate_driver_decomposition(
            primary_maps.observed.edge_table,
            primary_maps.null_persistence,
            primary_maps.edge_ids,
            climate_edges,
            minimum_supported_edges=minimum_edges,
        )
    except ValueError as exc:
        drivers = pd.DataFrame(
            [
                {
                    "status": "not_evaluable_secondary_driver_decomposition",
                    "reason": str(exc),
                }
            ]
        )
    drivers.to_csv(frozen.H2_DRIVERS, index=False, lineterminator="\n")

    sensitivity_rows: list[dict[str, Any]] = []
    for subset in ("within_biome", "within_realm"):
        try:
            result = frozen.climate_concordance_test(
                primary_maps.observed.edge_table,
                primary_maps.null_persistence,
                primary_maps.edge_ids,
                climate_edges,
                predictor="multivariate_climate_distance",
                subset=subset,
                minimum_supported_edges=minimum_edges,
            )
        except ValueError as exc:
            sensitivity_rows.append(
                {
                    "type": "subset",
                    "value": subset,
                    "status": f"not_evaluable:{str(exc)[:300]}",
                    "weighted_r": np.nan,
                    "p_upper": np.nan,
                    "supported_edges": 0,
                }
            )
        else:
            sensitivity_rows.append(
                {
                    "type": "subset",
                    "value": subset,
                    "status": "complete_evaluable",
                    "weighted_r": result.statistic,
                    "p_upper": result.p_upper,
                    "supported_edges": result.supported_edges,
                }
            )

    for scale_spec in h2_contract["environment_source"]["scale_sensitivities"]:
        try:
            source = frozen.verified_climate_source(
                scale_spec["path"], scale_spec["content_sha256"]
            )
            cells = frozen.aggregate_climate_to_h1_grid(source, grid=grid)
            edges = frozen.build_edge_climate_contrasts(cells, grid=grid)
            result = frozen.climate_concordance_test(
                primary_maps.observed.edge_table,
                primary_maps.null_persistence,
                primary_maps.edge_ids,
                edges,
                predictor="multivariate_climate_distance",
                subset="global",
                minimum_supported_edges=minimum_edges,
            )
        except ValueError as exc:
            sensitivity_rows.append(
                {
                    "type": "climate_scale_km",
                    "value": int(scale_spec["scale_km"]),
                    "status": f"not_evaluable:{str(exc)[:300]}",
                    "weighted_r": np.nan,
                    "p_upper": np.nan,
                    "supported_edges": 0,
                }
            )
        else:
            sensitivity_rows.append(
                {
                    "type": "climate_scale_km",
                    "value": int(scale_spec["scale_km"]),
                    "status": "complete_evaluable",
                    "weighted_r": result.statistic,
                    "p_upper": result.p_upper,
                    "supported_edges": result.supported_edges,
                }
            )

    sensitivity = pd.DataFrame(sensitivity_rows)
    sensitivity.to_csv(frozen.H2_SENSITIVITIES, index=False, lineterminator="\n")

    h1_pass = float(h1_result["primary"]["p_upper"]) < float(h1_result["alpha"])
    h2_pass = float(primary.p_upper) < float(h2_contract["primary_test"]["alpha"])
    if h1_pass and h2_pass:
        hierarchy = "support_hierarchical_recurrent_boundary_macroclimate_concordance"
    elif not h1_pass:
        hierarchy = "diagnostic_only_h1_not_supported_no_climate_mechanism_claim"
    else:
        hierarchy = "h1_supported_but_h2_macroclimate_concordance_not_supported"

    payload = {
        "protocol": h2_contract["protocol"],
        "status": "complete_h2_evaluable",
        "hierarchical_decision": hierarchy,
        "h1_primary_p_upper": float(h1_result["primary"]["p_upper"]),
        "h1_primary_supported": bool(h1_pass),
        "primary": {
            "climate_source_scale_km": int(source_spec["scale_km"]),
            "weighted_r": float(primary.statistic),
            "p_upper": float(primary.p_upper),
            "alpha": float(h2_contract["primary_test"]["alpha"]),
            "supported_edges": int(primary.supported_edges),
            "support": bool(h2_pass),
        },
        "driver_decomposition": drivers.to_dict(orient="records"),
        "sensitivities": sensitivity.to_dict(orient="records"),
        "files": {
            "climate_cells": str(frozen.H2_CLIMATE_CELLS.relative_to(frozen.ROOT)),
            "climate_edges": str(frozen.H2_CLIMATE_EDGES.relative_to(frozen.ROOT)),
            "drivers": str(frozen.H2_DRIVERS.relative_to(frozen.ROOT)),
            "sensitivities": str(frozen.H2_SENSITIVITIES.relative_to(frozen.ROOT)),
        },
        "claim_ceiling": h2_contract["claim_ceiling"],
    }
    frozen.write_json(frozen.H2_RESULT, payload)
    return payload


def main() -> int:
    frozen.run_h2 = robust_run_h2
    return int(frozen.main())


if __name__ == "__main__":
    raise SystemExit(main())

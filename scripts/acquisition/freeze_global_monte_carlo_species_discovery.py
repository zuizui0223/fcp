#!/usr/bin/env python3
"""Freeze a repeated metadata-only global species-discovery census.

This stage deliberately opens no candidate image pixels.  It repeats the same
outcome-blind equal-area iNaturalist metadata query for a fixed number of rounds
so the global discovery frame is not determined by one random page per cell.
Per-request failures are recorded and never replaced or retried.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pandas as pd

from fcp_pipeline.random_photo_pool import (
    InaturalistObservationClient,
    inat_query_for_cell,
    parse_candidate_observation,
)
from fcp_pipeline.shared_transition_surface import EqualAreaGrid

ROOT = Path(__file__).resolve().parents[2]
CONTRACT = ROOT / "docs/supporting/global_monte_carlo_species_discovery_contract_v1.json"
BASELINE = ROOT / "data/frozen/random_photo_first_candidate_pool_v1.csv"
OBS_INDEX = ROOT / "data/frozen/global_monte_carlo_species_discovery_observation_index_v1.csv.gz"
SPECIES = ROOT / "data/frozen/global_monte_carlo_species_discovery_species_v1.csv"
ROUND_AUDIT = ROOT / "data/frozen/global_monte_carlo_species_discovery_round_audit_v1.csv"
CELL_AUDIT = ROOT / "data/frozen/global_monte_carlo_species_discovery_cell_audit_v1.csv"
MANIFEST = ROOT / "docs/supporting/global_monte_carlo_species_discovery_manifest_v1.json"
OUTPUTS = (OBS_INDEX, SPECIES, ROUND_AUDIT, CELL_AUDIT, MANIFEST)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def jaccard(left: set[str], right: set[str]) -> float:
    union = left | right
    return float(len(left & right) / len(union)) if union else float("nan")


def compact_frame(frame: pd.DataFrame, *, round_id: int) -> pd.DataFrame:
    required = ["cell_id", "observation_id", "photo_id", "species", "inat_taxon_id"]
    missing = sorted(set(required) - set(frame.columns))
    if missing:
        raise RuntimeError(f"discovery frame missing required columns: {missing}")
    out = frame[required].copy()
    out.insert(0, "round_id", int(round_id))
    return out


def main() -> None:
    existing = [str(path.relative_to(ROOT)) for path in OUTPUTS if path.exists()]
    if existing:
        raise RuntimeError(
            "global species-discovery output already exists; refusing new random queries: "
            + ", ".join(existing)
        )
    contract = json.loads(CONTRACT.read_text(encoding="utf-8"))
    if contract["status"] != "frozen_before_any_new_global_discovery_query":
        raise RuntimeError("species-discovery contract is not in the frozen pre-query state")
    if contract["outcome_firewall"]["flower_colour_used"] is not False:
        raise RuntimeError("species-discovery contract unexpectedly permits colour")

    discovery = contract["fresh_discovery"]
    grid_spec = discovery["grid"]
    grid = EqualAreaGrid(n_lon=int(grid_spec["n_lon"]), n_sinlat=int(grid_spec["n_sinlat"]))
    if grid.n_cells != int(grid_spec["cells"]):
        raise RuntimeError("species-discovery grid drift")
    rounds = int(discovery["rounds"])
    expected_attempts = int(discovery["total_fresh_request_attempts"])
    if expected_attempts != rounds * grid.n_cells:
        raise RuntimeError("species-discovery request-count contract is inconsistent")

    baseline = pd.read_csv(BASELINE)
    if len(baseline) != int(contract["baseline_round_zero"]["observations"]):
        raise RuntimeError("baseline discovery observation count drift")
    if baseline["species"].nunique() != int(contract["baseline_round_zero"]["species"]):
        raise RuntimeError("baseline discovery species count drift")

    compact_parts: list[pd.DataFrame] = [compact_frame(baseline, round_id=0)]
    baseline_species = set(baseline["species"].astype(str))
    cumulative_species = set(baseline_species)
    round_sets: list[set[str]] = []
    round_rows: list[dict[str, object]] = [
        {
            "round_id": 0,
            "request_attempts": 0,
            "request_errors": 0,
            "api_returned": int(len(baseline)),
            "accepted": int(len(baseline)),
            "round_species": int(len(baseline_species)),
            "new_species": int(len(baseline_species)),
            "cumulative_species": int(len(cumulative_species)),
            "jaccard_previous_fresh_round": np.nan,
        }
    ]
    cell_species: dict[int, set[str]] = {
        cell: set(group["species"].astype(str))
        for cell, group in baseline.groupby("cell_id", sort=False)
    }
    for cell in range(grid.n_cells):
        cell_species.setdefault(cell, set())

    allowed = frozenset(str(x).casefold() for x in discovery["allowed_photo_licenses"])
    client = InaturalistObservationClient(
        request_interval_seconds=float(discovery["request_interval_seconds"]),
        timeout_seconds=45.0,
        max_retries=int(discovery["request_retries"]),
        user_agent="fcp-global-monte-carlo-species-discovery/1.0 (github.com/zuizui0223/fcp)",
    )

    total_attempts = 0
    total_errors = 0
    total_returned = 0
    total_accepted = 0
    previous_fresh: set[str] | None = None

    for round_id in range(1, rounds + 1):
        accepted_rows: list[dict[str, object]] = []
        round_errors = 0
        round_returned = 0
        round_accepted = 0
        for cell_id in range(grid.n_cells):
            total_attempts += 1
            params = inat_query_for_cell(
                grid,
                cell_id,
                per_page=int(discovery["per_cell_random_page_size"]),
                taxon_id=int(discovery["taxon_id"]),
                flowering_term_id=int(discovery["flowering_term_id"]),
                flowering_term_value_id=int(discovery["flowering_term_value_id"]),
                maximum_positional_accuracy_m=int(discovery["maximum_positional_accuracy_m"]),
                allowed_photo_licenses=tuple(discovery["allowed_photo_licenses"]),
            )
            try:
                payload = client.observations(params)
                raw = payload.get("results") or []
                if not isinstance(raw, Sequence) or isinstance(raw, (str, bytes)):
                    raise RuntimeError("iNaturalist results is not a sequence")
            except Exception:
                round_errors += 1
                total_errors += 1
                continue
            round_returned += len(raw)
            total_returned += len(raw)
            seen_obs_request: set[int] = set()
            seen_photo_request: set[int] = set()
            for observation in raw:
                if not isinstance(observation, Mapping):
                    continue
                parsed = parse_candidate_observation(
                    observation,
                    expected_cell_id=cell_id,
                    grid=grid,
                    maximum_positional_accuracy_m=float(discovery["maximum_positional_accuracy_m"]),
                    allowed_photo_licenses=allowed,
                )
                if parsed is None:
                    continue
                oid = int(parsed["observation_id"])
                pid = int(parsed["photo_id"])
                if oid in seen_obs_request or pid in seen_photo_request:
                    continue
                seen_obs_request.add(oid)
                seen_photo_request.add(pid)
                accepted_rows.append(parsed)
                round_accepted += 1
                total_accepted += 1
                cell_species[cell_id].add(str(parsed["species"]))

        fresh = pd.DataFrame(accepted_rows)
        if len(fresh):
            compact_parts.append(compact_frame(fresh, round_id=round_id))
            fresh_species = set(fresh["species"].astype(str))
        else:
            fresh_species = set()
        new_species = fresh_species - cumulative_species
        cumulative_species |= fresh_species
        round_rows.append(
            {
                "round_id": round_id,
                "request_attempts": int(grid.n_cells),
                "request_errors": int(round_errors),
                "api_returned": int(round_returned),
                "accepted": int(round_accepted),
                "round_species": int(len(fresh_species)),
                "new_species": int(len(new_species)),
                "cumulative_species": int(len(cumulative_species)),
                "jaccard_previous_fresh_round": (
                    np.nan if previous_fresh is None else jaccard(fresh_species, previous_fresh)
                ),
            }
        )
        round_sets.append(fresh_species)
        previous_fresh = fresh_species
        print(
            json.dumps(
                {
                    "round": round_id,
                    "accepted": round_accepted,
                    "round_species": len(fresh_species),
                    "new_species": len(new_species),
                    "cumulative_species": len(cumulative_species),
                    "errors": round_errors,
                },
                sort_keys=True,
            ),
            flush=True,
        )

    if total_attempts != expected_attempts:
        raise RuntimeError("species-discovery fresh request-attempt count drift")

    index = pd.concat(compact_parts, ignore_index=True)
    index = index.sort_values(["round_id", "cell_id", "observation_id", "photo_id"], kind="mergesort").reset_index(drop=True)
    # Species summary uses unique observations so repeated random rediscovery does not inflate information.
    unique_obs = index.drop_duplicates(["observation_id", "photo_id"], keep="first").copy()
    first_round = index.groupby(["species", "inat_taxon_id"], observed=True)["round_id"].min().rename("first_discovery_round")
    rounds_seen = index.groupby(["species", "inat_taxon_id"], observed=True)["round_id"].nunique().rename("n_rounds_seen")
    unique_obs_n = unique_obs.groupby(["species", "inat_taxon_id"], observed=True)["observation_id"].nunique().rename("unique_discovery_observations")
    cells_n = unique_obs.groupby(["species", "inat_taxon_id"], observed=True)["cell_id"].nunique().rename("occupied_discovery_cells")
    species_frame = pd.concat([first_round, rounds_seen, unique_obs_n, cells_n], axis=1).reset_index()
    species_frame = species_frame.sort_values(["inat_taxon_id", "species"], kind="mergesort").reset_index(drop=True)

    round_audit = pd.DataFrame(round_rows)
    cell_audit = pd.DataFrame(
        {
            "cell_id": list(range(grid.n_cells)),
            "unique_species_across_baseline_and_fresh_rounds": [len(cell_species[i]) for i in range(grid.n_cells)],
        }
    )
    fresh_union = set().union(*round_sets) if round_sets else set()
    odd = set().union(*(round_sets[i] for i in range(0, len(round_sets), 2))) if round_sets else set()
    even = set().union(*(round_sets[i] for i in range(1, len(round_sets), 2))) if len(round_sets) > 1 else set()
    error_fraction = float(total_errors / total_attempts) if total_attempts else 1.0
    maximum_error_fraction = float(contract["failure_policy"]["maximum_error_fraction_for_coverage_claim"])
    coverage_evaluable = bool(error_fraction <= maximum_error_fraction)

    OBS_INDEX.parent.mkdir(parents=True, exist_ok=True)
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    index.to_csv(OBS_INDEX, index=False, compression="gzip", lineterminator="\n")
    species_frame.to_csv(SPECIES, index=False, lineterminator="\n")
    round_audit.to_csv(ROUND_AUDIT, index=False, lineterminator="\n")
    cell_audit.to_csv(CELL_AUDIT, index=False, lineterminator="\n")

    manifest = {
        "protocol": contract["protocol"],
        "status": (
            "complete_metadata_only_global_species_discovery"
            if coverage_evaluable
            else str(contract["failure_policy"]["if_exceeded"])
        ),
        "candidate_image_pixels_opened": False,
        "flower_colour_used": False,
        "baseline": {
            "observations": int(len(baseline)),
            "species": int(len(baseline_species)),
        },
        "fresh": {
            "rounds": rounds,
            "request_attempts": int(total_attempts),
            "expected_request_attempts": expected_attempts,
            "request_errors": int(total_errors),
            "request_error_fraction": error_fraction,
            "coverage_error_fraction_ceiling": maximum_error_fraction,
            "coverage_evaluable": coverage_evaluable,
            "api_returned": int(total_returned),
            "accepted_query_rows": int(total_accepted),
            "fresh_unique_species": int(len(fresh_union)),
        },
        "combined": {
            "species": int(len(species_frame)),
            "unique_observations": int(len(unique_obs)),
            "new_species_beyond_baseline": int(len(set(species_frame["species"].astype(str)) - baseline_species)),
        },
        "stability": {
            "odd_even_fresh_round_species_jaccard": jaccard(odd, even),
            "baseline_vs_fresh_species_jaccard": jaccard(baseline_species, fresh_union),
            "final_round_new_species": int(round_audit.iloc[-1]["new_species"]),
            "final_round_cumulative_species": int(round_audit.iloc[-1]["cumulative_species"]),
        },
        "next_gate": contract["next_gate"],
        "lineage": {
            "contract_sha256": sha256_file(CONTRACT),
            "baseline_sha256": sha256_file(BASELINE),
        },
        "files": {
            "observation_index": {"path": str(OBS_INDEX.relative_to(ROOT)), "sha256": sha256_file(OBS_INDEX)},
            "species_frame": {"path": str(SPECIES.relative_to(ROOT)), "sha256": sha256_file(SPECIES)},
            "round_audit": {"path": str(ROUND_AUDIT.relative_to(ROOT)), "sha256": sha256_file(ROUND_AUDIT)},
            "cell_audit": {"path": str(CELL_AUDIT.relative_to(ROOT)), "sha256": sha256_file(CELL_AUDIT)},
        },
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2), flush=True)


if __name__ == "__main__":
    main()

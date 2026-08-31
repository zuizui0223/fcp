"""Coordinate firewall and completeness gate for atlas image measurement."""

from __future__ import annotations

from collections import Counter
import hashlib
from typing import Any, Mapping, Sequence


INFERENCE_PROTOCOL = "jbi-image-first-global-flower-colour-atlas-inference-v3"
WORKER_MANIFEST_FIELDS = {
    "measurement_id",
    "species_blind_id",
    "image_filename",
    "photo_license",
}
TERMINAL_MEASUREMENT_STATES = {
    "automated_colour_state_admitted",
    "automated_colour_state_not_evaluable",
    "image_acquisition_failed",
}


def validate_inference_contract(contract: Mapping[str, Any]) -> None:
    """Fail closed if the prospective post-simulation decision tree drifts."""

    if contract.get("protocol") != INFERENCE_PROTOCOL:
        raise ValueError("unexpected atlas inference protocol")
    if contract.get("frozen_before_scaleout_candidate_pixels") is not True:
        raise ValueError("inference contract was not frozen before scale-out pixels")
    evidence = contract.get("qualification_evidence", {})
    if evidence.get("geographic_shared_boundary", {}).get("decision") != (
        "not_evaluable_failed_preimage_signal_recovery"
    ):
        raise ValueError("v2 geographic STOP must remain not_evaluable")
    if evidence.get("environmental_pollinator_overlay_null", {}).get("decision") != (
        "pass_spatially_constrained_overlay_null"
    ):
        raise ValueError("overlay-null qualification evidence is missing")
    measurement = contract.get("scaleout_measurement_gate", {})
    expected = {
        "minimum_admitted_fraction_per_species": 0.7,
        "minimum_admitted_observations_per_species": 210,
        "minimum_background_control_fraction_among_admitted": 0.7,
        "minimum_evaluable_species_per_cohort": 20,
        "minimum_evaluable_species_total": 160,
        "all_eight_cohorts_must_finish": True,
        "blinding_salt": "fcp-atlas-v3-coordinate-firewall",
    }
    for key, value in expected.items():
        if measurement.get(key) != value:
            raise ValueError(f"scale-out measurement rule changed: {key}")
    branches = contract.get("ordered_inference_v3", {}).get("branches", [])
    if [row.get("branch") for row in branches] != [
        "shared_geographic_concentration",
        "environmental_concordance",
        "pollinator_biogeographic_concordance",
    ]:
        raise ValueError("v3 branch order changed")
    if (
        branches[0].get("fixed_outcome") != "not_evaluable"
        or branches[0].get("real_colour_test_permitted") is not False
    ):
        raise ValueError("stopped geographic branch reopened")


def _blind(salt: str, label: str, value: object) -> str:
    payload = f"{salt}\x1f{label}\x1f{value}".encode("utf-8")
    return hashlib.sha256(payload).hexdigest().upper()[:20]


def measurement_shard(measurement_id: str, shard_count: int) -> int:
    """Assign one blinded image to a stable zero-based compute shard."""

    if shard_count < 1:
        raise ValueError("shard_count must be positive")
    digest = hashlib.sha256(
        f"fcp-atlas-v3-worker-shard\x1f{measurement_id}".encode("utf-8")
    ).digest()
    return int.from_bytes(digest[:8], "big") % shard_count


def select_measurement_shard(
    rows: Sequence[Mapping[str, Any]], *, shard_index: int, shard_count: int
) -> list[dict[str, Any]]:
    """Validate the blind interface and return its deterministic shard."""

    if shard_index < 0 or shard_index >= shard_count:
        raise ValueError("shard_index must lie in [0, shard_count)")
    selected: list[dict[str, Any]] = []
    seen: set[str] = set()
    for raw in rows:
        if set(raw) != WORKER_MANIFEST_FIELDS:
            raise ValueError("measurement worker manifest fields changed or leaked")
        measurement_id = str(raw["measurement_id"])
        if not measurement_id or measurement_id in seen:
            raise ValueError("measurement worker IDs must be non-empty and unique")
        seen.add(measurement_id)
        expected_filename = f"{measurement_id}.jpg"
        if str(raw["image_filename"]) != expected_filename:
            raise ValueError("measurement image filename is not the blinded ID")
        if measurement_shard(measurement_id, shard_count) == shard_index:
            selected.append(dict(raw))
    return selected


def validate_measurement_result_rows(
    rows: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    """Validate location-free terminal output before shard aggregation."""

    forbidden_tokens = (
        "latitude",
        "longitude",
        "observed",
        "observer",
        "taxon",
        "species",
        "attribution",
        "photo_url",
        "environment",
        "pollinator",
    )
    validated: list[dict[str, Any]] = []
    seen: set[str] = set()
    for raw in rows:
        leaked = [
            str(key)
            for key in raw
            if any(token in str(key).casefold() for token in forbidden_tokens)
            and str(key) != "species_blind_id"
        ]
        if leaked:
            raise ValueError(f"measurement result leaked protected fields: {leaked}")
        measurement_id = str(raw.get("measurement_id") or "")
        if not measurement_id or measurement_id in seen:
            raise ValueError("measurement result IDs must be non-empty and unique")
        seen.add(measurement_id)
        status = str(raw.get("automated_colour_state_status") or "")
        if status not in TERMINAL_MEASUREMENT_STATES:
            raise ValueError(f"unknown measurement terminal state: {status!r}")
        background = raw.get("background_features_available")
        if not isinstance(background, bool):
            raise ValueError("background_features_available must be boolean")
        if status == "image_acquisition_failed" and background:
            raise ValueError("failed acquisition cannot have background features")
        validated.append(dict(raw))
    return validated


def build_measurement_firewall(
    rows: Sequence[Mapping[str, Any]],
    *,
    salt: str = "fcp-atlas-v3-coordinate-firewall",
) -> dict[str, list[dict[str, Any]]]:
    """Split frozen observations into location-blind measurement and sealed keys."""

    measurement_rows: list[dict[str, Any]] = []
    species_rows: list[dict[str, Any]] = []
    coordinate_rows: list[dict[str, Any]] = []
    seen_measurements: set[str] = set()
    for raw in rows:
        required = (
            "cohort_id",
            "species",
            "inat_taxon_id",
            "observation_id",
            "photo_id",
            "photo_url_large",
            "photo_license",
            "latitude",
            "longitude",
            "observed_month",
            "local_solar_quarter",
            "observer_id",
        )
        missing = [key for key in required if raw.get(key) in (None, "")]
        if missing:
            raise ValueError(f"scale-out firewall row is missing: {missing}")
        leaked = []
        for key, value in raw.items():
            folded = str(key).casefold()
            if folded == "candidate_image_pixels_opened" and value is False:
                continue
            if any(token in folded for token in ("colour", "color", "roi", "pixel")):
                leaked.append(str(key))
        if leaked:
            raise ValueError(f"pre-measurement rows contain image outcomes: {sorted(leaked)}")
        measurement_id = "FCPM-" + _blind(salt, "photo", raw["photo_id"])
        species_blind_id = "FCPS-" + _blind(salt, "species", raw["inat_taxon_id"])
        if measurement_id in seen_measurements:
            raise ValueError("measurement IDs must be unique")
        seen_measurements.add(measurement_id)
        measurement_rows.append(
            {
                "measurement_id": measurement_id,
                "species_blind_id": species_blind_id,
                "image_filename": f"{measurement_id}.jpg",
                "photo_license": str(raw["photo_license"]),
            }
        )
        species_rows.append(
            {
                "measurement_id": measurement_id,
                "species_blind_id": species_blind_id,
                "cohort_id": str(raw["cohort_id"]),
            }
        )
        coordinate_rows.append(
            {
                "measurement_id": measurement_id,
                "species_blind_id": species_blind_id,
                "cohort_id": str(raw["cohort_id"]),
                "species": str(raw["species"]),
                "inat_taxon_id": str(raw["inat_taxon_id"]),
                "observation_id": str(raw["observation_id"]),
                "photo_id": str(raw["photo_id"]),
                "photo_url_large": str(raw["photo_url_large"]),
                "photo_license": str(raw["photo_license"]),
                "attribution": str(raw.get("attribution") or ""),
                "latitude": float(raw["latitude"]),
                "longitude": float(raw["longitude"]),
                "positional_accuracy_m": float(raw.get("positional_accuracy_m") or 0),
                "observed_on": str(raw.get("observed_on") or ""),
                "observed_month": int(raw["observed_month"]),
                "local_solar_quarter": int(raw["local_solar_quarter"]),
                "observer_id": str(raw["observer_id"]),
                "observer": str(raw.get("observer") or ""),
                "primary_thinning_cell": str(raw.get("primary_thinning_cell") or ""),
                "sensitivity_thinning_cell": str(
                    raw.get("sensitivity_thinning_cell") or ""
                ),
            }
        )
    return {
        "measurement_manifest": measurement_rows,
        "sealed_species_key": species_rows,
        "sealed_coordinate_key": coordinate_rows,
    }


def evaluate_scaleout_measurement_gate(
    measurement_results: Sequence[Mapping[str, Any]],
    sealed_species_key: Sequence[Mapping[str, Any]],
    gate: Mapping[str, Any],
) -> dict[str, Any]:
    """Decide species/cohort completeness without opening the coordinate key."""

    keys = {str(row["measurement_id"]): dict(row) for row in sealed_species_key}
    if len(keys) != len(sealed_species_key):
        raise ValueError("sealed species key contains duplicate measurement IDs")
    results = {str(row["measurement_id"]): row for row in measurement_results}
    if len(results) != len(measurement_results) or set(results) != set(keys):
        raise ValueError("measurement results must exactly cover the sealed denominator")

    totals: Counter[str] = Counter()
    admitted: Counter[str] = Counter()
    background: Counter[str] = Counter()
    cohort_for_species: dict[str, str] = {}
    for measurement_id, key in keys.items():
        species_id = str(key["species_blind_id"])
        cohort = str(key["cohort_id"])
        previous = cohort_for_species.setdefault(species_id, cohort)
        if previous != cohort:
            raise ValueError("one blinded species cannot cross cohorts")
        totals[species_id] += 1
        result = results[measurement_id]
        status = str(result.get("automated_colour_state_status", ""))
        if status not in TERMINAL_MEASUREMENT_STATES:
            raise ValueError(f"unknown measurement terminal state: {status!r}")
        if status == "automated_colour_state_admitted":
            admitted[species_id] += 1
            if result.get("background_features_available") is True:
                background[species_id] += 1

    species_results: list[dict[str, Any]] = []
    evaluable_by_cohort: Counter[str] = Counter()
    for species_id in sorted(totals):
        total = totals[species_id]
        admitted_count = admitted[species_id]
        admitted_fraction = admitted_count / total
        background_fraction = (
            background[species_id] / admitted_count if admitted_count else 0.0
        )
        evaluable = (
            admitted_count >= int(gate["minimum_admitted_observations_per_species"])
            and admitted_fraction >= float(gate["minimum_admitted_fraction_per_species"])
            and background_fraction
            >= float(gate["minimum_background_control_fraction_among_admitted"])
        )
        cohort = cohort_for_species[species_id]
        if evaluable:
            evaluable_by_cohort[cohort] += 1
        species_results.append(
            {
                "species_blind_id": species_id,
                "cohort_id": cohort,
                "frozen_denominator": total,
                "admitted_observations": admitted_count,
                "admitted_fraction": admitted_fraction,
                "background_control_fraction": background_fraction,
                "status": "measurement_evaluable" if evaluable else "not_evaluable",
            }
        )

    cohorts = sorted({str(row["cohort_id"]) for row in sealed_species_key})
    cohort_results = [
        {
            "cohort_id": cohort,
            "evaluable_species": evaluable_by_cohort[cohort],
            "status": (
                "measurement_evaluable"
                if evaluable_by_cohort[cohort]
                >= int(gate["minimum_evaluable_species_per_cohort"])
                else "not_evaluable"
            ),
        }
        for cohort in cohorts
    ]
    total_evaluable = sum(row["status"] == "measurement_evaluable" for row in species_results)
    passed = (
        len(cohorts) == 8
        and all(row["status"] == "measurement_evaluable" for row in cohort_results)
        and total_evaluable >= int(gate["minimum_evaluable_species_total"])
    )
    return {
        "status": (
            "pass_scaleout_measurement_completeness"
            if passed
            else "not_evaluable_scaleout_measurement_completeness"
        ),
        "coordinates_opened": False,
        "frozen_measurements": len(keys),
        "evaluable_species": total_evaluable,
        "species_results": species_results,
        "cohort_results": cohort_results,
        "coordinate_join_permitted": passed,
    }

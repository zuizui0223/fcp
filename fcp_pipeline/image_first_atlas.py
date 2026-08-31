"""Prospective metadata and geometry freeze for the image-first FCP atlas.

The public interface is deliberately small:

``freeze_atlas_metadata(contract, adapter)``
    Select a metadata-only cohort and balanced observation manifest through an
    injected iNaturalist adapter.

``freeze_atlas_geometry(rows, contract)``
    Choose the finest predeclared spatial scale supported by species identities
    and coordinates alone.

Neither function accepts image pixels or colour values.  This keeps the cohort,
admission and scale decisions upstream of flower ROI and colour measurement.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import date
import hashlib
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Protocol, Sequence

import numpy as np

from .shared_transition_surface import (
    EqualAreaGrid,
    build_edge_cell_geometry,
    geometry_opportunity_summary,
)
from .spatial_graph import spherical_knn_edges


PROTOCOL = "jbi-image-first-global-flower-colour-atlas-v1"


class AtlasMetadataAdapter(Protocol):
    """Seam for the remote metadata source; tests use an in-memory adapter."""

    def species_counts(self, query: Mapping[str, Any]) -> Sequence[Mapping[str, Any]]:
        """Return ranked species-count records for the frozen metadata query."""

    def observations(
        self,
        taxon_id: int,
        query: Mapping[str, Any],
    ) -> Sequence[Mapping[str, Any]]:
        """Return candidate observation metadata for one taxon."""


@dataclass(frozen=True)
class AtlasMetadataFreeze:
    """Metadata-only output returned before any candidate image is opened."""

    cohort: tuple[dict[str, Any], ...]
    observations: tuple[dict[str, Any], ...]
    audit: dict[str, Any]


def sha256_text(*parts: object) -> str:
    payload = "\x1f".join(str(part) for part in parts).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def sha256_lf_canonical(path: Path) -> str:
    """Hash a frozen text artifact after canonicalising platform newlines to LF."""

    payload = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(payload).hexdigest()


def spatial_cell(latitude: float, longitude: float, degrees: float) -> str:
    if not math.isfinite(latitude) or not math.isfinite(longitude):
        raise ValueError("cell coordinates must be finite")
    if degrees <= 0:
        raise ValueError("cell size must be positive")
    lat_bin = math.floor(latitude / degrees) * degrees
    lon_bin = math.floor(longitude / degrees) * degrees
    return f"{lat_bin:.6f},{lon_bin:.6f}"


def local_solar_quarter(month: int, latitude: float) -> int:
    """Return a hemisphere-adjusted quarter for the frozen phenology sensitivity."""

    if month < 1 or month > 12:
        raise ValueError("month must lie in 1..12")
    shifted = (month - 1 + (6 if latitude < 0 else 0)) % 12
    return shifted // 3 + 1


def validate_atlas_contract(contract: Mapping[str, Any]) -> None:
    """Fail closed when a colour-dependent or literature-led admission rule appears."""

    if contract.get("protocol") != PROTOCOL:
        raise ValueError(f"unexpected atlas protocol: {contract.get('protocol')!r}")
    if contract.get("status") != "frozen_before_metadata_queries_and_image_pixels":
        raise ValueError("atlas contract is not frozen at the pre-metadata gate")

    firewall = contract.get("outcome_firewall", {})
    forbidden = [
        key
        for key in (
            "candidate_image_pixels_opened",
            "flower_roi_used",
            "continuous_colour_used",
            "literature_classification_used_for_admission",
            "stage_a_effects_used",
            "stage_b_surfaces_used",
            "environmental_layers_used",
        )
        if firewall.get(key) is not False
    ]
    if forbidden:
        raise ValueError(f"outcome firewall is open for: {forbidden}")

    admission = contract.get("admission", {})
    if int(admission.get("target_species", 0)) != 50:
        raise ValueError("atlas pilot must target exactly 50 prospectively admitted species")
    tiers = [int(value) for value in admission.get("sample_size_tiers_descending", [])]
    if tiers != [500, 400, 300]:
        raise ValueError("sample-size tiers must be frozen as 500, 400, 300")
    if int(admission.get("maximum_species_per_genus", 0)) < 1:
        raise ValueError("maximum_species_per_genus must be positive")

    source = contract.get("metadata_source", {})
    if int(source.get("angiosperm_taxon_id", 0)) != 47125:
        raise ValueError("atlas root must be the frozen iNaturalist Angiospermae taxon")
    if int(source.get("flowering_term_id", 0)) != 12 or int(
        source.get("flowering_term_value_id", 0)
    ) != 13:
        raise ValueError("flowering annotation term/value must be frozen as 12/13")
    if float(source.get("maximum_positional_accuracy_m", math.inf)) != 5000.0:
        raise ValueError("maximum positional accuracy must be exactly 5 km")

    thinning = contract.get("thinning", {})
    if float(thinning.get("primary_cell_degrees", 0)) != 0.25:
        raise ValueError("primary thinning must use 0.25-degree cells")
    if float(thinning.get("sensitivity_cell_degrees", 0)) != 0.5:
        raise ValueError("thinning sensitivity must use 0.5-degree cells")

    geometry = contract.get("geometry_only_scale_selection", {})
    scales = [int(row.get("scale_km", 0)) for row in geometry.get("candidates", [])]
    if scales != [100, 250, 500]:
        raise ValueError("geometry candidates must be ordered 100, 250, 500 km")
    if geometry.get("priority_rule") != "finest_passing_scale_before_colour_measurement":
        raise ValueError("geometry scale priority rule is not frozen")

    display = contract.get("display", {})
    if display.get("species_labels_visible") is not False:
        raise ValueError("atlas display must remain species-free")
    if display.get("inference_remains_species_conditioned") is not True:
        raise ValueError("species-free display cannot remove species from inference")


def _coordinates(observation: Mapping[str, Any]) -> tuple[float | None, float | None]:
    geojson = observation.get("geojson") or {}
    coords = geojson.get("coordinates") or []
    if isinstance(coords, list) and len(coords) >= 2:
        try:
            return float(coords[1]), float(coords[0])
        except (TypeError, ValueError):
            pass
    location = observation.get("location")
    if isinstance(location, str) and "," in location:
        try:
            latitude, longitude = location.split(",", 1)
            return float(latitude), float(longitude)
        except ValueError:
            pass
    if "latitude" in observation and "longitude" in observation:
        try:
            return float(observation["latitude"]), float(observation["longitude"])
        except (TypeError, ValueError):
            pass
    return None, None


def _month(observation: Mapping[str, Any]) -> int | None:
    raw = str(observation.get("observed_on") or "")
    if len(raw) < 10:
        return None
    try:
        parsed = date.fromisoformat(raw[:10])
    except ValueError:
        return None
    return parsed.month


def _large_photo_url(url: str | None) -> str | None:
    if not url:
        return None
    for size in ("square", "small", "medium"):
        token = f"/{size}."
        if token in url:
            return url.replace(token, "/large.")
    return url


def _prepare_observation(
    observation: Mapping[str, Any],
    taxon: Mapping[str, Any],
    contract: Mapping[str, Any],
) -> dict[str, Any] | None:
    source = contract["metadata_source"]
    allowed_licenses = {str(value).casefold() for value in source["allowed_photo_licenses"]}
    try:
        observation_id = int(observation["id"])
        accuracy = float(observation["positional_accuracy"])
    except (KeyError, TypeError, ValueError):
        return None
    if not math.isfinite(accuracy) or accuracy < 0 or accuracy > float(
        source["maximum_positional_accuracy_m"]
    ):
        return None
    if bool(observation.get("obscured")):
        return None
    if observation.get("geoprivacy") not in (None, "", "open"):
        return None

    latitude, longitude = _coordinates(observation)
    if latitude is None or longitude is None:
        return None
    if not (-90 <= latitude <= 90 and -180 <= longitude <= 180):
        return None
    month = _month(observation)
    if month is None:
        return None

    eligible_photos = []
    for photo in observation.get("photos") or []:
        license_code = str(photo.get("license_code") or "").casefold()
        try:
            photo_id = int(photo["id"])
        except (KeyError, TypeError, ValueError):
            continue
        if license_code in allowed_licenses and photo.get("url"):
            eligible_photos.append((photo_id, photo))
    if not eligible_photos:
        return None
    photo_id, photo = min(eligible_photos, key=lambda item: item[0])

    user = observation.get("user") or {}
    observer_id = user.get("id", observation.get("observer_id"))
    if observer_id in (None, ""):
        return None

    thinning = contract["thinning"]
    primary_degrees = float(thinning["primary_cell_degrees"])
    sensitivity_degrees = float(thinning["sensitivity_cell_degrees"])
    species = str(taxon["name"])
    salt = str(contract["admission"]["stable_hash_salt"])
    return {
        "species": species,
        "inat_taxon_id": int(taxon["id"]),
        "inat_genus_id": int(taxon["parent_id"]),
        "observation_id": observation_id,
        "photo_id": str(photo_id),
        "photo_url_large": _large_photo_url(str(photo.get("url"))),
        "photo_license": str(photo.get("license_code") or ""),
        "attribution": str(photo.get("attribution") or ""),
        "latitude": float(latitude),
        "longitude": float(longitude),
        "positional_accuracy_m": float(accuracy),
        "observed_on": str(observation.get("observed_on") or ""),
        "observed_month": int(month),
        "local_solar_quarter": local_solar_quarter(month, latitude),
        "observer_id": str(observer_id),
        "observer": str(user.get("login") or observation.get("observer") or ""),
        "primary_thinning_cell": spatial_cell(latitude, longitude, primary_degrees),
        "sensitivity_thinning_cell": spatial_cell(latitude, longitude, sensitivity_degrees),
        "selection_hash": sha256_text(
            salt,
            species,
            observation_id,
            photo_id,
        ),
    }


def _balanced_selection(
    rows: Sequence[Mapping[str, Any]],
    contract: Mapping[str, Any],
) -> list[dict[str, Any]]:
    caps = contract["thinning"]["hard_caps"]
    admission = contract["admission"]
    maximum = max(int(value) for value in admission["sample_size_tiers_descending"])
    species = str(rows[0]["species"]) if rows else ""
    salt = str(admission["stable_hash_salt"])

    by_cell: dict[str, list[dict[str, Any]]] = defaultdict(list)
    seen_observations: set[int] = set()
    seen_photos: set[str] = set()
    for raw in rows:
        row = dict(raw)
        observation_id = int(row["observation_id"])
        photo_id = str(row["photo_id"])
        if observation_id in seen_observations or photo_id in seen_photos:
            continue
        seen_observations.add(observation_id)
        seen_photos.add(photo_id)
        by_cell[str(row["primary_thinning_cell"])].append(row)
    for cell_rows in by_cell.values():
        cell_rows.sort(key=lambda row: (row["selection_hash"], int(row["observation_id"])))

    cell_order = sorted(
        by_cell,
        key=lambda cell: sha256_text(salt, species, "primary-cell", cell),
    )
    cursors = {cell: 0 for cell in cell_order}
    selected: list[dict[str, Any]] = []
    observer_counts: Counter[str] = Counter()
    primary_cell_counts: Counter[str] = Counter()
    sensitivity_cell_counts: Counter[str] = Counter()
    month_counts: Counter[int] = Counter()

    def can_take(row: Mapping[str, Any]) -> bool:
        return (
            observer_counts[str(row["observer_id"])]
            < int(caps["maximum_per_observer"])
            and primary_cell_counts[str(row["primary_thinning_cell"])]
            < int(caps["maximum_per_primary_cell"])
            and sensitivity_cell_counts[str(row["sensitivity_thinning_cell"])]
            < int(caps["maximum_per_sensitivity_cell"])
            and month_counts[int(row["observed_month"])]
            < int(caps["maximum_per_calendar_month"])
        )

    def take(row: dict[str, Any]) -> None:
        selected.append(row)
        observer_counts[str(row["observer_id"])] += 1
        primary_cell_counts[str(row["primary_thinning_cell"])] += 1
        sensitivity_cell_counts[str(row["sensitivity_thinning_cell"])] += 1
        month_counts[int(row["observed_month"])] += 1

    while len(selected) < maximum:
        progressed = False
        for cell in cell_order:
            cell_rows = by_cell[cell]
            while cursors[cell] < len(cell_rows):
                row = cell_rows[cursors[cell]]
                cursors[cell] += 1
                if not can_take(row):
                    continue
                take(row)
                progressed = True
                break
            if len(selected) >= maximum:
                break
        if not progressed:
            break

    for tier in [int(value) for value in admission["sample_size_tiers_descending"]]:
        if len(selected) >= tier:
            return selected[:tier]
    return []


def _selection_qc(
    selected: Sequence[Mapping[str, Any]],
    candidate_count: int,
    contract: Mapping[str, Any],
) -> dict[str, Any]:
    observers = Counter(str(row["observer_id"]) for row in selected)
    primary_cells = Counter(str(row["primary_thinning_cell"]) for row in selected)
    sensitivity_cells = Counter(str(row["sensitivity_thinning_cell"]) for row in selected)
    months = Counter(int(row["observed_month"]) for row in selected)
    quarters = Counter(int(row["local_solar_quarter"]) for row in selected)
    n = len(selected)
    gates = contract["thinning"]["admission_gates"]
    failures: list[str] = []
    checks = {
        "minimum_unique_observers": len(observers),
        "minimum_primary_spatial_cells": len(primary_cells),
        "minimum_sensitivity_spatial_cells": len(sensitivity_cells),
        "minimum_calendar_months": len(months),
        "minimum_local_solar_quarters": len(quarters),
    }
    for key, observed in checks.items():
        if observed < int(gates[key]):
            failures.append(key)

    fractions = {
        "maximum_observer_fraction": max(observers.values(), default=0) / n if n else 1.0,
        "maximum_primary_cell_fraction": max(primary_cells.values(), default=0) / n if n else 1.0,
        "maximum_sensitivity_cell_fraction": max(sensitivity_cells.values(), default=0) / n if n else 1.0,
    }
    for key, observed in fractions.items():
        if observed > float(gates[key]):
            failures.append(key)

    return {
        "candidate_count": int(candidate_count),
        "selected_count": int(n),
        **checks,
        **fractions,
        "calendar_month_counts": sorted(months.items()),
        "local_solar_quarter_counts": sorted(quarters.items()),
        "gate_failures": failures,
        "gate_pass": bool(n) and not failures,
    }


def _species_query(contract: Mapping[str, Any]) -> dict[str, Any]:
    source = contract["metadata_source"]
    return {
        "taxon_id": int(source["angiosperm_taxon_id"]),
        "rank": "species",
        "quality_grade": "research",
        "photos": True,
        "geo": True,
        "geoprivacy": "open",
        "captive": False,
        "term_id": int(source["flowering_term_id"]),
        "term_value_id": int(source["flowering_term_value_id"]),
        "acc_below": int(source["maximum_positional_accuracy_m"]),
        "created_d2": str(source["created_d2"]),
        "photo_license": list(source["allowed_photo_licenses"]),
        "limit": int(source["candidate_species_pool_size"]),
    }


def _observation_query(contract: Mapping[str, Any]) -> dict[str, Any]:
    source = contract["metadata_source"]
    return {
        **_species_query(contract),
        "limit": int(contract["admission"]["maximum_candidates_per_species"]),
        "order_by": "id",
        "order": "asc",
    }


def freeze_atlas_metadata(
    contract: Mapping[str, Any],
    adapter: AtlasMetadataAdapter,
) -> AtlasMetadataFreeze:
    """Freeze the 50-species metadata cohort without opening any image pixels."""

    validate_atlas_contract(contract)
    admission = contract["admission"]
    target_species = int(admission["target_species"])
    excluded_names = {str(value).casefold() for value in admission["excluded_frozen_species"]}
    maximum_per_genus = int(admission["maximum_species_per_genus"])

    raw_counts = list(adapter.species_counts(_species_query(contract)))
    candidates: list[dict[str, Any]] = []
    for record in raw_counts:
        taxon = record.get("taxon") or {}
        try:
            count = int(record["count"])
            taxon_id = int(taxon["id"])
            parent_id = int(taxon["parent_id"])
        except (KeyError, TypeError, ValueError):
            continue
        if str(taxon.get("rank") or "").casefold() != "species":
            continue
        if taxon.get("is_active") is not True:
            continue
        name = str(taxon.get("name") or "").strip()
        if not name or name.casefold() in excluded_names:
            continue
        candidates.append(
            {
                "count": count,
                "taxon": {
                    "id": taxon_id,
                    "name": name,
                    "rank": "species",
                    "parent_id": parent_id,
                },
            }
        )
    candidates.sort(
        key=lambda row: (-int(row["count"]), str(row["taxon"]["name"]), int(row["taxon"]["id"]))
    )

    cohort: list[dict[str, Any]] = []
    selected_rows: list[dict[str, Any]] = []
    audit_rows: list[dict[str, Any]] = []
    genus_counts: Counter[int] = Counter()
    used_photos: set[str] = set()

    for rank, record in enumerate(candidates, start=1):
        if len(cohort) >= target_species:
            break
        taxon = record["taxon"]
        genus_id = int(taxon["parent_id"])
        base = {
            "metadata_rank": rank,
            "species": str(taxon["name"]),
            "inat_taxon_id": int(taxon["id"]),
            "inat_genus_id": genus_id,
            "flowering_annotated_observation_count": int(record["count"]),
        }
        if genus_counts[genus_id] >= maximum_per_genus:
            audit_rows.append({**base, "status": "taxonomic_cap_skip", "gate_pass": False})
            continue

        raw_observations = adapter.observations(int(taxon["id"]), _observation_query(contract))
        prepared = [
            row
            for observation in raw_observations
            if (
                row := _prepare_observation(observation, taxon, contract)
            )
            is not None
            and str(row["photo_id"]) not in used_photos
        ]
        selected = _balanced_selection(prepared, contract)
        qc = _selection_qc(selected, len(prepared), contract)
        if not qc["gate_pass"]:
            audit_rows.append({**base, "status": "metadata_gate_failed", **qc})
            continue

        cohort_order = len(cohort) + 1
        for row in selected:
            row["cohort_order"] = cohort_order
            used_photos.add(str(row["photo_id"]))
        cohort.append(
            {
                "cohort_order": cohort_order,
                **base,
                "selected_photographs": len(selected),
                "candidate_image_pixels_opened": False,
                "flower_colour_used": False,
                "literature_classification_used_for_admission": False,
            }
        )
        selected_rows.extend(selected)
        genus_counts[genus_id] += 1
        audit_rows.append({**base, "status": "admitted", **qc})

    complete = len(cohort) == target_species
    sample_counts = Counter(int(row["selected_photographs"]) for row in cohort)
    audit = {
        "protocol": PROTOCOL,
        "status": "pass_50_species_metadata_only" if complete else "not_evaluable_insufficient_species",
        "candidate_image_pixels_opened": False,
        "flower_roi_used": False,
        "continuous_colour_used": False,
        "literature_classification_used_for_admission": False,
        "stage_a_effects_used": False,
        "stage_b_surfaces_used": False,
        "environmental_layers_used": False,
        "species_count_records_received": len(raw_counts),
        "species_candidates_after_static_filters": len(candidates),
        "species_audited": len(audit_rows),
        "species_admitted": len(cohort),
        "selected_observations": len(selected_rows),
        "sample_size_tier_counts": {str(key): value for key, value in sorted(sample_counts.items())},
        "species_results": audit_rows,
        "next_gate": (
            "freeze geometry-only 100/250/500-km scale before flower ROI or colour"
            if complete
            else "STOP: expand only the predeclared metadata candidate pool under a versioned amendment"
        ),
    }
    return AtlasMetadataFreeze(
        cohort=tuple(cohort),
        observations=tuple(selected_rows),
        audit=audit,
    )


def freeze_atlas_geometry(
    rows: Sequence[Mapping[str, Any]],
    contract: Mapping[str, Any],
) -> dict[str, Any]:
    """Select the finest passing 100/250/500-km scale from geometry only."""

    validate_atlas_contract(contract)
    geometry_contract = contract["geometry_only_scale_selection"]
    k = int(geometry_contract["knn_k"])
    min_edges_per_cell = int(geometry_contract["minimum_edges_per_species_cell"])
    criteria = geometry_contract["passing_criteria"]

    by_species: dict[str, list[Mapping[str, Any]]] = defaultdict(list)
    for row in rows:
        by_species[str(row["species"])].append(row)
    if len(by_species) != int(contract["admission"]["target_species"]):
        raise ValueError("geometry input must contain the complete 50-species cohort")

    base_graphs: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = {}
    for species, species_rows in sorted(by_species.items()):
        latitude = np.asarray([float(row["latitude"]) for row in species_rows], dtype=float)
        longitude = np.asarray([float(row["longitude"]) for row in species_rows], dtype=float)
        edges, distance = spherical_knn_edges(latitude, longitude, k=k)
        base_graphs[species] = (latitude, longitude, edges, distance)

    candidates: list[dict[str, Any]] = []
    selected_scale: int | None = None
    for scale in geometry_contract["candidates"]:
        scale_km = int(scale["scale_km"])
        grid = EqualAreaGrid(n_lon=int(scale["n_lon"]), n_sinlat=int(scale["n_sinlat"]))
        passing_geometries = []
        species_rows = []
        for species, (latitude, longitude, edges, distance) in base_graphs.items():
            try:
                geometry = build_edge_cell_geometry(
                    latitude,
                    longitude,
                    edges,
                    distance,
                    grid=grid,
                    max_edge_km=scale_km,
                    min_edges_per_cell=min_edges_per_cell,
                )
                retained_edges = int(len(geometry.retained_edges))
                detectable_cells = int(np.count_nonzero(geometry.detectable))
            except ValueError:
                geometry = None
                retained_edges = 0
                detectable_cells = 0
            evaluable = (
                geometry is not None
                and retained_edges >= int(criteria["minimum_retained_edges_per_species"])
                and detectable_cells >= int(criteria["minimum_detectable_cells_per_species"])
            )
            species_rows.append(
                {
                    "species": species,
                    "retained_edges": retained_edges,
                    "detectable_cells": detectable_cells,
                    "geometry_evaluable": evaluable,
                }
            )
            if evaluable and geometry is not None:
                passing_geometries.append(geometry)

        if passing_geometries:
            opportunity = geometry_opportunity_summary(
                passing_geometries,
                min_detectable_species=int(criteria["minimum_detectable_species_per_shared_cell"]),
            )
        else:
            opportunity = {
                "n_cells": grid.n_cells,
                "n_cells_A_ge_1": 0,
                "n_cells_A_ge_2": 0,
                "n_cells_A_ge_3": 0,
                "n_cells_A_ge_4": 0,
                "max_A": 0,
                "species_with_any_shared_opportunity": 0,
                "retained_edges_per_species": [],
                "detectable_cells_per_species": [],
            }
        evaluable_species = sum(row["geometry_evaluable"] for row in species_rows)
        passes = (
            evaluable_species >= int(criteria["minimum_evaluable_species"])
            and int(opportunity["n_cells_A_ge_3"]) >= int(criteria["minimum_cells_A_ge_3"])
            and int(opportunity["species_with_any_shared_opportunity"])
            >= int(criteria["minimum_species_with_shared_opportunity"])
        )
        candidate = {
            "scale_km": scale_km,
            "grid": {"n_lon": grid.n_lon, "n_sinlat": grid.n_sinlat, "n_cells": grid.n_cells},
            "approximate_cell_area_km2": grid.cell_area_km2,
            "evaluable_species": evaluable_species,
            "passes_geometry_only_criteria": passes,
            "opportunity": opportunity,
            "species_geometry": species_rows,
        }
        candidates.append(candidate)
        if selected_scale is None and passes:
            selected_scale = scale_km

    return {
        "protocol": PROTOCOL,
        "status": "geometry_scale_frozen" if selected_scale is not None else "not_evaluable_no_scale_passed",
        "selection_used_image_pixels": False,
        "selection_used_flower_roi": False,
        "selection_used_continuous_colour": False,
        "selection_used_species_identity": True,
        "selection_used_coordinates": True,
        "priority_rule": geometry_contract["priority_rule"],
        "selected_primary_scale_km": selected_scale,
        "sensitivity_scales_km": [
            int(row["scale_km"])
            for row in geometry_contract["candidates"]
            if int(row["scale_km"]) != selected_scale
        ],
        "candidate_scales": candidates,
        "next_gate": (
            "lock cohort, admission and scale hashes; only then permit flower ROI extraction"
            if selected_scale is not None
            else "STOP: no colour extraction; amend geometry support rules prospectively"
        ),
    }


DISPLAY_FORBIDDEN_FIELDS = frozenset(
    {
        "species",
        "taxon",
        "taxon_id",
        "inat_taxon_id",
        "inat_taxon_name",
        "inat_genus_id",
        "family",
    }
)


def species_free_display_rows(rows: Iterable[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Strip taxon identity from future map/photo-bar rows while retaining colour/ROI."""

    allowed = (
        "display_id",
        "photo_id",
        "latitude",
        "longitude",
        "roi_thumbnail_path",
        "colour_L",
        "colour_a",
        "colour_b",
        "colour_hex",
        "photo_bar_order",
    )
    output = []
    for raw in rows:
        forbidden_present = DISPLAY_FORBIDDEN_FIELDS.intersection(raw)
        if forbidden_present and raw.get("display_contains_species") is True:
            raise ValueError(f"species-free display was explicitly opened: {sorted(forbidden_present)}")
        display = {key: raw[key] for key in allowed if key in raw}
        if not {"display_id", "photo_id"}.issubset(display):
            raise ValueError("display rows require display_id and photo_id")
        output.append(display)
    return output

"""Outcome-blind occupancy and matching primitives for RGFCA G5.

This module contains no flower-colour code.  It implements only the prospective
metadata geometry and matched-control rules frozen for G5 sympatry colour
assembly.
"""
from __future__ import annotations

import hashlib
import math
from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd

N_LON = 320
N_SINLAT = 160
N_CELLS = N_LON * N_SINLAT


def equal_area_cell_id(latitude: float, longitude: float) -> int:
    """Return the frozen 320 x 160 longitude x sin(latitude) cell id."""
    lat = float(latitude)
    lon = float(longitude)
    if not math.isfinite(lat) or not math.isfinite(lon):
        raise ValueError("coordinates must be finite")
    if lat < -90.0 or lat > 90.0 or lon < -180.0 or lon > 180.0:
        raise ValueError("coordinates outside geographic bounds")
    wrapped = (lon + 180.0) % 360.0
    x = int(math.floor(wrapped / 360.0 * N_LON))
    s = math.sin(math.radians(lat))
    y = int(math.floor((s + 1.0) * 0.5 * N_SINLAT))
    x = min(max(x, 0), N_LON - 1)
    y = min(max(y, 0), N_SINLAT - 1)
    return y * N_LON + x


def cell_indices(cell_id: int) -> tuple[int, int]:
    cell = int(cell_id)
    if cell < 0 or cell >= N_CELLS:
        raise ValueError("cell_id outside frozen grid")
    return cell % N_LON, cell // N_LON


def cell_center(cell_id: int) -> tuple[float, float]:
    """Return latitude, longitude of the equal-area cell centre."""
    x, y = cell_indices(cell_id)
    lon = -180.0 + (x + 0.5) * 360.0 / N_LON
    sinlat = -1.0 + (y + 0.5) * 2.0 / N_SINLAT
    lat = math.degrees(math.asin(min(1.0, max(-1.0, sinlat))))
    return lat, lon


def dilate_cells(cells: Sequence[int] | set[int]) -> set[int]:
    """One-cell 8-neighbour Chebyshev dilation with longitude wrapping."""
    out: set[int] = set()
    for raw in cells:
        x, y = cell_indices(int(raw))
        for dy in (-1, 0, 1):
            yy = y + dy
            if yy < 0 or yy >= N_SINLAT:
                continue
            for dx in (-1, 0, 1):
                xx = (x + dx) % N_LON
                out.add(yy * N_LON + xx)
    return out


def occupancy_overlap(a: Sequence[int] | set[int], b: Sequence[int] | set[int]) -> tuple[int, float]:
    aa = set(int(x) for x in a)
    bb = set(int(x) for x in b)
    if not aa or not bb:
        return 0, 0.0
    shared = len(aa & bb)
    union = len(aa | bb)
    return shared, float(shared / union) if union else 0.0


def is_sympatric(a: Sequence[int] | set[int], b: Sequence[int] | set[int]) -> tuple[bool, int, float]:
    shared, jaccard = occupancy_overlap(a, b)
    return bool(shared >= 5 and jaccard >= 0.10), shared, jaccard


def is_allopatric_after_dilation(a: Sequence[int] | set[int], b: Sequence[int] | set[int]) -> bool:
    return len(dilate_cells(a) & dilate_cells(b)) == 0


def taxonomic_distance_class(
    focal_genus: str,
    focal_family: str,
    partner_genus: str,
    partner_family: str,
) -> str:
    fg = str(focal_genus).strip().casefold()
    ff = str(focal_family).strip().casefold()
    pg = str(partner_genus).strip().casefold()
    pf = str(partner_family).strip().casefold()
    if not fg or not ff or not pg or not pf:
        raise ValueError("genus and family must be non-empty")
    if fg == pg:
        return "congeneric"
    if ff == pf:
        return "same_family_noncongeneric"
    return "different_family"


def deterministic_tie_key(focal: str, partner: str, control: str) -> str:
    return hashlib.sha256(f"{focal}|{partner}|{control}".encode("utf-8")).hexdigest()


def build_matched_sets(
    species_metadata: pd.DataFrame,
    occupancy: Mapping[str, set[int]],
    *,
    maximum_distance: float = 1.5,
    controls_target: int = 5,
    controls_minimum: int = 2,
) -> pd.DataFrame:
    """Build directed focal x sympatric-partner matched sets without colour data."""
    required = {
        "species",
        "accepted_genus",
        "accepted_family",
        "dominant_realm",
        "z_bio01",
        "z_bio04",
        "z_bio12",
        "z_bio15",
        "z_log1p_occupied_cells",
        "z_log1p_gbif_occurrences",
        "matching_eligible",
    }
    missing = sorted(required - set(species_metadata.columns))
    if missing:
        raise ValueError(f"species metadata lacks required columns: {missing}")
    work = species_metadata.loc[species_metadata["matching_eligible"].astype(bool)].copy()
    work["species"] = work["species"].astype(str)
    work = work.sort_values("species", kind="mergesort").drop_duplicates("species", keep="first")
    by_species = {str(row.species): row for row in work.itertuples(index=False)}
    labels = sorted(set(by_species) & set(occupancy))
    rows: list[dict[str, object]] = []
    set_index = 0

    for focal in labels:
        f = by_species[focal]
        f_cells = occupancy[focal]
        for partner in labels:
            if partner == focal:
                continue
            p = by_species[partner]
            sympatric, shared, jaccard = is_sympatric(f_cells, occupancy[partner])
            if not sympatric:
                continue
            tax_class = taxonomic_distance_class(
                f.accepted_genus, f.accepted_family, p.accepted_genus, p.accepted_family
            )
            partner_realm = str(p.dominant_realm)
            partner_vector = np.asarray(
                [
                    float(p.z_bio01) - float(f.z_bio01),
                    float(p.z_bio04) - float(f.z_bio04),
                    float(p.z_bio12) - float(f.z_bio12),
                    float(p.z_bio15) - float(f.z_bio15),
                    float(p.z_log1p_occupied_cells),
                    float(p.z_log1p_gbif_occurrences),
                ],
                dtype=float,
            )
            if not np.isfinite(partner_vector).all():
                continue
            candidates: list[tuple[float, str, str]] = []
            for control in labels:
                if control in {focal, partner}:
                    continue
                c = by_species[control]
                if str(c.dominant_realm) != partner_realm:
                    continue
                control_tax = taxonomic_distance_class(
                    f.accepted_genus, f.accepted_family, c.accepted_genus, c.accepted_family
                )
                if control_tax != tax_class:
                    continue
                if not is_allopatric_after_dilation(f_cells, occupancy[control]):
                    continue
                control_vector = np.asarray(
                    [
                        float(c.z_bio01) - float(f.z_bio01),
                        float(c.z_bio04) - float(f.z_bio04),
                        float(c.z_bio12) - float(f.z_bio12),
                        float(c.z_bio15) - float(f.z_bio15),
                        float(c.z_log1p_occupied_cells),
                        float(c.z_log1p_gbif_occurrences),
                    ],
                    dtype=float,
                )
                if not np.isfinite(control_vector).all():
                    continue
                distance = float(np.linalg.norm(control_vector - partner_vector))
                if distance <= float(maximum_distance):
                    candidates.append((distance, deterministic_tie_key(focal, partner, control), control))
            candidates.sort(key=lambda item: (item[0], item[1], item[2]))
            selected = candidates[: int(controls_target)]
            if len(selected) < int(controls_minimum):
                continue
            set_index += 1
            set_id = f"g5set_{set_index:08d}"
            for rank, (distance, tie_key, control) in enumerate(selected, start=1):
                rows.append(
                    {
                        "set_id": set_id,
                        "focal_species": focal,
                        "sympatric_partner": partner,
                        "taxonomic_distance_class": tax_class,
                        "partner_dominant_realm": partner_realm,
                        "shared_cells": int(shared),
                        "occupancy_jaccard": float(jaccard),
                        "control_rank": int(rank),
                        "control_species": control,
                        "matching_distance": float(distance),
                        "tie_break_sha256": tie_key,
                    }
                )
    columns = [
        "set_id",
        "focal_species",
        "sympatric_partner",
        "taxonomic_distance_class",
        "partner_dominant_realm",
        "shared_cells",
        "occupancy_jaccard",
        "control_rank",
        "control_species",
        "matching_distance",
        "tie_break_sha256",
    ]
    return pd.DataFrame(rows, columns=columns)


__all__ = [
    "N_LON",
    "N_SINLAT",
    "N_CELLS",
    "equal_area_cell_id",
    "cell_indices",
    "cell_center",
    "dilate_cells",
    "occupancy_overlap",
    "is_sympatric",
    "is_allopatric_after_dilation",
    "taxonomic_distance_class",
    "deterministic_tie_key",
    "build_matched_sets",
]

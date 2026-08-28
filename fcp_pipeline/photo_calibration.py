"""Build and validate the blinded 480-photo measurement-calibration package."""

from __future__ import annotations

import hashlib
from typing import Iterable

import pandas as pd

CALIBRATION_PACKAGE_VERSION = "jbi-ch1-photo-calibration-v1"

VISIBILITY_VALUES = {"", "evaluable", "not_evaluable"}
VISIBILITY_FAILURE_VALUES = {
    "",
    "occlusion",
    "distance",
    "blur",
    "overexposure",
    "underexposure",
    "non_target_organ",
    "insufficient_flower_area",
    "other",
}
SEGMENTATION_VALUES = {"", "ok", "failed", "not_applicable"}
COLOUR_ASSIGNMENT_VALUES = {"", "resolved", "unresolved", "not_applicable"}

_LOCATION_OR_OUTCOME_HINTS = {
    "latitude",
    "longitude",
    "lat",
    "lon",
    "lng",
    "observed_on",
    "date",
    "observer",
    "user_login",
    "place_guess",
    "locality",
    "location",
    "flower_colour",
    "flower_color",
    "colour_state",
    "color_state",
    "boundary",
    "environment",
    "bioclim",
}

MEDIA_URL_CANDIDATES = (
    "photo_url",
    "image_url",
    "media_url",
    "medium_url",
    "large_url",
    "original_url",
)


def _norm(name: object) -> str:
    return str(name).strip().lower().replace("-", "_").replace(" ", "_")


def _blind_id(species: str, photo_id: str) -> str:
    payload = f"jbi-ch1-calibration-v1\x1f{species}\x1f{photo_id}".encode("utf-8")
    return hashlib.sha256(payload).hexdigest()[:16]


def default_display_columns(frame: pd.DataFrame, *, species_col: str, photo_id_col: str) -> list[str]:
    columns = [species_col, photo_id_col]
    mapping = {_norm(c): str(c) for c in frame.columns}
    for candidate in MEDIA_URL_CANDIDATES:
        if candidate in mapping:
            columns.append(mapping[candidate])
            break
    return columns


def build_calibration_sheet(
    frozen_split: pd.DataFrame,
    *,
    species_col: str = "species",
    photo_id_col: str = "photo_id",
    keep_columns: Iterable[str] | None = None,
) -> pd.DataFrame:
    required = {species_col, photo_id_col, "split"}
    if not required.issubset(frozen_split.columns):
        raise ValueError(f"split missing columns: {sorted(required - set(frozen_split.columns))}")

    calibration = frozen_split.loc[frozen_split["split"].astype(str) == "calibration"].copy()
    if len(calibration) != 480:
        raise ValueError(f"expected 480 calibration rows, found {len(calibration)}")
    counts = calibration.groupby(species_col).size()
    if len(counts) != 6 or not (counts == 80).all():
        raise ValueError("calibration sheet must contain exactly 80 rows for each of six species")

    if keep_columns is None:
        keep = default_display_columns(
            calibration, species_col=species_col, photo_id_col=photo_id_col
        )
    else:
        keep = list(keep_columns)
        missing = [c for c in keep if c not in calibration.columns]
        if missing:
            raise ValueError(f"requested display columns missing: {missing}")
        forbidden = [c for c in keep if _norm(c) in _LOCATION_OR_OUTCOME_HINTS]
        if forbidden:
            raise ValueError(
                "calibration display columns may not expose location/date/observer/outcome fields: "
                + ", ".join(forbidden)
            )
        for critical in (species_col, photo_id_col):
            if critical not in keep:
                keep.insert(0, critical)

    out = calibration[keep].copy()
    out.insert(
        0,
        "blind_id",
        [
            _blind_id(str(s).strip(), str(p).strip())
            for s, p in out[[species_col, photo_id_col]].itertuples(index=False, name=None)
        ],
    )
    out["flower_visibility"] = ""
    out["visibility_failure_code"] = ""
    out["segmentation_status"] = ""
    out["segmentation_failure_code"] = ""
    out["colour_assignment"] = ""
    out["colour_state"] = ""
    out["notes"] = ""

    # Hash order prevents source/geographic ordering from leaking into review order.
    out = out.sort_values([species_col, "blind_id"], kind="mergesort").reset_index(drop=True)
    return out


def validate_calibration_sheet(
    frame: pd.DataFrame,
    *,
    species_col: str = "species",
    photo_id_col: str = "photo_id",
    require_complete: bool = False,
) -> None:
    required = {
        "blind_id",
        species_col,
        photo_id_col,
        "flower_visibility",
        "visibility_failure_code",
        "segmentation_status",
        "segmentation_failure_code",
        "colour_assignment",
        "colour_state",
        "notes",
    }
    if not required.issubset(frame.columns):
        raise ValueError(f"calibration sheet missing fields: {sorted(required - set(frame.columns))}")
    if len(frame) != 480:
        raise ValueError("calibration sheet must have exactly 480 rows")
    if frame[photo_id_col].astype(str).duplicated().any():
        raise ValueError("calibration sheet contains duplicate photo IDs")
    counts = frame.groupby(species_col).size()
    if len(counts) != 6 or not (counts == 80).all():
        raise ValueError("calibration sheet must have six species with 80 rows each")

    visibility = frame["flower_visibility"].fillna("").astype(str).str.strip()
    bad = sorted(set(visibility) - VISIBILITY_VALUES)
    if bad:
        raise ValueError(f"invalid flower_visibility values: {bad}")

    failure = frame["visibility_failure_code"].fillna("").astype(str).str.strip()
    bad = sorted(set(failure) - VISIBILITY_FAILURE_VALUES)
    if bad:
        raise ValueError(f"invalid visibility_failure_code values: {bad}")

    segmentation = frame["segmentation_status"].fillna("").astype(str).str.strip()
    bad = sorted(set(segmentation) - SEGMENTATION_VALUES)
    if bad:
        raise ValueError(f"invalid segmentation_status values: {bad}")

    assignment = frame["colour_assignment"].fillna("").astype(str).str.strip()
    bad = sorted(set(assignment) - COLOUR_ASSIGNMENT_VALUES)
    if bad:
        raise ValueError(f"invalid colour_assignment values: {bad}")

    colour_state = frame["colour_state"].fillna("").astype(str).str.strip()

    if require_complete:
        if (visibility == "").any():
            raise ValueError("complete calibration requires flower_visibility for every row")
        not_eval = visibility == "not_evaluable"
        if (not_eval & (failure == "")).any():
            raise ValueError("not_evaluable rows require visibility_failure_code")
        if ((visibility == "evaluable") & (segmentation == "")).any():
            raise ValueError("evaluable rows require segmentation_status")
        if (not_eval & (segmentation != "not_applicable")).any():
            raise ValueError("not_evaluable rows must set segmentation_status=not_applicable")
        seg_ok = segmentation == "ok"
        if (seg_ok & (assignment == "")).any():
            raise ValueError("segmentation_status=ok requires colour_assignment")
        if ((assignment == "resolved") & (colour_state == "")).any():
            raise ValueError("resolved colour assignments require colour_state")
        if ((assignment != "resolved") & (colour_state != "")).any():
            raise ValueError("only resolved assignments may carry a colour_state")
        if ((segmentation == "failed") & (assignment != "not_applicable")).any():
            raise ValueError("failed segmentation requires colour_assignment=not_applicable")


def calibration_summary(frame: pd.DataFrame, *, species_col: str = "species") -> pd.DataFrame:
    validate_calibration_sheet(frame, species_col=species_col, require_complete=False)
    rows: list[dict[str, object]] = []
    for species, group in frame.groupby(species_col, sort=True):
        visibility = group["flower_visibility"].fillna("").astype(str).str.strip()
        segmentation = group["segmentation_status"].fillna("").astype(str).str.strip()
        assignment = group["colour_assignment"].fillna("").astype(str).str.strip()
        rows.append(
            {
                species_col: species,
                "n_total": int(len(group)),
                "n_visibility_completed": int((visibility != "").sum()),
                "n_evaluable": int((visibility == "evaluable").sum()),
                "n_not_evaluable": int((visibility == "not_evaluable").sum()),
                "n_segmentation_ok": int((segmentation == "ok").sum()),
                "n_segmentation_failed": int((segmentation == "failed").sum()),
                "n_colour_resolved": int((assignment == "resolved").sum()),
                "n_colour_unresolved": int((assignment == "unresolved").sum()),
            }
        )
    return pd.DataFrame(rows)

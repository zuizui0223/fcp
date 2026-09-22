"""FCP v2 measurement-validity helpers.

Pure technical utilities only. This module deliberately contains no species-level
biological outcome logic, H2 decision logic, spatial inference, or ecological
predictors. It exists to freeze and test technical measurement channels before
biological outcomes are opened.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import io
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

EV_LEVELS = (-1.0, -0.5, 0.0, 0.5, 1.0)
NEUTRAL_RGB = (128, 128, 128)

FORBIDDEN_TECHNICAL_COLUMN_TOKENS = (
    "species",
    "taxon",
    "morph",
    "q_white",
    "white_response",
    "spatial_outcome",
    "latitude",
    "longitude",
    "observer",
    "pollinator",
    "climate",
    "h3",
)


def srgb_to_linear(rgb01: np.ndarray) -> np.ndarray:
    x = np.asarray(rgb01, dtype=np.float64)
    if np.any(~np.isfinite(x)) or np.any((x < 0.0) | (x > 1.0)):
        raise ValueError("sRGB values must be finite and within [0,1]")
    return np.where(x <= 0.04045, x / 12.92, ((x + 0.055) / 1.055) ** 2.4)


def linear_to_srgb(linear01: np.ndarray) -> np.ndarray:
    x = np.asarray(linear01, dtype=np.float64)
    if np.any(~np.isfinite(x)):
        raise ValueError("linear-light values must be finite")
    x = np.clip(x, 0.0, 1.0)
    return np.where(
        x <= 0.0031308,
        12.92 * x,
        1.055 * np.power(x, 1.0 / 2.4) - 0.055,
    )


def apply_exposure_ev(rgb: np.ndarray, ev: float) -> np.ndarray:
    """Apply a deterministic exposure shift in linear light and return uint8 sRGB."""

    x = np.asarray(rgb, dtype=np.uint8)
    if x.ndim < 2 or x.shape[-1] != 3:
        raise ValueError("RGB input must end in three channels")
    if float(ev) not in EV_LEVELS:
        raise ValueError(f"EV level {ev!r} is not frozen in {EV_LEVELS}")
    lin = srgb_to_linear(x.astype(np.float64) / 255.0)
    shifted = np.clip(lin * (2.0 ** float(ev)), 0.0, 1.0)
    srgb = linear_to_srgb(shifted)
    return np.rint(np.clip(srgb * 255.0, 0.0, 255.0)).astype(np.uint8)


def relative_luminance(rgb: np.ndarray) -> np.ndarray:
    x = np.asarray(rgb, dtype=np.uint8)
    if x.ndim < 2 or x.shape[-1] != 3:
        raise ValueError("RGB input must end in three channels")
    lin = srgb_to_linear(x.astype(np.float64) / 255.0)
    return (
        0.2126729 * lin[..., 0]
        + 0.7151522 * lin[..., 1]
        + 0.0721750 * lin[..., 2]
    )


def exposure_metrics(rgb_pixels: np.ndarray) -> dict[str, float]:
    x = np.asarray(rgb_pixels, dtype=np.uint8)
    if x.ndim != 2 or x.shape[1] != 3 or len(x) == 0:
        raise ValueError("expected non-empty N x 3 uint8 RGB pixels")
    maxima = x.max(axis=1)
    lum = relative_luminance(x)
    qs = {
        q: float(np.quantile(lum, q))
        for q in (0.01, 0.10, 0.50, 0.90, 0.99)
    }
    return {
        "clip_fraction": float(np.mean(maxima == 255)),
        "near_clip_fraction": float(np.mean(maxima >= 250)),
        "clip_fraction_r": float(np.mean(x[:, 0] == 255)),
        "clip_fraction_g": float(np.mean(x[:, 1] == 255)),
        "clip_fraction_b": float(np.mean(x[:, 2] == 255)),
        "luminance_q01": qs[0.01],
        "luminance_q10": qs[0.10],
        "luminance_q50": qs[0.50],
        "luminance_q90": qs[0.90],
        "luminance_q99": qs[0.99],
        "dynamic_range_q99_minus_q01": float(qs[0.99] - qs[0.01]),
        "black_fraction": float(np.mean(np.max(x, axis=1) == 0)),
    }


def neutralize_background(
    rgb: np.ndarray,
    flower_mask: np.ndarray,
    *,
    neutral_rgb: Sequence[int] = NEUTRAL_RGB,
) -> np.ndarray:
    image = np.asarray(rgb, dtype=np.uint8)
    mask = np.asarray(flower_mask, dtype=bool)
    if image.ndim != 3 or image.shape[2] != 3:
        raise ValueError("RGB image must be H x W x 3")
    if mask.shape != image.shape[:2]:
        raise ValueError("flower mask shape does not match image")
    neutral = np.asarray(tuple(neutral_rgb), dtype=np.int64)
    if neutral.shape != (3,) or np.any((neutral < 0) | (neutral > 255)):
        raise ValueError("neutral RGB must contain three uint8-like values")
    if tuple(int(v) for v in neutral) != NEUTRAL_RGB:
        raise ValueError(f"neutral background is frozen at {NEUTRAL_RGB}")
    out = image.copy()
    out[~mask] = neutral.astype(np.uint8)
    return out


def mask_iou(a: np.ndarray, b: np.ndarray) -> float:
    ma = np.asarray(a, dtype=bool)
    mb = np.asarray(b, dtype=bool)
    if ma.shape != mb.shape:
        raise ValueError("mask shapes differ")
    union = int(np.count_nonzero(ma | mb))
    if union == 0:
        return 1.0
    return float(np.count_nonzero(ma & mb) / union)


def mask_boundary_fraction(mask: np.ndarray) -> float:
    """Fraction of positive mask pixels lying on the 4-neighbour boundary."""

    m = np.asarray(mask, dtype=bool)
    n = int(np.count_nonzero(m))
    if n == 0:
        return float("nan")
    padded = np.pad(m, 1, constant_values=False)
    core = padded[1:-1, 1:-1]
    interior = (
        core
        & padded[:-2, 1:-1]
        & padded[2:, 1:-1]
        & padded[1:-1, :-2]
        & padded[1:-1, 2:]
    )
    boundary = core & ~interior
    return float(np.count_nonzero(boundary) / n)


def jitter_box_xyxy(
    box: Sequence[float],
    *,
    dx_fraction: float = 0.0,
    dy_fraction: float = 0.0,
    scale_fraction: float = 0.0,
) -> tuple[float, float, float, float]:
    """Deterministically perturb one XYXY box without image-boundary clipping."""

    if dx_fraction not in (-0.05, 0.0, 0.05):
        raise ValueError("dx_fraction must be one of -0.05, 0, 0.05")
    if dy_fraction not in (-0.05, 0.0, 0.05):
        raise ValueError("dy_fraction must be one of -0.05, 0, 0.05")
    if scale_fraction not in (-0.10, 0.0, 0.10):
        raise ValueError("scale_fraction must be one of -0.10, 0, 0.10")
    x0, y0, x1, y1 = map(float, box)
    w, h = x1 - x0, y1 - y0
    if not (w > 0 and h > 0):
        raise ValueError("box must have positive width and height")
    cx = (x0 + x1) / 2.0 + dx_fraction * w
    cy = (y0 + y1) / 2.0 + dy_fraction * h
    scale = 1.0 + scale_fraction
    nw, nh = w * scale, h * scale
    return (cx - nw / 2.0, cy - nh / 2.0, cx + nw / 2.0, cy + nh / 2.0)


def frozen_prompt_jitter_set(box: Sequence[float]) -> tuple[tuple[float, float, float, float], ...]:
    """Return the exact deterministic v2 prompt-jitter set.

    We perturb one axis or one scale component at a time plus the unperturbed
    box, avoiding a large combinatorial grid that would be expensive and harder
    to audit.
    """

    params = (
        (0.0, 0.0, 0.0),
        (-0.05, 0.0, 0.0),
        (+0.05, 0.0, 0.0),
        (0.0, -0.05, 0.0),
        (0.0, +0.05, 0.0),
        (0.0, 0.0, -0.10),
        (0.0, 0.0, +0.10),
    )
    return tuple(
        jitter_box_xyxy(box, dx_fraction=dx, dy_fraction=dy, scale_fraction=scale)
        for dx, dy, scale in params
    )


def validate_response_blind_columns(columns: Iterable[str]) -> None:
    names = [str(c) for c in columns]
    leaked = [
        name
        for name in names
        if any(token in name.casefold() for token in FORBIDDEN_TECHNICAL_COLUMN_TOKENS)
    ]
    if leaked:
        raise ValueError(f"technical table leaked biological/outcome columns: {leaked}")


def canonical_technical_table_sha256(frame: pd.DataFrame) -> str:
    validate_response_blind_columns(frame.columns)
    if "measurement_id" not in frame:
        raise ValueError("technical table lacks measurement_id")
    if frame["measurement_id"].astype(str).nunique() != len(frame):
        raise ValueError("measurement_id must be unique")
    ordered = frame.sort_values("measurement_id", kind="mergesort").reset_index(drop=True)
    buffer = io.StringIO()
    ordered.to_csv(buffer, index=False, lineterminator="\n")
    return hashlib.sha256(buffer.getvalue().encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class TechnicalSeal:
    rows: int
    unique_measurement_ids: int
    table_sha256: str
    biological_outcomes_opened: bool = False


def seal_technical_table(frame: pd.DataFrame) -> TechnicalSeal:
    digest = canonical_technical_table_sha256(frame)
    return TechnicalSeal(
        rows=int(len(frame)),
        unique_measurement_ids=int(frame["measurement_id"].astype(str).nunique()),
        table_sha256=digest,
        biological_outcomes_opened=False,
    )


__all__ = [
    "EV_LEVELS",
    "NEUTRAL_RGB",
    "TechnicalSeal",
    "apply_exposure_ev",
    "canonical_technical_table_sha256",
    "exposure_metrics",
    "frozen_prompt_jitter_set",
    "jitter_box_xyxy",
    "mask_boundary_fraction",
    "mask_iou",
    "neutralize_background",
    "relative_luminance",
    "seal_technical_table",
    "srgb_to_linear",
    "linear_to_srgb",
    "validate_response_blind_columns",
]

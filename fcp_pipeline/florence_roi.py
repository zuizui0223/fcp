"""Pure geometry helpers for selecting one Florence flower ROI.

The selector is intentionally independent of colour pixels, species identity, candidate
states, geography, or prior labels. It uses only Florence-proposed bounding boxes and
image geometry, so changing colour outcomes cannot change which ROI is selected.
"""
from __future__ import annotations

import math
from typing import Iterable, Sequence


def box_area_fraction(box: Sequence[float], image_size: tuple[int, int]) -> float:
    """Return box area / image area after clipping negative width/height to zero."""
    if len(box) != 4:
        raise ValueError("box must contain four coordinates")
    width, height = image_size
    if width <= 0 or height <= 0:
        raise ValueError("image dimensions must be positive")
    x0, y0, x1, y1 = (float(v) for v in box)
    area = max(0.0, x1 - x0) * max(0.0, y1 - y0)
    return area / float(width * height)


def normalized_center_distance_sq(
    box: Sequence[float], image_size: tuple[int, int]
) -> float:
    """Squared box-centre distance from image centre, normalized by image dimensions."""
    width, height = image_size
    if width <= 0 or height <= 0:
        raise ValueError("image dimensions must be positive")
    x0, y0, x1, y1 = (float(v) for v in box)
    bx = (x0 + x1) / 2.0
    by = (y0 + y1) / 2.0
    return ((bx - width / 2.0) / width) ** 2 + ((by - height / 2.0) / height) ** 2


def choose_flower_box(
    boxes: Iterable[Sequence[float]], image_size: tuple[int, int]
) -> list[float] | None:
    """Select one plausible target-flower box from Florence proposals.

    Rules are frozen from geometry only:

    * no proposals -> None;
    * exactly one proposal -> keep it, even if it occupies much of the frame (important
      for close-up flowers such as the pilot Ipomoea image);
    * multiple proposals -> prefer individual-object boxes occupying 0.1%--20% of the
      image; if none exist, relax the upper bound to 35%; if still none exist, retain
      all proposals rather than silently fabricating a crop;
    * among retained proposals, balance object size against distance from image centre.

    The 20% filter exists specifically to prevent a scene/group box from beating many
    individual flower boxes merely because it has larger area. No colour/state result
    participates in this choice.
    """
    normalized: list[list[float]] = []
    for raw in boxes:
        if len(raw) != 4:
            continue
        box = [float(v) for v in raw]
        if box[2] <= box[0] or box[3] <= box[1]:
            continue
        normalized.append(box)

    if not normalized:
        return None
    if len(normalized) == 1:
        return normalized[0]

    primary = [
        box
        for box in normalized
        if 0.001 <= box_area_fraction(box, image_size) <= 0.20
    ]
    relaxed = [
        box
        for box in normalized
        if 0.001 <= box_area_fraction(box, image_size) <= 0.35
    ]
    candidates = primary or relaxed or normalized

    def score(box: Sequence[float]) -> float:
        area_fraction = max(box_area_fraction(box, image_size), 1e-12)
        centre_distance = normalized_center_distance_sq(box, image_size)
        return math.sqrt(area_fraction) / (1.0 + 6.0 * centre_distance)

    return list(max(candidates, key=score))

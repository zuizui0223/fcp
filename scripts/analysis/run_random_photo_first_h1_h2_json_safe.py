#!/usr/bin/env python3
"""JSON-safe launcher for the prospectively frozen H1 -> H2 runner.

This wrapper does not alter sampling, statistics, seeds, nulls, climate inputs, or
scientific decisions. It converts non-finite not-evaluable/QC scalar placeholders
to JSON null and invokes the fail-closed H2 orchestration layer.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np

import run_random_photo_first_h1_h2 as frozen
import run_random_photo_first_h1_h2_failclosed as failclosed


_ORIGINAL_DUMPS = json.dumps


def json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return [json_safe(item) for item in value.tolist()]
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        number = float(value)
        return number if math.isfinite(number) else None
    return value


def safe_dumps(value: Any, *args: Any, **kwargs: Any) -> str:
    return _ORIGINAL_DUMPS(json_safe(value), *args, **kwargs)


def safe_write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        _ORIGINAL_DUMPS(json_safe(payload), indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    frozen.write_json = safe_write_json
    # frozen imports the standard json module, so replace only its dumps entrypoint;
    # _ORIGINAL_DUMPS remains an immutable reference to the real serializer.
    frozen.json.dumps = safe_dumps
    return int(failclosed.main())


if __name__ == "__main__":
    raise SystemExit(main())

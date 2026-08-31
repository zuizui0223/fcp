#!/usr/bin/env python3
"""Validate the pre-image dated-source amendment."""

from __future__ import annotations

import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fcp_pipeline.atlas_dated_source import validate_dated_source_amendment


def main() -> None:
    path = ROOT / "docs/supporting/jbi_atlas_dated_source_amendment_v1.json"
    amendment = json.loads(path.read_text(encoding="utf-8"))
    validate_dated_source_amendment(amendment)
    print(json.dumps({"status": "pass", "protocol": amendment["protocol"]}))


if __name__ == "__main__":
    main()

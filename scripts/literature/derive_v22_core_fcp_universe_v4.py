#!/usr/bin/env python3
"""Final conservative display-core guard for the all-species FCP universe.

This wrapper hardens v3 against a specific false-positive path: papers whose focal trait
is stigma/gynoecium/anther/pollen colour can mention generic `flower colour polymorphism`
in background text. Such generic wording must not rescue them into the primary display
core. A non-display target is rescued only by explicit petal/corolla/perianth/tepal/
labellum/bract/inflorescence colour evidence. Expanded-universe membership is unchanged.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import derive_v22_core_fcp_universe_v3 as v3

base = v3.base

# Second independent guard. The inherited direct_reason() first checks NONDISPLAY_ONLY_RE
# against DISPLAY_STRUCTURE_RE. Restrict the latter to actual display structures, not the
# generic phrase `flower colour polymorphism`, which can appear only as background text.
base.DISPLAY_STRUCTURE_RE = re.compile(
    r"(?is)\b(?:petal|corolla|perianth|tepal|labellum|bract|inflorescence)\b.{0,35}\bcolou?r\b"
    r"|\bcolou?r\b.{0,35}\b(?:petal|corolla|perianth|tepal|labellum|bract|inflorescence)\b"
)

if __name__ == "__main__":
    base.main()
    # Version the realized core explicitly after the inherited writer finishes.
    # This changes metadata only, not membership or C/S evidence.
    import sys
    outdir = Path(sys.argv[sys.argv.index("--outdir") + 1])
    summary_path = outdir / "v22_core_fcp_universe_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["core_protocol_version"] = "display-core-v4-final"
    summary["nondisplay_guard"] = (
        "sexual-organ-target studies cannot enter the primary core through generic "
        "background flower-colour wording; explicit display-structure colour evidence is required"
    )
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

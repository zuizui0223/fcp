#!/usr/bin/env python3
"""Focal-taxon consistency guard for the primary display FCP core.

Taxonomic delimitation/circumscription papers can discuss multiple congeners. A direct
colour-polymorphism phrase belonging to one congener must not be transferred to the
species named as the focal taxon when that focal taxon is explicitly described as all
one flower colour. This wrapper keeps the v5 display rules and rejects that high-
specificity contradiction. C/S outcomes are not consulted.
"""
from __future__ import annotations

import json
import re
from pathlib import Path
import sys

import derive_v22_core_fcp_universe_v5 as v5

base = v5.base

TAXONOMIC_COMPARISON_TITLE_RE = re.compile(
    r"(?is)\b(?:species\s+delimitation|circumscription|subdivision|taxonom\w*|species\s+complex)\b"
)
FOCAL_MONOCHROME_RE = re.compile(
    r"(?is)\b(?:plants?\s+of\s+)?(?:[A-Z][a-z]+|[A-Z]\.)\s+[a-z][a-z-]+(?:\s+s\.\s*str\.)?"
    r".{0,180}\b(?:were|are|was|is)\s+all\s+"
    r"(?:white|yellow|pink|purple|blue|red|orange|green|cream|violet|magenta)[- ]flowered\b"
)

_original = base.direct_reason

def direct_reason(text: str, title: str) -> str:
    if TAXONOMIC_COMPARISON_TITLE_RE.search(title or "") and FOCAL_MONOCHROME_RE.search(text or ""):
        return ""
    return _original(text, title)

base.direct_reason = direct_reason

if __name__ == "__main__":
    base.main()
    outdir = Path(sys.argv[sys.argv.index("--outdir") + 1])
    p = outdir / "v22_core_fcp_universe_summary.json"
    s = json.loads(p.read_text(encoding="utf-8"))
    s["core_protocol_version"] = "display-core-v6-focal-consistent"
    s["focal_taxon_guard"] = (
        "taxonomic comparison sources cannot transfer another species' colour polymorphism "
        "to a focal taxon explicitly described as all one flower colour"
    )
    p.write_text(json.dumps(s, indent=2) + "\n", encoding="utf-8")

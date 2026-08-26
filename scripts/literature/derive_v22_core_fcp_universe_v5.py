#!/usr/bin/env python3
"""Final display-core rule after auditing false exclusions.

Use the target-aware v3 guard (which excludes titles explicitly targeting stigma/
gynoecium/anther/pollen colour) and add a direct display-colour variation rule so that
sources explicitly documenting natural/intraspecific flower-colour variation are not
lost merely because they use `variation` instead of `polymorphism` or `morph` in the
relevant wording. C/S positivity is never used for membership.
"""
from __future__ import annotations

import json
import re
from pathlib import Path
import sys

import derive_v22_core_fcp_universe_v3 as v3

base = v3.base

# General, outcome-independent display-FCP evidence. The expanded universe has already
# passed primary/natural/intraspecific/discrete eligibility. This rule only restores the
# explicit wording `flower colour variation` that the narrower morph/polymorphism rules
# missed (e.g. Aquilegia/Eruca-style wording).
base.DIRECT_PATTERNS = [
    ("explicit_flower_colour_variation", re.compile(
        r"(?is)\b(?:natural\s+|intraspecific\s+|within[- ]population\s+|geographic\w*\s+|spatial\s+)?"
        r"(?:variation|differentiation)\s+in\s+(?:flower|floral|petal|corolla|perianth|tepal|labellum|bract|inflorescence)\s+colou?r\b"
        r"|\b(?:flower|floral|petal|corolla|perianth|tepal|labellum|bract|inflorescence)\s+colou?r\s+(?:variation|differentiation)\b"
    )),
    *base.DIRECT_PATTERNS,
]

if __name__ == "__main__":
    base.main()
    outdir = Path(sys.argv[sys.argv.index("--outdir") + 1])
    p = outdir / "v22_core_fcp_universe_summary.json"
    s = json.loads(p.read_text(encoding="utf-8"))
    s["core_protocol_version"] = "display-core-v5-final"
    s["core_boundary"] = (
        "primary core excludes source sets supported only by sexual-organ colour targets "
        "or review/synthesis, but retains explicit natural discrete flower-display colour variation"
    )
    p.write_text(json.dumps(s, indent=2) + "\n", encoding="utf-8")
